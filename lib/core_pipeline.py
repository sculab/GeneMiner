#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:48
# @Author  : xiepulin
# @File    : core_pipeline.py
# @Software: PyCharm

import argparse
import sys
import subprocess
import datetime
import os
import logging
import traceback
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from concurrent import futures
import re
import multiprocessing
from Bio import pairwise2
from Bio.Seq import Seq
import pandas as pd
from openpyxl import load_workbook
from basic import *
from verify_parameters import *
from build_reference_database import *




#############################################
#############################################
'''
第四部分
挖掘基因核心流程

(1)filter         过滤reads
(2)assemble       组装reads
(3)check assembled_out  校验contigs
(4)输出结果信息
'''
################################################
################################################



class Pipeline():
    def __init__(self, configuration_information,type,out_dir_name,thread_number, kmer, wordsize, max_length,min_length,options):

        self.configuration_information=configuration_information
        self.type = type
        self.out_dir_name=out_dir_name
        self.thread_number = thread_number
        self.wordsize=wordsize
        self.kmer = kmer
        self.max_length=max_length
        self.min_length=min_length
        self.options=options



        self.my_software_name=self.configuration_information["my_software_name"]
        self.filter_software=self.configuration_information["filter_path"]    #filter_reads.pl / filter  绝对路径
        self.assemble_software=self.configuration_information["assemble_path"]  #minia                   绝对路径
        self.reference_database=self.configuration_information["reference_database"] #参考基因数据库
        self.filtered_out=self.configuration_information["filtered_out"]
        self.assembled_out=self.configuration_information["assembled_out"]
        self.assembled_log=self.configuration_information["assembled_log"]
        self.GM_results=self.configuration_information["GM_results"]
        self.system=self.configuration_information["system"]
        self.whole_log = self.configuration_information["whole_log"]
        self.results_information_excel=self.configuration_information["results_information_excel"] #"results_information.xlsx"





    '''
    1 多线程并行的方法过滤目标基因。产生两个过滤文件Filtered_reads__R1.fastq ， Filtered_reads__R2.fastq ,并且合并双端数据为filtered.fq
    '''

    def filter_reads(self, file):
        # file="ycf4.fasta""
        out_dir_name=self.out_dir_name
        type = str(self.type)  # tfa,tgb, cp ,target 中的某一种
        filter_software = self.filter_software  # lib中的filter  #绝对路径
        filtered_out=self.filtered_out  #filter_out  二级目录
        reference_database=self.reference_database
        gene_name = file.split(".fasta")[0] #ycf4
        system=self.system
        whole_log = self.whole_log
        wordsize=self.wordsize   #filter 的kmer


        ref = file
        dir = gene_name  # ycf4.fasta ---- ycf4
        path = os.path.join(out_dir_name, filtered_out, dir)  #out_dir/filtered_out/ycf4
        dir_make(path)

        data1_path = os.path.join(out_dir_name, "data1.fq")
        data2_path = os.path.join(out_dir_name, "data2.fq")
        reference_path = os.path.join(out_dir_name, reference_database, ref)
        whole_log_path = os.path.join(out_dir_name, whole_log)

        if is_exist(data1_path) and is_exist(data2_path):
            path_list = [path,filter_software,data1_path,data2_path,reference_path,whole_log_path]
            [path, filter_software, data1_path, data2_path, reference_path, whole_log_path] = get_absolute_and_map_path(path_list, system)  # 如果为windows环境，会批量映射路径


            cmd = "cd '{0}' && '{1}' -1 '{2}' -2 '{3}' -r '{4}' -k {5} >/dev/null  2>&1".format(path,
                  filter_software,data1_path,data2_path,reference_path,wordsize)
            runCommand(cmd,system)
            cmd = "cd '{0}' && cat {1} {2} > {3}".format(path, "Filtered_reads__R1.fastq", "Filtered_reads__R2.fastq",
                                                       "filtered.fq")
            runCommand(cmd,system)
            # 立刻将完成的基因打印
            # 立刻将完成的基因打印
            message = "{0}_gene: '{1}' has been filtered. DONE".format(type, gene_name)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)

        else:
            pass

    def wrap_filter_reads(self, file):
        self.filter_reads(file)

    def filter_reads_parallel(self):
        thread_number=self.thread_number
        type = str(self.type)  # mito , cp , tfa, tgb 中的某一种
        out_dir_name=self.out_dir_name  # 定位所有基因名，所有传参都是根据基因名
        reference_database=self.reference_database
        filtered_out=self.filtered_out

        reference_database_path= os.path.join(out_dir_name,reference_database)
        files = os.listdir(reference_database_path)  # ycf4.fasta
        filtered_out_path = os.path.join(out_dir_name,filtered_out)  #out/filter_out
        dir_make(filtered_out_path)

        executor = futures.ThreadPoolExecutor(max_workers=thread_number)
        task_pool = []
        result = []

        for file in files:
            task_pool.append(executor.submit(self.wrap_filter_reads, file))

        for task in tqdm(desc="{0:<22}".format("Filtering_reads"), iterable=futures.as_completed(task_pool),
                         total=len(task_pool)):
            result.append(task.result())

        executor.shutdown()  # 所有线程结束，再进行下一步


    '''
    2 多线程并行的方法，拼接目标基因  产生三个文件：assembled_out.contigs.fa  assembled_out.h5  assembled_out.unitigs.fa
    '''

    def assemble_reads(self, file, kmer):
        # file=matk
        # gene_name=matk

        type = str(self.type)
        assemble_software=self.assemble_software  #minia   #绝对路径
        out_dir_name=self.out_dir_name
        filterd_out=self.filtered_out            #二级目录
        assembled_out=self.assembled_out         #二级目录
        assembled_log=self.assembled_log         #二级目录
        system=self.system
        whole_log=self.whole_log                 #二级文件



        input = os.path.join(out_dir_name, filterd_out, file, "filtered.fq")
        whole_log_path = os.path.join(out_dir_name, whole_log)

        if is_exist(input) == 0:
            path_list=[whole_log_path]
            [whole_log_path]=get_absolute_and_map_path(path_list,system)# 如果为windows环境，会批量映射路径
            message = "{} :Assembling reads has failed.  It could be that the quality of the raw data is not good enough, the amount of input data is too small or the species  are not closely related".format(
                file)
            cmd = "echo {0}>>'{1}' ".format(message,whole_log_path)
            runCommand(cmd,system)
            return 0

        dir = os.path.join(out_dir_name, assembled_out, file)
        dir_make(dir)
        log = os.path.join(out_dir_name, assembled_log, file + "_log.txt")
        path_list= [dir,assemble_software, input, log,whole_log_path]
        [dir,assemble_software, input, log,whole_log_path]=get_absolute_and_map_path(path_list,system)  # 如果为windows环境，会批量映射路径
        try:
            cmd = "cd '{0}' && '{1}' -in '{2}' -out '{3}' -kmer-size {4} -nb-cores 1 >>'{5}' 2>&1".format(dir,assemble_software,input,
                 "assembled_out",kmer, log)
            runCommand(cmd,system)

            cmd="cd '{0}' && rm -rf assembled_out.h5".format(dir)    # assembled_out.h5文件很大，而且实际上没什么用，可以删除
            runCommand(cmd, system)

            message = "{0}_gene: '{1}' has been assembled successfully. DONE".format(type, file)
            cmd = "echo {0} >> '{1}' ".format(message,whole_log_path)
            runCommand(cmd,system)
        except:
            message = "{0}_gene: '{1}' assembly failed".format(type, file)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd,system)

        # assembled_out.unitigs.fa.glue1....  只保留三个文件 assembled_out.contigs.fa assembled_out.unitigs.fa assembled_out.h5
    def wrap_assemble_reads(self, file, kmer):
        self.assemble_reads(file, kmer)
        return file

    def assemble_reads_parallel(self):
        type=self.type
        thread_number = self.thread_number
        kmer = self.kmer
        out_dir_name=self.out_dir_name
        filtered_out=self.filtered_out
        assembled_out=self.assembled_out
        assembled_log=self.assembled_log

        filtered_out_path = os.path.join(out_dir_name, filtered_out)
        files = os.listdir(filtered_out_path)            #matK
        task_pool = []
        result = []
        executor = futures.ThreadPoolExecutor(max_workers=thread_number)

        assembled_log_path = os.path.join(out_dir_name, assembled_log)
        assembled_out_path= os.path.join(out_dir_name, assembled_out)
        dir_make(assembled_log_path)
        dir_make(assembled_out_path)

        for file in files:
            task_pool.append(executor.submit(self.wrap_assemble_reads, file, kmer))
        for task in tqdm(desc="{0:<22}".format("Assembling_reads"), iterable=futures.as_completed(task_pool),
                         total=len(task_pool)):
            result.append(task.result())


    '''
    3  检测assembled_out的结果，校正长度与方向，并且对齐剪切 
    '''

    def check_contigs(self, file):
        # file=psbA   pbsA
        configuration_information=self.configuration_information
        type = self.type  # "tfa"
        my_software_name = self.my_software_name  #GM
        out_dir_name=self.out_dir_name
        GM_results = self.GM_results
        reference_database=self.reference_database
        assembled_out = self.assembled_out
        system=self.system
        whole_log=self.whole_log
        max_length=self.max_length
        min_length=self.min_length
        options=self.options

        gene_name = file
        ref_fasta = file + ".fasta"  # psbA ---- pbsA.fasta
        ref_path = os.path.join(out_dir_name, reference_database, ref_fasta)
        assembled_path = os.path.join(out_dir_name, assembled_out, file, "assembled_out.contigs.fa")
        whole_log_path=os.path.join(out_dir_name,whole_log)

        path_list = [whole_log_path]
        [whole_log_path] = get_absolute_and_map_path(path_list, system)  # 如果为windows环境，会批量映射路径
        if is_exist(assembled_path) == 0:
            message = " '{0}': assembling reads has failed,because the data quality was poor or the amount of input data was too small".format(
                file)
            cmd = "echo {0} >> '{1}'".format(message,whole_log_path)
            runCommand(cmd,system)
            return 0

        GM_results_path = os.path.join(out_dir_name, GM_results)     #GM_results的路径
        GM_results_path_raw = os.path.join(GM_results_path, gene_name + "_raw.fasta")
        GM_results_path_raw_best=os.path.join(GM_results_path,gene_name+"_raw_best.fasta")
        GM_results_path_options=os.path.join(GM_results_path,gene_name+"_options.fasta")
        GM_results_path_no_trimmed=os.path.join(GM_results_path,gene_name+".fasta")
        GM_results_path_trimmed = os.path.join(GM_results_path,gene_name + "_trimmed.fasta")

        # 限定长度
        # 确定方向
        # 切齐 双序列比对(局部最优解)
        # 得到最佳剪切结果
        my_verify=Get_the_best_result(configuration_information,out_dir_name,assembled_path,ref_path,GM_results_path_raw, GM_results_path_raw_best,GM_results_path_options,GM_results_path_no_trimmed,GM_results_path_trimmed,gene_name,max_length,min_length,options)
        my_verify.my_makeblastdb_blastn()
        m8_information=my_verify.parse_blastn_m8()  #解析blastn结果
        m8_information=my_verify.add_query_reference_information_2_m8(m8_information) #向m8_information里面增添参考序列和查询序列的信息
        m8_information=my_verify.filter_blastn_m8(m8_information) #将一致度，覆盖度，以及不合符用户规定长度的序列剔除
        my_verify.m8_2_fasta(m8_information)                     #生成最原始的数据 xx.raw.fasta
        my_verify.cut_align_m8(m8_information)                   #剪切对齐       xx.fasta xx.trimmed.fasta //xx.raw_best.fasta
        my_verify.deal_temp_file()                               #将Blast的临时文件删除


    def wrap_check_contigs(self, file):
        self.check_contigs(file)
        return  file

    def check_contigs_parallel(self):
        type = self.type

        thread_number = self.thread_number
        out_dir_name=self.out_dir_name
        GM_results=self.GM_results
        assembled_out = self.assembled_out

        path=os.path.join(out_dir_name, GM_results)
        dir_make(path)

        files = os.listdir(os.path.join(out_dir_name, assembled_out))  #matk
        # executor = futures.ThreadPoolExecutor(max_workers=thread_number) #动态规划  cpu密集型计算，使用多进程不使用多线程
        executor = futures.ProcessPoolExecutor(max_workers=thread_number)
        task_pool = []
        results = []
        for file in files:
            task_pool.append(executor.submit(self.wrap_check_contigs, file))

        for task in tqdm(desc="{0:<22}".format("Verifying_contigs"), iterable=futures.as_completed(task_pool),total=len(task_pool)):
            results.append(task.result())
        executor.shutdown()

    '''
    4 记录重要信息
    '''

    def record_log(self, file):
        # file=#matK_log.txt

        temp = {}
        type = str(self.type)  #"tfa"
        my_software_name=self.my_software_name
        out_dir_name=self.out_dir_name
        reference_database=self.reference_database
        filterd_out=self.filtered_out
        assembled_out=self.assembled_out
        assemble_log=self.assembled_log
        GM_results=self.GM_results
        system=self.system

        name = file.split("_log.txt")[0]  # rps16_log.txt ----rps16,寻找其他文件夹的名字全靠基因名
        temp["{}_gene_name".format(type)] = name  # rps16，作为excel第一列

        '''
        第一部分 记录需要计算的信息 filter (1)生成的filter_reads条数，(2)丰度 richness
        如果filter 结果失败（超远缘） 将一条序列没有
        '''
        filtered_out_path = os.path.join(out_dir_name, filterd_out, name,
                                       "Filtered_reads__R1.fastq")  # out_dir_name/filtered_out/matK/Filtered_reads__R1.fastq
        ref_path = os.path.join(out_dir_name, reference_database, name + ".fasta")  # /out_dir/reference_database/matk.fasta



        if os.path.exists(filtered_out_path):
            path_list=[filtered_out_path]
            [filtered_out_path]=get_absolute_and_map_path(path_list,system)   # 如果为windows环境，会批量映射路径

            cmd6 = "awk 'NR==2' '{0}'|wc -c".format(filtered_out_path)  # 查看fastq格式第二行长度，即一个read的长度
            cmd7 = "cat '{0}' | wc -l".format(filtered_out_path)  # 查看read条数，双末端测序 正向和反向数据合在一起算一个计量单位.记得除以4

            cmd_list=[cmd6,cmd7]
            [cmd6,cmd7]=getCommand(cmd_list,system)  #迎合 subprocess.getoutput(cmd7)

            read_length = int(subprocess.getoutput(cmd6))  # 151
            filtered_reads_number = int(int(subprocess.getoutput(cmd7)) / 4)  # 1284行/4 = 321条
            length_seq = calculate_ref_size(ref_path)
            ref_length = length_seq[0]  # 参考基因组的大小
            richness = round((filtered_reads_number * read_length) / ref_length, 2)  # 34.00 保留两位小数
            temp["filtered_reads_number"] = filtered_reads_number  # 过滤后reads数量
            temp["richness"] = richness  # 测序深度，丰度
        else:
            temp["filtered_reads_number"] = "None"  # 过滤后reads数量
            temp["richness"] = "None"  # 测序深度，丰度

        '''
        第二部分，记录最简单的信息:minia拼接的 (1)组装百分比assembly (2)成图graph construction（3）最大长度max_length
        如果minia结果失败，assembled_log  EXCEPTION: error opening file: xx.contigs.fa (No such file or directory)
        '''

        assembled_out_path = os.path.join(out_dir_name, assembled_out, name, "assembled_out.contigs.fa")
        log_path=os.path.join(out_dir_name,assemble_log,file)

        if os.path.exists(assembled_out_path) and os.path.exists(log_path):
            path_list=[assembled_out_path,log_path]
            [assembled_out_path,log_path]=get_absolute_and_map_path(path_list,system)

            cmd3 = "cat '{0}'|grep  -e 'assembly.*:\s[0-9].[0-9]' ".format(log_path)  # assembly  : 0.123
            cmd4 = "cat '{0}'|grep  -e 'graph construction' ".format(log_path)
            cmd5 = "cat '{0}'|grep  -e 'max_length' ".format(log_path)


            cmd_list=[cmd3,cmd4,cmd5]
            [cmd3,cmd4,cmd5]=getCommand(cmd_list,system)  #迎合 subprocess.getoutput(cmd7)


            assembly_infomation = subprocess.getoutput(cmd3)  # assembly                                   : 0.134
            graph_construction_information = subprocess.getoutput(cmd4)  # graph construction              : 4.510
            max_length_infomation = subprocess.getoutput(cmd5)  # max_length                               : 1993

            assemble = re.findall(r"\d*\.\d*", assembly_infomation)[0]  # 0.134
            assemble = float(assemble)
            temp["assembled_percentage"] = assemble

            graph_construction = re.findall(r"\d*\.\d*", graph_construction_information)[0]  # 4.510
            graph_construction = float(graph_construction)
            temp["graph_construction"] = graph_construction

            max_length = re.findall(r"\d+", max_length_infomation)[0]  # 1993
            max_length = float(max_length)
            temp["assembled_max_length"] = max_length
        else:
            temp["assembled_percentage"] = "None"
            temp["graph_construction"] = "None"
            temp["assembled_max_length"] = "None"

        '''
        第三部分 记录是否生成的信息 (1)GM结果是否生成，(2)GM_trimmed结果是否生成  (3)GM最大长度，(4)trimmed_GM identity/coverage
        '''


        GM_exist_path = os.path.join(out_dir_name, GM_results,  name + ".fasta")  # out_dir_name/GM_results/matK.fasta
        GM_trimmed_exist_path = os.path.join(out_dir_name, GM_results, name+ "_trimmed.fasta")  # out_dir_name/GM_results/matK_trimmed.fasta

        if os.path.exists(GM_exist_path):
            GM_length_list = []  # 考虑存在但为空的情况
            for rec in SeqIO.parse(GM_exist_path, "fasta"):
                GM_length = len(rec.seq)
                GM_length_list.append(GM_length)

            if GM_length_list != []:
                GM_results_max_length = max(GM_length_list)
                temp["results_max_length"] = GM_results_max_length
                temp["gene_extraction"] = "successful"
            else:
                temp["results_max_length"] = "None"
                temp["gene_extraction"] = "failed"


        else:
            temp["results_max_length"] = "None"
            temp["gene_extraction"] = "failed"

        if os.path.exists(GM_trimmed_exist_path):
            temp["gene_trimmed"] = "successful"
            identity_and_coverage = get_identity_and_coverage_path(GM_trimmed_exist_path, ref_path)
            temp["identity_trimmed"] = str(identity_and_coverage[0]) + "%"
            temp["coverage_trimmed"] = str(identity_and_coverage[1]) + "%"
        else:
            temp["gene_trimmed"] = "failed"
            temp["identity_trimmed"] = "None"
            temp["coverage_trimmed"] = "None"

        return temp

    def wrap_record_log(self, file):
        temp = self.record_log(file)
        return temp

    def record_log_parallel(self):
        type = str(self.type)
        thread_number = self.thread_number
        my_software_name = self.my_software_name
        out_dir_name=self.out_dir_name
        assembled_log=self.assembled_log
        results_information_excel=self.results_information_excel  #results_information.xlsx"

        files = os.listdir(os.path.join(out_dir_name, assembled_log))               #matK_log.txt

        executor = futures.ProcessPoolExecutor(max_workers=thread_number)

        task_pool = []
        results = []
        result_information = []  # 记录所有核心输出结果
        for file in files:
            task_pool.append(executor.submit(self.wrap_record_log, file))

        for task in tqdm(desc="{0:<22}".format("Recording_information"), iterable=futures.as_completed(task_pool), total=len(task_pool)):
            result_information.append(task.result())

        # 调整顺序
        order = ["{}_gene_name".format(type), "filtered_reads_number", "richness", "graph_construction",
                 "assembled_percentage",
                 "assembled_max_length", "results_max_length", "identity_trimmed", "coverage_trimmed",
                 "gene_extraction", "gene_trimmed"]



        #如果filter一个结果都没有，会报错
        try:
            df_information = pd.DataFrame(result_information)
            df_information = df_information[order]

            #向excel追加sheet，防止被覆盖的写法
            writer = pd.ExcelWriter(os.path.join(out_dir_name, results_information_excel),
                                    engine='openpyxl')  # 可以向不同的sheet写入数据
            if is_exist(os.path.join(out_dir_name, results_information_excel)) == 0:
                df_information.to_excel(writer, sheet_name="{}_genes".format(type), header=True, index=False)
                writer.save()  # 保存
            else:
                book = load_workbook(os.path.join(out_dir_name, results_information_excel))
                writer.book = book
                df_information.to_excel(writer, sheet_name="{}_genes".format(type), header=True, index=False)
                writer.save()  # 保存
        except:
            pass






