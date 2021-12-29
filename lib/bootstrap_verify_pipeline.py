#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/9/14 16:19
# @Author  : xiepulin
# @File    : bootstrap_verify_pipeline.py
# @Software: PyCharm
import shutil
import time
import argparse
import sys
import subprocess
import datetime
import os
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
import random
import tempfile
from basic import *

################################################
#################################################

'''
自展检测  每次检测一个基因
'''
##################################################
#################################################

class Bootstrap_verify():
    def __init__(self, configuration_information, type, out_dir_name, gm_result,ref,thread_number, kmer, bootstrap_number,max_length,min_length,options):
        self.configuration_information = configuration_information
        self.type = type
        self.out_dir_name = out_dir_name
        self.gm_result = gm_result            #out_dir/GM_results/matk_trimmed.fasta  #绝对路径
        self.ref=ref                          #out_dir/reference_database/matk.fasta  绝对路径
        self.thread_number = thread_number
        self.kmer = kmer
        self.bootstrap_number = bootstrap_number
        self.max_length=max_length
        self.min_length=min_length
        self.options=options    #yes or no


        self.data1=self.configuration_information["data1"]
        self.data2=self.configuration_information["data2"]
        self.my_software_name = self.configuration_information["my_software_name"]
        self.filter_software = self.configuration_information["filter_path"]  # filter_reads.pl / filter  绝对路径
        self.assemble_software = self.configuration_information["assemble_path"]  # minia                   绝对路径
        self.muscle_software=self.configuration_information["muscle_path"]     #muscle

        self.reference_database = self.configuration_information["reference_database"]  # 参考基因数据库
        self.filtered_out = self.configuration_information["filtered_out"]
        self.assembled_out = self.configuration_information["assembled_out"]
        self.GM_results = self.configuration_information["GM_results"]
        self.bootstrap_out = self.configuration_information["bootstrap_out"]

        self.system = self.configuration_information["system"]
        self.results_information_excel = self.configuration_information[
            "results_information_excel"]  # "results_information.xlsx"
        self.whole_log = self.configuration_information["whole_log"]
        self.bootstrap_data_set=self.configuration_information["bootstrap_data_set"]


        file_name = os.path.basename(self.gm_result)

        if "_trimmed.fa" in file_name:
            file_name = file_name.split("_trimmed.fa")[0]
        elif "_trimmed.fas" in file_name:
            file_name = file_name.split("_trimmed.fas")[0]
        elif "_trimmed.fasta" in file_name:
            file_name = file_name.split("_trimmed.fasta")[0]
        elif ".fa" in file_name:
            file_name = file_name.split(".fa")[0]
        elif ".fas" in file_name:
            file_name = file_name.split(".fas")[0]
        elif ".fasta" in file_name:
            file_name = file_name.split(".fasta")

        else:
            file_name = file_name
        self.file_name = file_name






    '''
    1.1获得变异率
    由GM最优结果 和 参考基因组计算一致度而来
    '''
    def get_var_rate(self):
        #file matk.fasta
        gm_result=self.gm_result
        ref=self.ref
        out_dir_name = self.out_dir_name
        type = self.type
        reference_database = self.reference_database

        GM_record=SeqIO.read(gm_result,"fasta")
        GM_sequence=GM_record.seq

        var_rate_list = []
        for rec in SeqIO.parse(ref,"fasta"):
            ref_sequence = rec.seq
            alignments = pairwise2.align.globalxx(GM_sequence, ref_sequence)  # 全局比对，相同的残基就给1分，不同和gap不扣分
            matches = alignments[0][2]
            identity = (matches / len(GM_sequence))  # 86.89%
            var_rate_temp = 1 - identity
            var_rate_list.append(var_rate_temp)
        var_rate=max(var_rate_list)
        return var_rate


    '''
    1.2随机突变 要求输入变异率 GM_result的结果 ，pmsk结果变异后的结果
    '''
    def random_replacement(self,var_rate,file):
        # file bootstrap_0.fasta
        gm_result=self.gm_result
        out_dir_name=self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        bootstrap_reference_database=self.reference_database
        gene_name=self.file_name

        id_number = str(file)  # bootstrap_0
        my_records = []
        my_file_name = str(file) + ".fasta"
        for rec in SeqIO.parse(gm_result, "fasta"):
            id = id_number
            description = ""
            seq = rec.seq

            seq_to_list = list(seq)
            seq_length = len(seq_to_list)
            var_number = int(seq_length * var_rate)
            var_site_list = random.sample(range(0, seq_length), var_number)
            var_site_list = sorted(var_site_list)
            # 随机替换的对应表
            base_dict = {"A": ["T", "C", "G"],
                         "T": ["A", "C", "G"],
                         "C": ["A", "T", "G"],
                         "G": ["A", "T", "C"]}
            # 随机变异
            for i in var_site_list:
                var = random.randint(0, 2)
                base = seq_to_list[i]
                seq_to_list[i] = base_dict[base][var]
            seq_var = "".join(seq_to_list)
            # print(seq_var)#变异后的序列
            my_record = SeqRecord(seq=Seq(seq_var), id=id, description=description)
            my_records.append(my_record)

        if my_records != []:
            path = os.path.join(out_dir_name, sub_out_dir_name,bootstrap_reference_database, my_file_name)
            SeqIO.write(my_records, path, "fasta")

    def wrap_random_replacement(self,var_rate,file):
        self.random_replacement(var_rate,file)

    def random_replacement_parallel(self):

        out_dir_name = self.out_dir_name
        sub_out_dir_name=self.bootstrap_out  #bootstrap
        reference_database=self.reference_database
        thread_number=self.thread_number
        bootstrap_number=self.bootstrap_number
        gene_name=self.file_name  #matK

        var_rate=self.get_var_rate()
        files = ["bootstrap" + "_" + str(i) for i in range(bootstrap_number)]   #bootstrap_1
        result = []
        task_pool = []
        executor = futures.ThreadPoolExecutor(max_workers=thread_number)


        bootstrap_reference_database = os.path.join(out_dir_name, sub_out_dir_name,reference_database) #out_dir_name/bootstrap/matK
        dir_make(bootstrap_reference_database)
        for file in files:
            task_pool.append(executor.submit(self.wrap_random_replacement, var_rate, file))

        total = len(task_pool)
        number = 1
        sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Preparing_data", number, total))
        for task in futures.as_completed(task_pool):
            if number < total:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Preparing_data",number, total))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Preparing_data", number, total)+"\n")
                sys.stdout.flush()
            result.append(task.result())
            number = number + 1
        executor.shutdown()  # 所有线程结束，再进行下一步



    '''
    1.3 多线程并行的方法过滤基因      过滤reads
    产生两个过滤文件Filtered_reads__R1.fastq ， Filtered_reads__R2.fastq ,并且合并双端数据为filtered.fq
    '''

    def filter_reads(self, file):
        #file=bootstrap_0.fasta
        #file_name matk
        file_name=self.file_name
        data1=self.data1
        data2=self.data2
        out_dir_name = self.out_dir_name
        sub_out_dir_name = self.bootstrap_out
        filterd_out = self.filtered_out
        bootstrap_reference_database = self.reference_database
        type = str(self.type)  # tfa,tgb, cp ,target 中的某一种
        filter_software = self.filter_software  # lib中的filter
        system = self.system
        whole_log = self.whole_log


        gene_name = file.split(".fasta")[0]
        ref = file
        dir =  gene_name  # bootstrap_0
        filtered_out_path = os.path.join(out_dir_name, sub_out_dir_name, filterd_out, dir)
        dir_make(filtered_out_path)

        data1_path=os.path.join(out_dir_name,data1)
        data2_path=os.path.join(out_dir_name,data2)
        reference_path=os.path.join(out_dir_name,sub_out_dir_name,bootstrap_reference_database,ref)
        whole_log_path=os.path.join(out_dir_name,whole_log)

        if is_exist(data1_path) and is_exist(data2_path):
            path_list = [out_dir_name, filter_software, filtered_out_path, data1_path, data2_path, reference_path,
                         whole_log_path]
            [out_dir_name, filter_software, filtered_out_path, data1_path, data2_path, reference_path,
             whole_log_path] = get_absolute_and_map_path(path_list, system)  # 如果为windows环境，会批量映射路径
            cmd = "cd '{0}' && '{1}' -1 '{2}' -2 '{3}' -r '{4}'  >/dev/null  2>&1".format(filtered_out_path,
                                                                                          filter_software, data1_path,
                                                                                          data2_path, reference_path)
            runCommand(cmd, system)
            cmd = "cd '{0}' && cat {1} {2} > {3}".format(filtered_out_path, "Filtered_reads__R1.fastq",
                                                         "Filtered_reads__R2.fastq",
                                                         "filtered.fq")
            runCommand(cmd, system)

            message = "{0}_{1}_{2}:Reads filtering has been successfully completed. DONE".format(type, file_name,gene_name)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)
        else:
            pass

    def wrap_filter_reads(self, file):
        self.filter_reads(file)

    def filter_reads_parallel(self):
        thread_number=self.thread_number
        out_dir_name = self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        filterd_out=self.filtered_out


        bootstrap_reference_database = self.reference_database
        bootstrap_reference_database_path = os.path.join(out_dir_name,sub_out_dir_name,bootstrap_reference_database)

        files = os.listdir(bootstrap_reference_database_path)  # bootstrap_0.fasta
        filtered_out_path = os.path.join(out_dir_name, sub_out_dir_name,filterd_out)
        dir_make(filtered_out_path)
        executor = futures.ThreadPoolExecutor(max_workers=thread_number)
        task_pool = []
        result = []
        for file in files:
            task_pool.append(executor.submit(self.wrap_filter_reads, file))
        total = len(task_pool)
        number = 1
        sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Filtering_reads", number, total)) #直接打印，防止第一个运行太久，否则第一个完成才会打印出来
        for task in futures.as_completed(task_pool):
            if number < total:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Filtering_reads", number, total))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Filtering_reads", number, total) + "\n")
                sys.stdout.flush()
            result.append(task.result())
            number = number + 1

        executor.shutdown()  # 所有线程结束，再进行下一步


    '''
    1.4 多线程并行的方法，拼接目标基因  拼接reads
    产生三个文件：assembled_out.contigs.fa  assembled_out.h5  assembled_out.unitigs.fa
    '''

    def assemble_reads(self, file, kmer):
        # file=bootstrap_0
        #file_name  matk
        out_dir_name = self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        assemble_software=self.assemble_software  #minia
        filtered_out=self.filtered_out
        assembed_out=self.assembled_out
        system = self.system
        whole_log = self.whole_log
        file_name=self.file_name
        gene_name=file      #boootstrap_0
        type=self.type


        input = os.path.join(out_dir_name, sub_out_dir_name,filtered_out, gene_name, "filtered.fq")
        whole_log_path = os.path.join(out_dir_name, whole_log)

        if is_exist(input) == 0:
            message = "{0}_{1}_{2} :Assembling reads has failed.  It could be that the quality of the raw data is not good enough, the amount of input data is too small or the species  are not closely related".format(
                type,file_name,gene_name)
            cmd = "echo {0}>>'{1}'".format(message,whole_log)
            subprocess.call(cmd, shell=True)
            return 0

        dir = os.path.join(out_dir_name,sub_out_dir_name, assembed_out,gene_name)
        dir_make(dir)

        path_list = [dir, assemble_software, input, whole_log_path]
        [dir, assemble_software, input, whole_log_path] = get_absolute_and_map_path(path_list,
                                                                                    system)  # 如果为windows环境，会批量映射路径


        try:
            cmd = "cd '{0}' && '{1}' -in '{2}' -out '{3}' -kmer-size {4} >/dev/null 2>&1 ".format(dir,
                                                                                                  assemble_software,
                                                                                                  input,
                                                                                                  "assembled_out", kmer)
            runCommand(cmd, system)

            cmd = "cd '{0}' && rm -rf assembled_out.h5".format(dir)  # assembled_out.h5文件很大，而且实际上没什么用，可以删除
            runCommand(cmd, system)

            message = "{0}_{1}_{2}: The gene has been successfully assembled. DONE".format(type, file_name,gene_name)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)
        except:
            message = "{0}_{1}_{2}: Assembly failed".format(type,file_name,gene_name)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)

        # assembled_out.unitigs.fa.glue1....  只保留三个文件 assembled_out.contigs.fa assembled_out.unitigs.fa assembled_out.h5

    def wrap_assembled_reads(self, file, kmer):
        self.assemble_reads(file, kmer)
        return file

    def assembled_reads_parallel(self):
        out_dir_name = self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        filtered_out=self.filtered_out
        assembled_out=self.assembled_out
        thread_number = self.thread_number
        kmer = self.kmer

        filtered_out_path = os.path.join(out_dir_name, sub_out_dir_name,filtered_out)
        files = os.listdir(filtered_out_path)  # bootstrap_0
        task_pool = []
        result = []
        executor = futures.ThreadPoolExecutor(max_workers=thread_number)
        assembled_out_path= os.path.join(out_dir_name, sub_out_dir_name,assembled_out)
        dir_make(assembled_out_path)

        for file in files:
            task_pool.append(executor.submit(self.wrap_assembled_reads, file, kmer))

        total = len(task_pool)
        number = 1
        sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Assembling_reads", number, total))
        for task in futures.as_completed(task_pool):
            if number < total:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Assembling_reads", number, total))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Assembling_reads", number, total) + "\n")
                sys.stdout.flush()
            result.append(task.result())
            number = number + 1
        executor.shutdown()  # 所有线程结束，再进行下一步

        # for task in tqdm(desc="{0:<20}".format("Assembling_reads"), iterable=futures.as_completed(task_pool),
        #                  total=len(task_pool)):
        #     result.append(task.result())



    '''
    1.5  检测assembled_out的结果，校正长度与方向，并且对齐剪切 
    '''

    def check_contigs(self, file):
        # file=bootstrap_0
        #file_name    matk
        type = self.type  # "tfa"
        configuration_information = self.configuration_information
        file_name=self.file_name
        out_dir_name=self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        reference_database=self.reference_database
        assembled_out = self.assembled_out
        GM_results = self.GM_results
        system=self.system
        whole_log=self.whole_log
        max_length=self.max_length
        min_length=self.min_length
        options = self.options
        gene_name=file               #bootstrap_0



        reference = gene_name + ".fasta"  #bootstrap_0.fasta
        ref_path = os.path.join(out_dir_name, sub_out_dir_name,reference_database,reference)
        assembled_out_path = os.path.join(out_dir_name,sub_out_dir_name, assembled_out, file, "assembled_out.contigs.fa")
        whole_log_path=os.path.join(out_dir_name,whole_log)



        path_list = [whole_log_path]
        [whole_log_path] = get_absolute_and_map_path(path_list, system)  # 如果为windows环境，会批量映射路径
        if is_exist(assembled_out_path) == 0:
            message = "{0}_{1}_{2}:Failed to pass verification.".format(type,file_name,gene_name)
            cmd = "echo {0} >> '{1}'".format(message,whole_log_path)
            runCommand(cmd,system)
            return 0


        GM_results_path = os.path.join(out_dir_name, sub_out_dir_name,GM_results)     #GM_results的路径
        GM_results_path_raw = os.path.join(GM_results_path, gene_name + "_raw.fasta")
        GM_results_path_raw_best = os.path.join(GM_results_path, gene_name + "_raw_best.fasta")
        GM_results_path_options=os.path.join(GM_results_path,gene_name+"_.options.fasta")
        GM_results_path_no_trimmed = os.path.join(GM_results_path, gene_name + ".fasta")
        GM_results_path_trimmed = os.path.join(GM_results_path,gene_name + "_trimmed.fasta")

        my_verify = Get_the_best_result(configuration_information,out_dir_name,assembled_out_path, ref_path,GM_results_path_raw,
                                        GM_results_path_raw_best, GM_results_path_options, GM_results_path_no_trimmed,
                                        GM_results_path_trimmed, gene_name, max_length, min_length, options)

        my_verify.my_makeblastdb_blastn()
        m8_information = my_verify.parse_blastn_m8()  # 解析blastn结果
        m8_information = my_verify.add_query_reference_information_2_m8(
            m8_information)  # 向m8_information里面增添参考序列和查询序列的信息
        m8_information = my_verify.filter_blastn_m8(m8_information)  # 将一致度，覆盖度，以及不合符用户规定长度的序列剔除
        my_verify.m8_2_fasta(m8_information)  # 生成最原始的数据 xx.raw.fasta
        my_verify.cut_align_m8(m8_information)  # 剪切对齐       xx.fasta xx.trimmed.fasta //xx.raw_best.fasta
        my_verify.deal_temp_file()  # 将Blast的临时文件删除


    def wrap_check_contigs(self, file):
       self.check_contigs(file)

    def check_contigs_parallel(self):
        type = self.type

        thread_number = self.thread_number
        out_dir_name=self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        assembled_out = self.assembled_out
        GM_results = self.GM_results


        assembled_out_path=os.path.join(out_dir_name,sub_out_dir_name,assembled_out)
        GM_results_path=os.path.join(out_dir_name,sub_out_dir_name,GM_results)
        dir_make(GM_results_path)

        files = os.listdir(assembled_out_path)  #bootstrap_0
        # executor = futures.ThreadPoolExecutor(max_workers=thread_number) #动态规划  cpu密集型计算，使用多进程不使用多线程
        executor = futures.ProcessPoolExecutor(max_workers=thread_number)
        task_pool = []
        results = []
        for file in files:
            task_pool.append(executor.submit(self.wrap_check_contigs, file))


        total = len(task_pool)
        number = 1
        sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Verifying_contigs", number, total))
        for task in futures.as_completed(task_pool):
            if number < total:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Verifying_contigs", number, total))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r' + "{0:<22}:{1:>4}/{2}".format("Verifying_contigs", number, total) + "\n")
                sys.stdout.flush()

            results.append(task.result())   #没有return返回None 有的话返回函数的return值
            number = number + 1
        executor.shutdown()

    '''
    记录结果
    '''

    def record_bootstrap_information(self, file):
        # file   bootstrap_0.fasta
        name = re.findall(r"bootstrap_\d+", file)[0]  # bootstrap_0.fasta ---- bootstrap_0
        bootstrap_number = re.findall(r"\d+", file)[0]  # bootstrap_0 ----0



        out_dir_name = self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        GM_results = self.GM_results
        gm_result = self.gm_result  # 之前用软件Gene miner做出来的结果

        GM_results_path = os.path.join(out_dir_name, sub_out_dir_name, GM_results)  # out_dir_name/bootstrap_0/GM_results
        files = os.listdir(GM_results_path)  # bootstrap_0_trimmed.fasta, bootstrap_0.fasta

        '''
        用bootstrap做出来的结果：（1）可能一个结果也没有，（2）没能剪切bootstrap_0_raw.fasta,bootstrap_0_raw_bset.fasta （3）能够剪切 bootstrap_0_raw.fasta,bootstrap_0.fasta,bootstrap_0_trimmed.fasta
        '''

        raw_best = name + "_raw_best.fasta"  # 没有通过剪切规则
        no_trimmed = name + ".fasta"  # 通过剪切规则 剪切前
        trimmed = name + "_trimmed.fasta"  # 通过剪切规则 剪切后


        flag = 0
        for i in files:
            if raw_best in i:  # bootstrap_0_raw_best
                flag = 1
                break
            elif trimmed in i:  # bootstrap_0_trimmed.fasta
                flag = 2
                break
            else:
                flag = 0  # 啥也没有

        bootstrap_information = {}

        if flag == 0:  # 没有做出结果
            bootstrap_information["bootstrap_number"] = int(bootstrap_number)
            bootstrap_information["identity"] = "None"
            bootstrap_information["coverage"] = "None"

        elif flag == 1:  # bootstrap_0_raw_best

            bootstrap_gm_result = name + "_raw_best.fasta"
            ref_path = os.path.join(GM_results_path, bootstrap_gm_result)
            identity_and_coverage = get_identity_and_coverage_path(gm_result, ref_path)
            bootstrap_information["bootstrap_number"] = int(bootstrap_number)
            bootstrap_information["identity"] = identity_and_coverage[0]
            bootstrap_information["coverage"] = identity_and_coverage[1]

        else:  # bootstrap_0_trimmed.fasta
            bootstrap_gm_result = name + "_trimmed.fasta"
            ref_path = os.path.join(GM_results_path, bootstrap_gm_result)
            identity_and_coverage = get_identity_and_coverage_path(gm_result, ref_path)
            bootstrap_information["bootstrap_number"] = int(bootstrap_number)
            bootstrap_information["identity"] = identity_and_coverage[0]
            bootstrap_information["coverage"] = identity_and_coverage[1]

        return bootstrap_information

    def wrap_record_bootstrap_information(self, file):
        bootstrap_information = self.record_bootstrap_information(file)
        return bootstrap_information

    def record_bootstrap_information_parallel(self):
        out_dir_name = self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        reference_database=self.reference_database
        thread_number = self.thread_number
        file_name = self.file_name   #matk

        ref_path = os.path.join(out_dir_name, sub_out_dir_name,reference_database)  # out_dir_name/bootstrap_out/reference_database
        files = os.listdir(ref_path)  # bootstrap_0.fasta
        task_pool = []

        executor = futures.ProcessPoolExecutor(max_workers=thread_number)
        for file in files:
            task_pool.append(executor.submit(self.wrap_record_bootstrap_information, file))

        bootstrap_information_all = []
        message_all=[]



        for i in task_pool:
            bootstrap_information = i.result()
            bootstrap_information_all.append(bootstrap_information)
            nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')

            identity = bootstrap_information["identity"]
            coverage = bootstrap_information["coverage"]
            bootstrap_number = bootstrap_information["bootstrap_number"]

            if identity == "None":
                message = "{0}|{1}|bootstrap_{2:<4}:No results were obtained".format(nowTime, file_name,
                                                                                      bootstrap_number)
                print(message, flush=True)
                message_all.append(message)

            else:
                message = "{0}|{1}|bootstrap_{2:<4}:identity:{3}% coverage:{4}%".format(nowTime, file_name,
                                                                                        bootstrap_number, identity,
                                                                                        coverage)
                print(message, flush=True)
                message_all.append(message)

        print("")  # 所有的Bootstap打印结束后，分隔开

        return message_all

    '''
    将bootstrap结果归并在一起，写为temp文件 bootstrap_data_set.fasta
    '''

    def prepare_bootstrap_set(self):
        out_dir_name=self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        GM_results=self.GM_results
        muscle_path=self.muscle_software
        system=self.system


        bootstrap_data_set= self.bootstrap_data_set


        bootstrap_set_path=os.path.join(out_dir_name, sub_out_dir_name, GM_results, bootstrap_data_set)
        GM_results_path=os.path.join(out_dir_name,sub_out_dir_name,GM_results)

        bootstrap_path=[]
        files=os.listdir(GM_results_path)
        sign="_trimmed.fasta"
        for i in files:
            if sign in str(i):              #因为list下不仅有文件，还有文件夹。str(i) 让警告都不出现
                path=os.path.join(out_dir_name,sub_out_dir_name,GM_results,i)
                bootstrap_path.append(path)

        if bootstrap_path==[]:
            print("Low support and failed to generate consensus sequence")
            flag=0
            return  [flag,bootstrap_set_path]
        else:
            #写所有的fasta文件
            bootstrap_set_records = []
            for i in bootstrap_path:
                rec=SeqIO.read(i,"fasta")
                basename=os.path.basename(i).split(sign)[0]  #D:\Happy_life_and_work\scu\python\Gene_Miner\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\second\bootstrap_out\matK\GM_results\bootstrap_7_trimmed.fasta -----bootstrap_7

                id=basename
                seq=rec.seq
                description="length_"+str(len(seq))
                my_record=SeqRecord(id=id,seq=seq,description=description)
                bootstrap_set_records.append(my_record)
            SeqIO.write(bootstrap_set_records,bootstrap_set_path,"fasta")

            #多序列比对(覆盖式)
            input=bootstrap_set_path
            output=bootstrap_set_path
            path_list = [muscle_path, input, output]
            [muscle_path, input, output] = get_absolute_and_map_path(path_list, system)
            # cmd = "'{0}' -align '{1}' -output '{2}' >/dev/null 2>&1  ".format(muscle_path, input, output)

            cmd = "'{0}' -in '{1}' -out '{2}' >/dev/null 2>&1  ".format(muscle_path, input, output) #muscle3
            runCommand(cmd, system)

            flag=1   #但凡有一个Bootstrap结果 都会置为1
            return [flag,bootstrap_set_path]



    '''
    用bootstrap序列生成的一致序列
    并计算一致序列每个位点的支持度
    '''
    def bootstrap_support_rate(self):
        out_dir_name = self.out_dir_name
        sub_out_dir_name = self.bootstrap_out
        GM_results=self.GM_results
        file_name=self.file_name #matk

        #获得bootstrap打印结果
        message_all = self.record_bootstrap_information_parallel()

        #获得支持率结果
        bootstrap_set_information=self.prepare_bootstrap_set()
        flag=bootstrap_set_information[0]
        if flag==0:
            return 0
        else:
            bootstrap_set=bootstrap_set_information[1]

        # (1)统计bootstrap长度，选择其中的一条即可，此时的bootstrap_set_data是多序列比对之后的
        number = 0
        length = 0
        for rec in SeqIO.parse(bootstrap_set, "fasta"):
            if number == 0:
                length = len(rec.seq)
                break

        # 生成位点信息记录字典，记录每个位点出现碱基及其对应概率
        all_site = []
        for i in range(length):
            temp = {}
            temp["{}".format(str(i))] = []
            all_site.append(temp)
        for rec in SeqIO.parse(bootstrap_set, "fasta"):
            for i in range(length):
                all_site[i]["{}".format(str(i))].append(rec.seq[i])

        consensus_information = []
        for i in range(length):
            temp = {}
            number_A = all_site[i]["{}".format(str(i))].count("A")
            number_C = all_site[i]["{}".format(str(i))].count("C")
            number_T = all_site[i]["{}".format(str(i))].count("T")
            number_G = all_site[i]["{}".format(str(i))].count("G")
            number_N = all_site[i]["{}".format(str(i))].count("N")
            number_gap = all_site[i]["{}".format(str(i))].count("-")

            number_list = [number_A, number_T, number_C, number_G, number_N,number_gap]
            base_list = all_site[i]["{}".format(str(i))]  # ['A', 'T', 'A', 'A', 'A', 'A']
            most_likely_base = max(base_list, key=base_list.count)  # 当两种情况出现的时候，会默认出现第一个
            temp["base"] = most_likely_base
            temp["support_rate"] = round(max(number_list) / len(base_list), 4)
            consensus_information.append(temp)

        # 记录输出信息
        consensus_sequence_list = []
        consensus_support_rate_list = []
        for i in consensus_information:
            consensus_sequence_list.append(i["base"])
            consensus_support_rate_list.append(i["support_rate"])

        consensus_sequence = "".join(consensus_sequence_list)
        consensus_support_rate = ""
        for i in consensus_support_rate_list:
            consensus_support_rate = consensus_support_rate + str("{0:<8}".format(i))


        consensus_file=os.path.join(out_dir_name,sub_out_dir_name,GM_results,"{}_bootstrap.txt".format(file_name))
        with open(consensus_file,"a") as f:
            f.write("bootstrap information:"+"\n")
            for i in message_all:
                f.write(i+"\n")
            f.write("\n")
            f.write("consensus sequenece:\n")
            f.write(consensus_sequence+"\n")
            f.write("\n")
            # f.write(split_line+"\n")
            f.write("support rate:"+"\n")
            f.write(consensus_support_rate+"\n")


        #生成bootstrap一致序列 fasta 记得去除空位"-"
        consensus_sequence_list=consensus_sequence.split("-")
        consensus_sequence_ultimate=""
        for i in consensus_sequence_list:
            consensus_sequence_ultimate=consensus_sequence_ultimate+i

        id = file_name + "_" + "bootstrap_consensus"
        description="lengtn_"+str(len(consensus_sequence_ultimate))
        seq = Seq(consensus_sequence_ultimate)
        my_consensus_fasta_file=SeqRecord(id=id,seq=seq,description=description)
        my_consensus_fasta_file_path=os.path.join(out_dir_name,sub_out_dir_name,GM_results,"{}_bootstrap_consensus.fasta".format(file_name))
        SeqIO.write(my_consensus_fasta_file,my_consensus_fasta_file_path,"fasta")
        #删除临时文件
        # os.remove(bootstrap_set)  #temp文件



    '''
    打包一个基因的bootstrap结果
    '''
    def pack_bootstrap_one(self):
        out_dir_name=self.out_dir_name
        sub_out_dir_name=self.bootstrap_out
        reference_database=self.reference_database
        filtered_out=self.filtered_out
        assembled_out=self.assembled_out
        GM_results=self.GM_results
        file_name=self.file_name   #matk

        pack_dir=os.path.join(out_dir_name,sub_out_dir_name,file_name)
        dir_make(pack_dir)

        reference_database_path=os.path.join(out_dir_name,sub_out_dir_name,reference_database)
        filtered_out_path=os.path.join(out_dir_name,sub_out_dir_name,filtered_out)
        assembled_out_path=os.path.join(out_dir_name,sub_out_dir_name,assembled_out)
        GM_results_path=os.path.join(out_dir_name,sub_out_dir_name,GM_results)
        folder_needed_pack = [reference_database_path,filtered_out_path,assembled_out_path, GM_results_path]

        try:
            for i in folder_needed_pack:
                shutil.move(i,pack_dir)
        except:
            pass


















################################
#################################
'''
批量自展检测
'''
##################################
###################################

class Bootstrap_pipeline():
    def __init__(self,configuration_information, type, out_dir_name, thread_number, kmer,bootstrap_number,max_length,min_length,options):
        self.configuration_information = configuration_information
        self.type = type
        self.out_dir_name = out_dir_name
        self.thread_number = thread_number
        self.kmer = kmer
        self.bootstrap_number = bootstrap_number

        self.my_software_name = self.configuration_information["my_software_name"]
        self.filter_software = self.configuration_information["filter_path"]  # filter_reads.pl / filter  绝对路径
        self.assemble_software = self.configuration_information["assemble_path"]  # minia                   绝对路径
        self.reference_database = self.configuration_information["reference_database"]  # 参考基因数据库
        self.filtered_out = self.configuration_information["filtered_out"]
        self.assembled_out = self.configuration_information["assembled_out"]
        self.GM_results = self.configuration_information["GM_results"]
        self.bootstrap_out = self.configuration_information["bootstrap_out"]

        self.system = self.configuration_information["system"]
        self.results_information_excel = self.configuration_information[
            "results_information_excel"]  # "results_information.xlsx"
        self.whole_log = self.configuration_information["whole_log"]
        self.max_length = max_length
        self.min_length = min_length
        self.options = options  # yes or no

    '''
    1.获得可以用来Bootstrap的基因和ref
    '''

    def prepare_bootstrap_data(self):
        out_dir_name = self.out_dir_name
        bootstrap_out = self.bootstrap_out
        type = self.type
        reference_database = self.reference_database
        GM_results = self.GM_results

        reference_database_path = os.path.join(out_dir_name, reference_database)  # out_dir_name/reference_database
        GM_results_path = os.path.join(out_dir_name, GM_results)
        ref_files_list = os.listdir(reference_database_path)


        GM_results_list = []  # 存储挖掘出来的基因，trimmed优先于raw
        ref_files_list_ultimate = []  # 存储对应的参考序列
        failed_gene = []  # 存储未能挖掘出的基因
        for i in ref_files_list:
            i = str(i)
            gene_name = i.split(".fasta")[0]  # trnI-GAU.fasta ---- trnI-GAU

            GM_result_trimmed_name = gene_name + "_trimmed.fasta"  # trnI-GAU----trnI-GAU_trimmed.fasta
            GM_result_raw_name = gene_name + ".fasta"  # trnI-GAU----trnI-GAU.fasta

            GM_result_trimmed_path = os.path.join(GM_results_path, GM_result_trimmed_name)
            GM_result_raw_path = os.path.join(GM_results_path, GM_result_raw_name)

            ref_path = os.path.join(reference_database_path, i)
            if os.path.exists(GM_result_trimmed_path):
                GM_results_list.append(GM_result_trimmed_path)
                ref_files_list_ultimate.append(ref_path)
            elif os.path.exists(GM_result_raw_path):
                pass
            else:
                failed_gene.append(gene_name)

        prepare_data_information=[GM_results_list,ref_files_list_ultimate,failed_gene]
        return  prepare_data_information

    def run_bootstrap_pipeline(self):
        prepare_data_information=self.prepare_bootstrap_data()
        gm_results_list=prepare_data_information[0]   #out_dir/GM_result/matk_trimmed.fasta psbA
        system=self.system
        whole_log=self.whole_log
        out_dir_name=self.out_dir_name
        max_length=self.max_length
        min_length=self.min_length
        options=self.options

        whole_log_path=os.path.join(out_dir_name,whole_log)
        if gm_results_list==[]:
            path_list=[whole_log_path]
            [whole_log_path]=get_absolute_and_map_path(path_list,system)
            message="GeneMiner can't do bootstrap verification, because it didn't generate the good quality trimmed sequences."
            #打印
            print(message)
            #记录
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)
            return 0

        configuration_information=self.configuration_information
        type=self.type
        out_dir=self.out_dir_name
        ref=prepare_data_information[1]             #matk 对应的参考序列 绝对路径
        thread_number=self.thread_number
        kmer=self.kmer
        bootstrap_number=self.bootstrap_number



        index=0
        task_number=len(gm_results_list)

        for i in range(task_number):
            my_bootstrap=Bootstrap_verify(configuration_information,type, out_dir,gm_results_list[i],ref[i],thread_number, kmer, bootstrap_number,max_length,min_length,options)

            index=index+1
            message="{0:<22}:{1:>4}/{2}".format("Bootstrap",index,task_number)
            print(message,flush=True)

            my_bootstrap.random_replacement_parallel()
            my_bootstrap.filter_reads_parallel()
            my_bootstrap.assembled_reads_parallel()
            my_bootstrap.check_contigs_parallel()
            my_bootstrap.bootstrap_support_rate()
            my_bootstrap.pack_bootstrap_one()




# '''
# windows
# '''
#
if __name__ == '__main__':

    '''
    单个基因
    '''

    # system=platform.system().lower()
    # configuration_information = {'mito_dir': 'mito_genes', 'cp_dir': 'cp_genes', 'tfa_dir': 'target_genes_from_fa', 'tgb_dir': 'target_genes_from_gb', 'data1': 'data1.fq', 'data2': 'data2.fq', 'results_information_excel': 'results_information.xlsx', 'reference_database': 'reference_database', 'filtered_out': 'filtered_out', 'assembled_out': 'assembled_out', 'assembled_log': 'assembled_log', 'callback_out': 'callback_out', 'bootstrap_out': 'bootstrap_out', 'GM_results': 'GM_results', 'filter_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\\lib\\filter', 'assemble_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\\lib\\minia', 'muscle_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\\lib\\muscle3', 'my_software_name': 'GM', 'system': 'windows', 'whole_log': 'log.txt', 'bootstrap_data_set': 'bootstrap_data_set.fasta', 'bootstrap_concensusu': 'bootstrap_concensus.fasta'}
    #
    # out_dir=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\problem_test"
    # thread_number=8
    # kmer=43
    # bootstrap_number=30
    # max_length=5000
    # min_length=300
    # options="no"
    #
    # gm_result=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\problem_test\GM_results\Acronema_trimmed.fasta"
    # ref=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\problem_test\reference_database\Acronema.fasta"
    #
    # my_bootstrap=Bootstrap_verify(configuration_information,"rtfa", out_dir, gm_result,ref,thread_number, kmer, bootstrap_number,max_length,min_length,options )
    # my_bootstrap.random_replacement_parallel()
    # my_bootstrap.filter_reads_parallel()
    # my_bootstrap.assembled_reads_parallel()
    # my_bootstrap.check_contigs_parallel()
    # my_bootstrap.bootstrap_support_rate()
    # my_bootstrap.pack_bootstrap_one()
    #


    '''
    批量基因
    '''


    system=platform.system().lower()
    configuration_information = {'mito_dir': 'mito_genes', 'cp_dir': 'cp_genes', 'tfa_dir': 'target_genes_from_fa', 'tgb_dir': 'target_genes_from_gb', 'data1': 'data1.fq', 'data2': 'data2.fq', 'results_information_excel': 'results_information.xlsx', 'reference_database': 'reference_database', 'filtered_out': 'filtered_out', 'assembled_out': 'assembled_out', 'assembled_log': 'assembled_log', 'callback_out': 'callback_out', 'bootstrap_out': 'bootstrap_out', 'GM_results': 'GM_results', 'filter_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\\lib\\filter', 'assemble_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\\lib\\minia', 'muscle_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\\lib\\muscle3', 'my_software_name': 'GM', 'system': 'windows', 'whole_log': 'log.txt', 'bootstrap_data_set': 'bootstrap_data_set.fasta', 'bootstrap_concensusu': 'bootstrap_concensus.fasta'}

    # out_dir=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\second"
    out_dir=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\weilingcai"
    thread_number=8
    kmer=43
    bootstrap_number=5
    max_length=5000
    min_length=300
    options="no"
    my_bootstrap_pipeline=Bootstrap_pipeline(configuration_information,"tfa",out_dir,thread_number,kmer,bootstrap_number,max_length,min_length,options)
    my_bootstrap_pipeline.run_bootstrap_pipeline()