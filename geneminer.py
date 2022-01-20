#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:58
# @Author  : xiepulin
# @File    : geneminer.py
# @Software: PyCharm



import argparse
import sys
import subprocess
import datetime
import os
import logging
import re
import multiprocessing
import signal
import threading


'''
导入第三方库（非标准库）
'''
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import pairwise2
from openpyxl import  load_workbook
from tqdm import tqdm
from concurrent import futures
# import PySimpleGUI as sg

cur_path=os.path.realpath(sys.argv[0]) # 脚本当前路径
father_path=os.path.dirname(cur_path) # 脚本的父目录
sys.path.append(os.path.join(father_path, "lib"))  # filter ，assemble的软件存放在脚本同级目录下的lib




from lib.basic import *
from lib.verify_parameters import *
from lib.extract_data import *
from lib.build_reference_database import *
from lib.core_pipeline  import *
from lib.iterative_pipeline import *
from lib.bootstrap_verify_pipeline import *
from lib.pack_results import pack_the_results
gv._init_()#在basic 中已经申明过了

def get_filter_assemble_software_path():

    cur_path = os.path.realpath(sys.argv[0])  # 脚本当前路径
    cur_path = os.path.dirname(cur_path)  # 脚本的父目录,father_path 覆盖

    filter_path = os.path.join(cur_path, "lib", "filter")
    minia_path = os.path.join(cur_path, "lib", "minia")
    # muscle_path=os.path.join(cur_path,"lib","muscle")  #muscle5
    muscle_path = os.path.join(cur_path, "lib", "muscle3") #muscle3
    makeblastdb_path=os.path.join(cur_path,"lib","makeblastdb")
    blastn_path=os.path.join(cur_path,"lib","blastn")

    return  [filter_path,minia_path,muscle_path,makeblastdb_path,blastn_path]






def get_configuration_information():
    configuration_information={}
    #一级目录
    mito_dir="mito_genes"
    cp_dir="cp_genes"
    tfa_dir="target_genes_from_fa"
    tgb_dir="target_genes_from_gb"
    data1 = "data1.fq"
    data2 = "data2.fq"
    results_information_excel = "results_information.xlsx"
    #二级目录
    reference_database="reference_database"
    filtered_out= "filtered_out"
    assembled_out="assembled_out"
    assembled_log="assembled_log"
    iterated_out = "iterated_out"
    bootstrap_out="bootstrap_out"
    GM_results="GM_results"
    blastn_out="blastn_out"

    #位置信息
    path=get_filter_assemble_software_path()
    filter_path=path[0]    #filter绝对路径
    assemble_path=path[1]  #minia绝对路径
    muscle_path=path[2]    #muscle绝对路径
    makeblastdb_path=path[3] #makeblastdb_path 绝对路径
    blastn_path=path[4]      #blastn 绝对路径
    #其他信息
    my_software_name="GM"
    system=get_platform()
    whole_log="log.txt"
    bootstrap_data_set="bootstrap_data_set.fasta"    #bootstarp  将所有Bootstrap结果写在一起的fasta
    bootstrap_concensus="bootstrap_concensus.fasta"  #导出共识序列
    #构建配置信息
    configuration_information["mito_dir"] = mito_dir
    configuration_information["cp_dir"] =cp_dir
    configuration_information["tfa_dir"] = tfa_dir
    configuration_information["tgb_dir"] = tgb_dir
    configuration_information["data1"] = data1
    configuration_information["data2"] = data2
    configuration_information["results_information_excel"] = results_information_excel
    configuration_information["reference_database"] =  reference_database
    configuration_information["filtered_out"] =  filtered_out
    configuration_information["assembled_out"] = assembled_out
    configuration_information["assembled_log"] =assembled_log
    configuration_information["iterated_out"] = iterated_out
    configuration_information["bootstrap_out"] = bootstrap_out
    configuration_information["GM_results"] = GM_results
    configuration_information["blastn_out"]=blastn_out
    configuration_information["filter_path"] = filter_path
    configuration_information["assemble_path"] = assemble_path
    configuration_information["muscle_path"]=muscle_path
    configuration_information["makeblastdb_path"]=makeblastdb_path
    configuration_information["blastn_path"]=blastn_path
    configuration_information["my_software_name"]=my_software_name
    configuration_information["system"]=system
    configuration_information["whole_log"] = whole_log
    configuration_information["bootstrap_data_set"]=bootstrap_data_set
    configuration_information["bootstrap_concensusu"]=bootstrap_concensus



    return  configuration_information

# print(get_configuration_information())


def main(args):
    # signal.signal(signal.SIGINT, signal_handler)  #检测函数终止退出（ctrl+c）
    gv.set_value("my_gui_flag",1)  #用于判定GUI是否处于运行状态，1代表运行，0代表没有运行

    configuration_information=get_configuration_information()
    '''
    从程序外获得的参数信息
    '''
    data1 = args.data1
    data2 = args.data2
    paired = args.paired
    single= args.single
    target_reference_fa = args.target_reference_fa
    target_reference_gb=args.target_reference_gb
    cp_reference = args.cp_reference
    mito_reference = args.mito_reference
    out_dir = args.out



    max_length = args.max
    min_length = args.min
    soft_boundary = args.soft_boundary
    data_size = args.number
    kmer = args.kmer
    thread = args.thread
    bootstrap_number=args.bootstrap_number
    iterative_number=args.iterative_number
    sf=args.reference_filtering


    data1=get_absolute(data1)
    data2=get_absolute(data2)
    paired=get_absolute(paired)
    paired=get_absolute(paired)
    single=get_absolute(single)
    target_reference_fa=get_absolute(target_reference_fa)
    target_reference_gb=get_absolute(target_reference_gb)
    cp_reference=get_absolute(cp_reference)
    mito_reference=get_absolute(mito_reference)





    #校验信息
    check_python_version()
    system=get_platform().lower()
    check_out_dir_GUI(out_dir,system)
    out_dir=check_out_dir(out_dir)                    #输出文件夹检测
    out_dir = get_absolute(out_dir)
    gv.set_value("out_dir", out_dir)
    thread_number = check_threads_number(thread)  # 确定线程数量
    kmer=check_kmer(kmer)    #检测kmer
    soft_boundary=check_soft_boundary(soft_boundary)  #检测软边界
    check_max_min_length(max_length,min_length)      #检测基因最大最小长度
    check_ref_format(target_reference_fa,target_reference_gb,cp_reference,mito_reference)  #检测参考序列格式
    sf= check_reference_filtering_strategy(sf)  #检测参考序列过滤策略
    data_size=check_datasize(data_size)                        #检测数据量 大小
    raw_data_dict=check_raw_data_type(data1,data2,paired,single)  #第一次出现：检测原始数据类型，方便判断是否运行 （有可能是相对路径）

    bootstrap_information=check_bootstrap_parameter(bootstrap_number)  #检测自展检测参数["Yes",100]
    bootstrap_number=bootstrap_information[1]

    iterative_information=check_iterative_parameter(iterative_number)     #检测梯度逼近参数 ["Yes",10]
    iterative_number=iterative_information[1]


    task_pool=check_whole_task(target_reference_fa,target_reference_gb,cp_reference,mito_reference)



    parameter_information_dict = {"project name": out_dir,
                                  "data1": data1, "data2": data2, "paired": paired, "single": single,
                                  "reference (fa)": target_reference_fa, "reference (gb)": target_reference_gb,
                                  "reference (cp)": cp_reference, "reference (mito)": mito_reference,
                                  "max length": max_length, "min length": min_length,
                                  "soft boundary": soft_boundary, "data size":data_size,
                                  "k-mer": kmer, "threads": thread_number,
                                  "bootstrap": bootstrap_information[0],
                                  "bootstrap number": bootstrap_information[1],
                                  "iterative":iterative_information[0],
                                  "iterative number":iterative_information[1]}

    print_parameter_information(parameter_information_dict)
    extract_raw_data(out_dir,raw_data_dict,data1,data2,paired,single,data_size,configuration_information) # 检查原始数据格式是否正确，提取原始数据,产生out_dir文件夹



    if "target_gene_from_fa" in task_pool:
        my_fasta=Extract_reference_fasta(configuration_information,"tfa", out_dir, target_reference_fa, soft_boundary, sf, max_length, min_length)
        my_fasta.extract_reference_from_fasta()
        my_fasta.filter_refrence()

        my_core_pipeline=Pipeline(configuration_information,"tfa",out_dir,thread_number,kmer,max_length,min_length,"yes")
        my_core_pipeline.filter_reads_parallel()
        my_core_pipeline.assemble_reads_parallel()
        my_core_pipeline.check_contigs_parallel()
        my_core_pipeline.record_log_parallel()


        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(configuration_information, "tfa", out_dir, thread_number, kmer, iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)

        if bootstrap_number!="None":
           my_bootstrap_pipeline=Bootstrap_pipeline(configuration_information,"tfa",out_dir,thread_number,kmer,bootstrap_number,max_length,min_length,"no")
           my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("tfa", configuration_information, out_dir, bootstrap_number, iterative_number)




    #Genbank格式目标基因流程
    if "target_gene_from_gb" in task_pool:
        my_fasta = Extract_reference_fasta(configuration_information, "tgb", out_dir, target_reference_gb, soft_boundary, sf, max_length, min_length)
        my_fasta.extract_reference_from_gb()
        my_fasta.filter_refrence()

        my_core_pipeline = Pipeline(configuration_information, "tgb", out_dir, thread_number, kmer,max_length,min_length,"yes")
        my_core_pipeline.filter_reads_parallel()
        my_core_pipeline.assemble_reads_parallel()
        my_core_pipeline.check_contigs_parallel()
        my_core_pipeline.record_log_parallel()

        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(configuration_information, "tgb", out_dir, thread_number, kmer, iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)




        if bootstrap_number != "None":
            my_bootstrap_pipeline = Bootstrap_pipeline(configuration_information, "tgb", out_dir, thread_number, kmer,bootstrap_number,max_length,min_length,"no")
            my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("tgb", configuration_information, out_dir, bootstrap_number, iterative_number)

    #叶绿体基因流程
    if "cp_gene" in task_pool:
        my_fasta = Extract_reference_fasta(configuration_information, "cp", out_dir, cp_reference, soft_boundary, sf, max_length, min_length)
        my_fasta.extract_reference_from_gb()

        my_core_pipeline = Pipeline(configuration_information, "cp", out_dir, thread_number, kmer,max_length,min_length,"yes")
        my_core_pipeline.filter_reads_parallel()
        my_core_pipeline.assemble_reads_parallel()
        my_core_pipeline.check_contigs_parallel()
        my_core_pipeline.record_log_parallel()

        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(configuration_information, "cp", out_dir, thread_number, kmer,
                                                         iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)


        if bootstrap_number != "None":
            my_bootstrap_pipeline = Bootstrap_pipeline(configuration_information, "cp", out_dir, thread_number, kmer,bootstrap_number,max_length,min_length,"no")
            my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("cp", configuration_information, out_dir, bootstrap_number, iterative_number)


    #线粒体基因流程
    if "mito_gene" in task_pool:
        mito_fasta = Extract_reference_fasta(configuration_information, "mito", out_dir, mito_reference, soft_boundary, sf, max_length, min_length)
        mito_fasta.extract_reference_from_gb()
        mito_fasta.filter_refrence()

        mito_pipeline = Pipeline(configuration_information, "mito", out_dir, thread_number, kmer,max_length,min_length,"yes")
        mito_pipeline.filter_reads_parallel()
        mito_pipeline.assemble_reads_parallel()
        mito_pipeline.check_contigs_parallel()
        mito_pipeline.record_log_parallel()

        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(configuration_information, "mito", out_dir, thread_number, kmer,
                                                         iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)

        if bootstrap_number != "None":
            my_bootstrap_pipeline = Bootstrap_pipeline(configuration_information, "mito", out_dir, thread_number, kmer,bootstrap_number,max_length,min_length,"no")
            my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("mito", configuration_information, out_dir, bootstrap_number, iterative_number)
    my_gui_flag=0  #运行结束，释放



if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)  #检测函数终止退出（ctrl+c） #必须放在主程序中
    multiprocessing.freeze_support()  # windows上Pyinstaller打包多进程程序需要添加特殊指令
    gv.set_value("my_gui_flag", 0)  #用于判定脚本是否跑完，还可以防止run双击覆盖事件
    parser = argparse.ArgumentParser(usage="%(prog)s <-1 -2|-s|-12>  <-rn|rcp|-rmito|rt>  [options]",
                                     description="GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data\n"
                                                 "Version: 1.0\n"
                                                 "Copyright (C) 2021 Pulin Xie\n"
                                                 "Please contact <xiepulin@stu.edu.scu.cn> if you have any questions",

                                     formatter_class=argparse.RawTextHelpFormatter,
                                     # help信息中会自动取消掉换行符和空格，argparse.RawTextHelpFormatter
                                     # 可以将help信息恢复为原始的文本信

                                     # epilog="xiepulin", #参数说明之后，显示程序的其他说明
                                     # add_help=False   #禁用帮助信息
                                     )
    #原始数据输入部分
    basic_option_group = parser.add_argument_group(title="Basic option")
    basic_option_group.add_argument("-1", dest="data1",
                                    help="One end of the paired-end reads, support fastq/fastq.gz/fastq.bz2", metavar="")
    basic_option_group.add_argument("-2", dest="data2",
                                    help="Another end of the paired-end reads, support fastq/fastq.gz/fastq.bz2", metavar="")
    basic_option_group.add_argument("-12", dest="paired",
                                    help="Paired-end reads, support fastq/fastq.gz/fastq.bz2",
                                    metavar="")
    basic_option_group.add_argument("-s", "--single", dest="single",
                                    help="Single-read, support fastq/fastq.gz/fastq.bz2", metavar="")
    basic_option_group.add_argument("-o", "--out", dest="out", help="Specify the result folder [default='auto']",
                                    default="auto", metavar="")

    basic_option_group.add_argument("-rcp", dest="cp_reference",
                                    help="Chloroplast reference genome, only support GenBank-format",
                                    metavar="<file|dir>")
    basic_option_group.add_argument("-rmito", dest="mito_reference",
                                    help="Mitochondrial reference genome, only support GenBank-format",
                                    metavar="<file|dir>")
    basic_option_group.add_argument("-rtfa", dest="target_reference_fa",
                                    help="References of target genes, only support fasta format", metavar="<file|dir>")
    basic_option_group.add_argument("-rtgb", dest="target_reference_gb",
                                    help="References of target genes, only support GenBank format", metavar="<file|dir>")


    #高级参数部分
    advanced_option_group = parser.add_argument_group(title="Advanced option")
    advanced_option_group.add_argument("-n", "--number",
                                       help="The number of rows of raw data.\n"
                                            "If you want to use all the data, you can set as ['all']",
                                       default=10000000, metavar="")
    advanced_option_group.add_argument("-k", "--kmer", dest="kmer", help="Size of kmer  [default = 31]", default=31,
                                       type=int, metavar="")
    advanced_option_group.add_argument("-max", dest="max", help="The maximum length of genes [default=5000]", default=5000,
                                       type=int, metavar="")
    advanced_option_group.add_argument("-min", dest="min", help="The minimum length of genes [default = 300]", default=300,
                                       type=int, metavar="")
    advanced_option_group.add_argument("-t", "--thread",
                                       help="The number of threads [default = 'auto']",
                                       default="auto", metavar="")
    advanced_option_group.add_argument("-b", "--boundary", dest="soft_boundary",
                                       help="Extend the length to both sides of the gene while extracting  genes from  Genbank file [default=75]",
                                       default=75, type=int, metavar="")

    advanced_option_group.add_argument("-sf",dest="reference_filtering",help="Select the reference sequences to reduce the computation.\n"
                                                                          "s1: do nothing;\n"
                                                                          "s2,3,4: only use the reference sequence with the shortest/median/longest length;\n"
                                                                          "s5: remove sequences with abnormal length."
                                                                          "[default = 's1']",

                                       default="s1",metavar="")



    #梯度逼近部分
    gradient_approximation_group=parser.add_argument_group(title="Gradient approximation option")


    gradient_approximation_group.add_argument("-in","--iterative_number", dest="iterative_number",
                                        help="Specify the number of iterative loop to gradually approximate the best results",type=int,
                                              metavar="")

    #自展参数部分
    bootstrap_option_group = parser.add_argument_group(title="Bootstrap option")
    bootstrap_option_group.add_argument("-bn","--bootstrap_number",dest="bootstrap_number",type=int,
                                        help="Specify the bootstrap number.Evaluate the results based on the bootstrap method",
                                        metavar="")
    args = parser.parse_args()
    main(args)

















