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
from openpyxl import load_workbook
from tqdm import tqdm
from concurrent import futures
import PySimpleGUI as sg

my_version = 'Version 1.0b build 20220122'
my_cite = 'Cite: https://github.com/sculab/Geneminer '
# current_path = os.path.abspath(__file__)  # 脚本当前路径
# father_path = os.path.abspath(os.path.dirname(current_path))  # 脚本的父目录
# sys.path.append(os.path.join(father_path, "lib"))  # filter ，assemble的软件存放在脚本同级目录下的lib

# cur_path = os.getcwd()
# sys.path.append(os.path.join(cur_path,"lib"))  # filter ，assemble的软件存放在脚本同级目录下的lib


cur_path = os.path.realpath(sys.argv[0])  # 脚本当前路径
father_path = os.path.dirname(cur_path)  # 脚本的父目录
# filter ，assemble的软件存放在脚本同级目录下的lib
sys.path.append(os.path.join(father_path, "lib"))

from lib.pack_results import pack_the_results
from lib.bootstrap_verify_pipeline import *
from lib.iterative_pipeline import *
from lib.core_pipeline import *
from lib.build_reference_database import *
from lib.extract_data import *
from lib.verify_parameters import *
from lib.basic import *
gv._init_()  # 在basic 中已经申明过了


def get_filter_assemble_software_path():
    cur_path = os.path.realpath(sys.argv[0])  # 脚本当前路径
    cur_path = os.path.dirname(cur_path)  # 脚本的父目录,father_path 覆盖

    filter_path = os.path.join(cur_path, "lib", "filter")
    minia_path = os.path.join(cur_path, "lib", "minia")
    # muscle_path=os.path.join(cur_path,"lib","muscle")  #muscle5
    muscle_path = os.path.join(cur_path, "lib", "muscle3")  # muscle3
    makeblastdb_path = os.path.join(cur_path, "lib", "makeblastdb")
    blastn_path = os.path.join(cur_path, "lib", "blastn")

    return [filter_path, minia_path, muscle_path, makeblastdb_path, blastn_path]


def get_configuration_information():
    configuration_information = {}
    # 一级目录
    mito_dir = "mito_genes"
    cp_dir = "cp_genes"
    tfa_dir = "target_genes_from_fa"
    tgb_dir = "target_genes_from_gb"
    data1 = "data1.fq"
    data2 = "data2.fq"
    results_information_excel = "results_information.xlsx"
    # 二级目录
    reference_database = "reference_database"
    filtered_out = "filtered_out"
    assembled_out = "assembled_out"
    assembled_log = "assembled_log"
    iterated_out = "iterated_out"
    bootstrap_out = "bootstrap_out"
    GM_results = "GM_results"
    blastn_out = "blastn_out"

    # 位置信息
    path = get_filter_assemble_software_path()
    filter_path = path[0]  # filter绝对路径
    assemble_path = path[1]  # minia绝对路径
    muscle_path = path[2]  # muscle绝对路径
    makeblastdb_path = path[3]  # makeblastdb_path 绝对路径
    blastn_path = path[4]  # blastn 绝对路径
    # 其他信息
    my_software_name = "GM"
    system = get_platform()
    whole_log = "log.txt"
    # bootstarp  将所有Bootstrap结果写在一起的fasta
    bootstrap_data_set = "bootstrap_data_set.fasta"
    bootstrap_concensus = "bootstrap_concensus.fasta"  # 导出共识序列
    # 构建配置信息
    configuration_information["mito_dir"] = mito_dir
    configuration_information["cp_dir"] = cp_dir
    configuration_information["tfa_dir"] = tfa_dir
    configuration_information["tgb_dir"] = tgb_dir
    configuration_information["data1"] = data1
    configuration_information["data2"] = data2
    configuration_information["results_information_excel"] = results_information_excel
    configuration_information["reference_database"] = reference_database
    configuration_information["filtered_out"] = filtered_out
    configuration_information["assembled_out"] = assembled_out
    configuration_information["assembled_log"] = assembled_log
    configuration_information["iterated_out"] = iterated_out
    configuration_information["bootstrap_out"] = bootstrap_out
    configuration_information["GM_results"] = GM_results
    configuration_information["blastn_out"] = blastn_out
    configuration_information["filter_path"] = filter_path
    configuration_information["assemble_path"] = assemble_path
    configuration_information["muscle_path"] = muscle_path
    configuration_information["makeblastdb_path"] = makeblastdb_path
    configuration_information["blastn_path"] = blastn_path
    configuration_information["my_software_name"] = my_software_name
    configuration_information["system"] = system
    configuration_information["whole_log"] = whole_log
    configuration_information["bootstrap_data_set"] = bootstrap_data_set
    configuration_information["bootstrap_concensusu"] = bootstrap_concensus

    return configuration_information
# print(get_configuration_information())


def main(args):
    # signal.signal(signal.SIGINT, signal_handler)  #检测函数终止退出（ctrl+c）
    gv.set_value("my_gui_flag", 1)  # 用于判定GUI是否处于运行状态，1代表运行，0代表没有运行

    configuration_information = get_configuration_information()
    '''
    从程序外获得的参数信息
    '''
    data1 = args.data1
    data2 = args.data2
    paired = args.paired
    single = args.single
    target_reference_fa = args.target_reference_fa
    target_reference_gb = args.target_reference_gb
    cp_reference = args.cp_reference
    mito_reference = args.mito_reference
    out_dir = args.out
    max_length = args.max
    min_length = args.min
    soft_boundary = args.soft_boundary
    data_size = args.number
    kmer = args.kmer
    wordsize = args.wordsize
    thread = args.thread
    bootstrap_number = args.bootstrap_number
    iterative_number = args.iterative_number
    sf = args.reference_filtering

    data1 = get_absolute(data1)
    data2 = get_absolute(data2)
    paired = get_absolute(paired)
    paired = get_absolute(paired)
    single = get_absolute(single)
    target_reference_fa = get_absolute(target_reference_fa)
    target_reference_gb = get_absolute(target_reference_gb)
    cp_reference = get_absolute(cp_reference)
    mito_reference = get_absolute(mito_reference)

    # 校验信息
    check_python_version()
    system = get_platform().lower()
    check_out_dir_GUI(out_dir, system)
    out_dir = check_out_dir(out_dir)  # 输出文件夹检测
    out_dir = get_absolute(out_dir)
    gv.set_value("out_dir", out_dir)
    thread_number = check_threads_number(thread)  # 确定线程数量
    kmer = check_kmer(kmer)  # 检测kmer
    wordsize = check_wordsize(wordsize)  # 检测wordsize
    soft_boundary = check_soft_boundary(soft_boundary)  # 检测软边界
    check_max_min_length(max_length, min_length)  # 检测基因最大最小长度
    check_ref_format(target_reference_fa, target_reference_gb,
                     cp_reference, mito_reference)  # 检测参考序列格式
    sf = check_reference_filtering_strategy(sf)  # 检测参考序列过滤策略
    data_size = check_datasize(data_size)  # 检测数据量 大小
    # 第一次出现：检测原始数据类型，方便判断是否运行 （有可能是相对路径）
    raw_data_dict = check_raw_data_type(data1, data2, paired, single)

    bootstrap_information = check_bootstrap_parameter(
        bootstrap_number)  # 检测自展检测参数["Yes",100]
    bootstrap_number = bootstrap_information[1]

    iterative_information = check_iterative_parameter(
        iterative_number)  # 检测梯度逼近参数 ["Yes",10]
    iterative_number = iterative_information[1]

    task_pool = check_whole_task(
        target_reference_fa, target_reference_gb, cp_reference, mito_reference)

    parameter_information_dict = {"project name": out_dir,
                                  "data1": data1, "data2": data2, "paired": paired, "single": single,
                                  "reference (fa)": target_reference_fa, "reference (gb)": target_reference_gb,
                                  "reference (cp)": cp_reference, "reference (mito)": mito_reference,
                                  "max length": max_length, "min length": min_length,
                                  "soft boundary": soft_boundary, "data size": data_size,
                                  "k-mer": kmer, "threads": thread_number,
                                  "word-size": wordsize,
                                  "bootstrap": bootstrap_information[0],
                                  "bootstrap number": bootstrap_information[1],
                                  "iterative": iterative_information[0],
                                  "iterative number": iterative_information[1]}

    print_parameter_information(parameter_information_dict)
    extract_raw_data(out_dir, raw_data_dict, data1, data2, paired, single,
                     data_size, configuration_information)  # 检查原始数据格式是否正确，提取原始数据,产生out_dir文件夹

    if "target_gene_from_fa" in task_pool:
        my_fasta = Extract_reference_fasta(
            configuration_information, "tfa", out_dir, target_reference_fa, soft_boundary, sf, max_length, min_length)
        my_fasta.extract_reference_from_fasta()
        my_fasta.filter_refrence()

        my_core_pipeline = Pipeline(configuration_information, "tfa", out_dir,
                                    thread_number, kmer, wordsize, max_length, min_length, "yes")
        my_core_pipeline.filter_reads_parallel()
        my_core_pipeline.assemble_reads_parallel()
        my_core_pipeline.check_contigs_parallel()
        my_core_pipeline.record_log_parallel()

        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(
                configuration_information, "tfa", out_dir, thread_number, kmer, wordsize, iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)

        if bootstrap_number != "None":
            my_bootstrap_pipeline = Bootstrap_pipeline(
                configuration_information, "tfa", out_dir, thread_number, kmer, bootstrap_number, max_length, min_length, "no")
            my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("tfa", configuration_information,
                         out_dir, bootstrap_number, iterative_number)

    # Genbank格式目标基因流程
    if "target_gene_from_gb" in task_pool:
        my_fasta = Extract_reference_fasta(
            configuration_information, "tgb", out_dir, target_reference_gb, soft_boundary, sf, max_length, min_length)
        my_fasta.extract_reference_from_gb()
        my_fasta.filter_refrence()

        my_core_pipeline = Pipeline(configuration_information, "tgb", out_dir,
                                    thread_number, kmer, wordsize, max_length, min_length, "yes")
        my_core_pipeline.filter_reads_parallel()
        my_core_pipeline.assemble_reads_parallel()
        my_core_pipeline.check_contigs_parallel()
        my_core_pipeline.record_log_parallel()

        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(
                configuration_information, "tgb", out_dir, thread_number, kmer, wordsize, iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)

        if bootstrap_number != "None":
            my_bootstrap_pipeline = Bootstrap_pipeline(
                configuration_information, "tgb", out_dir, thread_number, kmer, bootstrap_number, max_length, min_length, "no")
            my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("tgb", configuration_information,
                         out_dir, bootstrap_number, iterative_number)

    # 叶绿体基因流程
    if "cp_gene" in task_pool:
        my_fasta = Extract_reference_fasta(
            configuration_information, "cp", out_dir, cp_reference, soft_boundary, sf, max_length, min_length)
        my_fasta.extract_reference_from_gb()

        my_core_pipeline = Pipeline(configuration_information, "cp", out_dir,
                                    thread_number, kmer, wordsize, max_length, min_length, "yes")
        my_core_pipeline.filter_reads_parallel()
        my_core_pipeline.assemble_reads_parallel()
        my_core_pipeline.check_contigs_parallel()
        my_core_pipeline.record_log_parallel()

        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(configuration_information, "cp", out_dir, thread_number, kmer, wordsize,
                                                           iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)

        if bootstrap_number != "None":
            my_bootstrap_pipeline = Bootstrap_pipeline(
                configuration_information, "cp", out_dir, thread_number, kmer, bootstrap_number, max_length, min_length, "no")
            my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("cp", configuration_information,
                         out_dir, bootstrap_number, iterative_number)

    # 线粒体基因流程
    if "mito_gene" in task_pool:
        mito_fasta = Extract_reference_fasta(
            configuration_information, "mito", out_dir, mito_reference, soft_boundary, sf, max_length, min_length)
        mito_fasta.extract_reference_from_gb()
        mito_fasta.filter_refrence()

        mito_pipeline = Pipeline(configuration_information, "mito", out_dir,
                                 thread_number, kmer, wordsize, max_length, min_length, "yes")
        mito_pipeline.filter_reads_parallel()
        mito_pipeline.assemble_reads_parallel()
        mito_pipeline.check_contigs_parallel()
        mito_pipeline.record_log_parallel()

        if iterative_number != "None":
            my_iterative_pipeline = Run_iterative_pipeline(configuration_information, "mito", out_dir, thread_number, kmer, wordsize,
                                                           iterative_number, max_length, min_length, "no")
            iterative_information = my_iterative_pipeline.prepare_iterative_data()
            my_iterative_pipeline.run_iterative_parallel(iterative_information)

        if bootstrap_number != "None":
            my_bootstrap_pipeline = Bootstrap_pipeline(
                configuration_information, "mito", out_dir, thread_number, kmer, bootstrap_number, max_length, min_length, "no")
            my_bootstrap_pipeline.run_bootstrap_pipeline()
        pack_the_results("mito", configuration_information,
                         out_dir, bootstrap_number, iterative_number)

    gv.set_value("my_gui_flag", 0)  # 运行结束，释放


def geneminer_GUI():
    sg.theme('Default 1')
    # 图标
    icon = b'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'
    icon_48 = b'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'
    input_frame = [
        # -1
        [sg.Text('Data1:', size=(11, 1), justification='right', font=("Arial", 12)),
         sg.Input(key='-1-', size=(20, 1), font=("Arial", 12), expand_x=True,
                  tooltip='One end of the paired-end reads, support fastq/fastq.gz/fastq.bz2',
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Arial", 12), target='-1-', file_types=(("-1", "*.fq"), ("-1", "*.fastq"), ("-1", "*.fq.gz"), ("-1", "*.fastq.gz"), ("-1", "*.fq.bz2"), ("-1", "*.fastq.bz2"), ("-1", "*.Fq")))],

        # -2
        [sg.Text('Data2:', size=(11, 1), justification='right', font=("Arial", 12)),
         sg.Input(key='-2-', size=(20, 1), font=("Arial", 12), expand_x=True,
                  tooltip="Another end of the  paired-end reads, support fastq/fastq.gz/fastq.bz2",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Arial", 12), target='-2-', file_types=(("-2", "*.fq"), ("-2", "*.fastq"), ("-2", "*.fq.gz"), ("-2", "*.fastq.gz"), ("-2", "*.fq.bz2"), ("-2", "*.fastq.bz2")))],

        #-Paired -12
        [sg.Text('Paired reads:', size=(11, 1), justification='right', font=("Arial", 12)),
         sg.Input(key='-12-', size=(20, 1), font=("Arial", 12), expand_x=True,
                  tooltip='Paired-end reads,support fastq/fastq.gz/fastq.bz2',
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Arial", 12), target='-12-', file_types=(("-12", "*.fq"), ("-12", "*.fastq"), ("-12", "*.fq.gz"), ("-12", "*.fastq.gz"), ("-12", "*.fq.bz2"), ("-12", "*.fastq.bz2")))],

        #-Single -s
        [sg.Text('Single reads:', size=(11, 1), justification='right', font=("Arial", 12)),
         sg.Input(key='-s-', size=(20, 1), font=("Arial", 12), expand_x=True,
                  tooltip="Single-reads, support fastq/fastq.gz/fastq.bz2",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Arial", 12), target='-s-', file_types=(("-s", "*.fq"), ("-s", "*.fastq"), ("-s", "*.fq.gz"), ("-s", "*.fastq.gz"), ("-s", "*.fq.bz2"), ("-s", "*.fastq.bz2")))]
    ]
    reference_frame = [
        # -rtfa
        [sg.Text('Ref. (fasta):', size=(11, 1), justification='right', font=("Arial", 12)),
         sg.Input(key='-rtfa-', size=(20, 1), font=("Arial", 12), expand_x=True,
                  tooltip="Reference sequences, only support fasta format",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Arial", 12), target='-rtfa-', file_types=(("rtfa", "*.fasta"), ("rtfa", "*.fa"), ("rtfa", "*.fas"), ("rtfa", "*.FASTA"), ("rtfa", "*.FAS"), ("rtfa", "*.FA"))), sg.FolderBrowse("Folder", font=("Arial", 12), target='-rtfa-')],

        # -rtgb
        [sg.Text('Ref. (gb):', size=(11, 1), justification='right', font=("Arial", 12)),
         sg.Input(key='-rtgb-', size=(20, 1), font=("Arial", 12), expand_x=True,
                  tooltip="Reference sequences,only support GenBank format",
                  readonly=False, background_color='#BAD9BC'),
         sg.FileBrowse("File", font=("Arial", 12), target='-rtgb-', file_types=(("rtgb", "*.gb"),)), sg.FolderBrowse("Folder", font=("Arial", 12), target='-rtgb-')],

        # # -rcp
        # [sg.Text('Chlo. genome:', size=(11, 1), justification='right', font=("Arial", 12)),
        #  sg.Input(key='-rcp-', size=(20, 1), font=("Arial", 12), expand_x=True,
        #           tooltip="Chloroplast reference genome, only support GenBank format",
        #           readonly=False, background_color='#BAD9BC'),
        #  sg.FileBrowse("File", font=("Arial", 12), target='-rcp-', file_types=(("rcp", "*.gb"),)), sg.FolderBrowse("Folder", font=("Arial", 12), target='-rcp-')],

        # # rmito
        # [sg.Text('Mito. genome:', size=(11, 1), justification='right', font=("Arial", 12)),
        #  sg.Input(key='-rmito-', size=(20, 1), font=("Arial", 12), expand_x=True,
        #           tooltip="Mitochondrial reference genome, only support GenBank format",
        #           readonly=False, background_color='#BAD9BC'),
        #  sg.FileBrowse("File", font=("Arial", 12), target='-rmito-', file_types=(("rmito", "*.gb"),)), sg.FolderBrowse("Folder", font=("Arial", 12), target='-rmito-')]

    ]
    out_frame = [
        #-OUT -o
        [sg.Text('Output Folder:', size=(11, 1), justification='right', font=("Arial", 12)),
         sg.Input(key='-o-', size=(20, 1), font=("Arial", 12), expand_x=True,
                  tooltip="Specify the result folder ", readonly=False, background_color='#BAD9BC'),
         sg.FolderBrowse("Folder", font=("Arial", 12), target="-o-")]]
    filter_frame = [[sg.Text('Rows:', size=(11, 1), justification='right', font=("Arial", 12)),
                     sg.Input("10000000", key='-n-', size=(11, 1), font=("Arial", 12), expand_x=False, background_color = '#BFE9FF',
                              tooltip="The number of rows of raw data from NGS data.\nSet to 'all' if you want to use all the data."),
                     # -wordsize -kmer
                     sg.Text('Word-size:', size=(8, 1),
                             justification='right', font=("Arial", 12)),
                     sg.Input(16, key='-w-', size=(11, 1), font=("Arial", 12), expand_x=False, background_color = '#BFE9FF',
                              tooltip="The length of a word-size  [default=16]")]]
    assemble_frame = [[sg.Text('K-mer:', size=(11, 1), justification='right', font=("Arial", 12)),
                       sg.Input(31, key='-k-', size=(11, 1), font=("Arial", 12), expand_x=False, background_color = '#BFE9FF',
                                tooltip="The size of a k-mer [default=31]"),
                       # Boundaty threads
                       sg.Text('Boundary:', size=(8, 1),
                               justification='right', font=("Arial", 12)),
                       sg.Input(75, key='-b-', size=(11, 1), font=("Arial", 12), expand_x=False, background_color = '#BFE9FF',
                                tooltip="Extend the length to both sides of the gene while extracting genes from Genbank file [default=75]")],
                      # min max
                      [sg.Text('Min:', size=(11, 1), justification='right', font=("Arial", 12)),
                       sg.Input(300, key='-min-', size=(11, 1), font=("Arial", 12), expand_x=False, background_color = '#BFE9FF',
                                tooltip="The minimum length of a gene  [default=300]"),
                       sg.Text('Max:', size=(8, 1),
                               justification='right', font=("Arial", 12)),
                       sg.Input(5000, key='-max-', size=(11, 1), font=("Arial", 12), expand_x=False, background_color = '#BFE9FF',
                                tooltip="The maximum length of a gene [default=5000]"),
                       ],


                      [sg.Text('Threads:', size=(11, 1), justification='right', font=("Arial", 12)),
                       sg.Input("auto", key='-t-', size=(11, 1), font=("Arial", 12), expand_x=False, background_color = '#BFE9FF',
                                tooltip="Specify the number of threads you want to run [default=auto]"),
                       ]]
    verify_frame = [
        [sg.Text('Iterative:', size=(12, 1), justification='right', font=("Arial", 12)),
         sg.Input(2, key='-in-', size=(5, 1), font=("Arial", 12), expand_x=False, readonly=False, background_color = '#BFE9FF',
                  tooltip="Specify the number of iterative loop to filter reads correctly around indels"),
         sg.CB(' ', font=("Arial", 12), size=(3, 1),
               key='-cbin-', enable_events=True),


         sg.Text('Bootstrap:', size=(8, 1),
                 justification='right', font=("Arial", 12)),
         sg.Input(10, key='-bn-', size=(6, 1), font=("Arial", 12), expand_x=False, readonly=False, background_color = '#BFE9FF',
                  tooltip="Specify the bootstrap number. Evaluate the results based on the bootstrap method"),
         sg.CB(' ', font=("Arial", 12), size=(3, 1),
               key='-cbbn-', enable_events=True)
         ],
         # Reference filtering
        [sg.Text("Computation:", size=(12, 1), justification='right', font=("Arial", 12)),
         sg.Input("s1", key='-sf-', size=(20, 1), font=("Arial", 12), expand_x=True, background_color = '#BFE9FF',
                  tooltip="Selecting the method to reduce the computation\n"
                          "s1: do nothing;\n"
                          "s2/s3/s4:only use the reference sequence with the shortest/median/longest length;\n"
                          "s5: remove the reference sequences with abnormal length;\n"
                          "file: specify the name of reference sequences to use.\n"
                          "[default='s1']"), sg.FileBrowse("File", target="-sf-", font=("Arial", 12), file_types=(("txt", "*.txt"),))],

         ]
    left_col = [
        # 基本设定
        [sg.Text('Basic Option', size=(12, 1),
                 justification='left', font=("Arial", 12, 'bold'))],
        [sg.Frame('Data', layout=input_frame,
                  expand_x=True, font=("Arial", 12))],
        [sg.Frame('Reference', layout=reference_frame,
                  expand_x=True, font=("Arial", 12))],
        [sg.Frame('Output', layout=out_frame,
                  expand_x=True, font=("Arial", 12))],

        # 高级设定
        [sg.Text('Advanced Option', size=(20, 1),
                 justification='left', font=("Arial", 12, 'bold'))],
        [sg.Frame('Filter', layout=filter_frame,
                  expand_x=True, font=("Arial", 12))],
        [sg.Frame('Assembler', layout=assemble_frame,
                  expand_x=True, font=("Arial", 12))],
        [sg.Frame('Verifier', layout=verify_frame,
                  expand_x=True, font=("Arial", 12))],

        



        
        # 空行备用
        # [sg.T(' ', size=(40,1), justification='left', font=("Arial",12,'bold'))],
    ]

    right_col = [
        # [sg.Output(size=(100, 25), font=("Arial", 12, 'bold'), text_color="black", key='-OUTPUT-')],
        [sg.MLine(key='-ML1-', size=(85, 10), enter_submits=True, reroute_stdout=True, reroute_stderr=True,
                  enable_events=True, visible=False)],
        [sg.MLine(key='-ML2-', size=(85, 20), autoscroll=True,
                  font='TkFixedFont', background_color="#DEDBBD")],  # 等宽字体

        # 按钮
        [sg.Button('Run', size=(8, 1), font=("Arial", 12)),
         sg.Button('Results', size=(8, 1), font=("Arial", 12)),
         sg.Button('Clear', size=(8, 1), font=("Arial", 12)),
         sg.Button('Reset', size=(8, 1), font=("Arial", 12)),
         sg.Button('Close', size=(8, 1), font=("Arial", 12))]
    ]
    layout = [
        [sg.Image(icon_48),
         sg.Text('GeneMiner', justification="left", font=("Arial", 24, "bold"))
         ],
        [sg.T('Authors: Xie Pulin, Guo Yongling, Zhang Zhen, Yu Yan', font=("Arial", 12)),
         sg.T('Email: xiepulin@stu.scu.edu.cn', font=("Arial", 12))],
        [sg.Pane([sg.Column(left_col, element_justification='l', expand_x=True, expand_y=True),
                  sg.Column(right_col, element_justification='c', expand_x=True, expand_y=True)], orientation='h',
                 relief=sg.RELIEF_SUNKEN, key='-PANE-')],
        [
            sg.T(my_cite, font=("Arial", 12)),
            sg.Text(my_version, expand_x=True,
                    justification="right", font=("Arial", 12)),
        ],
    ]
    window = sg.Window('GeneMiner', layout, finalize=True,
                       icon=icon, resizable=True, use_default_focus=False)
    # window['-OUTPUT-'].expand(True, True, True)
    window["-ML2-"].expand(True, True, True)
    window['-PANE-'].expand(True, True, True)
    window['-bn-'].update(disabled=True)
    window['-in-'].update(disabled=True)

    '''
    ML1_length 统计原始长度
    count_clock  类似外部计时器
    ML1_ML1NEW_length 统计上一次统计长度和本次统计长度 若发生变化才会刷新
    '''
    ML1_length = -1
    # count_clock = 0
    # count_clock_force=0
    # ML1_ML1NEW_length = [0, 0]
    while True:
        event, values = window.read(timeout=250)
        if event == sg.WIN_CLOSED or event == 'Close':  # if user closes window or clicks cancel
            break

        if values['-cbbn-']:
            window['-bn-'].update(disabled=False)
        else:
            window['-bn-'].update(disabled=True)

        if values['-cbin-']:
            window['-in-'].update(disabled=False)
        else:
            window['-in-'].update(disabled=True)

        if event == 'Results':
            file_to_show = values['-o-']
            if file_to_show:
                try:
                    file_to_show = get_absolute(file_to_show)
                    if get_platform() == "windows":
                        os.startfile(file_to_show)
                    else:
                        subprocess.call(["open", file_to_show])
                except:
                    pass
            else:
                try:
                    file_to_show = gv.get_value("out_dir")
                    if get_platform() == "windows":
                        os.startfile(file_to_show)
                    else:
                        subprocess.call(["open", file_to_show])
                except:
                    pass

        if event == "Clear":
            window['-ML1-'].update("")
            window['-ML2-'].update("")
            values['-ML1-'] = ""
            ML1_length = -1

        # 1代表运行 0代表运行结束。
        if event == 'Run' and gv.get_value("my_gui_flag") == 0:
            # 基本参数部分
            if values["-1-"]:
                args.data1 = values["-1-"]
            else:
                args.data1 = None
            if values["-2-"]:
                args.data2 = values["-2-"]
            else:
                args.data2 = None
            if values["-12-"]:
                args.paired = values["-12-"]
            else:
                args.paired = None
            if values["-s-"]:
                args.single = values["-s-"]
            else:
                args.single = None
            if values["-rtfa-"]:
                args.target_reference_fa = values["-rtfa-"]
            else:
                args.target_reference_fa = None
            if values["-rtgb-"]:
                args.target_reference_gb = values["-rtgb-"]
            else:
                args.target_reference_gb = None
            # if values["-rcp-"]:
            #     args.cp_reference = values["-rcp-"]
            # else:
            #     args.cp_reference = None
            # if values["-rmito-"]:
            #     args.mito_reference = values["-rmito-"]
            # else:
            #     args.mito_reference = None
            if values["-o-"]:
                if not os.path.exists(values["-o-"]):
                    os.makedirs(values["-o-"])
                args.out = values["-o-"]
            else:
                sg.Popup("You must specify an output folder!",
                         title='Info', keep_on_top=True, font=("Arial", 12))
                continue

            # 高级参数部分
            if values["-k-"]:
                args.kmer = int(values["-k-"])
            else:
                args.kmer = 31
            if values["-w-"]:
                args.wordsize = int(values["-w-"])
            else:
                args.wordsize = 16
            if values["-min-"]:
                args.min = int(values["-min-"])
            else:
                args.min = 300
            if values["-max-"]:
                args.max = int(values["-max-"])
            else:
                args.max = 5000
            if values["-b-"]:
                args.soft_boundary = int(values["-b-"])
            if values["-n-"]:
                if values["-n-"].upper() == "ALL":
                    args.number = "all"
                elif values["-n-"].isdigit():
                    args.number = int(values["-n-"])
            else:
                args.number = 10000000
            if values["-t-"]:
                if values["-t-"].upper() == "AUTO":
                    args.thread = "auto"
                elif values["-t-"].isdigit():
                    args.thread = int(values["-t-"])
            else:
                args.thread = "auto"

            if values["-sf-"]:
                args.reference_filtering = values["-sf-"]
            else:
                args.reference_filtering = 's1'

            # 自展检测部分
            if values["-cbbn-"]:
                if values["-bn-"]:
                    args.bootstrap_number = int(values["-bn-"])
            else:
                args.bootstrap_number = None

            # 梯度逼近部分
            if values["-cbin-"]:
                if values["-in-"]:
                    args.iterative_number = int(values["-in-"])
            else:
                args.iterative_number = None

            threading.Thread(target=main, args=(args,), daemon=True).start()

        if event == "Reset":
            gv.set_value("my_gui_flag", 0)
            values['-ML1-'] = ""
            ML1_length = -1  # 初始值赋值为-1，直接打印
            # 基础参数部分
            window['-1-'].update("")
            window['-2-'].update("")
            window['-12-'].update("")
            window['-s-'].update("")
            window['-rtfa-'].update("")
            window['-rtgb-'].update("")
            # window['-rcp-'].update("")
            # window['-rmito-'].update("")
            window['-o-'].update("")
            # 高级参数部分
            window['-k-'].update("31")
            window['-w-'].update("16")
            window['-min-'].update("300")
            window['-max-'].update("5000")
            window['-b-'].update("75")
            window['-n-'].update("10000000")
            window['-t-'].update("auto")
            window['-sf-'].update("s1")
            window['-bn-'].update("10")
            window['-in-'].update("2")
            # 输出部分
            # 清空输出界面 ML2复制的ML1，所以清空ML1就等效于清空ML2(实时刷新的时候)
            # 清空输出界面 由于现在加了判断，ML2虽然复制ML1，但最后的时候不再进行实时复制，所以需要单独删除ML2
            window['-ML1-'].update("")
            window["-ML2-"].update("")

            '''
            清空参数空间
            '''
            args.data1 = None
            args.data2 = None
            args.paired = None
            args.single = None
            args.target_reference_fa = None
            args.target_reference_gb = None
            args.cp_reference = None
            args.mito_reference = None
            args.out = None
            args.kmer = 31
            args.wordsize = 16
            args.min = 300
            args.max = 5000
            args.soft_boundary = 75
            args.number = 10000000
            args.thread = "auto"
            args.reference_filtering = "s1"
            args.bootstrap_number = 10
            args.iterative_number = 2

        '''
        打印实时输出，需要判断后一次的输出是否比前一次的输出长，如果是才实时刷新
        '''
        if ML1_length != len(values['-ML1-']):
            window['-ML2-'].update(ML_str_re(values['-ML1-']))
            ML1_length = len(values['-ML1-'])

    window.close()


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal_handler)  # 检测函数终止退出（ctrl+c） #必须放在主程序中
    multiprocessing.freeze_support()  # windows上Pyinstaller打包多进程程序需要添加特殊指令
    gv.set_value("my_gui_flag", 0)  # 用于判定脚本是否跑完，还可以防止run双击覆盖事件
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
    # 原始数据输入部分
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

    # 高级参数部分
    advanced_option_group = parser.add_argument_group(title="Advanced option")
    advanced_option_group.add_argument("-n", "--number",
                                       help="The number of rows of raw data.\n"
                                            "If you want to use all the data, you can set as ['all']",
                                       default=10000000, metavar="")
    advanced_option_group.add_argument("-k", "--kmer", dest="kmer", help="Size of kmer  [default = 31]", default=31,
                                       type=int, metavar="")
    advanced_option_group.add_argument("-w", "--wordsize", dest="wordsize", help="length of a word-size [default=16]", default=16,
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

    advanced_option_group.add_argument("-sf", dest="reference_filtering", help="Select the reference sequences to reduce the computation.\n"
                                       "s1: do nothing;\n"
                                       "s2,3,4: only use the reference sequence with the shortest/median/longest length;\n"
                                       "s5: remove sequences with abnormal length."
                                       "[default = 's1']",

                                       default="s1", metavar="")

    # 梯度逼近部分
    gradient_approximation_group = parser.add_argument_group(
        title="Gradient approximation option")

    gradient_approximation_group.add_argument("-in", "--iterative_number", dest="iterative_number",
                                              help="Specify the number of iterative loop to gradually approximate the best results", type=int,
                                              metavar="")

    # 自展参数部分
    bootstrap_option_group = parser.add_argument_group(
        title="Bootstrap option")
    bootstrap_option_group.add_argument("-bn", "--bootstrap_number", dest="bootstrap_number", type=int,
                                        help="Specify the bootstrap number.Evaluate the results based on the bootstrap method",
                                        metavar="")
    args = parser.parse_args()
    # main(args)
    geneminer_GUI()
