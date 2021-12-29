#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:53
# @Author  : xiepulin
# @File    : verify_parameters.py
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


###############################################
'''
检测系统环境部分
(1)python版本
'''


##############################################



def check_python_version():
    this_python = sys.version_info[:2]
    min_version = (3, 6)
    if this_python < min_version:
        message_parts = [
            "GeneMiner does not work on Python {}.{}.".format(*this_python),            #*可以用来接受元组
            "The minimum supported Python version is {}.{}.".format(*min_version),
            "Please download the eligible version from https://www.python.org/.".format(*this_python)]
        print("ERROR: " + " ".join(message_parts))
        gv.set_value("my_gui_flag", 0)
        sys.exit()


##############################################################
'''
第一部分 检测各项参数是否正确,正确后打印参数
(1)线程/进程数量的选择数量 -t
(2)kmer 大小的检测 -k  
2.2wordsize 大小的检测       过滤的kmer
(3)软边界大小检测 -b
(4)检测基因max_length,min_length 
(5)参考基因组格式是否对应  2.1 -rtfa fasta;  -rtgb/-rcp/-rmito genebank格式  2.2选择了哪些参考基因组  2.3参考基因组是否有重名文件
(6)数据量大小的校验 -n
(7) 原始数据校验 即输入的序列是否是真实存在的 -1 -2 /-12 / -s 
(8) 所有任务检测
(9)自展检测参数检验
(10)梯度逼近参数检验
(11)输出文件夹检测
(12)参考序列过滤策略 

'''
#############################################################


'''
1.检测线程数量默认为所有线程数量的1/4 默认不大于8
'''
def check_threads_number(thread):
    thread_number = 1  # 初始化，至少一个
    thread_number_all = multiprocessing.cpu_count()
    if thread == "auto":
        thread_number = int(thread_number_all * 0.25)
        if thread_number == 0:
            thread_number = thread_number + 1
        elif thread_number>=8:
            thread_number=8
        else:
            thread_number=thread_number
    else:
        thread=str(thread)
        if thread.isdigit():         #判断是否为纯数字
            thread = int(thread)
            if thread > thread_number_all:
                print("Number of threads exceed the  maximum, please check the -t parameter")
                gv.set_value("my_gui_flag", 0)
                sys.exit()
            elif thread <= 0:
                print("Number of threads shuold exceed  0, please check the -t parameter")
            else:
                thread_number = thread

        else:  #防止输入乱七八糟的东西
            print("The number of threads must be an integer greater than 0, please check the -t parameter")
            gv.set_value("my_gui_flag", 0)
            sys.exit()

    return thread_number

'''
2限定kmer的大小，[15，127]
'''
def check_kmer(kmer):
    if kmer==31:
        kmer=31
    elif kmer<=0 or kmer>=127:
        print("Kmer sizes range from 15 to 127, please check the -k parameter")
        gv.set_value("my_gui_flag", 0)
        sys.exit()
    else:
        kmer=kmer
    return  kmer

'''
2.2限定wordsize的大小[16，32]
'''
def check_wordsize(wordsize):
    max=32
    min=16
    if wordsize<min or wordsize>max:
        print("Wordsize should be between 16 and 32, please check the -w parameter")
        gv.set_value("my_gui_flag", 0)
        sys.exit()
    else:
        wordsize=wordsize
    return wordsize



'''
3 检测软边界 默认为75 范围为0~200之间
'''
def check_soft_boundary(soft_boundary):
    soft_boundary = int(soft_boundary)
    max = 200
    min = 0
    if soft_boundary > max or soft_boundary < min:
        print(
            "The length of the soft boundary is limited to 0 to 200, and the recommended length is 0.5 * reads_length, please check the -b parameter")
        gv.set_value("my_gui_flag", 0)
        sys.exit()
    else:
        soft_boundary = soft_boundary
    return soft_boundary

'''
4 检测长度
'''
def check_max_min_length(max_length,min_length):
    if max_length<= min_length:
        print("The maximum gene length should be greater than the minimum gene length, please check the -max parameter ")
        gv.set_value("my_gui_flag", 0)
        sys.exit()

    if max_length<=0:
        print("The maximum gene length should be greater than zero, please check the -max parameter")
        gv.set_value("my_gui_flag", 0)
        sys.exit()

    if min_length < 0:
        print("The minimum gene length should not be less than zero, please check the -min parameter")
        gv.set_value("my_gui_flag", 0)
        sys.exit()

'''
5参考基因组格式是否对应  
（1） -rtfa fasta;  -rtgb/-rcp/-rmito genebank格式  
（2）选择了哪些参考基因组  
（3）参考基因组是否有重名文件
'''
def check_ref_format(target_reference_fa, target_reference_gb,cp_reference, mito_reference):
    ref_list = []
    if target_reference_fa != None:
        ref_list.append(target_reference_fa)
    if target_reference_gb != None:
        ref_list.append(target_reference_gb)
    if cp_reference != None:
        ref_list.append(cp_reference)
    if mito_reference != None:
        ref_list.append((mito_reference))

    # 0代表文件夹 1代表文件
    if target_reference_fa in ref_list:
        flag = file_or_directory(target_reference_fa)
        if flag == 0:
            files = get_files(target_reference_fa)
            flag = 0
            Nonconforming_file = []  #记录不合格的文件
            for file in files:
                answer = is_fasta(file)  # answer 返回False/True
                if answer == False:
                    Nonconforming_file.append(file)
                    flag = flag + 1

            if flag != 0:
                print("references should be in fasta format,please check -rtfa parameter")
                Nonconforming_file = [os.path.basename(file) for file in Nonconforming_file]
                if len(Nonconforming_file) == 1:
                    Nonconforming_file = "".join(Nonconforming_file)
                    print("{} is not in fasta format".format(Nonconforming_file))
                else:
                    Nonconforming_file = ",".join(Nonconforming_file)
                    print("{} are not in fasta format".format(Nonconforming_file))
                gv.set_value("my_gui_flag", 0)
                sys.exit()

            file_name_list = [os.path.basename(file) for file in files]
            for i in file_name_list:
                flag = 0
                if file_name_list.count(i) > 1:
                    print("{0} appears repeatedly. Ensure that the file name is unique".format(i))
                    file_name_list.remove(i)
                    flag = flag + 1
            if flag != 0:
                gv.set_value("my_gui_flag", 0)
                sys.exit()

        if flag == 1:
            answer = is_fasta(target_reference_fa)
            if answer == False:
                print("references should be in fasta format,please check -rtfa parameter")
                print("{0} is not in fasta format".format(target_reference_fa))
                gv.set_value("my_gui_flag", 0)
                sys.exit()

    if target_reference_gb in ref_list:
        flag = file_or_directory(target_reference_gb)
        if flag == 0:
            files = get_files(target_reference_gb)
            flag = 0
            Nonconforming_file = []
            for file in files:
                answer = is_gb(file)  # answer 返回False/True
                if answer == False:
                    Nonconforming_file.append(file)
                    flag = flag + 1

            if flag != 0:
                print("references should be in GenBank-format, please check -rtgb parameter")
                Nonconforming_file = [os.path.basename(file) for file in Nonconforming_file]
                if len(Nonconforming_file) == 1:
                    Nonconforming_file = "".join(Nonconforming_file)
                    print("{} is not in GenBank-format".format(Nonconforming_file))
                else:
                    Nonconforming_file = ",".join(Nonconforming_file)
                    print("{} are not in GenBank-format".format(Nonconforming_file))
                gv.set_value("my_gui_flag", 0)
                sys.exit()

            file_name_list = [os.path.basename(file) for file in files]
            for i in file_name_list:
                flag = 0
                if file_name_list.count(i) > 1:
                    print("{0} appears repeatedly. Ensure that the file name is unique".format(i))
                    file_name_list.remove(i)
                    flag = flag + 1
            if flag != 0:
                gv.set_value("my_gui_flag", 0)
                sys.exit()

        if flag == 1:
            answer = is_gb(target_reference_gb)
            if answer == False:
                print("references should be in GenBank-format, please check -rtgb parameter")
                print("{0} is not in GenBank-format".format(target_reference_gb))
                gv.set_value("my_gui_flag", 0)
                sys.exit()

    if cp_reference in ref_list:
        flag = file_or_directory(cp_reference)
        if flag == 0:
            files = get_files(cp_reference)
            flag = 0
            Nonconforming_file = []
            for file in files:
                answer = is_gb(file)  # answer 返回False/True
                if answer == False:
                    Nonconforming_file.append(file)
                    flag = flag + 1

            if flag != 0:
                print("cp_reference should be in GenBank-format, please check -rcp parameter")
                Nonconforming_file = [os.path.basename(file) for file in Nonconforming_file]
                if len(Nonconforming_file) == 1:
                    Nonconforming_file = "".join(Nonconforming_file)
                    print("{} is not in GenBank-format ".format(Nonconforming_file))
                else:
                    Nonconforming_file = ",".join(Nonconforming_file)
                    print("{} are not in GenBank-format".format(Nonconforming_file))
                gv.set_value("my_gui_flag", 0)
                sys.exit()

            file_name_list = [os.path.basename(file) for file in files]
            for i in file_name_list:
                flag = 0
                if file_name_list.count(i) > 1:
                    print("{0} appears repeatedly. Ensure that the file name is unique".format(i))
                    file_name_list.remove(i)
                    flag = flag + 1
            if flag != 0:
                gv.set_value("my_gui_flag", 0)
                sys.exit()

        if flag == 1:
            answer = is_gb(cp_reference)
            if answer == False:
                print("cp_reference should be in GenBank-format, please check -rcp parameter")
                print("{0} is not in  GenBank-format".format(cp_reference))
                gv.set_value("my_gui_flag", 0)
                sys.exit()

    if mito_reference in ref_list:
        flag = file_or_directory(mito_reference)
        if flag == 0:
            files = get_files(mito_reference)
            flag = 0
            Nonconforming_file = []     #记录不合格的文件
            for file in files:
                answer = is_gb(file)  # answer 返回False/True
                if answer == False:
                    Nonconforming_file.append(file)
                    flag = flag + 1

            if flag != 0:
                print("mito_references should be in GenBank-format, please check -rmito parameter")
                Nonconforming_file = [os.path.basename(file) for file in Nonconforming_file]
                if len(Nonconforming_file) == 1:
                    Nonconforming_file = "".join(Nonconforming_file)
                    print("{} is not in GenBank-format".format(Nonconforming_file))
                else:
                    Nonconforming_file = ",".join(Nonconforming_file)
                    print("{} are not in GenBank-format".format(Nonconforming_file))
                gv.set_value("my_gui_flag", 0)
                sys.exit()

            file_name_list = [os.path.basename(file) for file in files]
            for i in file_name_list:
                flag = 0
                if file_name_list.count(i) > 1:
                    print("{0} appears repeatedly. Ensure that the file name is unique".format(i))
                    file_name_list.remove(i)
                    flag = flag + 1
            if flag != 0:
                gv.set_value("my_gui_flag", 0)
                sys.exit()

        if flag == 1:
            answer = is_gb(mito_reference)
            if answer == False:
                print("mito_references should be in GenBank-format, please check -rmito parameter")
                print("{0} is not in  GenBank-format".format(mito_reference))
                gv.set_value("my_gui_flag", 0)
                sys.exit()


'''
6 检测数据量大小 -n
根据后缀确定文件类型，根据用户输入行数，提取reads。产生data1.fq，data2.fq
目前做的简单处理， 将-12(paired_end) 和-s(single_read) -1 -2 都使用相同的算法处理
'''
def check_datasize(data_size):
    data_size=str(data_size) # 先把数字或者"all" 都转为字符串处理   防止'int' object has no attribute 'isdigit'


    min_data_size = 1000000
    if data_size.lower()=="all":
        ultimate_data_size=data_size.lower()
        return ultimate_data_size

    elif data_size.isdigit():
        data_size=int(data_size)
        if data_size < min_data_size:
            print(
                "Please check -n parameter. for better results, input data should not be less than {0} lines.".format(
                    min_data_size))
            gv.set_value("my_gui_flag", 0)
            sys.exit()
        else:
            ultimate_data_size = data_size - data_size % 100000  # 10w起步，保证是4的倍数
        return ultimate_data_size
    else:
        print("Please check -n parameter, it must be an integer greater than {} or 'all' ".format(min_data_size))
        gv.set_value("my_gui_flag", 0)
        sys.exit()


'''
7 检测 -1 -2  -12  -single 是否真实存在，是否正确选择
'''
def check_raw_data_type(data1, data2, paired, single):
    '''
    必须三选一，且只能三选一
    -1 -2 必须同时存在
    data1, data2, paired, single 必须真实存在
    '''
    if (data1 == None or data2 == None) and paired == None and single == None:
        print("Please select one of the parameters from [-1 -2], [-s] and  [-12]")
        gv.set_value("my_gui_flag", 0)
        sys.exit()
    if (data1 == None and data2 != None) or (data1 != None and data2 == None):
        print("You must select both the -1 and -2 parameters")
        gv.set_value("my_gui_flag", 0)
        sys.exit()
    flag = 0
    raw_data_list = []  # 记录data1,data2,paired,single的参数选择
    if data1 != None and data2 != None:
        raw_data_list.append(data1)
        raw_data_list.append(data2)
        flag = flag + 1
    if paired != None:
        flag = flag + 1
        raw_data_list.append(paired)
    if single != None:
        flag = flag + 1
        raw_data_list.append(single)
    if flag != 1:
        print("Only one set of parameter can be selected from [-1 -2], [-s] and [-12]")
        gv.set_value("my_gui_flag", 0)
        sys.exit()

    # print(raw_data_list)
    flag = 0
    for i in raw_data_list:
        if is_exist(i) == 0:
            print("{}:the file does not exist".format(i))
            flag = flag + 1
    if flag != 0:
        gv.set_value("my_gui_flag", 0)
        sys.exit()


    '''
    检测原始数据的格式类型 fq .gz .bz2
    '''
    if data1 in raw_data_list and data2 in raw_data_list:
        if "fastq.gz" in data1 or "fq.gz" in data1 or ".gz" in data1:  # fastq.gz/fastq.bz2 包含了fastq所以要进一步判断
            type1 = 1
        elif "fastq.bz2" in data1 or "fq.bz2" in data1 or "bz2" in data1:
            type1 = 2
        elif "fastq" in data1 or ".fq" in data1:
            type1 = 3
        else:
            print("-1  need .fq/.fastq/.fq.gz /.fq.bz2 /.fastq.gz /.fastq.bz2 as the suffix")
            gv.set_value("my_gui_flag", 0)
            sys.exit(1)

        if "fastq.gz" in data2 or "fq.gz" in data2 or ".gz" in data2:  # fastq.gz/fastq.bz2 包含了fastq所以要进一步判断
            type2 = 1
        elif "fastq.bz2" in data2 or "fq.bz2" in data2 or "bz2" in data2:
            type2 = 2
        elif "fastq" in data2 or ".fq" in data2:
            type2 = 3
        else:
            print(" -2 need .fq/.fastq/.fq.gz /.fq.bz2 /.fastq.gz /.fastq.bz2 as the suffix")
            gv.set_value("my_gui_flag", 0)
            sys.exit(1)

        if type1 != type2:
            print("data1 reads and data2 reads should be the same format ")
            gv.set_value("my_gui_flag", 0)
            sys.exit(1)
        else:
            type = type1

    elif paired in raw_data_list:
        if "fastq.gz" in paired or "fq.gz" in paired or ".gz" in paired:  # fastq.gz/fastq.bz2 包含了fastq所以要进一步判断
            type = 1
        elif "fastq.bz2" in paired or "fq.bz2" in paired or "bz2" in paired:
            type = 2
        elif "fastq" in paired or ".fq" in paired:
            type = 3
        else:
            print("-12  need .fq/.fastq/.fq.gz /.fq.bz2 /.fastq.gz /.fastq.bz2 as the suffix")
            gv.set_value("my_gui_flag", 0)
            sys.exit(1)

    else:
        if "fastq.gz" in single or "fq.gz" in single or ".gz" in single:  # fastq.gz/fastq.bz2 包含了fastq所以要进一步判断
            type = 1
        elif "fastq.bz2" in single or "fq.bz2" in single or "bz2" in single:
            type = 2
        elif "fastq" in single or ".fq" in single:
            type = 3
        else:
            print("-s  need .fq/.fastq/.fq.gz /.fq.bz2 /.fastq.gz /.fastq.bz2 as the suffix")
            gv.set_value("my_gui_flag", 0)
            sys.exit(1)

    raw_data_dict = {type: raw_data_list}

    return raw_data_dict                              ##raw_data_dict={3:['gw1.fq']}}



'''
8 检测所有任务
'''

def check_whole_task(target_reference_fa, target_reference_gb, cp_reference, mito_reference):
    task_pool = []  # 任务池
    task_pool_exist = []  # 任务是否存在的池子
    if target_reference_fa == None and target_reference_gb == None and cp_reference == None and mito_reference == None:
        print("Please select at least one from -rtfa, -rtgb, -rcp and -rmito")
        gv.set_value("my_gui_flag", 0)
        sys.exit()
    if target_reference_fa != None:
        task = "target_gene_from_fa"
        task_pool.append(task)
        task_pool_exist.append(target_reference_fa)
    if target_reference_gb != None:
        task = "target_gene_from_gb"
        task_pool.append(task)
        task_pool_exist.append(target_reference_gb)
    if cp_reference != None:
        task = "cp_gene"
        task_pool.append(task)
        task_pool_exist.append(cp_reference)
    if mito_reference != None:
        task = "mito_gene"
        task_pool.append(task)
        task_pool_exist.append(mito_reference)

    flag = 0
    # print(task_pool)
    # print(task_pool_exist)
    for i in task_pool_exist:
        if is_exist(i) == 0:
            print("{} : does not exist".format(i))
            flag = flag + 1
    if flag != 0:
        gv.set_value("my_gui_flag", 0)
        sys.exit()


    return task_pool




'''
9.检测自展参数
'''
def check_bootstrap_parameter(bootstrap_number):
    if bootstrap_number==None or bootstrap_number=="None":
        flag = "No"
        bootstrap_number = "None"
        bootstrap_information = [flag, bootstrap_number]

    else:
        if  bootstrap_number <= 0:
            print("The number must be greater than 0, please check the -bn parameter")
            gv.set_value("my_gui_flag", 0)
            sys.exit()
        else:
            flag = "Yes"
            bootstrap_number = bootstrap_number
            bootstrap_information = [flag, bootstrap_number]

    return bootstrap_information


'''
10.检测梯度逼近参数
'''
def check_iterative_parameter(iterative_number):
    if iterative_number==None or iterative_number=="None":
        flag = "No"
        iterative_number = "None"
        iterative_information = [flag, iterative_number]

    else:
        if  iterative_number <= 0:
            print("The number must be greater than 0, please check the -in parameter")
            gv.set_value("my_gui_flag", 0)
            sys.exit()
        else:
            flag = "Yes"
            iterative_number = iterative_number
            iterative_information = [flag, iterative_number]

    return iterative_information





'''
11.检测输出文件夹
'''


'''
windows下额外的输出文件夹检测
'''
def check_out_dir_GUI(out,system):
    if system=="windows":
        #pysimplegui 会把输入数据转换为str格式
        if out=='' or out=="auto"  or out==None:  #out==None 针对图形界面用户reset之后，out_dir忘记输入的情况
            print("You should specify an output folder, please check the -o parameter")
            gv.set_value("my_gui_flag", 0)
            sys.exit()



def check_out_dir(out):
    nowTime = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
    result_name = "GM" + "_" + nowTime

    if out != "auto":
        if os.path.isdir(out):
            if len(os.listdir(out))==0:  #文件夹存在但里面为空是能够使用的
                out_dir_name=out
            else:
                print("{} already exists and there are files under the folder, please check the -o parameter".format(out))
                gv.set_value("my_gui_flag", 0)
                sys.exit(0)
        else:
            out_dir_name=out

    else:
        out_dir_name = result_name

    return out_dir_name



'''
12 参考序列过滤策略检测
'''

'''
两类数据，一类是txt，一类是s1,s2,s3,s4,s5
'''
def check_reference_filtering_strategy(sf):
    choose_data=["s1","s2","s3","s4","s5"]
    if sf ==None:
        sf= "s1"
        return sf

    file_flag = os.path.isfile(sf)
    if file_flag==True:
        flag=is_txt_file(sf)
        if flag==0:
            message = "You can choose one from ['s1', 's2', 's3', 's4', 's5'] or specify a plain text file with 'txt' as the suffix, please check the -sf parameter"
            print(message, flush=True)
            gv.set_value("my_gui_flag", 0)
            sys.exit()
        else:
            sf=sf

    elif file_flag==False:
        if sf in choose_data:
            sf=sf
        else:
            message = "You can choose one from ['s1', 's2', 's3', 's4', 's5'] or specify a plain text file with 'txt' as the suffix， please check the -sf parameter "
            print(message, flush=True)
            gv.set_value("my_gui_flag", 0)
            sys.exit()
    else:
        pass
    return  sf







# s=None
# m="s1"
# n="m2"
# ll="s2"
# lsl=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\1.txt"
# oo=check_reference_filtering_strategy(lsl)
# print(oo)




'''
参数校验合格之后，打印参数表
在extract_raw_data之前
'''
def print_parameter_information(parameter_information_dict):
    temp={}
    for key,value in parameter_information_dict.items():
        if parameter_information_dict[key] != None:
            temp[key] = value
    message_all = []
    for key, value in temp.items():
        message = "{0:<22}:{1:<}".format(key, value) #左对齐 22位宽  target_reference最长 16位宽
        message_all.append(message)

    symbol="-"
    print(symbol*22,flush=True)
    header="GeneMiner: a software for extracting phylogenetic markers from next generation sequencing data\n" \
           "Version: 1.0\n" \
           "Copyright (C) 2021 Pulin Xie\n" \
           "Please contact <xiepulin@stu.edu.scu.cn>, if you have any bugs or questions"

    print(header,flush=True)
    for i in message_all:
        # i=i.replace("_"," ")
        i=i.capitalize()#首字母大写
        print(i, flush=True)
    print(symbol * 22, flush=True)






