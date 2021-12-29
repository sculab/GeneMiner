#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:47
# @Author  : xiepulin
# @File    : basic.py
# @Software: PyCharm


import sys
import subprocess
import os
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqRecord import  SeqRecord
from Bio.Seq import  Seq
from pathlib import Path
import platform
import numpy as np
import shutil
import re
import global_var as gv
################################################################
##################################################################
'''
流程部分
'''
###################################################################
###################################################################

'''
检测终止信号，优雅退出
'''
def signal_handler(signal, frame):
    print("GeneMiner has been terminated")
    gv.set_value("my_gui_flag", 0)
    sys.exit(0)
'''
判断文件是否为fasta格式
'''
def is_fasta(filename):
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            flag=any(fasta)  # any:如果都不是fasta 会提示false。但凡有一个是fasta,会提示True。
            # 但每次都是对一个文件进行判定，所以功能得以实现
    except:
        flag=False
    return flag


'''
判断文件是否为Genbank格式
'''
def is_gb(filename):
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "gb")
            flag= any(fasta)
    except:
        flag=False
    return flag


'''
判断文件是否是txt文本文件，布尔值
'''
def is_txt_file(file):
    try:
        flag=1
        suffix=os.path.splitext(file)[-1]
        if suffix!=".txt":
            flag=0
    except:
        flag=0
    return flag






'''
为了防止数据质量太差，导致某一步骤做不出来结果，无法接续运行
为了提高程序程序鲁棒性，判断文件是否存在，如果不存在或者文件大小为0，认定为无效返回0  有效返回1 
'''

def is_exist(file):
    if os.path.isfile(file):
        if os.path.getsize(file) > 0:
            flag = 1
        else:
            flag = 0

    elif os.path.isdir(file):
        files = get_files(file)
        if files==[]:
            flag=0
        else:
            flag = 1
            for i in files:
                if os.path.getsize(i) > 0:
                    continue
                else:
                    flag = 0
                    break
    else:
        flag=0

    return flag


'''
获取目录下所有文件
'''
def get_files(ref):
    file_path_list = []
    for root, dirs, files in os.walk(ref):
        for file in files:
            file_path = os.path.join(root, file)
            file_path_list.append(file_path)
    return file_path_list

'''
新建文件夹
'''
def dir_make(out_dir):
    if os.path.exists(out_dir) and len(os.listdir(out_dir)) > 0:
        print("\nerror! {}：the folder has already existed,and there are files under the folder.".format(out_dir))
        sys.exit()
    #针对于用户先创建文件夹，再使用
    elif os.path.exists(out_dir) and len(os.listdir(out_dir))==0:
        pass
    else:
        os.makedirs(out_dir)



def  cutting_line(message=""):
    element = ["-"]
    cutting_line_number = 85
    element_all = element * cutting_line_number
    length = len(message)
    start = int((cutting_line_number - length) / 2)

    for i in range(len(message)):
        element_all[i + start] = message[i]
    element_all = "".join(element_all)
    print(element_all)
    return element_all





'''
判断输入参考基因组是文件还是文件夹，0代表文件夹，1代表文件
'''
def file_or_directory(ref):
    if os.path.isdir(ref):
        files = os.listdir(ref)
        if len(files) == 0:
            print("the input  reference genome folder is empty")
            sys.exit()
        else:
            flag = 0
            # print("this is a dir")
    elif os.path.isfile(ref):
        size = os.path.getsize(ref)
        if size == 0:
            print("The input reference genome file is empty")
            sys.exit()
        else:
            flag = 1
            # print("this is a file")
    else:
        sys.exit()
    return flag

#########################################################################
#########################################################################
'''
reference 过滤 
(build_reference_database)
策略一：不做处理
策略二：仅保留最长序列
策略三：仅保留中位数长度的序列
策略四：箱线图剔除异常长度的序列
策略五：根据一致度聚类
'''
#########################################################################
#########################################################################






'''
策略2：最短序列
'''
def get_shortest_sequence(sequence_path):
    sequence_list = []

    # 取最长的seq （结果序列）
    for rec in SeqIO.parse(sequence_path, "fasta"):
        temp_sequence = rec.seq
        sequence_list.append(temp_sequence)
    if len(sequence_list) == 1:
        seq_sequence = sequence_list[0]
    else:
        length_list = [len(i) for i in sequence_list]
        index = length_list.index(min(length_list))
        seq_sequence = sequence_list[index]

    my_record=[]
    for rec in SeqIO.parse(sequence_path,"fasta"):
        if len(rec.seq)>len(seq_sequence):
            continue
        else:
            id = rec.id
            description=rec.description
            seq=rec.seq
            record=SeqRecord(id=id,seq=seq,description=description)
            my_record.append(record)
            break

    SeqIO.write(my_record,sequence_path,"fasta")

'''
策略3：中位数序列
'''
def get_median_sequence(sequence_path):
    sequence_list = []

    # 取中位数长度的seq （结果序列）。
    # 若为奇数，必能取中位数长度的序列;若为偶数，先排序，删除掉最短的序列，再找出中位数长度的序列,即倾向于越长越好
    for rec in SeqIO.parse(sequence_path, "fasta"):
        temp_sequence = rec.seq
        sequence_list.append(temp_sequence)


    if len(sequence_list) == 1:
        seq_sequence = sequence_list[0]
    else:
        length_list = [len(i) for i in sequence_list]

        sequence_list=sorted(sequence_list,key=lambda i:len(i),reverse=False)
        length_list=sorted(length_list,reverse=False)

        if len(length_list)%2!=0:
            median=np.median(length_list)
        else:
            sequence_list.remove(sequence_list[0])
            length_list.remove(length_list[0])
            median = np.median(length_list)

        median=int(median)           #前面的数据结构是Numpy的
        index = length_list.index(median)
        seq_sequence = sequence_list[index]



    my_record = []
    for rec in SeqIO.parse(sequence_path, "fasta"):
        if len(rec.seq) < len(seq_sequence):
            continue
        else:
            id = rec.id
            description = rec.description
            seq = rec.seq
            record = SeqRecord(id=id, seq=seq, description=description)
            my_record.append(record)
            break

    SeqIO.write(my_record, sequence_path, "fasta")

'''
策略4：最长序列
'''
def get_longest_sequence(sequence_path):
    sequence_list = []

    # 取最长的seq （结果序列）
    for rec in SeqIO.parse(sequence_path, "fasta"):
        temp_sequence = rec.seq
        sequence_list.append(temp_sequence)
    if len(sequence_list) == 1:
        seq_sequence = sequence_list[0]
    else:
        length_list = [len(i) for i in sequence_list]
        index = length_list.index(max(length_list))
        seq_sequence = sequence_list[index]

    my_record=[]
    for rec in SeqIO.parse(sequence_path,"fasta"):
        if len(rec.seq)<len(seq_sequence):
            continue
        else:
            id = rec.id
            description=rec.description
            seq=rec.seq
            record=SeqRecord(id=id,seq=seq,description=description)
            my_record.append(record)
            break

    SeqIO.write(my_record,sequence_path,"fasta")

def get_longest_sequence_only(sequence_path):
    sequence_list = []

    # 取最长的seq （结果序列）
    for rec in SeqIO.parse(sequence_path, "fasta"):
        temp_sequence = str(rec.seq)
        sequence_list.append(temp_sequence)
    if len(sequence_list) == 1:
        seq_sequence = sequence_list[0]
    else:
        length_list = [len(i) for i in sequence_list]
        index = length_list.index(max(length_list))
        seq_sequence = sequence_list[index]
    return  seq_sequence


'''
策略5：
利用箱线图，得到序列的的上界长度和下界长度
输入参考序列的路径，返回上界长度和下界长度
'''
def box_plots_bottom_and_top(reference_path):
    length_list=[]
    for rec in SeqIO.parse(reference_path, "fasta"):
        sequence = rec.seq
        length = len(sequence)
        length_list.append(length)
    q1=np.percentile(length_list,(25))  #lower_quartile
    q3=np.percentile(length_list,(75))   #upper_quartile
    iqr=q3-q1    #四分位间距
    top=q3+1.5*iqr
    bottom=q1-1.5*iqr
    return  [bottom,top]

'''
利用箱线图，将异常长度的reference剔除.
要求输入参考序列的路径
'''
def box_plots_filter_reference(reference_path):
    [bottom,top]=box_plots_bottom_and_top(reference_path)

    my_records=[]
    for rec in SeqIO.parse(reference_path, "fasta"):
        id = rec.name
        sequence = rec.seq
        description = rec.description
        if description=="":
            description=" "

        length=len(sequence)
        if length>=bottom and length<=top:  #保留正常值==过滤异常值
            record=SeqRecord(id=id,seq=sequence,description=description)
            my_records.append(record)
    if my_records!=[]:
        SeqIO.write(my_records,reference_path,"fasta")


'''
根据用户提供的文本文件中的名字  找到对应的reference
'''
def get_specify_txt_reference(reference_path, txt_file):
    information = []
    with open(txt_file, "r") as f:
        for line in f.readlines():
            line = line.rstrip()
            information.append(line)

    my_records = []
    for rec in SeqIO.parse(reference_path, "fasta"):
        name = rec.name
        if name in information:
            sequence = rec.seq
            description = rec.description
            name = rec.name
            record = SeqRecord(id=name, seq=sequence, description=description)
            my_records.append(record)

    if my_records != []:
        SeqIO.write(my_records, reference_path, "fasta")
















'''
计算参考基因组的大小
返回平均长度，用于过滤太长或者太短的序列;返回参考基因组的第一条序列，用于判定方向是否和参考基因组一致
[average_length, sequence_result]
'''
def calculate_ref_size(ref):
    flag = file_or_directory(ref)
    # 文件
    if flag == 1:
        length = []
        number = 0
        sum = 0
        sequence_result = ""
        for rec in SeqIO.parse(ref, "fasta"):
            sequence = rec.seq
            length.append(len(sequence))

            if number == 0:  # 只记录第一条
                sequence_result = rec.seq
                sequence_result = str(sequence_result)
            number = number + 1

        for i in length:
            sum = sum + i
        average_length = int(sum / number)
        return [average_length, sequence_result]
    # 文件夹的形式
    if flag == 0:
        files = get_files(ref)
        length = []
        number = 0
        sum = 0
        sequence_result = ""
        for file in files:
            for rec in SeqIO.parse(file, "fasta"):
                sequence = rec.seq
                length.append(len(sequence))

                if number == 0:  # 只记录第一条
                    sequence_result = rec.seq
                    sequence_result = str(sequence_result)
                number = number + 1
        for i in length:
            sum = sum + i
        average_length = int(sum / number)

        # print(length)
        return [average_length, sequence_result]



'''
双序列全局比对，采用动态规划算法。返回双序列全局比对得分，可以用来判断 序列方向和参考序列方向 是否一致
'''
def pairwise(seq1, seq2):
    seq1=str(seq1)
    seq2=str(seq2)

    seq1 = str.upper(seq1)
    seq2 = str.upper(seq2)
    alignment = pairwise2.align.globalms(seq1, seq2, 1, -1, -3, -2) #同blast准则 75%保守度
    # alignment = pairwise2.align.localms(seq1, seq2, 2, -1, -2, -0.5)
    score = alignment[0][2]
    return score


def cut_align(seq, ref):
    seq=str(seq) #防止输入的是Seq格式（biopython）
    ref=str(ref)
    seq = str.upper(seq)
    ref = str.upper(ref)
    alignment = pairwise2.align.localms(seq, ref, 1, -1, -3, -2)
    seq_align = alignment[0][
        0]  # 第一条双序列比对结果  PMSK_align ,此时和ref_align相同长度  Alignment(seqA='-TCGAAAAAAAAAAAAAAA', seqB='ATCGATCG-----------', score=8.0, start=1, end=5)
    start_site = alignment[0][3]
    end_site = alignment[0][4]
    cut_align_sequence = seq_align[start_site:end_site]

    cut_align_sequence=cut_align_sequence.replace("-","")
    ref_length = len(ref)
    if len(cut_align_sequence) < 0.8 * ref_length:
        trimmed_contigs = ""
    elif len(cut_align_sequence) >= 0.8 * ref_length and len(cut_align_sequence) <= 1.2 * ref_length:  # 剪切的序列既不能太短了，又不能在剪切的时候引入太多空位
        trimmed_contigs = cut_align_sequence
        # trimmed_contigs = trimmed_contigs.replace("-", "")
    else:
        trimmed_contigs = ""
    return trimmed_contigs






###############################################
'''
计算一致度和覆盖度     平均时间 0.42290181159973145s
'''
#################################333


def get_identity_number(seq_sequence, ref_sequence):
    alignments = pairwise2.align.globalxx(seq_sequence, ref_sequence)  # 全局比对，相同的残基就给1分，不同和gap不扣分
    matches = alignments[0][2] #比对得分，即相同个数
    return  matches

'''
计算连续空位数量 将连续gap视为一个base
Gap-compressed identity
'''
def get_gap_compressed_number(sequence):
    s=re.findall(r"-{2,}",sequence)  #统计连续空位
    # print(s)
    number=0
    if s==[]:
        number=number
    else:
        for i in s:
            number=number+len(i)-1
    return number

def get_identity_and_coverage_path(seq_path, ref_path):  # query 为seq结果序列,target为目标参考序列（一个fasta文件，可能为多条）
    '''
      以往的策略是选择trimmed中最长的一条作为 seq_path
      经过更改之后，trimmed中有且仅只有一条序列.
      考虑到bootstrap如果未能做出trimmed,只好使用较次的raw_results，所以选择最长的一条，也能最大程度的弥补
      '''
    sequence_list = []
    reference_list = []
    # 取最长的seq （结果序列）
    for rec in SeqIO.parse(seq_path, "fasta"):
        temp_sequence = rec.seq
        sequence_list.append(temp_sequence)
    if len(sequence_list) == 1:
        seq_sequence = sequence_list[0]
    else:
        length_list = [len(i) for i in sequence_list]
        index = length_list.index(max(length_list))
        seq_sequence = sequence_list[index]

    # 取最长的参考序列 （参考序列）
    for rec in SeqIO.parse(ref_path, "fasta"):
        temp_sequence = rec.seq
        reference_list.append(temp_sequence)
    if len(reference_list) == 1:
        ref_sequence = reference_list[0]
    else:
        length_list = [len(i) for i in reference_list]
        index = length_list.index(max(length_list))
        ref_sequence = reference_list[index]

    # 用全局比对计算一致度
    alignments = pairwise2.align.globalms(seq_sequence, ref_sequence, 5, -4, -1, -0.1)  # 全局比对,blast替换矩阵
    matches1 = alignments[0][0]
    matches2 = alignments[0][1]
    identity_number = get_identity_number(matches1, matches2)
    gap_compressed_number = get_gap_compressed_number(matches1)  # 将连续gap视作一个base

    identity = (identity_number / (len(matches1) - gap_compressed_number)) * 100

    # 同样使用全局比对计算覆盖度
    gap = alignments[0][0].count("-")
    length = len(alignments[0][0])
    coverage = ((length - gap) / length) * 100

    identity_and_coverage = [round(identity, 2), round(coverage, 2)]
    return identity_and_coverage


def get_identity_and_coverage_sequence(seq_sequence, ref_sequence):  # query 为seq结果序列,target为目标参考序列（一个fasta文件，可能为多条）
    # 用全局比对计算一致度
    alignments = pairwise2.align.globalms(seq_sequence, ref_sequence,5, -4, -1, -0.1)  # 全局比对,blast替换矩阵
    matches1 = alignments[0][0]
    matches2=  alignments[0][1]
    identity_number=get_identity_number(matches1,matches2)
    gap_compressed_number=get_gap_compressed_number(matches1)  #将连续gap视作一个base

    identity= (identity_number / (len(matches1) - gap_compressed_number)) * 100

    #同样使用全局比对计算覆盖度
    gap = alignments[0][0].count("-")
    length = len(alignments[0][0])
    coverage = ((length - gap) / length) * 100

    identity_and_coverage = [round(identity, 2), round(coverage, 2)]
    return identity_and_coverage





###########################################################################################
############################################################################################
'''
windows linux 切换部分
'''
#####################################################################
####################################################################

'''
将windows的路径映射为wsl可以用的路径  wsl的路径相较于linux的路径，只在于加上了/mnt (等价于windows的盘符)
'''
def windows2wsl_path(path):
    if path==None:
        return 0

    temp_path=Path(path)           #将路径作为pathlib库下面的Path对象
    path_parts=temp_path.parts   #返回元组
    path_parts=list(path_parts)  #列表便于修改  ['D:\\', 'Happy_life_and_work', 'scu', 'python', 'Gene_Miner', 'eeee4 重构版本  大修']
    drive=path_parts[0] #盘符 D:
    drive=drive.split(":")[0].lower() #d
    path_parts[0]=drive          #将windows下的盘符映射为为wsl下的路径
    result_path="/mnt"
    for i in path_parts:
        result_path=os.path.join(result_path,i)
    result_path=Path(result_path).as_posix()  #转为linux 风格

    return result_path

'''
得到系统类型
'''
def get_platform():
    current_os = platform.system().lower()   #linux windows darwin
    return  current_os





'''
按照某种规则，获得某目录下文件的绝对路径

raw_best 获得绝对路径
trimmed 获得剪切前的序列
'''
def get_absolute_path_by_name(file_dir,sign):
    L = []
    for root, dirs, files in os.walk(file_dir):  # 获取所有文件
        for file in files:  # 遍历所有文件名
            if sign in file and sign in "_trimmed":
                file=file.split("_trimmed.fasta")[0]+".fasta"
                L.append(os.path.join(root, file))  # 拼接处绝对路径并放入列表
            elif sign in file:
                L.append(os.path.join(root, file))
    return L







'''
得到绝对路径
'''
def get_absolute(path):
    if path==None:
        return None
    else:
        if os.path.isabs(path):
            return path
        else:
            path = os.path.abspath(path)
            return path


'''
得到绝对路径并作映射
允许批量
'''
def get_absolute_and_map_path(path_list,system):
    path_list=list(path_list)
    if system == "windows":
        for i in range(len(path_list)):
            if path_list[i] == None:
                path_list[i]=None
            else:
                if os.path.isabs(path_list[i]):
                    path_list[i] = windows2wsl_path(path_list[i])
                else:
                    path_list[i] = os.path.abspath(path_list[i])
                    path_list[i] = windows2wsl_path(path_list[i])

    elif system == "linux" or system == "darwin":
        for i in range(len(path_list)):
            if path_list[i] == None:
                path_list[i]=None
            else:
                if os.path.isabs(path_list[i]):
                    path_list[i]=path_list[i]
                else:
                    path_list[i] = os.path.abspath(path_list[i])
                    path_list[i] = path_list[i]
    else:
        pass
    return  path_list


'''
运行cmd命令
'''
def runCommand(cmd,system):
    if system == "linux":
        subprocess.call(cmd, shell=True)

    elif system == "darwin":
        subprocess.call(cmd, shell=True)

    elif system == "windows":
        cmd = "wsl /bin/bash -c {} ".format('"'+cmd+'"')
        subprocess.call(cmd, shell=True)
    else:
        pass

'''
批量得到cmd命令，但不运行
'''
def getCommand(cmd_list,system):
    cmd_list=list(cmd_list)
    if system == "linux":
        for i in range(len(cmd_list)):
            cmd_list[i]=cmd_list[i]
    elif system == "darwin":
        for i in range(len(cmd_list)):
            cmd_list[i] = cmd_list[i]
    elif system == "windows":
        for i in range(len(cmd_list)):
            cmd_list[i] = "wsl /bin/bash -c {} ".format('"'+cmd_list[i]+'"')
    else:
        pass
    return  cmd_list









############################################
############################################

'''
bootstrap_pipeline /iterative_pipelin
'''

###########################################
###########################################



'''
利用muscle多序列比对
'''
def multiple_sequence_alignment(muscle_path,input,output,system):
    path_list = [muscle_path,input,output]
    [muscle_path,input,output] = get_absolute_and_map_path(path_list, system)

    cmd="'{0}' -align '{1}' -output '{2}' ".format(muscle_path,input,output)
    runCommand(cmd,system)




'''
校正序列
（1）长度校正： [0.8:3] and [min:max]
（2）方向校正 
（3）剪切     有多个剪切就有多个option
（4）最佳结果 （一致度最高且最长）
'''
class Get_the_best_result():
    def __init__(self, configuration_information,out_dir_name, seq_path, ref_path, output_raw, output_raw_best, output_options,
                 output_no_trimmed, output_trimmed, gene_name, max_length, min_length, options):
        self.configuration_information = configuration_information
        self.out_dir_name=out_dir_name
        self.seq_path = seq_path
        self.ref_path = ref_path
        self.gene_name = gene_name
        self.max_length = max_length
        self.min_length = min_length
        self.length_seq = calculate_ref_size(self.ref_path)  # 返回平均长度和第一条seq
        self.options = options  # 确定是否要输出options

        '''
        blast相关参数
        '''
        assembled_path_father = os.path.dirname(seq_path)  # minia_out的位置
        reference_basename = "soft_link_reference.fa"
        seq_basename = os.path.basename(seq_path)
        blastn_result = "result.m8"
        index_db = "index"
        self.assembled_path_father = assembled_path_father
        self.reference_basename = reference_basename
        self.seq_basename = seq_basename
        self.blastn_result = blastn_result
        self.index_db = index_db
        self.index_db_path = os.path.join(assembled_path_father, index_db)  # index
        self.ref_path_soft_link_path = os.path.join(assembled_path_father, reference_basename)
        self.blastn_result_path = os.path.join(assembled_path_father, blastn_result)  # result.m8

        self.GM_results = self.configuration_information["GM_results"]
        self.system = self.configuration_information["system"]
        self.whole_log = self.configuration_information["whole_log"]
        self.whole_log_path = os.path.join(self.out_dir_name, self.whole_log)
        self.makeblastdb_path = self.configuration_information["makeblastdb_path"]
        self.blastn_path = self.configuration_information["blastn_path"]
        self.blastn_out = self.configuration_information["blastn_out"]

        '''
        GM_result部分，后续复制到指定位置
        '''

        self.output_raw = output_raw  # 长度+方向校正后的序列   xx.raw.fasta
        self.output_raw_best = output_raw_best  # 长度+方向校正后的序列，但没有通过剪切规则  xx.raw_best.fasta
        self.output_options = output_options  # 通过剪切规则后的序列    xx.options.fasta
        self.output_no_trimmed = output_no_trimmed  # 剪切前                xx.fasta
        self.output_trimmed = output_trimmed  # 剪切后                xx.trimmed.fasta

    '''
    makeblastdb 路径中不能带有空格，即使用引号包裹也不可以 解决思路：（1）用临时文件取代它 （2）切换到该工作目录后，再将对应文件复制到需要使用的地方
    50条序列 3.5091826915740967s
    '''

    def my_makeblastdb_blastn(self):
        system = self.system
        gene_name = self.gene_name
        reference_path = self.ref_path
        seq_path = self.seq_path
        makeblastdb_path = self.makeblastdb_path
        blastn_path = self.blastn_path
        assembled_path_father = self.assembled_path_father

        reference_basename = self.reference_basename
        seq_basename = self.seq_basename
        blastn_result = self.blastn_result  # result.m8
        index_db = self.index_db  # index

        # print(reference_basename,seq_basename)

        path_list = [assembled_path_father, reference_path, seq_path, makeblastdb_path, blastn_path]
        [assembled_path_father, reference_path, seq_path, makeblastdb_path, blastn_path] = get_absolute_and_map_path(
            path_list, system)

        # 先准备好ref 和query序列
        cmd = "cd '{0}' &&ln -s '{1}' {2}".format(assembled_path_father, reference_path, reference_basename)
        runCommand(cmd, system)

        # t1=time.time()
        # 建库
        cmd = "cd '{0}' && '{1}' -in '{2}' -out index/index  -dbtype nucl >/dev/null  2>&1".format(
            assembled_path_father, makeblastdb_path, reference_basename)
        runCommand(cmd, system)
        # 比对 (1,-1,3,2)打分矩阵适用于75%保守度
        reward = 1
        penalty = -1
        gapopen = 3
        gapextend = 2
        number = 1
        threads = 1
        blast_word_size=25
        cmd = "cd '{0}' && '{1}' -query '{2}' -out {3} -db {4}/index -outfmt 6  -reward {5} -penalty {6} -gapopen {7} -gapextend {8} -num_alignments {9} -num_threads {10} -word_size {11}   >/dev/null  2>&1".format(
            assembled_path_father, blastn_path, seq_basename, blastn_result, index_db, reward, penalty, gapopen,
            gapextend, number, threads,blast_word_size)
        runCommand(cmd, system)
        # t2=time.time()
        # print(t2-t1)

    def parse_blastn_m8(self):
        blastn_result_path = self.blastn_result_path
        m8_information = []
        if is_exist(blastn_result_path):
            with open(blastn_result_path, "r") as f:
                for line in f.readlines():
                    information = {}
                    line = line.strip()
                    line = line.split("\t")
                    information["query_id"] = line[0]
                    information["subject_id"] = line[1]
                    information["identity_local"] = line[2]
                    information["alignment_length"] = line[3]
                    information["mismatches"] = line[4]
                    information["gap"] = line[5]
                    information["query_start"] = line[6]
                    information["query_end"] = line[7]
                    information["subject_start"] = line[8]
                    information["subject_end"] = line[9]
                    information["e_value"] = line[10]
                    information["score"] = line[11]

                    if int(information["query_end"]) >= int(information["query_start"]):
                        flag1 = 1
                    else:
                        flag1 = -1

                    if int(information["subject_end"]) >= int(information["subject_start"]):
                        flag2 = 1
                    else:
                        flag2 = -1

                    strand = flag1 * flag2
                    information["strand"] = strand  # 利用strand判断是否需要反向互补

                    m8_information.append(information)
        return m8_information

    '''
    将reference和seq的相关信息，（sequence,name,description）添加到m8信息中
    '''

    def add_query_reference_information_2_m8(self, m8_information):
        reference_path = self.ref_path  # 这里不要使用软连接，用真实的路径
        seq_path = self.seq_path

        if m8_information == []:
            return []
        reference_id_list = []
        query_id_list = []
        for i in m8_information:
            if i["subject_id"] not in reference_id_list:  # 去重了
                reference_id_list.append(i["subject_id"])
            if i["query_id"] not in query_id_list:
                query_id_list.append(i["query_id"])

        subject_information = []
        query_information = []

        # 将reference信息添加到m8信息中
        for rec in SeqIO.parse(reference_path, "fasta"):
            name = rec.name
            sequence = rec.seq
            description = rec.description
            if name in reference_id_list:
                length = len(str(rec.seq))
                temp = {}
                temp["subject_id"] = name
                temp["subject_length"] = length
                temp["sequence"] = sequence
                temp["description"] = description
                subject_information.append(temp)
        for i in m8_information:
            for j in subject_information:
                if i["subject_id"] == j["subject_id"]:
                    i["subject_length"] = j["subject_length"]
                    i["subject_sequence"] = j["sequence"]
                    i["subject_description"] = j["description"]

        # 将assembled_contigs信息添加到m8信息中
        for rec in SeqIO.parse(seq_path, "fasta"):
            name = rec.name
            sequence = rec.seq
            description = rec.description
            if name in query_id_list:
                length = len(str(rec.seq))
                temp = {}
                temp["query_id"] = name
                temp["query_length"] = length
                temp["sequence"] = sequence
                temp["description"] = description
                query_information.append(temp)
        for i in m8_information:
            for j in query_information:
                if i["query_id"] == j["query_id"]:
                    i["query_length"] = j["query_length"]
                    i["query_sequence"] = j["sequence"]
                    i["query_description"] = j["description"]

        return m8_information

    '''
    根据blast的结果，进行粗过滤。identity,coverage.原始序列的长度

    这里有个隐藏的bug，如果序列中间有专门挖的空洞，将一条序列拆成质量很高的两部分，由于一致度不够的原因，就没有保存
    '''

    def filter_blastn_m8(self, m8_information):
        if m8_information == []:
            return []
        max_length = self.max_length
        min_length = self.min_length

        subject_coverage_threshold = 0.4  # 除去过短的就ok了
        identity_local_threshold = 50
        # 根据局部比对的一致度和覆盖度进行过滤
        m8_information_filtered = []
        for i in m8_information:
            start = int(i["subject_start"])
            end = int(i["subject_end"])
            subject_length = int(i["subject_length"])
            subject_coverage = abs(end - start + 1) / subject_length
            identity_local = float(i["identity_local"])
            query_length = int(i["query_length"])

            # print(identity_local,subject_coverage)

            if subject_coverage > subject_coverage_threshold and identity_local > identity_local_threshold and query_length <= max_length and query_length >= min_length:  # 一致度覆盖度有一定要求，最长最短看的是全长，而不是部分序列
                m8_information_filtered.append(i)

        return m8_information_filtered

    def m8_2_fasta(self, m8_information):
        out_path = self.output_raw
        if m8_information == []:
            return 0

        my_records = []
        for i in m8_information:
            if i["strand"] == 1:
                sequence = i["query_sequence"]
            else:
                sequence = i["query_sequence"].reverse_complement()
            description = i["query_description"]
            id = i["query_id"]
            record = SeqRecord(seq=sequence, id=id, description=description)
            my_records.append(record)
        if my_records != []:
            SeqIO.write(my_records, out_path, "fasta")

    def cut_align_m8(self, m8_information):

        GM_results_path_raw = self.output_raw
        GM_results_path_raw_best = self.output_raw_best
        GM_results_path_options = self.output_options
        GM_results_path_no_trimmed = self.output_no_trimmed
        GM_results_path_trimmed = self.output_trimmed
        options = self.options
        system = self.system
        whole_log_path = self.whole_log_path
        gene_name = self.gene_name

        if m8_information == []:
            return 0

        my_records_trimmed = []
        my_records_no_trimmed = []
        my_records_raw_best = []
        for i in m8_information:
            if i["strand"] == 1:
                sequence_query = i["query_sequence"]
                sequence_ref = i["subject_sequence"]
            else:
                sequence_query = i["query_sequence"].reverse_complement()
                sequence_ref = i["subject_sequence"]

            # print(sequence_query)
            # print(sequence_ref)
            trimmed_contigs = cut_align(sequence_query, sequence_ref)

            if trimmed_contigs != "":  # 能够剪切
                description_trimmed = "trimmed_length_" + str(len(trimmed_contigs))
                description_no_trimmed = "lenght_" + str(i["query_length"])
                id = i["query_id"]
                record_trimmed = SeqRecord(seq=Seq(trimmed_contigs), id=id, description=description_trimmed)
                record_no_trimmed = SeqRecord(seq=sequence_query, id=id, description=description_no_trimmed)
                my_records_trimmed.append(record_trimmed)
                my_records_no_trimmed.append(record_no_trimmed)

                identity_coverage = get_identity_and_coverage_sequence(sequence_query, sequence_ref)
                i["identity_global"] = identity_coverage[0]
                i["coverage_global"] = identity_coverage[1]

            else:  # 不能够剪切
                description_raw_best = "lenght_" + str(i["query_length"])
                id = i["query_id"]
                record_raw_best = SeqRecord(seq=sequence_query, id=id, description=description_raw_best)
                my_records_raw_best.append(record_raw_best)
                identity_coverage = get_identity_and_coverage_sequence(sequence_query, sequence_ref)
                i["identity_global"] = identity_coverage[0]
                i["coverage_global"] = identity_coverage[1]

        # print(m8_information)  #此时拿到最完全版本的m8_information

        path_list = [whole_log_path]
        [whole_log_path] = get_absolute_and_map_path(path_list, system)

        # 没能通过剪切规则，但需要在里面挑选一条最优越的。根据覆盖度
        if my_records_trimmed == []:
            SeqIO.write(my_records_raw_best, GM_results_path_raw_best, "fasta")  # xx.raw_best.fasta 后续选择最佳的一条
            message = " '{0}':Satisfactory results have been obtained, " \
                      "but after the pairwise sequence alignment, " \
                      "the program can not cut and align the sequences," \
                      "requiring users to cut or align according" \
                      " to their needs".format(gene_name)
            cmd = "echo {0} >>'{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)

            self.get_geneminer_best_result(m8_information, GM_results_path_raw_best)


        # 通过剪切规则，但需要在里面挑选一条最优越的。根据覆盖度
        elif my_records_trimmed != []:
            SeqIO.write(my_records_no_trimmed, GM_results_path_no_trimmed, "fasta")  # xx.fasta  剪切前  (可能一条或者多条)
            SeqIO.write(my_records_trimmed, GM_results_path_trimmed, "fasta")  # xx.trimmed.fasta （可能一条或者多条）

            if len(my_records_trimmed) > 1 and options == "yes":
                SeqIO.write(my_records_no_trimmed, GM_results_path_options, "fasta")

            message = "'{0}': check succeeded, aligned_and_cut sequence has also been generated".format(gene_name)
            cmd = "echo {0} >>'{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)

            self.get_geneminer_best_result(m8_information, GM_results_path_no_trimmed)
            self.get_geneminer_best_result(m8_information, GM_results_path_trimmed)

        else:
            pass

        return m8_information

    '''
    根据上一步向m8_information（已经经过过滤了）中添加的identity_global，coverage_global 计算最佳比对结果

    input 需要筛选最优解的序列
    '''

    def get_geneminer_best_result(self, m8_information, input_fasta):
        if m8_information == []:
            return 0

        number = len(list(SeqIO.parse(input_fasta, "fasta")))
        if number == 1:
            return 0

        index = 0
        for i in range(len(m8_information)):
            if i == 0:
                index = 0
            else:  # 优先一致度，然后覆盖度
                if m8_information[i]["identity_global"] > m8_information[i]["identity_global"]:
                    index = index
                elif m8_information[i]["identity_global"] > m8_information[i]["identity_global"] and \
                        m8_information[i]["coverage_global"] > m8_information[i]["coverage_global"]:
                    index = i

        flag = 0
        best_record = []
        for rec in SeqIO.parse(input_fasta, "fasta"):
            if flag == index:
                id = rec.name
                description = rec.description
                seq = rec.seq

                record = SeqRecord(seq=seq, id=id, description=description)
                best_record.append(record)
                break
            else:
                flag = flag + 1
                continue

        if best_record != []:
            SeqIO.write(best_record, input_fasta, "fasta")

    '''
    删掉临时文件 index, result.m8, soft_link_reference.fa
    '''
    def deal_temp_file(self):
        system = self.system
        assembled_path_father = self.assembled_path_father
        index_db = self.index_db  # index
        blastn_result = self.blastn_result  # result.m8
        reference_basename = self.reference_basename  # soft_link_matk.fasta

        index_db_path = os.path.join(assembled_path_father, index_db)
        blastn_result_path = os.path.join(assembled_path_father, blastn_result)
        reference_basename_path = os.path.join(assembled_path_father, reference_basename)

        path_list = [assembled_path_father, index_db_path, blastn_result_path, reference_basename_path]
        [assembled_path_father, index_db_path, blastn_result_path, reference_basename_path] = get_absolute_and_map_path(
            path_list, system)
        cmd = "cd '{0}' && rm -rf '{1}' '{2}' '{3}' >/dev/null  2>&1".format(assembled_path_father, index_db_path,
                                                                             blastn_result_path,
                                                                             reference_basename_path)
        runCommand(cmd, system)
















class Simplified_pipeline_bootstrap_iterative():
    def __init__(self, configuration_information, type, gm_result,filter_reference,cut_align_reference, out_dir_name, sub_out_dir_name, thread_number, kmer,wordsize,max_length,min_length,options):

        self.configuration_information = configuration_information
        self.type = type
        self.gm_result=gm_result                        #主要作用用于确定基因名
        self.filter_reference=filter_reference          #过滤参考序列，绝对路径
        self.cut_align_reference=cut_align_reference    #剪切参考序列，绝对路径
        self.out_dir_name = out_dir_name  # 总文件目录
        self.sub_out_dir_name = sub_out_dir_name  # 总文件目录下的一个目录
        self.thread_number = thread_number
        self.kmer = kmer
        self.wordsize=wordsize
        self.max_length = max_length
        self.min_length = min_length
        self.options = options

        self.my_software_name = self.configuration_information["my_software_name"]
        self.filter_software = self.configuration_information["filter_path"]  # filter_reads.pl / filter  绝对路径
        self.assemble_software = self.configuration_information["assemble_path"]  # minia                   绝对路径
        self.reference_database = self.configuration_information["reference_database"]  # 参考基因数据库
        self.filtered_out = self.configuration_information["filtered_out"]
        self.assembled_out = self.configuration_information["assembled_out"]
        self.GM_results = self.configuration_information["GM_results"]
        self.system = self.configuration_information["system"]
        self.results_information_excel = self.configuration_information[
            "results_information_excel"]  # "results_information.xlsx"
        self.whole_log = self.configuration_information["whole_log"]

        gene_name = os.path.basename(self.gm_result)
        if "_trimmed.fasta" in gene_name:
            gene_name = gene_name.split("_trimmed.fasta")[0]  # 方便后续扩展
            file_name = gene_name + ".fasta"
        elif "_raw_best.fasta" in gene_name:
            gene_name = gene_name.split("_raw_best.fasta")[0]
            file_name = gene_name + ".fasta"
        else:
            gene_name = gene_name.split(".fasta")[0]
            file_name = gene_name +".fasta"
        self.gene_name = gene_name
        self.file_name = file_name

    '''
    1产生两个过滤文件Filtered_reads__R1.fastq ， Filtered_reads__R2.fastq ,并且合并双端数据为filtered.fq
    '''

    def filter_reads(self):
        # file  matk.fasta

        out_dir_name = self.out_dir_name  # 定位所有基因名，所有传参都是根据基因名
        sub_out_dir_name = self.sub_out_dir_name  # bootstrap / iterative
        filtered_out = self.filtered_out  # 二级目录
        reference_database = self.reference_database
        type = str(self.type)  # tfa,tgb, cp ,target 中的某一种
        filter_software = self.filter_software  # lib中的filter  #绝对路径
        system = self.system
        whole_log=self.whole_log
        gene_name = self.gene_name  # ycf4
        reference_path = self.filter_reference
        wordsize=self.wordsize

        filtered_out_path = os.path.join(out_dir_name, sub_out_dir_name,filtered_out)
        dir_make(filtered_out_path)

        data1_path = os.path.join(out_dir_name, "data1.fq")
        data2_path = os.path.join(out_dir_name, "data2.fq")
        whole_log_path=os.path.join(out_dir_name,whole_log)


        # print(filter_software)
        if is_exist(data1_path) and is_exist(data2_path):
            path_list = [out_dir_name, filter_software, filtered_out_path, data1_path, data2_path, reference_path,whole_log_path]
            [out_dir_name, filter_software, filtered_out_path, data1_path, data2_path,reference_path,whole_log_path] = get_absolute_and_map_path(path_list, system)  # 如果为windows环境，会批量映射路径

            # print([out_dir_name, filter_software, filtered_out_path, data1_path, data2_path,reference_path,whole_log_path])
            cmd = "cd '{0}' && '{1}' -1 '{2}' -2 '{3}' -r '{4}' -k {5} >/dev/null  2>&1".format(filtered_out_path,filter_software, data1_path,data2_path, reference_path,wordsize)
            # print(cmd)
            runCommand(cmd, system)
            cmd = "cd '{0}' && cat {1} {2} > {3}".format(filtered_out_path, "Filtered_reads__R1.fastq",
                                                         "Filtered_reads__R2.fastq",
                                                         "filtered.fq")
            runCommand(cmd, system)
            # 立刻将完成的基因打印
            message = "{0}_gene: '{1}' has been filtered. DONE".format(type, gene_name)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)
        else:
            pass

    '''
    2 拼接目标基因  产生三个文件：assembled_out.contigs.fa  assembled_out.h5  assembled_out.unitigs.fa
    '''
    def assemble_reads(self):
        # file=matk.fasta
        type = str(self.type)
        kmer = self.kmer
        assemble_software = self.assemble_software  # minia   #绝对路径
        out_dir_name = self.out_dir_name
        sub_out_dir_name = self.sub_out_dir_name
        filterd_out = self.filtered_out
        assembled_out = self.assembled_out
        system = self.system
        whole_log = self.whole_log
        gene_name = self.gene_name  #matk

        input = os.path.join(out_dir_name, sub_out_dir_name,filterd_out, "filtered.fq")
        whole_log_path = os.path.join(out_dir_name, whole_log)

        if is_exist(input) == 0:
            path_list = [whole_log_path]
            [whole_log_path] = get_absolute_and_map_path(path_list, system)  # 如果为windows环境，会批量映射路径
            message = "'{0}' :Assembling reads has failed.  It could be that the quality of the raw data is not good enough, the amount of input data is too small or the species  are not closely related".format(
                gene_name)
            cmd = "echo {0}>>'{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)
            return 0

        dir = os.path.join(out_dir_name, sub_out_dir_name,assembled_out)  # assemble reads输出的文件夹
        dir_make(dir)

        path_list = [dir, assemble_software, input, whole_log_path]
        [dir, assemble_software, input, whole_log_path] = get_absolute_and_map_path(path_list,
                                                                                    system)  # 如果为windows环境，会批量映射路径
        try:
            cmd = "cd '{0}' && '{1}' -in '{2}' -out '{3}' -kmer-size {4} -nb-cores 1 >/dev/null 2>&1 ".format(dir,
                                                                                                  assemble_software,
                                                                                                  input,
                                                                                                  "assembled_out", kmer)
            runCommand(cmd, system)

            cmd = "cd '{0}' && rm -rf assembled_out.h5".format(dir)  # assembled_out.h5文件很大，而且实际上没什么用，可以删除
            runCommand(cmd, system)

            message = "{0}_gene: '{1}' has been assembled successfully. DONE".format(type, gene_name)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)
        except:
            message = "{0}_gene: '{1}' assembly failed".format(type, gene_name)
            cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
            runCommand(cmd, system)

        # assembled_out.unitigs.fa.glue1....  只保留三个文件 assembled_out.contigs.fa assembled_out.unitigs.fa assembled_out.h5

    '''
    3  检测assembled_out的结果，校正长度与方向，并且对齐剪切 
    '''

    def check_contigs(self):
        # file=matk.fasta
        configuration_information=self.configuration_information
        out_dir_name = self.out_dir_name
        sub_out_dir_name = self.sub_out_dir_name
        GM_results = self.GM_results
        reference_database = self.reference_database
        assembled_out = self.assembled_out
        system = self.system
        whole_log = self.whole_log
        gene_name = self.gene_name # matk.fasta ---- matk
        max_length = self.max_length
        min_length = self.min_length
        options = self.options
        ref_path = self.cut_align_reference

        assembled_path = os.path.join(out_dir_name, sub_out_dir_name, assembled_out,
                                      "assembled_out.contigs.fa")

        whole_log_path = os.path.join(out_dir_name, whole_log)

        path_list = [whole_log_path]
        [whole_log_path] = get_absolute_and_map_path(path_list, system)  # 如果为windows环境，会批量映射路径
        if is_exist(assembled_path) == 0:
            message = "'{0}': assembling reads has failed,because the data quality was poor or the amount of input data was too small".format(
                gene_name)
            cmd = "echo {0} >> '{1}'".format(message, whole_log_path)
            runCommand(cmd, system)
            return 0
        GM_results_path = os.path.join(out_dir_name, sub_out_dir_name, GM_results)  # bootstrap或者iterative下 GM_results的路径
        dir_make(GM_results_path)

        GM_results_path_raw = os.path.join(GM_results_path, gene_name + "_raw.fasta")
        GM_results_path_raw_best=os.path.join(GM_results_path,gene_name+"_raw_best.fasta")
        GM_results_path_options=os.path.join(GM_results_path,gene_name+"_options.fasta")
        GM_results_path_no_trimmed=os.path.join(GM_results_path,gene_name+".fasta")
        GM_results_path_trimmed = os.path.join(GM_results_path,gene_name + "_trimmed.fasta")




        # 限定长度
        # 确定方向
        # 切齐 双序列比对(局部最优解)
        # 得到最佳剪切结果
        # 限定长度
        # 确定方向
        # 切齐 双序列比对(局部最优解)
        # 得到最佳剪切结果
        my_verify = Get_the_best_result(configuration_information,out_dir_name, assembled_path, ref_path, GM_results_path_raw,
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










##################################
'''
GUI换行输出
'''

def ML_str_re(in_str=""):
    in_str = in_str.replace("\r", "\r\n")
    global del_count
    in_str_list = in_str.split("\n")
    out_str = ""
    for i in in_str_list:
        if i.endswith('\r') == False:
            out_str += i + "\n"
    return out_str

#################################3































