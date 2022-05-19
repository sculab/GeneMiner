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
from concurrent import futures
from Bio import pairwise2
import random
import re
from  collections import  defaultdict
from basic import get_basename,get_fasta_file,run_command,get_file_list,is_exist,mylog,cutting_line,get_files


################################################
#################################################

def get_bootstrap_hashdict(reference, merSize):
    kmer_dict = defaultdict(list)
    infile = open(reference, 'r', encoding='utf-8', errors='ignore')
    name=""
    seq=""
    my_list=[]
    while True:
        line = infile.readline()
        line = line.strip()
        if (line.startswith('>') or not line) and name:  # 保证最后一条序列能能在保存后退出
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
        if line.startswith('>'):
            name = line[1:]
            seq = ''
        else:
            seq += line
        if not line:
            break

    for i in my_list:
        gene_name=list(i.keys())[0]
        refseq=list(i.values())[0]
        for j in range(0,len(refseq)-merSize+1):
            temp_list, kmer = [], refseq[j:j + merSize]
            # print(kmer,j)
            if kmer in kmer_dict:
                if gene_name not in kmer_dict[kmer]:
                    kmer_dict[kmer].append(gene_name)
            else:
                kmer_dict[kmer]=[gene_name]
    return  kmer_dict

#根据名字获得对应序列
def get_seq_from_name(file,name_list):
    infile = open(file, 'r', encoding='utf-8', errors='ignore')
    name = ""
    seq = ""
    my_list=[]

    while True:
        line = infile.readline()
        line = line.strip()
        if (line.startswith('>') or not line) and (name in name_list):  # 保证最后一条序列能能在保存后退出,保证第一条是有效的
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
        if line.startswith('>'):
            name = line[1:]
            seq = ''
        else:
            seq += line
        if not line:
            break
    return  my_list  #[{gene:seq}]


'''
从fasta序列中获取seq，可指定条目
'''
def get_seq(fasta_file,max_seq_number=100,seq_count_limit=False):
    infile = open(fasta_file, 'r', encoding='utf-8', errors='ignore')
    seq, name = "", ""
    my_list = []
    seq_number = 0
    while True:
        line = infile.readline()
        line = line.strip()
        line=line.replace("N","")  #特定对于scaffold
        if (line.startswith('>') or not line) and name:  # 保证最后一条序列能能在保存后退出
            temp = {}
            temp = {name: seq}
            my_list.append(temp)
            seq_number += 1
        if line.startswith('>'):
            name = line[1:]
            seq = ''
        else:
            seq += line
        if not line:
            break
        if seq_count_limit and max_seq_number:
            if seq_number >= max_seq_number:
                break
    infile.close()
    return my_list


def get_gap_information(alignment_seq):
    gap_number_list = re.findall(r"-", alignment_seq)  # 统计空位数量
    continuous_gap_list = re.findall(r"-{2,}", alignment_seq)  # 统计连续空位
    gap_number = len(gap_number_list)
    gap_extend = 0  # 将连续gap记为一个gap,被压缩gap数量=gap长度 - 1
    if continuous_gap_list == []:
        gap_extend = gap_extend
    else:
        for i in continuous_gap_list:
            gap_extend += len(i) - 1

    gap_open = gap_number - gap_extend

    continuous_gap_inter_list = re.findall(r"[ACGT]-{2,}[ACGT]", alignment_seq)  # 记录中间 连续gap数量
    single_gap_inter_list = re.findall(r"[ACGT]-[ACGT]", alignment_seq)  # 记录中间 单一gap数量
    continuous_gap_inter_list_new = []
    single_gap_inter_list_new = []
    if continuous_gap_inter_list != []:
        for i in continuous_gap_inter_list:
            temp = re.sub(r"[ACGT]", "", i)
            continuous_gap_inter_list_new.append(temp)
    if single_gap_inter_list != []:
        for i in single_gap_inter_list:
            temp = re.sub(r"[ACGT]", "", i)
            single_gap_inter_list_new.append(temp)
    gap_extend_inter = 0
    if continuous_gap_inter_list_new == []:
        gap_extend_inter = 0
    else:
        for i in continuous_gap_inter_list_new:
            gap_extend_inter += len(i) - 1

    continuous_gap_inter = "".join(continuous_gap_inter_list_new)
    single_gap_inter = "".join((single_gap_inter_list_new))
    continuous_gap_inter_number = len(continuous_gap_inter)
    single_gap_inter_number = len(single_gap_inter)

    start_end_gap_number = gap_number - continuous_gap_inter_number - single_gap_inter_number  # 首尾gap = 总gap - 中间gap(单一+连续)

    # print("con:{}".format(continuous_gap_inter_list))
    # print("con_new:{}".format(continuous_gap_inter_list_new))
    # print("sin:{}".format(single_gap_inter_list))
    # print("sin_new:{}".format(single_gap_inter_list_new))

    return [gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number]


def get_identity(seq1, seq2):  # seq1 seq2(ref)
    #打分矩阵 Global alignment with free end gaps   5/-4/-12/-3 (65%)   1/-1/-3/-2 (75%)
    match_socre = 5
    mismatich_score = -4.0
    gap_open_socre = -12
    gap_extend_score = -3
    alignments = pairwise2.align.globalms(seq1, seq2, match_socre, mismatich_score, gap_open_socre,
                                          gap_extend_score)  # 全局比对，相同的残基就给1分，不同和gap不扣分


    alignment_seq1 = alignments[0][0]
    alignment_seq2 = alignments[0][1]
    score = alignments[0][2]  # 比对得分
    length = len(alignment_seq1)
    gap_information1 = get_gap_information(alignment_seq1)
    gap_information2 = get_gap_information(alignment_seq2)
    gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number = gap_information2 if gap_information1[0] <= \
                                                                                 gap_information2[
                                                                                     0] else gap_information1
    # print( gap_number, gap_open, gap_extend, gap_extend_inter, start_end_gap_number)
    # print("alignment1:{}".format(alignment_seq1))
    # print("alignment2:{}".format(alignment_seq2))

    # x = Symbol("x")  # match_number
    # y = Symbol("y")  # mismatch_number
    # expression = [x + y + gap_number - length,
    #               x * match_socre + y * mismatich_score + gap_open * gap_open_socre + gap_extend * gap_extend_score - score]
    # result = solve(expression, [x, y])
    # match_number = int(result[x])
    match_number=  ( (score -gap_open * gap_open_socre - gap_extend * gap_extend_score) - (mismatich_score * (length-gap_number))  )/(match_socre-mismatich_score)
    mismatch_number = ( match_socre * (length-gap_number) + gap_open * gap_open_socre + gap_extend * gap_extend_score - score ) / (match_socre-mismatich_score)
    match_number= int(match_number)
    mismatch_number=int(mismatch_number)
    # print(match_number,mismatch_number)
    identity = round((match_number / (length - gap_extend_inter-start_end_gap_number)), 4)
    return identity





def random_replacement(file,var_rate,output,bootstrap_number):
    file_name=get_basename(file)
    # 随机替换的对应表
    base_dict = {"A": ["T", "C", "G"],
                 "T": ["A", "C", "G"],
                 "C": ["A", "T", "G"],
                 "G": ["A", "T", "C"]}

    seq=get_seq(file)
    sequence=list( seq[0].values() )[0]
    gene_name=list( seq[0].keys())  [0]


    sequence_length=len(sequence)
    var_number=int (sequence_length*var_rate)



    #将自展基因写在用一个文件下
    for i in range(bootstrap_number):
        #随机变异，如果变异数为0,相当于直接输出该序列
        var_site_list = random.sample(range(0, sequence_length), var_number) #每一次的随机种子都不一样
        seq_to_list=list(sequence)
        for j in var_site_list:
            var=random.randint(0,2)
            base=seq_to_list[j]
            seq_to_list[j]=base_dict[base][var]
        seq_var="".join(seq_to_list)

        bootstrap_out=file_name+".fasta"
        bootstrap_out_path=os.path.join(output,bootstrap_out)
        with open(bootstrap_out_path,"a") as f:
            f.write('>'+gene_name+"_bootstrap_"+str(i)+"\n")
            f.write(seq_var+"\n")





'''
自展检测  每次检测一个基因
'''
##################################################
#################################################

class BootstrapPipeLine():
    def __init__(self,configuration_information):
        self.configuration_information=configuration_information
        self.data1=configuration_information["data1"]
        self.data2=configuration_information["data2"]
        self.single=configuration_information["single"]
        self.out_dir = configuration_information["out_dir"]

        self.k1=configuration_information["k1"]
        self.k2=configuration_information["k2"]
        self.step_length=configuration_information["step_length"]
        self.data_size = configuration_information["data_size"]
        self.limit_count = configuration_information["limit_count"]
        self.limit_min_length = configuration_information["limit_min_length"]
        self.limit_max_length = configuration_information["limit_max_length"]
        self.change_seed = configuration_information["change_seed"]
        self.scaffold_or_not = configuration_information["scaffold_or_not"]
        self.max_length=configuration_information["max_length"]
        self.min_length=configuration_information["min_length"]
        self.thread_number=configuration_information["thread_number"]
        self.bootstrap_number=configuration_information["bootstrap_number"]


        self.reference_database=configuration_information["reference_database"]
        self.filtered_out=configuration_information["filtered_out"]
        self.assembled_out=configuration_information["assembled_out"]
        self.GM_results=configuration_information["GM_results"]
        self.boostrap_out=configuration_information["bootstrap_out"]
        self.results_log=configuration_information["results_log"]
        self.my_software_name=configuration_information["my_software_name"]


        '''
        路径
        '''
        #软件路径
        self.filter_path = configuration_information["filter_path"]
        self.assemble_path = configuration_information["assemble_path"]
        #一级路径
        self.reference_database_path= os.path.join(self.out_dir,self.reference_database)
        self.filtered_out_path=os.path.join(self.out_dir,self.filtered_out)
        self.assembled_out_path=os.path.join(self.out_dir,self.assembled_out)
        self.GM_results_path=os.path.join(self.out_dir,self.GM_results)
        self.boostrap_out_path=os.path.join(self.out_dir,self.boostrap_out)
        #二级路径
        self.boostrap_out_reference_database_path=os.path.join(self.out_dir,self.boostrap_out,self.reference_database) #out_dir\bootstrap_out\reference_database
        self.boostrap_out_filtered_out_path=os.path.join(self.out_dir,self.boostrap_out,self.filtered_out)
        self.boostrap_out_assembled_out_path = os.path.join(self.out_dir, self.boostrap_out, self.assembled_out)
        self.boostrap_out_GM_results_path = os.path.join(self.out_dir, self.boostrap_out, self.GM_results)




    def get_mutated_sequence(self,gm_result_path, ref_path):
        kmer=self.k2   #拼接的kmer
        output=self.boostrap_out_reference_database_path
        bootstrap_number=self.bootstrap_number

        # hash 建库
        ref_kmer_dict = get_bootstrap_hashdict(ref_path, kmer)
        gm_kmer_dict = get_bootstrap_hashdict(gm_result_path, kmer)
        kmer_count = defaultdict(int)
        for i in gm_kmer_dict:
            if i in ref_kmer_dict:
                for z in ref_kmer_dict[i]:
                    kmer_count[z] += 1

        # 获得中位数kmercount 参考序列作为平均变异度的计算
        sorted_list = sorted(kmer_count.items(), key=lambda x: x[1], reverse=True)
        name_list = []
        length = len(sorted_list)  # 144
        name_list.append(sorted_list[int(length/2)][0])  # 取kmercount中位数

        if not kmer_count:
            limit_kmer = 17
            kmer = limit_kmer
            ref_kmer_dict_limit = get_bootstrap_hashdict(ref_path, kmer)
            gm_kmer_dict_limit = get_bootstrap_hashdict(gm_result_path, kmer)
            for i in gm_kmer_dict_limit:
                if i in ref_kmer_dict_limit:
                    for z in ref_kmer_dict_limit[i]:
                        kmer_count[z] += 1

        if not kmer_count:
            return 0


        # 获得一致度，变异度
        ref_seq = get_seq_from_name(ref_path, name_list)
        ref_seq = list(ref_seq[0].values())[0]
        gm_seq = get_seq(gm_result_path)
        gm_seq = list(gm_seq[0].values())[0]
        identity = get_identity(gm_seq, ref_seq)
        var_rate = round((1 - identity), 4)
        # 根据变异度 随机变异
        random_replacement(gm_result_path, var_rate, output, bootstrap_number)

    # 如果GM_results不存在，就没有后续了
    def get_mutated_sequence_parallel(self):
        GM_results_path=self.GM_results_path
        reference_database_path = self.reference_database_path
        thread_number=self.thread_number
        bootstrap_out_path=self.boostrap_out_path
        bootstrap_out_reference_database_path=self.boostrap_out_reference_database_path
        files=get_fasta_file(GM_results_path)
        if files==[]:
            return 0
        if not os.path.isdir(bootstrap_out_path):
            os.mkdir( bootstrap_out_path)
        if not os.path.isdir(bootstrap_out_reference_database_path):
            os.mkdir(bootstrap_out_reference_database_path)
        task_pool = []
        results=[]
        executor = futures.ProcessPoolExecutor(max_workers=thread_number)  #24s
        # executor = futures.ThreadPoolExecutor(max_workers=thread_number) #40s
        for i in files:
            name=get_basename(i)
            ref_path=os.path.join(reference_database_path,name+".fasta")
            task_pool.append(executor.submit(self.get_mutated_sequence,i,ref_path))
        total = len(task_pool)
        number = 1
        for i in task_pool:
            if number < total:
                sys.stdout.write('\r' + "{0:}:{1:>4}/{2}".format("Preparing bootstrap data",number, total))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r' + "{0:}:{1:>4}/{2}".format("Preparing bootstrap data", number, total) + "\n")
                sys.stdout.flush()
            results.append(i.result())
            number = number + 1

        executor.shutdown()
        return  1


    def bootstrap_filter(self):
        data1 = self.data1
        data2 = self.data2
        single = self.single
        boostrap_out_reference_database_path = self.boostrap_out_reference_database_path
        boostrap_out_filtered_out_path=self.boostrap_out_filtered_out_path
        k1 = self.k1
        step_length = self.step_length
        thread_number = self.thread_number
        data_size = self.data_size
        filter_path = self.filter_path

        files = get_files(boostrap_out_reference_database_path)
        if files == []:
            return 0

        if data1 and data2:
            cmd = 'python "{0}" -1 "{1}" -2 "{2}" -o "{3}"  -r "{4}" -k1 {5} -step_length {6} -t {7} -d {8} '.format(
                filter_path,data1, data2,boostrap_out_filtered_out_path,boostrap_out_reference_database_path, k1, step_length, thread_number, data_size)

            # print(cmd)
            run_command(cmd)
        elif single:
            cmd = 'python "{0}" -s "{1}" -o "{2}"  -r "{3}" -k1 {4} -step_length {5} -t {6} -d {7}'.format(filter_path,
                            single,boostrap_out_filtered_out_path,boostrap_out_reference_database_path,k1,step_length,thread_number,
                                                                                                        data_size)
            run_command(cmd)
        else:
            pass


    # def bootstrap_filter(self):
    #     data1 = self.data1
    #     data2 = self.data2
    #     single = self.single
    #     boostrap_out_reference_database_path = self.boostrap_out_reference_database_path
    #     boostrap_out_filtered_out_path=self.boostrap_out_filtered_out_path
    #     k1 = self.k1
    #     step_length = self.step_length
    #     thread_number = self.thread_number
    #     data_size = self.data_size
    #     filter_path = self.filter_path
    #
    #     filter_configuration_information={}
    #     filter_configuration_information["data1"]=data1
    #     filter_configuration_information["data2"] = data2
    #     filter_configuration_information["single"] = single
    #     filter_configuration_information["thread_number"] = thread_number
    #     filter_configuration_information["k1"] = k1
    #     filter_configuration_information["out_dir"]=boostrap_out_filtered_out_path
    #     filter_configuration_information["step_length"] = step_length
    #     filter_configuration_information["reference"] = boostrap_out_reference_database_path
    #     filter_configuration_information["data_size"] = data_size
    #     my_filter_main_internal(filter_configuration_information)




    def bootstrap_assemble(self):
        boostrap_out_reference_database_path = self.boostrap_out_reference_database_path
        boostrap_out_filtered_out_path = self.boostrap_out_filtered_out_path
        boostrap_out_assembled_out_path = self.boostrap_out_assembled_out_path
        thread_number = self.thread_number
        limit_count = self.limit_count
        limit_min_length = self.limit_min_length
        limit_max_length = self.limit_max_length
        scaffold_or_not = self.scaffold_or_not

        change_seed = self.change_seed
        k2 = self.k2
        assemble_path = self.assemble_path  # assemble.py

        files = get_file_list(boostrap_out_filtered_out_path)
        if files == []:
            return 0
        cmd = 'python "{0}" -i "{1}" -r "{2}" -o "{3}" -k2 {4} -limit_count {5} -limit_min_length {6}  -limit_max_length {7}  -change_seed {8}  -scaffold {9}  -t {10}  '.format(
            assemble_path,
            boostrap_out_filtered_out_path, boostrap_out_reference_database_path, boostrap_out_assembled_out_path, k2, limit_count, limit_min_length,
            limit_max_length, change_seed, scaffold_or_not, thread_number)
        run_command(cmd)

    def bootstrap_get_results_contig(self):
        path1 = self.boostrap_out_assembled_out_path
        path2 = os.path.join(path1, "short_contig")
        path3 = os.path.join(path1, "contig")
        path4 = os.path.join(path1, "scaffold")

        boostrap_out_GM_results_path = self.boostrap_out_GM_results_path
        if not os.path.isdir(boostrap_out_GM_results_path):
            os.mkdir(boostrap_out_GM_results_path)
        path_list = [path2, path3, path4]  # short contig scaffold

        results = []  # 冗余
        GM_results_list = []
        gene_name_list = []
        for i in path_list:
            if is_exist(i):
                fasta_file = get_fasta_file(i)
                results.extend(fasta_file)
        if results == []:
            return 0
        for i in results:
            name = get_basename(i)
            if name not in gene_name_list:
                gene_name_list.append(name)
                GM_results_list.append(i)

        for i in GM_results_list:
            new_path = os.path.join(boostrap_out_GM_results_path, get_basename(i) + ".fasta")
            shutil.copy(i, new_path)

    #默认剔除N后计算一致度
    def get_bootstrap_information(self,gm_result, bootstrap_result, log_path, bootstrap_number, gene_name):
        gene = gene_name
        bootstrap_number = bootstrap_number
        if is_exist(gm_result) and is_exist(bootstrap_result):
            seq1 = get_seq(gm_result, max_seq_number=1, seq_count_limit=True)
            gm_result_seq = list(seq1[0].values())[0]
            seq2 = get_seq(bootstrap_result, max_seq_number=1, seq_count_limit=True)
            bootstrap_result_seq = list(seq2[0].values())[0]
            score = format(get_identity(gm_result_seq, bootstrap_result_seq)*100,".2f") #format控制位数比round好

        else:
            score = "failed"
        sth = [gene, bootstrap_number, str(score)]

        mylog(log_path, sth)

    def get_bootstrap_information_paralle(self):
        GM_results_path=self.GM_results_path
        thread_number=self.thread_number
        boostrap_out_path=self.boostrap_out_path
        boostrap_out_GM_results_path=self.boostrap_out_GM_results_path
        bootstrap_csv_path=os.path.join(boostrap_out_path,"bootstrap.csv")
        bootstrap_number=self.bootstrap_number

        files=get_fasta_file(GM_results_path)
        if files==[]:
            return 0
        task_pool=[]
        results=[]
        executor=futures.ProcessPoolExecutor(thread_number)

        header=["gene","bootstrap_number","score"]
        mylog(bootstrap_csv_path,header)
        for i in files:
            name=get_basename(i)
            bootstrap_result=os.path.join(boostrap_out_GM_results_path,name+".fasta")
            task_pool.append(executor.submit(self.get_bootstrap_information,i,bootstrap_result,bootstrap_csv_path,bootstrap_number,name))
        for i in task_pool:
            results.append(i.result())
        executor.shutdown()




def my_bootstrap_pipeline_main(configuration_information):
    t1=time.time()
    print("")
    cutting_line(" Bootstrap ")
    print("Using GeneMiner...")
    my_bootstrap_pipeline = BootstrapPipeLine(configuration_information)
    flag=my_bootstrap_pipeline.get_mutated_sequence_parallel()
    if flag:
        my_bootstrap_pipeline.bootstrap_filter()
        my_bootstrap_pipeline.bootstrap_assemble()
        my_bootstrap_pipeline.bootstrap_get_results_contig()
        my_bootstrap_pipeline.get_bootstrap_information_paralle()
        t2=time.time()
        print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
        used_temp_time = format((t2 - t1), ".2f")
        print("Bootstrap time used: {}s".format(used_temp_time))
    else:
        t2 = time.time()
        print(" " * 50, flush=True, end='\r')  # 之前的信息比较长，冲刷掉
        used_temp_time = format((t2 - t1), ".2f")
        print("Bootstrap Failed: {}s".format(used_temp_time))




if __name__ == '__main__':
    data1=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\data1_2000w.fq"
    data2 =r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\data2_2000w.fq"
    single=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\data1.fq"

    out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\demo"
    rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\cp_gene"
    rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"

    # out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\matk_bootstrap"
    # rtfa = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\matk_ref"
    # rtgb = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\example\ref_gb\chuanxiong.gb"

    k1=17
    k2=31
    data_size='all'
    step_length=4
    limit_count="auto"
    limit_min_length=0.5
    limit_max_length=2
    change_seed=32
    scaffold_or_not=True
    max_length = 50000
    min_length = 0
    thread_number=4
    soft_boundary = 0
    bootstrap_information=[True,10]
    bootstrap=bootstrap_information[0]
    bootstrap_number=bootstrap_information[1]



    reference_database = "reference_database"
    filtered_out = "filtered_out"
    assembled_out = "assembled_out"
    bootstrap_out = "bootstrap_out"
    GM_results = "GM_results"
    results_log = "results.log"

    filter_path=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\lib\my_filter.py"
    assemble_path=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\lib\my_assemble.py"
    muscle_path=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeeeeee9 重构filter\lib\muscle3"

    #其他信息
    my_software_name = "GM"
    configuration_information = {"out_dir": out_dir,
                                 "data1": data1, "data2": data2, "single": single,
                                 "rtfa": rtfa, "rtgb": rtgb,
                                 "k1": k1, "k2": k2, "thread_number": thread_number,
                                 "step_length": step_length,
                                 "limit_count": limit_count,
                                 "limit_min_length": limit_min_length,
                                 "limit_max_length": limit_max_length,
                                 "change_seed": change_seed,
                                 "scaffold_or_not":scaffold_or_not,
                                 "max_length": max_length, "min_length": min_length,
                                 "soft_boundary": soft_boundary, "data_size": data_size,
                                 "bootstrap": bootstrap_information[0], "bootstrap_number": bootstrap_information[1],
                                 "reference_database": reference_database,
                                 "filtered_out": filtered_out, "assembled_out": assembled_out,
                                 "bootstrap_out": bootstrap_out,
                                 "GM_results": GM_results,
                                 "results_log": results_log,
                                 "my_software_name": my_software_name,
                                 "filter_path": filter_path, "assemble_path": assemble_path, "muscle_path": muscle_path
                                 }

    # my_bootstrap_pipeline=BootstrapPipeLine(configuration_information)
    # my_bootstrap_pipeline.get_mutated_sequence_parallel()
    # my_bootstrap_pipeline.bootstrap_filter()
    # my_bootstrap_pipeline.bootstrap_assemble()
    # my_bootstrap_pipeline.bootstrap_get_results_contig()
    # my_bootstrap_pipeline.get_bootstrap_information_paralle()
    my_bootstrap_pipeline_main(configuration_information)





