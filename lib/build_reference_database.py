#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:48
# @Author  : xiepulin
# @File    : build_reference_database.py
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



##########################################################
##########################################################
'''
第三部分
提取出参考基因的fasta序列，用于构建参考基因库。主要分为两类：一类从gb格式中根据基因名提取，一类从fasta格式中提取
包括 核基因 叶绿体基因 线粒体基因 目标基因
'''
###########################################################
###########################################################



class Extract_reference_fasta():
    def __init__(self,configuration_information,type,out_dir_name,ref,soft_boundary,sf,gene_max_length=5000, gene_min_length=300):

        self.configuration_information=configuration_information  #包含各级文件名字信息
        self.type=type                  #前缀 tfa tgb cp mito
        self.out_dir_name=out_dir_name  #最大一级的输出文件夹
        self.ref=ref                    #参考基因组
        self.soft_boundary=soft_boundary #软边界
        self.sf=sf
        self.gene_max_length=gene_max_length  #基因最大长度
        self.gene_min_length=gene_min_length   #基因最小长度

        self.reference_database=self.configuration_information["reference_database"]  #"reference_database"




    '''
    extract_location_from_gb() identify_location() collect_fasta_information_from_gb() 三个函数服务于extract_fasta_from_gb（）
    '''

    '''
    提取出具体的基因位置
    非复合基因  [297:2890](-)  复合基因 join{[69127:69541](-), [97170:97602](-), [96808:98834](-)}
    '''

    def extract_location_from_gb(self, feature, gene_location):
        if "join" not in gene_location:
            # 软边界
            start = int(feature.location.start)
            end = int(feature.location.end)

            location = [start, end]
            return location
        else:
            location = re.findall(r"\d+:\d+",
                                  gene_location)  # ['69427:69541', '97370:97602', '96808:96834']
            location = ",".join(location).replace(":", ",")  # 69427,69541,97370,97602,96808,96834
            location = location.split(",")  # ['69427', '69541', '97370', '97602', '96808', '96834']
            return location

    '''
    增添软边界 增添长度限制
    '''

    def identify_location(self, gene_location, soft_boundary, start_all, end_all, gene_max_length,
                          gene_min_length):
        # 过滤起始位点 和终止位点。
        # (1)加上和减去软边界，不得超过0 或者148512 （起点和终点） (2)长度在要求范围中间
        new_location = []
        for i in range(0, len(gene_location), 2):
            start = int(gene_location[i])
            end = int(gene_location[i + 1])
            soft_start = start - soft_boundary
            soft_end = end + soft_boundary
            if soft_start <= start_all:  # 不能超出最左的边界 0
                start = start
            else:
                start = soft_start
            if soft_end >= end_all:  # 不能超出右边界 全长148512
                end = end
            else:
                end = soft_end
            if (end - start >= gene_min_length) and (end - start <= gene_max_length):
                new_location.append(start)
                new_location.append(end)
        return new_location

    '''
    提取出基因的所有信息  
    '''

    def collect_reference_information_from_gb(self, feature, gene_location, rec, strand_all, identifier,organism,accession):
        gene_information_all = []

        strand = feature.strand
        sequence = rec.seq

        if len(gene_location) == 0:
            return gene_information_all
        elif len(gene_location) == 2:
            try:
                if strand_all == strand:
                    gene_information = {}
                    gene_information["gene_name"] = feature.qualifiers["gene"][0]
                    gene_information["gene_sequence"] = sequence[gene_location[0]:gene_location[-1]]
                    gene_information["identifier"] = identifier
                    gene_information["organism"] = organism
                    gene_information["accession"] =accession

                else:
                    gene_information = {}
                    gene_information["gene_name"] = feature.qualifiers["gene"][0]
                    gene_information["gene_sequence"] = sequence[
                                                        gene_location[0]:gene_location[-1]].reverse_complement()
                    gene_information["identifier"] = identifier
                    gene_information["organism"] = organism
                    gene_information["accession"] = accession
                gene_information_all = [gene_information]
            except:
                pass
        else:
            index_nmuber = 1
            for i in range(0, len(gene_location), 2):
                if strand_all == strand:
                    gene_information = {}
                    gene_information["gene_name"] = feature.qualifiers["gene"][0] + "_" + str(index_nmuber)
                    gene_information["gene_sequence"] = sequence[gene_location[i]:gene_location[i + 1]]
                    gene_information["identifier"] = identifier
                    gene_information["organism"] = organism
                    gene_information["accession"] = accession
                else:
                    gene_information = {}
                    gene_information["gene_name"] = feature.qualifiers["gene"][0] + "_" + str(index_nmuber)
                    gene_information["gene_sequence"] = sequence[
                                                        gene_location[i]:gene_location[i + 1]].reverse_complement()
                    gene_information["identifier"] = identifier
                    gene_information["organism"] = organism
                    gene_information["accession"] = accession

                gene_information_all.append(gene_information)
                index_nmuber = index_nmuber + 1

        return gene_information_all

    def extract_reference_from_gb(self):
        type = str(self.type)
        out_dir_name=self.out_dir_name
        ref = self.ref
        soft_boundary = self.soft_boundary
        gene_max_length = self.gene_max_length
        gene_min_length = self.gene_min_length
        reference_database=self.reference_database

        out_dir = os.path.join(out_dir_name, reference_database)
        dir_make(out_dir)
        flag = file_or_directory(ref)  # 0代表文件夹，1代表文件
        if flag == 1:
            gene_name_list = []
            Genes_information = []
            for rec in SeqIO.parse(ref, "gb"):
                sequence = rec.seq
                strand_all = 1  # genbank默认为正义链
                start_all = 0
                end_all = len(rec.seq)
                organism = rec.annotations["organism"].replace(" ", "_")  # 物种名  Ligusticum_chuanxiong
                accession = rec.name.replace(" ", "_")  # 样本名/genbank id
                identifier = organism + "_" + accession

                for feature in rec.features:
                    if feature.type == "source":
                        strand_all = feature.strand
                        start_all = int(feature.location.start)
                        end_all = int(feature.location.end)
                        # print(start_all,end_all)
                        # print(strand_all)

                    # OrderedDict([('gene', ['matK']), ('locus_tag', ['KQ413_pgp084']), ('db_xref', ['GeneID:65316243'])])  鲁棒性
                    elif "gene" not in feature.qualifiers:
                        continue

                    # 两个可能跨越原点的基因：trnH-GUG 和 psbA
                    elif feature.type == "gene" and feature.qualifiers["gene"][0] == "psbA":
                        gene_information = {}
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name not in gene_name_list:
                            gene_name_list.append(gene_name)
                        else:
                            continue
                        gene_information["gene_name"] = feature.qualifiers["gene"][0]
                        gene_information["gene_sequence"] = feature.location.extract(sequence)
                        gene_information["identifier"] = identifier
                        gene_information["organism"]=organism
                        gene_information["accession"] = accession
                        Genes_information.append(gene_information)

                    elif feature.type == "gene" and feature.qualifiers["gene"][0] == "trnH-GUG":
                        gene_information = {}
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name not in gene_name_list:
                            gene_name_list.append(gene_name)
                        else:
                            continue
                        gene_information["gene_name"] = feature.qualifiers["gene"][0]
                        gene_information["gene_sequence"] = feature.location.extract(sequence)
                        gene_information["identifier"] = identifier
                        gene_information["organism"] = organism
                        gene_information["accession"] = accession
                        Genes_information.append(gene_information)

                    elif feature.type == "gene" and feature.qualifiers["gene"][0] == "rps12":
                        continue
                    elif feature.type == "gene":
                        gene_name = feature.qualifiers["gene"][0]
                        if gene_name not in gene_name_list:
                            gene_name_list.append(gene_name)
                        else:
                            continue
                        gene_location = str(feature.location)
                        # print(gene_location)
                        gene_location = self.extract_location_from_gb(feature,
                                                                      gene_location)  # 分为join和非join两类，提取出具体的基因位置
                        # print(gene_location)
                        gene_location = self.identify_location(gene_location, soft_boundary, start_all, end_all,
                                                               gene_max_length, gene_min_length)
                        # print(gene_location)
                        gene_information_all = self.collect_reference_information_from_gb(feature, gene_location, rec,
                                                                                          strand_all, identifier,organism,accession)
                        if gene_information_all != []:

                            for i in gene_information_all:
                                Genes_information.append(i)
                    else:
                        pass

            # print(Genes_information)
            '''
            写fasta格式的文件
            '''
            # gene_name_list 只记录基因名，没有考虑join情况下 (gene_1 gene_2的情况)
            new_gene_name_list = []
            for i in Genes_information:
                gene_name = i["gene_name"]
                new_gene_name_list.append(gene_name)
            # print(new_gene_name_list)
            for i in new_gene_name_list:
                my_records = []
                for j in Genes_information:
                    if j["gene_name"] == i and len(j["gene_sequence"]) >= gene_min_length and len(
                            j["gene_sequence"]) <= gene_max_length:
                        my_record = SeqRecord(seq=j["gene_sequence"], id=j["organism"],                                  #文件名为基因名  >后面的名字为物种名  description为accession号码
                                              description=j["accession"])
                        my_records.append(my_record)
                        Genes_information.remove(j)
                # print(my_records)
                # print(len(my_records))
                if my_records != []:
                    SeqIO.write(my_records, os.path.join(out_dir, i + ".fasta"), "fasta")










        if flag == 0:
            files = get_files(ref)
            Genes_information = []
            # 多一层循环，与文件的唯一区别
            for file in files:
                gene_name_list = []   #针对每一个物种来说，重复基因我们就不再统计了
                for rec in SeqIO.parse(file, "gb"):  # ref 换成file
                    sequence = rec.seq
                    strand_all = 1  # genbank默认为正义链
                    start_all = 0
                    end_all = len(rec.seq)
                    organism = rec.annotations["organism"].replace(" ", "_")  # 物种名
                    accession = rec.name.replace(" ", "_")  # 样本名/genbank id  NC_038088
                    identifier = organism + "_" + accession

                    for feature in rec.features:
                        if feature.type == "source":
                            strand_all = feature.strand
                            start_all = int(feature.location.start)
                            end_all = int(feature.location.end)
                            # print(start_all,end_all)
                            # print(strand_all)
                        # OrderedDict([('gene', ['matK']), ('locus_tag', ['KQ413_pgp084']), ('db_xref', ['GeneID:65316243'])])  鲁棒性
                        elif "gene" not in feature.qualifiers:
                            continue
                        # 两个可能跨越原点的基因：trnH-GUG 和 psbA
                        elif feature.type == "gene" and feature.qualifiers["gene"][0] == "psbA":
                            gene_information = {}
                            gene_name = feature.qualifiers["gene"][0]
                            if gene_name not in gene_name_list:
                                gene_name_list.append(gene_name)
                            else:
                                continue
                            gene_information["gene_name"] = feature.qualifiers["gene"][0]
                            gene_information["gene_sequence"] = feature.location.extract(sequence)
                            gene_information["identifier"] = identifier
                            gene_information["organism"] = organism
                            gene_information["accession"] = accession
                            Genes_information.append(gene_information)

                        elif feature.type == "gene" and feature.qualifiers["gene"][0] == "trnH-GUG":
                            gene_information = {}
                            gene_name = feature.qualifiers["gene"][0]
                            if gene_name not in gene_name_list:
                                gene_name_list.append(gene_name)
                            else:
                                continue
                            gene_information["gene_name"] = feature.qualifiers["gene"][0]
                            gene_information["gene_sequence"] = feature.location.extract(sequence)
                            gene_information["identifier"] = identifier
                            gene_information["organism"] = organism
                            gene_information["accession"] = accession
                            Genes_information.append(gene_information)

                        elif feature.type == "gene" and feature.qualifiers["gene"][0] == "rps12":
                            continue
                        elif feature.type == "gene":
                            gene_name = feature.qualifiers["gene"][0]
                            if gene_name not in gene_name_list:
                                gene_name_list.append(gene_name)
                            else:
                                continue
                            gene_location = str(feature.location)
                            # print(gene_location)
                            gene_location = self.extract_location_from_gb(feature,
                                                                          gene_location)  # 分为join和非join两类，提取出具体的基因位置
                            # print(gene_location)
                            gene_location = self.identify_location(gene_location, soft_boundary, start_all, end_all,
                                                                   gene_max_length, gene_min_length)
                            # print(gene_location)
                            gene_information_all = self.collect_reference_information_from_gb(feature, gene_location,
                                                                                              rec, strand_all,
                                                                                              identifier,organism,accession)
                            if gene_information_all != []:
                                #print(gene_information_all) #[{'gene_name': 'trnK-UUU', 'gene_sequence': Seq('GAAGAGGTGGTCTATAAATTTTGACATTTATCTATCTTTTTTTTTTGATTGTAT...AGA'), 'identifier': 'Ligusticum_chuanxiong_NC_038088', 'organism': 'Ligusticum_chuanxiong', 'accession': 'NC_038088'}]
                                for i in gene_information_all:
                                    Genes_information.append(i)
                                    # print(Genes_information)
                        else:
                            pass

                # print(Genes_information)

            # print(Genes_information)
            '''
            写fasta格式的文件
            '''
            # gene_name_list 只记录基因名，没有考虑join情况下 (gene_1 gene_2的情况)
            new_gene_name_list = []
            for i in Genes_information:
                gene_name = i["gene_name"]
                new_gene_name_list.append(gene_name)
            # print(new_gene_name_list)
            for i in new_gene_name_list:
                my_records = []
                for j in Genes_information:
                    if j["gene_name"] == i and len(j["gene_sequence"]) >= gene_min_length and len(
                            j["gene_sequence"]) <= gene_max_length:
                        my_record = SeqRecord(seq=j["gene_sequence"], id=j["organism"],
                                              description=j["accession"])
                        my_records.append(my_record)
                        Genes_information.remove(j)
                # print(my_records)
                # print(len(my_records))
                if my_records != []:
                    SeqIO.write(my_records, os.path.join(out_dir, i + ".fasta"), "fasta")



    def extract_reference_from_fasta(self):
        ref=self.ref
        type=self.type
        out_dir_name=self.out_dir_name
        reference_databese=self.reference_database
        out_dir = os.path.join(out_dir_name, reference_databese)

        dir_make(out_dir)

        flag = file_or_directory(ref)
        if flag == 0:
            files = get_files(ref)
            for file in files:
                file_name = str(os.path.basename(file))  # root/out_dir_name/PMSK.fa ---- PMSK.fa
                #统一后缀
                if ".fa" in file_name:
                    file_name=file_name.split(".fa")[0] + ".fasta"

                elif  ".fas" in file_name:
                    file_name=file_name.split(".fas")[0]+".fasta"

                elif ".fasta" in file_name:
                    file_name=file_name

                else:
                    file_name=file_name+".fasta"

                path = os.path.join(out_dir, file_name)
                self.get_pure_fasta_format_sequence(file,path)

        if flag == 1:
            file_name = str(os.path.basename(ref))  # PMSK.fa
            #统一后缀
            if ".fa" in file_name:
                file_name = file_name.split(".fa")[0] + ".fasta"
            elif ".fas" in file_name:
                file_name = file_name.split(".fas")[0] + ".fasta"
            elif ".fasta" in file_name:
                file_name = file_name
            else:
                file_name = file_name + ".fasta"

            path = os.path.join(out_dir, file_name)
            self.get_pure_fasta_format_sequence(ref,path)



    #剔除ACGTU之外的序列，如？-等等
    def get_pure_fasta_format_sequence(self,file, output):
        my_records = []
        for rec in SeqIO.parse(file, "fasta"):
            id = rec.id
            seq = str(rec.seq).upper()
            description = rec.description

            seq = list(seq)

            dataset = ["A", "C", "G", "T", "U"]

            new_seq = []
            for i in range(len(seq)):
                if seq[i] in dataset:
                    new_seq.append(seq[i])

            if new_seq != []:
                new_seq = "".join(new_seq)
                record = SeqRecord(id=id, seq=Seq(new_seq), description=description)
                my_records.append(record)

        if my_records != []:
            SeqIO.write(my_records, output, "fasta")



    def filter_refrence(self):
        type = self.type
        sf=self.sf
        out_dir_name = self.out_dir_name
        reference_databese = self.reference_database
        out_dir = os.path.join(out_dir_name, reference_databese)
        ref=out_dir   #将之前生成的reference_database当作新的参考序列的路径

        #["s1", "s2", "s3", "s4","s5"] s1不处理,s2最短，s3中位数,s4最长，s5箱线图

        if sf=="s1":
            return 0
        elif sf=="s2":
            files = get_files(ref)
            for file in files:
               get_shortest_sequence(file)

        elif sf=="s3":
            files = get_files(ref)
            for file in files:
               get_median_sequence(file)

        elif sf == "s4":
            files = get_files(ref)
            for file in files:
                get_longest_sequence(file)

        elif sf=="s5":
            files = get_files(ref)
            for file in files:
                box_plots_filter_reference(file)

        elif is_txt_file(sf):
            files=get_files(ref)
            for file in files:
                get_specify_txt_reference(file,sf)
        else:
            pass

















if __name__ == '__main__':
    configuration_information={'mito_dir': 'mito_genes', 'cp_dir': 'cp_genes', 'tfa_dir': 'target_genes_from_fa', 'tgb_dir': 'target_genes_from_gb', 'data1': 'data1.fq', 'data2': 'data2.fq', 'results_information_excel': 'results_information.xlsx', 'reference_database': 'reference_database', 'filtered_out': 'filtered_out', 'assembled_out': 'assembled_out', 'assembled_log': 'assembled_log', 'callback_out': 'callback_out', 'bootstrap_out': 'bootstrap_out', 'GM_results': 'GM_results', 'blastn_out': 'blastn_out', 'filter_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\filter', 'assemble_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\minia', 'muscle_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\muscle3', 'makeblastdb_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\makeblastdb', 'blastn_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\blastn', 'my_software_name': 'GM', 'system': 'windows', 'whole_log': 'log.txt', 'bootstrap_data_set': 'bootstrap_data_set.fasta', 'bootstrap_concensusu': 'bootstrap_concensus.fasta'}

    # sf=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeeee7 xie2yu版本\reference.txt"
    sf="s1"
    target=Extract_reference_fasta(configuration_information,"tgb",r"..\new4",r"..\example\ref_gb",75,sf,5000,0)
    # target = Extract_reference_fasta(configuration_information, "tgb", r"..\new4", r"..\example\ref_gb\chuanxiong.gb", 0, sf, 100000,
    #                                  0)
    target.extract_reference_from_gb()
    target.filter_refrence()

    # target=Extract_reference_fasta(configuration_information,"tfa","out",r"../ref.fa",75,sf,5000,300)
    # target.extract_reference_from_fasta()
    # target.filter_refrence()

    # target=Extract_reference_fasta(configuration_information,"tfa","out",r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee4 重构版本  大修\test_fasta\ORTHOMCL4460.fas",75,5000,300)
    # target.extract_reference_from_fasta()


    # target=Extract_reference_fasta(configuration_information,"tfa","out",r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee4 重构版本  大修\test_fasta",75,5000,300)
    # target.extract_reference_from_fasta()

    # target=Extract_reference_fasta(configuration_information,"tfa","out",r"D:\Happy_life_and_work\scu\library-data\委陵菜验证数据\Rosaceae\Rosaceae",75,5000,300)
    # target.extract_reference_from_fasta()
    #








    #
    # target=Extract_reference_fasta(configuration_information,"tfa","out2","../ref_gb",75,5000,300)
    # target.extract_reference_from_gb()
    # #
    # target=Extract_reference_fasta(configuration_information,"tfa","out3","../ref.fa",75,5000,300)
    # target.extract_reference_from_fasta()

    # target=Extract_reference_fasta(configuration_information,"tfa","out4","../mito.gb",75,5000,300)
    # target.extract_reference_from_gb()

































