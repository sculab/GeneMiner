#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/12/5 17:38
# @Author  : xiepulin
# @File    : zzzzzzzzzzzz callback.py
# @Software: PyCharm
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:52
# @Author  : xiepulin
# @File    :  callback_pipeline.py
# @Software: PyCharm
import os
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent import futures
from tqdm import tqdm
from basic import *


class Callback_pipeline():
    def __init__(self, configuration_information, type, gm_result, filter_reference, cut_align_reference, out_dir_name,
                 thread_number, kmer, callback_number, max_length,
                 min_length, options):

        self.configuration_information = configuration_information
        self.type = type
        self.gm_result = gm_result  # 主要作用用于确定基因名
        self.filter_reference = filter_reference  # 过滤参考序列，绝对路径
        self.cut_align_reference = cut_align_reference  # 剪切参考序列，绝对路径
        self.out_dir_name = out_dir_name
        self.thread_number = thread_number
        self.kmer = kmer
        self.callback_number = callback_number
        self.max_length = max_length
        self.min_length = min_length
        self.options = options  # yes or no

        self.my_software_name = self.configuration_information["my_software_name"]
        self.filter_software = self.configuration_information["filter_path"]  # filter_reads.pl / filter  绝对路径
        self.assemble_software = self.configuration_information["assemble_path"]  # minia                   绝对路径
        self.reference_database = self.configuration_information["reference_database"]  # 参考基因数据库
        self.filtered_out = self.configuration_information["filtered_out"]
        self.assembled_out = self.configuration_information["assembled_out"]
        self.GM_results = self.configuration_information["GM_results"]
        self.callback_out = self.configuration_information["callback_out"]

        self.system = self.configuration_information["system"]
        self.results_information_excel = self.configuration_information[
            "results_information_excel"]  # "results_information.xlsx"
        self.whole_log = self.configuration_information["whole_log"]

        gene_name = os.path.basename(self.gm_result)
        if "_trimmed.fasta" in gene_name:
            gene_name = gene_name.split("_trimmed.fasta")[0]  # 方便后续扩展
        elif "_raw_best.fasta" in gene_name:
            gene_name = gene_name.split("_raw_best.fasta")[0]
        else:
            gene_name = gene_name.split(".fasta")[0]

        self.gene_name = gene_name

    def callback(self):
        configuration_information=self.configuration_information
        type = self.type
        out_dir_name = self.out_dir_name
        sub_out_dir_name = self.callback_out
        thread_number = 1
        callback_number = self.callback_number
        kmer = self.kmer
        max_length = self.max_length
        min_length = self.min_length
        options = self.options
        gm_result_first=self.gm_result
        gene_name=self.gene_name
        whole_log=self.whole_log
        system=self.system

        whole_log_path=os.path.join(out_dir_name,whole_log)


        for i in range(callback_number):
            gene_name = self.gene_name
            cut_align_reference = self.cut_align_reference
            real_sub_out_dir_name = os.path.join(sub_out_dir_name, gene_name, "callback_{}".format(str(i)))

            information = self.get_gm_result(i)
            if information == []:
                return 0
            else:
                gm_result = information[0]
                filter_reference = information[1]
                # print(information)
                my_pipeline = Simplified_pipeline_bootstrap_callback(configuration_information, type, gm_result,
                                                                     filter_reference, cut_align_reference,
                                                                     out_dir_name, real_sub_out_dir_name, thread_number,
                                                                     kmer, max_length, min_length, options)
                my_pipeline.filter_reads()
                my_pipeline.assemble_reads()
                my_pipeline.check_contigs()


                information_next=self.get_gm_result(i+1)
                gm_result_next=information_next[0]

                identity_and_coverage=get_identity_and_coverage_path(gm_result_next, gm_result_first)
                identity=identity_and_coverage[0]
                coverage=identity_and_coverage[1]
                # print(identity)
                # print(coverage)
                if identity<=60 or coverage<=50:
                    print("okkk")
                    path_list=[whole_log_path]
                    [whole_log_path]=get_absolute_and_map_path(path_list,system)
                    message="{0}_{1}:Warning, the raw data does not support further iterations.The callback failed".format(type,gene_name)
                    cmd = "echo {0} >> '{1}' ".format(message, whole_log_path)
                    runCommand(cmd, system)
                    break







    '''
    根据提供的索引 返回gm_result,和reference_filter
    如从 callback_0  callback_1   callback_2  callback_3   callback_4 中找对应的gm_result，此时的gm_result就是下一次的filter_reference
    information= [gm_result, filter_reference]
    '''

    def get_gm_result(self, index):
        out_dir_name = self.out_dir_name
        sub_out_dir_name = self.callback_out
        gene_name = self.gene_name
        GM_results = self.GM_results

        information = []

        '''
        第一次： gm_result 软件直接生成的结果,此时的index=0
        第二次： 从callback_0中找，此时的index=1
        ....   :....
        第n次： 从callback(n-1)中找，此时的index=n  刚好未用户设定的次数

        '''
        if index == 0:
            gm_result = self.gm_result
            filter_reference = self.filter_reference
            information = [gm_result, filter_reference]
        else:
            trimmed = "_trimmed"
            raw_best = "_raw_best"
            callback_dir = "callback_{}".format(str(index - 1))  # 从上一次的结果中找
            path = os.path.join(out_dir_name, sub_out_dir_name, gene_name, callback_dir, GM_results)

            if is_exist(path) == 1:
                temp1 = get_absolute_path_by_name(path, "_trimmed")
                temp2 = get_absolute_path_by_name(path, "_raw_best")
                if temp1 != []:
                    gm_result = temp1[0]
                    filter_reference = temp1[0]
                    information = [gm_result, filter_reference]
                elif temp2 != []:
                    gm_result = temp2[0]
                    filter_reference = temp2[0]
                    information = [gm_result, filter_reference]
                else:
                    pass

        return information


class Run_callback_pipeline():
    def __init__(self, configuration_information, type, out_dir_name, thread_number, kmer, callback_number, max_length,
                 min_length, options):

        self.configuration_information = configuration_information
        self.type = type
        self.out_dir_name = out_dir_name
        self.thread_number = thread_number
        self.kmer = kmer
        self.callback_number = callback_number
        self.max_length = max_length
        self.min_length = min_length
        self.options = options  # yes or no

        self.my_software_name = self.configuration_information["my_software_name"]
        self.filter_software = self.configuration_information["filter_path"]  # filter_reads.pl / filter  绝对路径
        self.assemble_software = self.configuration_information["assemble_path"]  # minia                   绝对路径
        self.reference_database = self.configuration_information["reference_database"]  # 参考基因数据库
        self.filtered_out = self.configuration_information["filtered_out"]
        self.assembled_out = self.configuration_information["assembled_out"]
        self.GM_results = self.configuration_information["GM_results"]
        self.callback_out = self.configuration_information["callback_out"]

        self.system = self.configuration_information["system"]
        self.results_information_excel = self.configuration_information[
            "results_information_excel"]  # "results_information.xlsx"
        self.whole_log = self.configuration_information["whole_log"]

        # out_dir\GM_results\ORTHOMCL16390_raw_best.fasta  -----gene_name

    def prepare_callback_data(self):
        out_dir_name = self.out_dir_name
        reference_database = self.reference_database
        GM_results = self.GM_results

        reference_database_path = os.path.join(out_dir_name, reference_database)  # out_dir_name/reference_database
        GM_results_path = os.path.join(out_dir_name, GM_results)
        ref_files_list = os.listdir(reference_database_path)

        call_back_task_pool_raw_best = []  # raw_best 剪切未能成功
        call_back_task_pool_trimmed = []  # trimmed 剪切成功
        call_back_cut_align_ref = []

        sign_raw_best = "_raw_best.fasta"
        sign_trimmed = "_trimmed.fasta"

        for i in ref_files_list:
            i = str(i)
            gene_name = i.split(".fasta")[0]  # trnI-GAU.fasta ---- trnI-GAU

            GM_result_trimmed_name = gene_name + "_trimmed.fasta"  # trnI-GAU----trnI-GAU_trimmed.fasta
            GM_result_best_raw_name = gene_name + "_raw_best.fasta"  # trnI-GAU----trnI-GAU_raw_best.fasta
            GM_result_trimmed_path = os.path.join(GM_results_path, GM_result_trimmed_name)
            GM_result_raw_best_path = os.path.join(GM_results_path, GM_result_best_raw_name)

            ref_path = os.path.join(reference_database_path, i)
            if os.path.exists(GM_result_trimmed_path):  # 对于剪切成功的序列，再做梯度逼近没啥意义
                pass
            elif os.path.exists(GM_result_raw_best_path):
                call_back_task_pool_raw_best.append(GM_result_raw_best_path)
                call_back_cut_align_ref.append(ref_path)
            else:
                pass

        callback_information = [call_back_task_pool_raw_best, call_back_cut_align_ref]
        # [['D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\weilingcai\\GM_results\\ORTHOMCL16390_raw_best.fasta'],
        # ['D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\weilingcai\\reference_database\\ORTHOMCL16390.fasta']]
        return callback_information

    def run_callback(self, raw_best, ref_filter):
        configuration_information = self.configuration_information
        type = self.type
        out_dir = self.out_dir_name
        gm_result = raw_best
        filter_reference = raw_best
        cut_align_reference = ref_filter
        thread_number = 1
        kmer = self.kmer
        callback_number = self.callback_number
        max_length = self.max_length
        min_length = self.min_length
        options = self.options

        tfa_callback_pipeline = Callback_pipeline(configuration_information, type, gm_result, filter_reference,
                                                  cut_align_reference, out_dir, thread_number, kmer, callback_number,
                                                  max_length, min_length, options)
        tfa_callback_pipeline.callback()

    def wrap_callback(self, raw_best, ref_filter):
        self.run_callback(raw_best, ref_filter)
        return raw_best

    def run_call_back_parallel(self, callback_information):
        thread_number=self.thread_number
        raw_best_list = callback_information[0]
        ref_filter = callback_information[1]

        if callback_information[0]==[]:    #如果都剪切好了，就没必要做梯度逼近了
            return 0

        task_pool = []
        result = []
        executor = futures.ThreadPoolExecutor(max_workers=thread_number)

        for i in range(len(raw_best_list)):
            task_pool.append(executor.submit(self.wrap_callback, raw_best_list[i], ref_filter[i]))

        for i in tqdm(desc="{0:<22}".format("Callback"), iterable=futures.as_completed(task_pool),
                      total=len(task_pool)):
            result.append(i.result())

        executor.shutdown()







if __name__ == '__main__':
    '''
    windows
    '''

    # '''
    # 单个基因
    # '''
    # system = platform.system().lower()
    # configuration_information = {'mito_dir': 'mito_genes', 'cp_dir': 'cp_genes', 'tfa_dir': 'target_genes_from_fa',
    #                              'tgb_dir': 'target_genes_from_gb', 'data1': 'data1.fq', 'data2': 'data2.fq',
    #                              'results_information_excel': 'results_information.xlsx',
    #                              'reference_database': 'reference_database', 'filtered_out': 'filtered_out',
    #                              'assembled_out': 'assembled_out', 'assembled_log': 'assembled_log',
    #                              'callback_out': 'callback_out', 'bootstrap_out': 'bootstrap_out',
    #                              'GM_results': 'GM_results',
    #                              'filter_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\filter',
    #                              'assemble_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\minia',
    #                              'muscle_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\muscle3',
    #                              'my_software_name': 'GM', 'system': 'windows', 'whole_log': 'log.txt',
    #                              'bootstrap_data_set': 'bootstrap_data_set.fasta',
    #                              'bootstrap_concensusu': 'bootstrap_concensus.fasta'}
    #
    #
    # # out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\third"
    # # gm_result=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\third\GM_results\matK_raw_best.fasta"
    # # filter_reference=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\third\GM_results\matK_raw_best.fasta"
    # # cut_align_reference=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\third\reference_database\matk.fasta"
    #
    # out_dir=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\weilingcai"
    # gm_result=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\weilingcai\GM_results\ORTHOMCL16390_raw_best.fasta"
    # filter_reference=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\weilingcai\GM_results\ORTHOMCL16390_raw_best.fasta"
    # cut_align_reference=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\weilingcai\reference_database\ORTHOMCL16390.fasta"
    #
    #
    # thread_number = 8
    # kmer = 43
    # callback_number = 5
    # max_length = 5000
    # min_length = 300
    # options = "no"
    #
    # tfa_callback_pipeline=Callback_pipeline(configuration_information, "tfa",gm_result,filter_reference,cut_align_reference,out_dir, thread_number, kmer, callback_number, max_length,min_length, options)
    # tfa_callback_pipeline.callback()
    #
    #

    '''
    批量基因
    '''
    system = platform.system().lower()
    configuration_information={'mito_dir': 'mito_genes', 'cp_dir': 'cp_genes', 'tfa_dir': 'target_genes_from_fa', 'tgb_dir': 'target_genes_from_gb', 'data1': 'data1.fq', 'data2': 'data2.fq', 'results_information_excel': 'results_information.xlsx', 'reference_database': 'reference_database', 'filtered_out': 'filtered_out', 'assembled_out': 'assembled_out', 'assembled_log': 'assembled_log', 'callback_out': 'callback_out', 'bootstrap_out': 'bootstrap_out', 'GM_results': 'GM_results', 'blastn_out': 'blastn_out', 'filter_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\filter', 'assemble_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\minia', 'muscle_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\muscle3', 'makeblastdb_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\makeblastdb', 'blastn_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee6 重构版本  大修 identity\\lib\\blastn', 'my_software_name': 'GM', 'system': 'windows', 'whole_log': 'log.txt', 'bootstrap_data_set': 'bootstrap_data_set.fasta', 'bootstrap_concensusu': 'bootstrap_concensus.fasta'}


    out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\third"
    # out_dir = r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee6 重构版本  大修 identity\weilingcai"
    thread_number = 8
    kmer = 43
    callback_number = 5
    max_length = 5000
    min_length = 300
    options = "no"

    tfa_callback_pipeline = Run_callback_pipeline(configuration_information, "tfa", out_dir, thread_number, kmer,
                                                  callback_number, max_length, min_length, options)
    callback_information = tfa_callback_pipeline.prepare_callback_data()
    tfa_callback_pipeline.run_call_back_parallel(callback_information)


























