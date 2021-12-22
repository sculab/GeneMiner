#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/28 11:39
# @Author  : xiepulin
# @File    : pack_results.py
# @Software: PyCharm

import os
import shutil
import platform
from basic import dir_make,cutting_line


def pack_the_results(type,configuration,out_dir_name,bootstrap_number,callback_number):
    message="Packing_results...."
    print(message,end="",flush=True)
    out_dir_name=out_dir_name
    reference_database=configuration["reference_database"]
    filtered_out=configuration["filtered_out"]
    assembled_out=configuration["assembled_out"]
    assembled_log=configuration["assembled_log"]
    GM_results=configuration["GM_results"]
    bootstrap_out=configuration["bootstrap_out"]
    callback_out=configuration["callback_out"]

    reference_database_path=os.path.join(out_dir_name,reference_database)
    filtered_out_path=os.path.join(out_dir_name,filtered_out)
    assembled_out_path=os.path.join(out_dir_name,assembled_out)
    assembled_log_path=os.path.join(out_dir_name,assembled_log)
    GM_results_path=os.path.join(out_dir_name,GM_results)
    bootstrap_out_path=os.path.join(out_dir_name,bootstrap_out)
    callback_out_path=os.path.join(out_dir_name,callback_out)

    pack_list=[reference_database_path,filtered_out_path,assembled_out_path,GM_results_path]
    rm_list=[assembled_log_path]

    try:
        if bootstrap_number!=None and bootstrap_number !="None":        #考虑虽然加了bootstrap参数，但是没有做出结果的情况
            pack_list.append(bootstrap_out_path)
        if callback_number!=None and callback_number !="None":
            pack_list.append(callback_out_path)

        dir=os.path.join(out_dir_name,type+"_genes")
        dir_make(dir)
        for i in rm_list:
            shutil.rmtree(i)
        for i in pack_list:
            target=os.path.join(dir,os.path.basename(i))
            shutil.move(i,target)


        print("....ok",flush=True)


        #判断是否正确挖掘出基因(一个都没有做出来的情况)
        GM_results_path = os.path.join(dir,GM_results)       #out_dir/sub_out_dir/GM_results
        GM_results_list = os.listdir(GM_results_path)

        if len(GM_results_list) == 0:
            message1 = "Failed to extract genes"
            cutting_line(message1)
            # print("select a closer reference or increase the amount of datasize")
        else:
            message = "Finished successfully, the genes have been extracted"
            cutting_line(message)
    except:
        print("....ok", flush=True)


    print("")  # 打印空行，为了好看




# system=platform.system().lower()
# configuration_information = {'mito_dir': 'mito_genes', 'cp_dir': 'cp_genes', 'tfa_dir': 'target_genes_from_fa',
#                              'tgb_dir': 'target_genes_from_gb', 'data1': 'data1.fq', 'data2': 'data2.fq',
#                              'results_information_excel': 'results_information.xlsx',
#                              'reference_database': 'reference_database', 'filtered_out': 'filtered_out',
#                              'assembled_out': 'assembled_out', 'assembled_log': 'assembled_log',
#                              'callback_out': 'callback_out', 'bootstrap_out': 'bootstrap_out',
#                              'GM_results': 'GM_results',
#                              'filter_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee4 重构版本  大修\\lib\\filter',
#                              'assemble_path': 'D:\\Happy_life_and_work\\scu\\python\\Gene_Miner\\eeee4 重构版本  大修\\lib\\minia',
#                              'my_software_name': 'GM', 'system': system, 'whole_log': 'log.txt'}
#
# # out_dir=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee4 重构版本  大修\first"
# out_dir=r"D:\Happy_life_and_work\scu\python\Gene_Miner\eeee5 重构版本  大修 更新bootstrap 更新梯度逼近\new6"
# thread_number=8
# kmer=43
# bootstrap_number=10
# callback_number ="None"
#
# pack_the_results("tfa",configuration_information,out_dir,bootstrap_number,callback_number)
#
#
#
