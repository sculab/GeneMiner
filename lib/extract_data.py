#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/20 16:47
# @Author  : xiepulin
# @File    : extract_data.py
# @Software: PyCharm

import subprocess
import os
import sys
# from basic import dir_make,runCommand,get_absolute_and_map_path
from basic import *




##########################################################
##########################################################
'''
第二部分
提取原始数据
'''
###########################################################
###########################################################

def extract_raw_data(out_dir,raw_data_dict, data1, data2, paired, single, data_size,configuration_information):
    #raw_data_dict={3:['gw1.fq']}}
    system=configuration_information["system"]

    type = list(raw_data_dict.keys())[0]
    raw_data_list = list(raw_data_dict.values())[0]
    if type == [] or raw_data_list == []:
        gv.set_value("my_gui_flag", 0)
        sys.exit()
    print("Collecting_raw_data....", end="", flush=True)
    dir_make(out_dir)

    path1=os.path.join(out_dir,"data1.fq")
    path2=os.path.join(out_dir,"data2.fq")

    path_list=[data1,data2,paired,single,path1,path2]
    if data_size=="all":
        if data1 in raw_data_list and data2 in raw_data_list:
            [data1,data2,paired,single,path1,path2]=get_absolute_and_map_path(path_list,system)  #如果为windows环境，会批量映射路径
            if type == 1 or type == 2:
                cmd1 = "zcat <'{0}' >'{1}' ".format(data1, path1)
                runCommand(cmd1,system)
                cmd2 = "zcat <'{0}' >'{1}' ".format(data2, path2)
                runCommand(cmd2,system)
            elif type == 3:
                cmd1 = "cat '{0}'  >'{1}' ".format(data1, path1)
                runCommand(cmd1,system)
                cmd2 = "cat '{0}'  >'{1}' ".format(data2, path2)
                runCommand(cmd2,system)
            else:
                pass
        elif paired in raw_data_list:
            [data1, data2, paired, single, path1, path2] = get_absolute_and_map_path(path_list,system)  # 如果为windows环境，会批量映射路径
            if type == 1 or type == 2:
                cmd1 = "zcat <'{0}'  >'{1}' ".format(paired, path1)
                runCommand(cmd1, system)
                cmd2 = "zcat <'{0}' >'{1}' ".format(paired, path2)
                runCommand(cmd2, system)
            elif type == 3:
                cmd1 = "cat '{0}'  >'{1}' ".format(paired, path1)
                runCommand(cmd1, system)
                cmd2 = "cat '{0}'  >'{1}' ".format(paired, path2)
                runCommand(cmd2, system)
            else:
                pass
        elif single in raw_data_list:
            [data1, data2, paired, single, path1, path2] = get_absolute_and_map_path(path_list,
                                                                                     system)  # 如果为windows环境，会批量映射路径
            if type == 1 or type == 2:
                cmd1 = "zcat <'{0}'  >'{1}' ".format(single,path1)
                runCommand(cmd1, system)
                cmd2 = "zcat <'{0}' >'{1}' ".format(single,path2)
                runCommand(cmd2, system)
            elif type == 3:
                cmd1 = "cat '{0}' >'{1}' ".format(single,path1)
                runCommand(cmd1, system)
                cmd2 = "cat '{0}' >'{1}' ".format(single,path2)
                runCommand(cmd2, system)
            else:
                pass
        else:
            pass

    else:
        if data1 in raw_data_list and data2 in raw_data_list:
            [data1, data2, paired, single, path1, path2] = get_absolute_and_map_path(path_list,
                                                                                     system)  # 如果为windows环境，会批量映射路径

            if type == 1 or type == 2:
                cmd1 = "zcat <'{0}' |head -{1} >'{2}' ".format(data1, data_size,path1)
                runCommand(cmd1, system)
                cmd2 = "zcat <'{0}' |head -{1} >'{2}' ".format(data2, data_size,path2)
                runCommand(cmd2, system)
            elif type == 3:
                cmd1 = "cat '{0}' |head -{1} >'{2}' ".format(data1, data_size,path1)
                runCommand(cmd1, system)
                cmd2 = "cat '{0}' |head -{1} >'{2}' ".format(data2, data_size,path2)
                runCommand(cmd2, system)
            else:
                pass
        elif paired in raw_data_list:
            [data1, data2, paired, single, path1, path2] = get_absolute_and_map_path(path_list,
                                                                                     system)  # 如果为windows环境，会批量映射路径
            if type == 1 or type == 2:
                cmd1 = "zcat <'{0}' |head -{1} >'{2}' ".format(paired, data_size, path1)
                runCommand(cmd1, system)
                cmd2 = "zcat <'{0}' |head -{1} >'{2}' ".format(paired, data_size, path2)
                runCommand(cmd2, system)
            elif type == 3:
                cmd1 = "cat '{0}' |head -{1} >'{2}' ".format(paired, data_size, path1)
                runCommand(cmd1, system)
                cmd2 = "cat '{0}' |head -{1} >'{2}' ".format(paired, data_size, path2)
                runCommand(cmd2, system)
            else:
                pass
        elif single in raw_data_list:
            [data1, data2, paired, single, path1, path2] = get_absolute_and_map_path(path_list,
                                                                                     system)  # 如果为windows环境，会批量映射路径
            if type == 1 or type == 2:
                cmd1 = "zcat <'{0}' |head -{1} >'{2}' ".format(single, data_size, path1)
                runCommand(cmd1, system)
                cmd2 = "zcat <'{0}' |head -{1} >'{2}' ".format(single, data_size, path2)
                runCommand(cmd2, system)
            elif type == 3:
                cmd1 = "cat '{0}' |head -{1} >'{2}' ".format(single, data_size, path1)
                runCommand(cmd1, system)
                cmd2 = "cat '{0}' |head -{1} >'{2}' ".format(single, data_size, path2)
                runCommand(cmd2, system)
            else:
                pass
        else:
            pass
    print("....ok", flush=True)








