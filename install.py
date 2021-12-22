#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/10/12 16:17
# @Author  : xiepulin
# @File    : setup.py.py
# @Software: PyCharm

import sys
import subprocess
import tempfile
import re
import os
import platform

#检测python版本
def check_python_version():
    this_python = sys.version_info[:2]
    min_version = (3, 6)
    if this_python < min_version:
        message_parts = [
            "GeneMiner does not work on Python {}.{}.".format(*this_python),            #*可以用来接受元组
            "The minimum supported Python version is {}.{}.".format(*min_version),
            "Please download the eligible version from https://www.python.org/.".format(*this_python)]
        print("ERROR: " + " ".join(message_parts))
        sys.exit(1)

#检测pip是否存在
def check_pip():
    cmd="pip3 --version >/dev/null "
    flag=subprocess.call(cmd,shell=True)
    if flag!=0:
        message_parts=["pip3 does not detect in the system.",
                       "You can download the appropriate version from https://pip.pypa.io/en/stable/installation/ ."
                       ]
        print("ERROR: " + " ".join(message_parts))
        sys.exit(1)


#将字符串转换为元组
def version_str2tuple(version):
    # version 1.79    pandas :1.1.15
    version = version.split(".")  # 1.1.15 --- [1,1,15]
    version = [int(i) for i in version]  # str--int
    version = tuple(version)  # [1,1,15] --- (1,1,15)
    return version


#检测依赖并安装
def check_dependencies():
    pip_temp_file = tempfile.NamedTemporaryFile()
    pip_temp_file_name = pip_temp_file.name

    requirements_temp_file = tempfile.NamedTemporaryFile()
    requirements_temp_file_name = requirements_temp_file.name

    requirements = {"biopython": "1.79", "pandas": "1.1.5", "tqdm": "4.62.3", "openpyxl": "3.0.9"}
    undownloaded_requirements = []
    key_list = list(requirements.keys())

    # DEPRECATION: The default format will switch to columns in the future. You can use --format=(legacy|columns)
    # (or define a format=(legacy|columns) in your pip.conf under the [list] section) to disable this warning.
    # 本质原因，9.0.1之前的pip 需要加入--format才能正常输出
    cmd = "pip3 list --format=columns >{}".format(pip_temp_file_name)
    subprocess.call(cmd, shell=True)
    for i in key_list:
        if i != "biopython":
            cmd = "cat {0}|grep -e {1} ".format(pip_temp_file_name, i)
            information = subprocess.getoutput(cmd)
            version = re.findall(r"\d+\.\d+\.\d+", information)
            if version == []:
                undownloaded_requirements.append(i)
            else:
                version = version[0]
                version = version_str2tuple(version)
                if version < version_str2tuple(requirements[i]):
                    undownloaded_requirements.append(i)
        else:
            cmd = "cat {0}|grep -e {1} ".format(pip_temp_file_name, i)
            information = subprocess.getoutput(cmd)
            version = re.findall(r"\d+\.\d+", information)
            if version == []:
                undownloaded_requirements.append(i)
            else:
                version = version[0]
                version = version_str2tuple(version)
                if version < version_str2tuple(requirements[i]):
                    undownloaded_requirements.append(i)

    if undownloaded_requirements == []:
        print("OK, All dependencies have been satisfied")
        sys.exit(0)

    uninstalled_requirements_file = tempfile.NamedTemporaryFile(mode="w+t")  # 读写txt的形式
    uninstalled_requirements_file_name = uninstalled_requirements_file.name
    for i in undownloaded_requirements:
        message = "{0}=={1}\n".format(i, requirements[i])
        uninstalled_requirements_file.write(message)

    uninstalled_requirements_file.seek(0)
    print(uninstalled_requirements_file.read())

    try:
        cmd = "pip3 install -r {} --user".format(uninstalled_requirements_file_name)
        subprocess.call(cmd, shell=True)
    except:
        pass
    pip_temp_file.close()
    requirements_temp_file.close()

#检测路径并写入
def check_path():
    py_path = os.path.dirname(os.path.abspath(__file__))
    current_os = platform.system().lower()
    has_path = False
    if current_os == "darwin": #如果是macos
        for i in os.environ.get("PATH").split(":"):
            if i == py_path:
                has_path = True
                break
        if has_path == False:
            write_path = input("Write GeneMiner's path to user's environment variable? (y/n)") or "y"
            if write_path.lower() == "y" or write_path.lower() == "yes":
                # 兼容macOS Catalina之前的系统
                with open(os.environ['HOME'] + '/.bash_profile', 'a') as f:
                    f.write("\nexport PATH=$PATH:" + py_path)
                # 兼容macOS Catalina之后的系统
                with open(os.environ['HOME'] + '/.zshrc', 'a') as f:
                    f.write("\nexport PATH=$PATH:" + py_path)
                os.environ["PATH"] += ":" + py_path
                os.system('zsh')

    elif current_os == "linux":
        pass

if __name__ == '__main__':
    check_python_version()
    check_pip()
    check_dependencies()
    check_path()