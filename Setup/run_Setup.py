# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:01:54 2022

@author: ZR48SA
"""

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd()))
print(os.getcwd())

#%% Run Setup

import subprocess


Scripts=[
"1_Download_diamond.py",
"2_Create_taxdf.py",
"3_Download_Uniref50.py",
"4_Download_unclustered_database.py",
"5_Download_SMSNet_model.py",
"6_Setup_SMSNet_environment.py"]

#ADD Pepnet Setup

command="cd "+basedir+" && "
command+=" && ".join(["python "+ str(Path(basedir,script)) for script in Scripts])


stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()