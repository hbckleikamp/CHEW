# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 16:22:09 2022

@author: ZR48SA
"""

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd()).parents[0]) #change base directory to HybridCycler
os.chdir(basedir)
print(os.getcwd())

#%%

import subprocess


command=" && ".join(["conda env remove -n smsnet ",  
                     "conda create -n smsnet python=3.5.2=5 ",
                     "conda activate smsnet ",
                     "pip install -r "+str(Path(basedir, "SMSNet-master","requirements.txt"))])



stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

