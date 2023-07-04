# -*- coding: utf-8 -*-


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:34:40 2021

@author: hugokleikamp
"""

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd()).parents[0]) #change base directory to HybridCycler
os.chdir(basedir)
print(os.getcwd())

#%% import 
import requests
import ftputil, urllib, gzip, zipfile, shutil, tarfile
import subprocess
import time





#%% Download single file

def download_extract(urls,path):
    
    if type(urls)==str:
        urls=[urls]
    
    for url in urls:
        print(url)
        
        if not os.path.exists(path): os.mkdir(path)
        filename  = str(Path(path,url.split("/")[-1]))
        
        #download
        retry_counter=0
        if not os.path.exists(filename): #check if file is already downloaded
            while True:
                try:
                    urllib.request.urlretrieve(url,filename)
                    break
                except:
                    time.sleep(2)
                    retry_counter+=1
                    print("retry: "+str(retry_counter))
                    if retry_counter>5:
                        print("downloading failed, try manual download using url specified in script")
                        break
                    
        
        if not os.path.exists(Path(filename).stem): #check if file is already there extracted
            #recursive extraction            
            while any([f.endswith((".zip",".gz",".tar")) for f in os.listdir(path)]):
                
                for f in os.listdir(path):
                    
                    i=str(Path(path,Path(f)))
                    o=str(Path(path,Path(f).stem))
                    
                    if f[0].isalnum() and f.endswith(".zip"):
                        print("extracting "+f)
                        with zipfile.ZipFile(i, 'r') as zip_ref:
                            zip_ref.extractall(path)
                        if os.path.exists(i): os.remove(i)
        
                    if f[0].isalnum() and f.endswith(".gz"):
                        print("extracting "+f)
                        with gzip.open(i,'rb') as f_in:
                            with open(o,'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                        if os.path.exists(i): os.remove(i)
                        
                    if f[0].isalnum() and f.endswith(".tar"):
                        print("extracting "+f)
                        tar = tarfile.open(i, "r:")
                        tar.extractall(path)
                        tar.close()
                        if os.path.exists(i): os.remove(i)
                

def merge(outpath,files):
    with open(outpath,'wb') as wfd:
        for f in files:
            print(f)
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return


#%% Download & extract

#check if system is windows or linux

from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    url="https://github.com/bbuchfink/diamond/releases/download/v2.0.12/diamond-linux64.tar.gz"
elif platform == "win32":
    # Windows
    url="https://github.com/bbuchfink/diamond/releases/download/v2.0.12/diamond-windows.zip"    
    
path=basedir
download_extract(url,path)

