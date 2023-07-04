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

basedir=str(Path(os.getcwd()).parents[0]) #change base directory to CHEW
os.chdir(basedir)
print(os.getcwd())


#%% import 
import requests
import ftputil, urllib, gzip, zipfile, shutil, tarfile
import subprocess
import time

import pandas as pd
import numpy as np






#%%

def download(urls,path):
    
    if type(urls)==str:
        urls=[urls]
    
    for url in urls:
        print(url)
        
        if not os.path.exists(path): os.makedirs(path)
        filename  = str(Path(path,url.split("/")[-1]))
    
        
        #download
        if not os.path.exists(filename): #check if file is already downloaded
            c=0
            while True:
                try:
                    urllib.request.urlretrieve(url,filename)
                    break
                except:
                    c+=1
                    time.sleep(2)
                    print("retry")
                    if c>5: #max of 5 retries
                        break
        
        return filename
        
        

def extract(path): #path to folder

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
                    
        
def extract_subfolders(folder):

    for subdir, dirs, files in os.walk(path):
        for d in dirs:
            if os.path.isdir(d):
                extract(d)

        if os.path.isdir(subdir):
            extract(subdir)
        

   


def merge(outpath,files):
    with open(outpath,'wb') as wfd:
        for f in files:
            print(f)
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return

def download_extract(urls,path):
    
    download(urls,path)
    extract_subfolders(path)
    


#%% Download & extract

#Uniref 50
url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
path=basedir
download_extract(url,path)


