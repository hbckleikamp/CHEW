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
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()

#%% import 
import requests
import ftputil, urllib, gzip, zipfile, shutil, tarfile
import Bio
from Bio import SeqIO

import urllib
import ftputil
import requests
import zipfile
from openpyxl import Workbook, load_workbook 
import datetime
import pandas as pd
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
        if not os.path.exists(filename): #check if file is already downloaded
            while True:
                try:
                    urllib.request.urlretrieve(url,filename)
                    break
                except:
                    time.sleep(2)
                    print("retry")
        
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


#%% Protein database construction

#Diamond
url="https://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-windows.zip"    
path=str(Path(basedir,"Diamond"))
download_extract(url,path)

##Uniprot based
#Swissprot
url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
path=str(Path(basedir,"swissprot"))
download_extract(url,path)

#Trembl
url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
path=str(Path(basedir,"Trembl"))
download_extract(url,path)

# #%% merge down databases (swissprot+Trembl, refseq/GTDB)

#Uniprot KB
dirs=["swissprot","Trembl"]
path="UniprotKB"
outfile="UniprotKB.fasta"
outpath=str(Path(basedir,path,outfile))
if not os.path.exists(path): os.mkdir(path)
files=[]
for d in dirs:
    [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]
merge(outpath,files)

#Cleanup
for f in files:
    os.remove(f)

# #Uniref 100
# url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
# path=str(Path(basedir,"Uniref100"))
# download_extract(url,path)

# #Uniref 90
# url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
# path=str(Path(basedir,"Uniref90"))
# download_extract(url,path)

# #Uniref 50
# url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
# path=str(Path(basedir,"Uniref50"))
# download_extract(url,path)

# #Refseq (protein and nonredundant protein files)
# folder="Refseq_protein_db"
# ftpbase='ftp.ncbi.nlm.nih.gov'
# ftpdir='/refseq/release/complete/'
# host = ftputil.FTPHost(ftpbase, 'anonymous', 'password')
# host.chdir(ftpdir)
# dir_list = host.listdir(host.curdir)
# path=str(Path(basedir,"Refseq"))
# for link in dir_list:  
#     #if link.endswith(".faa.gz"): #all
#     if link.endswith(".faa.gz") and Path(link).stem.startswith("complete.nonredundant_protein") : #only nonredundant
        
#         print(link)
#         url="https://"+ftpbase+ftpdir+link
#         download_extract(url,path)
        

# #Refseq 
# dirs=["Refseq"]
# files=[]
# for d in dirs:
#     [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]

# #NR
# path="Refseq_NR_Protein"
# outfile="Refseq_NR_Protein.fasta"
# outpath=str(Path(basedir,path,outfile))
# if not os.path.exists(path): os.mkdir(path)
# nrfiles=[f for f in files if "nonredundant" in f]
# merge(outpath,nrfiles)

# #Protein
# path="Refseq_Protein"
# outfile="Refseq_Protein.fasta"
# outpath=str(Path(basedir,path,outfile))
# if not os.path.exists(path): os.mkdir(path)
# merge(outpath,files)

# #Cleanup
# for f in files:
#     os.remove(f)


# # # GTDB 207
# url="https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"
# url="https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"

# path=str(Path(basedir,"GTDB"))
# download_extract(url,path)

# #rename GTDB files, equate IL to J, remove ambiguous amino acids
# dirs=[str(Path("GTDB","protein_faa_reps",d)) for d in ["bacteria","archaea"]]
# path=Path("GTDB","renamed")
# outfile="GTDB.fasta"
# outpath=str(Path(basedir,path,outfile))
# if not os.path.exists(path): os.mkdir(path)
# files=[]
# for d in dirs:
#     [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]

# for file in files:
#     print(file)
#     rs=[]
#     for record in SeqIO.parse(file,format="fasta"):
#         org=Path(file).stem.split("_protein")[0]
#         record.id=record.id+"_"+org
#         record.name=record.name+"_"+org
#         rs.append(record)
#     SeqIO.write(rs,
#                 str(Path(basedir,path,Path(file).stem+".fasta")),
#                 "fasta")



# # #GTDB (r207)
# dirs=[str(Path("GTDB","renamed"))]
# path="GTDB_merged"
# outfile="GTDB_merged.fasta"
# outpath=str(Path(basedir,path,outfile))
# if not os.path.exists(path): os.mkdir(path)
# files=[]
# for d in dirs:
#     [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]
# merge(outpath,files)


# Metadata


# #Uniprot KB mapper
# url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab+.gz"
# path=str(Path(basedir,"Metadata","UniprotKB"))
# download_extract(url,path)

# #GTDB 

# #taxonomy
# urls=["https://data.gtdb.ecogenomic.org/releases/latest/ar53_taxonomy.tsv.gz",
#       "https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz"]

# path=str(Path(basedir,"Metadata"))
# download_extract(urls,path)


# #metadata
# urls=["https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tar.gz",
#       "https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz"]
# path=str(Path(basedir,"Metadata"))
# download_extract(urls,path)

