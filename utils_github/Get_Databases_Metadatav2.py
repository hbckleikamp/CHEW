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

def download(urls,path):
    
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
        
        return filename
        
# def extract(filename):
#     path=Path(filename).parents[0]        
#     if not os.path.exists(Path(filename).stem): #check if file is already there extracted
#         #recursive extraction            
#         while any([f.endswith((".zip",".gz",".tar")) for f in os.listdir(path)]):
            
#             for f in os.listdir(path):
                
#                 i=str(Path(path,Path(f)))
#                 o=str(Path(path,Path(f).stem))
                
#                 if f[0].isalnum() and f.endswith(".zip"):
#                     print("extracting "+f)
#                     with zipfile.ZipFile(i, 'r') as zip_ref:
#                         zip_ref.extractall(path)
#                     if os.path.exists(i): os.remove(i)
    
#                 if f[0].isalnum() and f.endswith(".gz"):
#                     print("extracting "+f)
#                     with gzip.open(i,'rb') as f_in:
#                         with open(o,'wb') as f_out:
#                                 shutil.copyfileobj(f_in, f_out)
#                     if os.path.exists(i): os.remove(i)
                    
#                 if f[0].isalnum() and f.endswith(".tar"):
#                     print("extracting "+f)
#                     tar = tarfile.open(i, "r:")
#                     tar.extractall(path)
#                     tar.close()
#                     if os.path.exists(i): os.remove(i)
                        
        

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

    for subdir, dirs, files in os.walk(folder):
        for s in subdir():

            extract(s)
   


def merge(outpath,files):
    with open(outpath,'wb') as wfd:
        for f in files:
            print(f)
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return

#%%




    #%%
# #%% Protein database construction

# #Diamond
# url="https://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-windows.zip"    
# path=str(Path(basedir,"Diamond"))
# download_extract(url,path)

##Uniprot based
#Swissprot
# url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
# path=str(Path(basedir,"swissprot"))
# download_extract(url,path)

# #Trembl
# url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
# path=str(Path(basedir,"Trembl"))
# download_extract(url,path)

# # #%% merge down databases (swissprot+Trembl, refseq/GTDB)

# #Uniprot KB
# dirs=["swissprot","Trembl"]
# path="UniprotKB"
# outfile="UniprotKB.fasta"
# outpath=str(Path(basedir,path,outfile))
# if not os.path.exists(path): os.mkdir(path)
# files=[]
# for d in dirs:
#     [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]
# merge(outpath,files)

# #Cleanup
# for f in files:
#     os.remove(f)

# # #Uniref 100
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
#%%
#Refseq (protein and nonredundant protein files)
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

# # #NR
# # path="Refseq_NR_Protein"
# # outfile="Refseq_NR_Protein.fasta"
# # outpath=str(Path(basedir,path,outfile))
# # if not os.path.exists(path): os.mkdir(path)
# # nrfiles=[f for f in files if "nonredundant" in f]
# # merge(outpath,nrfiles)

# #Protein
# path="Refseq_Protein"
# outfile="Refseq_Protein.fasta"
# outpath=str(Path(basedir,path,outfile))
# if not os.path.exists(path): os.mkdir(path)
# merge(outpath,files)

# #Cleanup
# for f in files:
#     os.remove(f)




# # GTDB 207
# url="https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"
# url="https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"

# path=str(Path(basedir,"GTDB"))
# download_extract(url,path)

#rename GTDB files, equate IL to J, remove ambiguous amino acids
dirs=[str(Path("GTDB","protein_faa_reps",d)) for d in ["bacteria","archaea"]]
path=Path("GTDB","renamed")
outfile="GTDB.fasta"
outpath=str(Path(basedir,path,outfile))
if not os.path.exists(path): os.mkdir(path)
files=[]
for d in dirs:
    [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]

for file in files:
    print(file)
    rs=[]
    for record in SeqIO.parse(file,format="fasta"):
        org=Path(file).stem.split("_protein")[0]
        record.id=record.id+"|"+org
        record.name=record.name+"|"+org
        rs.append(record)
    SeqIO.write(rs,
                str(Path(basedir,path,Path(file).stem+".fasta")),
                "fasta")

#%%

# # #GTDB (r207)
dirs=[str(Path("GTDB","renamed"))]
path="GTDB_merged"
outfile="GTDB_merged.fasta"
outpath=str(Path(basedir,path,outfile))
if not os.path.exists(path): os.mkdir(path)
files=[]
for d in dirs:
    [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]
merge(outpath,files)


#%%
#%% Metadata


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

#%% Genneral Diamond db prepping

# #Variables
# Path_to_db=outpath         # path to protein database that uses NCBI taxonomy


# #for Bacterial_only
# Path_to_taxonomy="" #unused
# Taxid_delimiter="" #unused

# #for Remove_ambiguous
# Ambiguous_AAs=["B","X","Z","[","("] #"O","U"

# #Add decoy
# decoy_delimiter="" #unused
# decoy_method="" #unused


# import Bio
# from Bio import SeqIO
# import pandas as pd
# import itertools
# import random
# import subprocess

# # Options
# Bacterial_only=False  # retain only Bacteria (and Archaea(!)) in database  
# Equate_IL=True        # change I and J into L 
# Remove_ambiguous=True # remove ambiguous amino acids "B","O","U","X","Z","[","(" , and J in case IL is not equated
# Add_decoy=False        # append decoy of reversed peptides

# # Functions

# def chunk_gen(it,size=10**6):
#     c=itertools.count()
#     for _,g in itertools.groupby(it,lambda _:next(c)//size):
#         yield g

# #parse output_path
# Output_path=Path_to_db
# if Bacterial_only:   Output_path=Output_path.replace(".fasta","_BacArch.fasta")
# if Remove_ambiguous: Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
# if Equate_IL:        Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
# if Add_decoy:        Output_path=Output_path.replace(".fasta","_Decoy.fasta")

# # tax database and files
# if Bacterial_only:
#     ranks=["superkingdom","phylum","class","order","family","genus","species"] 
#     taxdf=pd.read_csv(Path_to_taxonomy,sep="\t")
#     taxdf.columns=['OX']+ranks+["OS"]
#     taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)

# if not Equate_IL: Ambiguous_AAs+=["J"]

# #read    
# recs=SeqIO.parse(Path_to_db,format="fasta")
# chunks=chunk_gen(recs)

# #write IL datbase
# print("writing "+Path(Output_path).stem)
# with open(Output_path,"w+") as f:

#     for ic,c in enumerate(chunks):
#         print("chunk "+ str(ic))
        
#         taxs=[]
#         rs=[]
#         for r in c:
            
#             s=str(r.seq)
#             d=r.description
            
#             if Bacterial_only:
#                 if Taxid_delimiter in d:
#                     tax=r.description.split(Taxid_delimiter)[1].split()[0] #get taxonomic identifiers
#                     taxs.append(tax)
#                 else: print("taxonomy delimiter not found!")
            
#             if Equate_IL:
#                 s=s.replace("I","L").replace("J","L")
            
#             if Remove_ambiguous:
#                 if sum([i in s for i in Ambiguous_AAs])!=0: 
#                     continue
                
#             if Add_decoy:
#                 if decoy_method=="reverse":
#                     rs.append([decoy_delimiter+d,s[::-1]])
#                 if decoy_method=="scramble":
#                     rs.append([decoy_delimiter+d,''.join(random.sample(s, len(s)))])
#                 if Bacterial_only:
#                     taxs.append(tax)
            
#             rs.append([d,s])
            
        
#         if Bacterial_only:
#             rs=pd.Series(rs)[pd.Series(taxs).isin(taxdf.OX)].tolist() #select only those that have superkingdom Archaea or Bacteria

#         s="\n".join([">"+"\n".join(r) for r in rs])
#         f.write(s+"\n")
    

# #Cleanup
# os.remove(Path_to_db)


# #%% construct diamond database

# command="cd" +' "'+basedir+'" && '
# #./diamond makedb --in reference.fasta -d reference
# command+="diamond makedb --in "+str(Path(Output_path)) + " -d "+Path(Output_path).stem
# stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

# print(stderr)

# # #%% Cleanup

# # #if database construction was sucessful:
# # if "No space left on device" in str(stderr): #retry with smaller file batch
# #     print("out of disk space, delete Unmerged Trembl/Swissprot Folders and UniprotKB file without IL equation to free up space,")
# #     print("then try to run diamond manually according to: https://github.com/bbuchfink/diamond/wiki")

# # if "returned non-zero exit status" not in str(stderr):
# #     print("no error found, proceeding with cleanup")


# #     if os.path.exists(str(Path(basedir,"UniprotKB.dmnd"))):
# #         if os.path.getsize(str(Path(basedir,"UniprotKB.dmnd")))>os.path.getsize(str(Path(Output_path)))*0.9: #check if corresponds to  size of IL equated UniprotKB
            
# #             paths=[str(Path(basedir,path)) for path in ["swissprot","Trembl","UniprotKB"]]
# #             for path in paths:
# #                 try:
# #                     shutil.rmtree(path)
# #                 except:
# #                     pass
# #         else:
# #             print("wrong filesize, database construction might be incomplete, cleanup halted!")
# #     else:
# #         print("diamond database not found, cleanup halted!")
