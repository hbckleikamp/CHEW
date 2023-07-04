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


import Bio
from Bio import SeqIO
import pandas as pd
import itertools
import random


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


import argparse



svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])
#%% parameters



DB="UniprotKB"
rm_merge=True               # remove database after merging 
output_folder=basedir       # where to put databases
Path_to_taxonomy=str(Path(basedir,"parsed_taxonomy.tsv"))

prepdb=True #after downloading prep DB with following arguments
make_dmnd=True


Bacterial_only=True    # retain only Bacteria (and Archaea(!)) in database  
Equate_IL=True         # change I and J into L 
Remove_ambiguous=True  # remove ambiguous amino acids "B","X","Z","[","(" , and J in case IL is not equated
No_Fragments=False     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
No_Dump=True           # remove dump taxa (unspecific filler names of NCBI with bloated annotations, like "uncultured" or "bacterium")
Add_decoy=False        # append decoy of reversed or scrambled peptides
Add_taxid=True         # add taxonomy id to header id, only write id of header, not description (smaller output files, needed for CHEW)

rm_prep=False #remove database after prepping


### can also be used from CLI:
    #Example syntax: 4_Download_unclustered_databse -DB "Swiss-Prot" -make_dmnd 0 -prepb 0

#%%


parser = argparse.ArgumentParser(description="Download database Input arguments",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()

### Define keyword dictionary 
cvars=set([str(k)+":#|%"+str(v)  for k,v in locals().copy().items() if k!="svars"])
kws={i.split(":#|%")[0]:i.split(":#|%")[1]   for i in list(cvars-svars)}
for k in kws.keys(): parser.add_argument("-"+k)
args = {k:v for k,v in vars(parser.parse_args()).items() if v is not None}
kws.update(args)
locals().update(kws)


ranks=["superkingdom","phylum","class","order","family","genus","species"] 
Ambiguous_AAs=["B","O","U","X","Z","[","("]
decoy_delimiter="decoy_"
decoy_method="reverse" #or "scramble"



#%%

def make_diamond_database(input_file,make_dmnd=make_dmnd):
    
    if make_dmnd:
        output_path=str(Path(output_folder,Path(input_file).stem))
        command="cd "+'"'+basedir +'"'+ " && "
        command+="diamond makedb --in "+'"'+input_file+'"' + " -d "+'"'+output_path+'"'
        print(command)
        stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        
        return output_path+".dmnd"

def is_fasta(input_file):
    fasta=SeqIO.parse(input_file,"fasta")
    return any(fasta)

def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g

def prep_db(Path_to_db,Ambiguous_AAs=Ambiguous_AAs):
    
    
    output_paths=[]
    if prepdb:
        
        if os.path.isdir(Path_to_db):
            input_paths=[str(Path(Path_to_db,i)) for i in os.listdir(Path_to_db)]
        else:
            input_paths=[input_paths]
         
        #parse output_path
        for input_path in input_paths:
            
            if is_fasta(input_path):
                
                Output_path=input_path
                if Bacterial_only:   Output_path=Output_path.replace(".fasta","_BacArch.fasta")
                if Remove_ambiguous: Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
                if No_Dump:          Output_path=Output_path.replace(".fasta","_NoDump.fasta")
                if No_Fragments:     Output_path=Output_path.replace(".fasta","_NoFrag.fasta")
                if Equate_IL:        Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
                if Add_decoy:        Output_path=Output_path.replace(".fasta","_Decoy.fasta")
                if Add_taxid:        Output_path=Output_path.replace(".fasta","_taxid.fasta")
        
                output_paths.append(Output_path)
                
                # tax database and files
                if Bacterial_only or No_Dump:
                    taxdf=pd.read_csv(Path_to_taxonomy,sep="\t")
                
                if Bacterial_only:
                    taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)
                
                if No_Dump:
                    taxdf=taxdf[taxdf["Dump_taxid"].astype(str)=="False"] 
                
                if not Equate_IL: Ambiguous_AAs+=["J"]
                
                once=True
                Taxid_delimiter="GTDB"
                recs=SeqIO.parse(Path_to_db,format="fasta")
                chunks=chunk_gen(recs)
                
                #write IL datbase
                print("writing "+Path(Output_path).stem)
                with open(Output_path,"w+") as f:
                
                    for ic,c in enumerate(chunks):
                        print("chunk "+ str(ic))
                
        
                        chunk_df=pd.DataFrame([[r.id,str(r.seq),r.description] for r in c],columns=["id","seq","description"])
                        
                        #Check taxid delimiter
                        if once:
                            if chunk_df.head(10).description.str.contains("OX=").any():     Taxid_delimiter="OX="
                            if chunk_df.head(10).description.str.contains("TaxID=").any():  Taxid_delimiter="TaxID=" #uniref style
                            once=False
                        
                        
                        if Bacterial_only or No_Dump: chunk_df=chunk_df[chunk_df.description.str.split(Taxid_delimiter).apply(lambda x:x[-1]).str.split(" ").apply(lambda x: x[0]).isin(taxdf["OX"])]
                
                        if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
                        if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Ambiguous_AAs],axis=1).any(axis=1)]
                        if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)]
                        if No_Dump:          chunk_df=chunk_df[chunk_df.description.str.split(Taxid_delimiter).apply(lambda x:x[-1]).str.split(" ").apply(lambda x: x[0]).isin(taxdf["OX"])]
                        
                        
                        if Add_taxid:        
                            if Taxid_delimiter=="GTDB": #only used for renaming GTDB folders 
                                chunk_df["id"]=chunk_df["id"]+"|"+Path(Path_to_db).parents[0].name
                            else:
                                chunk_df["id"]=chunk_df["id"]+"|"+chunk_df.description.str.split(Taxid_delimiter).apply(lambda x: x[-1]).str.split(" ").apply(lambda x:x[0])
                    
                
                
                        if Add_decoy:
                            decoy=chunk_df.copy()
                            if decoy_method=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
                            if decoy_method=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
                            decoy["id"]=decoy_delimiter+decoy["id"]
                            chunk_df=pd.concat([chunk_df,decoy])
                            
                        f.write("\n"+"\n".join(">"+chunk_df["id"]+"\n"+chunk_df["seq"]))
                        
                if rm_prep:
                    shutil.rmtree(input_path)
            
        return output_paths    
            
            

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
        

           


def merge(outpath,files,rm=rm_merge):
    with open(outpath,'wb') as wfd:
        for f in files:
            print(f)
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
                
            if rm:
                shutil.rmtree(f)
    return outpath

def download_extract(urls,path):
    
    download(urls,path)
    extract_subfolders(path)


#%% UniprotKB based

#Swissprot
if DB=="Swiss-Prot" or DB=="UniProtKB":
    url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    path=str(Path(basedir,"Swiss-Prot"))
    download_extract(url,path)
    if DB=="Swiss-Prot" :
        db=prep_db(path)


#Trembl
if DB=="TrEMBL" or DB=="UniProtKB":
    url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
    path=str(Path(basedir,"TrEMBL"))
    download_extract(url,path)
    if DB=="TrEMBL":
        db=prep_db(path)

if DB=="UniProtKB":

    #Uniprot KB
    dirs=["Swiss-Prot","TrEMBL"]
    path="UniProtKB"
    outfile="UniProtKB.fasta"
    outpath=str(Path(basedir,path,outfile))
    if not os.path.exists(path): os.mkdir(path)
    files=[]
    for d in dirs:
        [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]
    merged=merge(outpath,files)
    db=prep_db(merged)
    

if DB=="UniRef100": #Uniref 100
    url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
    path=str(Path(basedir,"UniRef100"))
    download_extract(url,path)
    db=prep_db(path)
    

if DB=="UniRef90": #Uniref 90
    url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
    path=str(Path(basedir,"UniRef90"))
    download_extract(url,path)
    db=prep_db(path)
    

if DB=="UniRef50": #Uniref 50
    url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
    path=str(Path(basedir,"UniRef50"))
    download_extract(url,path)
    db=prep_db(path)
    

#%% NCBI



if DB=="RefSeq": #Uniref 90 #Refseq (protein and nonredundant protein files)
    folder="RefSeq_protein_db"
    ftpbase='ftp.ncbi.nlm.nih.gov'
    ftpdir='/refseq/release/complete/'
    host = ftputil.FTPHost(ftpbase, 'anonymous', 'password')
    host.chdir(ftpdir)
    dir_list = host.listdir(host.curdir)
    path=str(Path(basedir,"RefSeq"))
    for link in dir_list:  
        if link.endswith(".faa.gz"): #all
        #if link.endswith(".faa.gz") and Path(link).stem.startswith("complete.nonredundant_protein") : #only nonredundant
            
            print(link)
            url="https://"+ftpbase+ftpdir+link
            download_extract(url,path)
            
    
    #Refseq 
    dirs=["Refseq"]
    files=[]
    for d in dirs:
        [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]

    path="Refseq_Protein"
    outfile="Refseq_Protein.fasta"
    outpath=str(Path(basedir,path,outfile))
    if not os.path.exists(path): os.mkdir(path)
    merged=merge(outpath,files)
    db=prep_db(merged)
    

if DB=="NCBI_NR":
    url="https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
    path=str(Path(basedir,"NCBI_NR"))
    download_extract(url,path)
    db=prep_db(path)
    
    

#%% GTDB

if DB=="GTDB":

    url="https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"
    
    path=str(Path(basedir,"GTDB"))
    download_extract(url,path)
    
    
    
    # #GTDB (r207)
    dirs=[str(Path("GTDB","renamed"))]
    path="GTDB_merged"
    outfile="GTDB_merged.fasta"
    outpath=str(Path(basedir,path,outfile))
    if not os.path.exists(path): os.mkdir(path)
    files=[]
    for d in dirs:
        [files.append(str(Path(basedir,d,i))) for i in os.listdir(str(Path(basedir,d)))]
        
    files=sum([prep_db(file) for file in files],[])
    db=merge(outpath,files)
    
    
#%% Make diamond db

make_diamond_database(db)


