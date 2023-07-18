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
basedir=str(Path(os.getcwd()))#.parents[0]) #change base directory to HybridCycler
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
import argparse


svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])
#%% parameters



DB="RefSeq" #UniProt, Swiss-Prot, TrEMBL, UniRef50, UniRef90, UniRef100, GTDB, NCBI_NR 
rm_merge=True               # remove database after merging 
output_folder=basedir       # where to put databases
diamond_output_folder=""

Path_to_taxonomy=str(Path(basedir,"parsed_taxonomy.tsv"))
diamond_path=str(Path(Path(basedir).parents[0],"diamond"))

prepdb=True #after downloading prep DB with following arguments
make_dmnd=True


BacArch_only=True      # retain only Bacteria (and Archaea(!)) in database  
Equate_IL=True         # change I and J into L 
Remove_ambiguous=True  # remove ambiguous amino acids "B","X","Z","[","(" , and J in case IL is not equated
No_Fragments=False     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
No_Dump=True           # remove dump taxa (unspecific filler names of NCBI with bloated annotations, like "uncultured" or "bacterium")
No_Desc=True           # only write id(accession) to fasta header instead of description (saves space)
Add_decoy=False        # append decoy of reversed or scrambled peptides
Add_taxid=True         # add taxonomy id to header id, only write id of header, not description (smaller output files, needed for CHEW)
Taxid_delimiter=""     #custom taxid delimiter 

rm_prep=False #remove database after prepping
rm_merge=False #remove unmerged database after merging

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

if not os.path.exists(output_folder): os.makedirs(output_folder) 

#%%

def make_diamond_database(input_files,make_dmnd=make_dmnd,output_folder=diamond_output_folder):
    
    
    if make_dmnd:

        if type(input_files)==str:
            if os.path.isdir(input_files):
                input_files=[str(Path(input_files,i)) for i in os.listdir(input_files)]
            else:
                input_files=input_files.split()
         
        for input_file in input_files:
    
            if not len(output_folder):
                
                output_folder=Path(input_file).parents[0]
            output_path=str(Path(output_folder,Path(input_file).stem))
   
            
            if not os.path.exists(output_folder): os.makedirs(output_folder) #check if file is already downloaded
    
            output_path=str(Path(output_folder,Path(input_file).stem))
            command='"'+diamond_path+'"'+" makedb --in "+'"'+input_file+'"' + " -d "+'"'+output_path+'"'
            print(command)
            stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

def is_fasta(input_file):
    fasta=SeqIO.parse(input_file,"fasta")
    return any(fasta)

def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g

def prep_db(Path_to_db,Ambiguous_AAs=Ambiguous_AAs,Taxid_delimiter=Taxid_delimiter):


    Path_to_db="H:/Databases/Swiss-Prot/Swiss-Prot/uniprot_sprot.fasta"

    # tax database and files
    if BacArch_only or No_Dump:
        taxdf=pd.read_csv(Path_to_taxonomy,sep="\t")


    if BacArch_only:
        taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)
    
    if No_Dump:
        taxdf=taxdf[taxdf["Dump_taxid"].astype(str)=="False"] 
    
    output_paths=[]
    if prepdb:
        if type(Path_to_db)==str:
            if os.path.isdir(Path_to_db):
                Path_to_db=[str(Path(Path_to_db,i)) for i in os.listdir(Path_to_db)]
            else:
                Path_to_db=Path_to_db.split()
         
        #parse output_path
        for input_path in Path_to_db:
            
            if is_fasta(input_path):
                
                Output_path=str(Path(Path(input_path).parents[0],Path(input_path).stem))
                
                if BacArch_only:     Output_path+="_BacArch" 
                if Remove_ambiguous: Output_path+="_NoAmb"
                if No_Dump:          Output_path+="_NoDump"
                if No_Desc:          Output_path+="_NoDesc"
                if No_Fragments:     Output_path+="_NoFrag"
                if Equate_IL:        Output_path+="_IJeqL"
                if Add_decoy:        Output_path+="_Decoy"
                if Add_taxid:        Output_path+="_taxid"
                Output_path+=".fa"
        
                output_paths.append(Output_path)
                

        
                if not Equate_IL: Ambiguous_AAs=Ambiguous_AAs+["J"]
                
                once=True
                
                recs=SeqIO.parse(input_path ,format="fasta")
                chunks=chunk_gen(recs)
                
                #write IL datbase
                print("writing "+Path(Output_path).stem)
                with open(Output_path,"w+") as f:
                
                    for ic,c in enumerate(chunks):
                        print("chunk "+ str(ic))
                
                        chunk_df=pd.DataFrame([[r.id,str(r.seq),r.description] for r in c],columns=["id","seq","description"])
                        
                        #Check taxid delimiter
                        if not Taxid_delimiter:
                            if once:
                                if chunk_df.head(10).description.str.contains("OX=").any():     Taxid_delimiter="OX="
                                elif chunk_df.head(10).description.str.contains("TaxID=").any():  Taxid_delimiter="TaxID=" #uniref style
                                elif chunk_df.head(10).description.str.contains("[",regex=False).any():  Taxid_delimiter="RefSeq" #uniref style
                                else: Taxid_delimiter="GTDB"
                                once=False
                        
                        if Add_taxid or BacArch_only or No_Dump:
                            
                            if Taxid_delimiter=="GTDB":
                                chunk_df["OX"]=Path(Path_to_db).parents[0].name
                            elif Taxid_delimiter=="RefSeq":
                
                                chunk_df["OX"]=taxdf.merge(chunk_df.description.str.split("[",regex=False).apply(lambda x: x[-1]).str.strip(" ]").rename("OS"),on="OS",how="right")["OX"].fillna("")
                            else:
                                chunk_df["OX"]=chunk_df.description.str.split(Taxid_delimiter).apply(lambda x:x[-1]).str.split(" ").apply(lambda x: x[0])
                            
                    
                        if BacArch_only or No_Dump: 
                            chunk_df=chunk_df[chunk_df["OX"].isin(taxdf["OX"])]
                
                        if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
                        if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Ambiguous_AAs],axis=1).any(axis=1)]
                        if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)]
          
                        if Add_taxid:        
                            if No_Desc:
                                chunk_df["id"]=chunk_df["id"]+"|"+chunk_df["OX"]
                    
                            else:
                                chunk_df["description"]=chunk_df["id"]+"|"+chunk_df["OX"]+chunk_df["description"].str.split(" ",1).apply(lambda x: x[2])
                
                
                        if Add_decoy:
                            decoy=chunk_df.copy()
                            if decoy_method=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
                            if decoy_method=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
                            decoy["id"]=decoy_delimiter+decoy["id"]
                            chunk_df=pd.concat([chunk_df,decoy])
                        
        
                        
                        ##### this only writes id!!!!
                        if No_Desc:
                            f.write("\n"+"\n".join(">"+chunk_df["id"]+"\n"+chunk_df["seq"])+"\n")
                        else:  
                            f.write("\n"+"\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"])+"\n")
                   
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
    
    if os.path.isdir(files):
        files=[str(Path(files,i)) for i in os.listdir(files)]
     
    with open(outpath,'wb') as wfd:
        for f in files:
            #if is_fasta(f):
            print(f)
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
                
            if rm:
                shutil.rmtree(f)
    return outpath

def download_extract(urls,path):
    
    download(urls,path)
    extract_subfolders(path)


# #%% Testing
# import numpy as np
# from collections import Counter
# tdf=pd.read_excel("C:/MP-CHEW/Datasets/PXD005776_MIX24/mix24names.xlsx",engine="openpyxl")

# taxa=tdf["NCBI OX"].astype(str).tolist()# tdf["Mix24 species name"].tolist()+tdf["NCBI synonym"].dropna().tolist()
# taxcounts=np.array([0]*len(taxa))

# tax_dict=dict()
# for t in taxa:
#     tax_dict.update({t:0})


# input_path="H:/Databases/Refseq_NR/RefSeq_merged_BacArch_NoAmb_NoDump_NoDesc_IJeqL_taxid.fa" #"H:/Databases/Refseq_NR/RefSeq_merged.fasta"
# recs=SeqIO.parse(input_path ,format="fasta")
# chunks=chunk_gen(recs)

# #write IL datbase



# for ic,c in enumerate(chunks):
#     print("chunk "+ str(ic))

#     ids=pd.Series([r.id for r in c])
#     ts=ids.str.split("|").apply(lambda x: x[-1])
#     counts=Counter(ts[ts.isin(taxa)])

#     for k,v in counts.items():
#         tax_dict.update({k:tax_dict.get(k)+v})

                                
                                

#     #for r in c: 
#         # for ixt,t in enumerate(taxa):
#         #     if t in r.description:
#         #         taxcounts[ixt]+=1



#%% UniprotKB based

#Swissprot
if DB=="Swiss-Prot" or DB=="UniProtKB":
    url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    path=str(Path(output_folder,"Swiss-Prot"))
    download_extract(url,path)
    if DB=="Swiss-Prot" :
        db=prep_db(path)


#Trembl
if DB=="TrEMBL" or DB=="UniProtKB":
    url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
    path=str(Path(output_folder,"TrEMBL"))
    download_extract(url,path)
    if DB=="TrEMBL":
        db=prep_db(path)

if DB=="UniProtKB":

    #Uniprot KB
    dirs=["Swiss-Prot","TrEMBL"]
    path="UniProtKB"
    outfile="UniProtKB.fasta"
    outpath=str(Path(output_folder,path,outfile))
    if not os.path.exists(path): os.mkdir(path)
    files=[]
    for d in dirs:
        [files.append(str(Path(output_folder,d,i))) for i in os.listdir(str(Path(output_folder,d)))]
    merged=merge(outpath,files)
    db=prep_db(merged)
    

if DB=="UniRef100": #Uniref 100
    url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
    path=str(Path(output_folder,"UniRef100"))
    download_extract(url,path)
    db=prep_db(path)
    

if DB=="UniRef90": #Uniref 90
    url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
    path=str(Path(output_folder,"UniRef90"))
    download_extract(url,path)
    db=prep_db(path)
    

if DB=="UniRef50": #Uniref 50
    url="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
    path=str(Path(output_folder,"UniRef50"))
    download_extract(url,path)
    db=prep_db(path)
    

#%% NCBI



if DB=="RefSeq": #Uniref 90 #Refseq (protein and nonredundant protein files)

    ftpbase='ftp.ncbi.nlm.nih.gov'
    ftpdir='/refseq/release/complete/'
    host = ftputil.FTPHost(ftpbase, 'anonymous', 'password')
    host.chdir(ftpdir)
    dir_list = host.listdir(host.curdir)
    path=str(Path(output_folder,"RefSeq"))
    for link in dir_list:  
        #if link.endswith(".faa.gz"): #all
        if link.endswith(".faa.gz") and Path(link).stem.startswith("complete.nonredundant_protein") : #only nonredundant
            
            print(link)
            url="https://"+ftpbase+ftpdir+link
            download_extract(url,path)
            

    outpath=str(Path(Path(path).parents[0],"RefSeq_merged.fasta"))
    merged=merge(outpath,files=path)
    

    db=prep_db(merged)
    

if DB=="NCBI_NR":
    url="https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
    path=str(Path(output_folder,"NCBI_NR"))
    download_extract(url,path)
    db=prep_db(path)
    
    

#%% GTDB

if DB=="GTDB":

    url="https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_proteins_aa_reps.tar.gz"
    
    path=str(Path(output_folder,"GTDB"))
    download_extract(url,path)
    
    
    
    # #GTDB (r207)
    dirs=[str(Path("GTDB","renamed"))]
    path="GTDB_merged"
    outfile="GTDB_merged.fasta"
    outpath=str(Path(output_folder,path,outfile))
    if not os.path.exists(path): os.mkdir(path)
    files=[]
    for d in dirs:
        [files.append(str(Path(output_folder,d,i))) for i in os.listdir(str(Path(output_folder,d)))]
        
    files=sum([prep_db(file) for file in files],[])
    db=merge(outpath,files)
    
    
#%% Make diamond db
make_diamond_database(db)


