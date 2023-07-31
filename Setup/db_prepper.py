# -*- coding: utf-8 -*-


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:34:40 2021

@author: hugokleikamp
"""

#%% clear variables and console, stor current variables

try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass


#% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd())) #change base directory 
os.chdir(basedir)
print(os.getcwd())

#% import 
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
from load_vars import *

svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])
#% parameters


Path_to_db=""
Path_to_taxonomy=str(Path(basedir,"parsed_taxonomy.tsv"))

BacArch_only=False     # retain only Bacteria (and Archaea(!)) in database  
Equate_IL=True         # change I and J into L 
Remove_ambiguous=True  # remove ambiguous amino acids "B","X","Z","[","(" , and J in case IL is not equated
No_Fragments=False     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
No_Dump=False           # remove dump taxa (unspecific filler names of NCBI with bloated annotations, like "uncultured" or "bacterium")
No_Desc=False          # only write id(accession) to fasta header instead of description (saves space)
Add_decoy=False        # append decoy of reversed or scrambled peptides
Add_taxid=True         # add taxonomy id to header id, only write id of header, not description (smaller output files, needed for CHEW)
Taxid_delimiter=""     #custom taxid delimiter 
output_folder=""
#rm_prep=False #remove database after prepping


    





### Define keyword dictionary 

cvars=set([str(k)+":#|%"+str(v)  for k,v in locals().copy().items() if k!="svars"])# "here it makes it into strs!"
kws=parse_kws(cvars,svars) #use this everywhere

#kws.update(args)
locals().update(kws)

#print(kws)

ranks=["superkingdom","phylum","class","order","family","genus","species"] 
Ambiguous_AAs=["B","O","U","X","Z","[","("]
decoy_delimiter="decoy_"
decoy_method="reverse" #or "scramble"



#%%

def is_fasta(input_file):
    fasta=SeqIO.parse(input_file,"fasta")
    return any(fasta)

def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g

#%%
def prep_db(Path_to_db,
            
            Taxid_delimiter=Taxid_delimiter,
            Ambiguous_AAs=Ambiguous_AAs,
            Path_to_taxonomy=Path_to_taxonomy,
            BacArch_only=BacArch_only,      # retain only Bacteria (and Archaea(!)) in database  
            Equate_IL=Equate_IL,        # change I and J into L 
            Remove_ambiguous=Remove_ambiguous,  # remove ambiguous amino acids "B","X","Z","[","(" , and J in case IL is not equated
            No_Fragments=No_Fragments,     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
            No_Dump=No_Dump,         # remove dump taxa (unspecific filler names of NCBI with bloated annotations, like "uncultured" or "bacterium")
            No_Desc=No_Desc,           # only write id(accession) to fasta header instead of description (saves space)
            Add_decoy=Add_decoy,        # append decoy of reversed or scrambled peptides
            Add_taxid=Add_taxid         # add taxonomy id to header id, only write id of header, not description (smaller output files, needed for CHEW)

            
            
            ):
    
    output_paths=[]

    if os.path.isdir(Path_to_db):
        Path_to_db=[str(Path(Path_to_db,i)) for i in os.listdir(Path_to_db)]
    else:
        Path_to_db=Path_to_db.split()
     


   
    # tax database and files
    if BacArch_only or No_Dump:
        taxdf=pd.read_csv(Path_to_taxonomy,sep="\t",dtype=str)


    if BacArch_only:
        taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)
    
    if No_Dump:
        taxdf=taxdf[taxdf["Dump_taxid"].astype(str)=="False"] 
    
    output_paths=[]
   
    if type(Path_to_db)==str:
        if os.path.isdir(Path_to_db):
            Path_to_db=[str(Path(Path_to_db,i)) for i in os.listdir(Path_to_db)]
        else:
            Path_to_db=Path_to_db.split()
     
    #parse output_path
    for input_path in Path_to_db:
        
        if is_fasta(input_path):
            
            
            Output_path=str(Path(Path(input_path).parents[0],Path(input_path).stem))
            if len(output_folder):
                if os.isdir(output_folder):
                   
                    if not os.path.exists(output_folder): 
                        os.makedir(output_folder)
                
                    Output_path=str(Path(output_folder,Path(input_path).stem))
            
            
            
            suf=Path(input_path).suffix
            Output_path=input_path
            if BacArch_only:   Output_path=Output_path.replace(suf,"_BacArch"+suf)
            if Remove_ambiguous: Output_path=Output_path.replace(suf,"_NoAmb"+suf)
            if No_Dump:          Output_path=Output_path.replace(suf,"_NoDump"+suf)
            if No_Fragments:     Output_path=Output_path.replace(suf,"_NoFrag"+suf)
            if Equate_IL:        Output_path=Output_path.replace(suf,"_IJeqL"+suf)
            if Add_decoy:        Output_path=Output_path.replace(suf,"_Decoy"+suf)
            if Add_taxid:        Output_path=Output_path.replace(suf,"_taxid"+suf)
    

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
                            chunk_df["description"]=chunk_df["id"]+"|"+chunk_df["OX"]+" "+chunk_df["description"].str.split(" ",n=1).apply(lambda x: x[-1])
            
            
                    if Add_decoy:
                        decoy=chunk_df.copy()
                        if decoy_method=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
                        if decoy_method=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
                        decoy["id"]=decoy_delimiter+decoy["id"]
                        decoy["description"]=decoy_delimiter+decoy["description"]
                        chunk_df=pd.concat([chunk_df,decoy])
                    
    
                    
                    ##### this only writes id!!!!
                    if No_Desc:
                        f.write("\n"+"\n".join(">"+chunk_df["id"]+"\n"+chunk_df["seq"])+"\n")
                    else:  
                        f.write("\n"+"\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"])+"\n")
                
            # if rm_prep:
            #     shutil.rmtree(input_path)

    return output_paths    
        

prep_db(Path_to_db)

#%%


