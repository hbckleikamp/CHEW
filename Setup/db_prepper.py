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


Path_to_db=""
Path_to_taxonomy=str(Path(basedir,"parsed_taxonomy.tsv"))

prepdb=True #after downloading prep DB with following arguments

Bacterial_only=True    # retain only Bacteria (and Archaea(!)) in database  
Equate_IL=True         # change I and J into L 
Remove_ambiguous=True  # remove ambiguous amino acids "B","X","Z","[","(" , and J in case IL is not equated
No_Fragments=False     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
No_Dump=True           # remove dump taxa (unspecific filler names of NCBI with bloated annotations, like "uncultured" or "bacterium")
Add_decoy=False        # append decoy of reversed or scrambled peptides
Add_taxid=True         # add taxonomy id to header id, only write id of header, not description (smaller output files, needed for CHEW)

rm_prep=False #remove database after prepping

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
            Path_to_db=[str(Path(Path_to_db,i)) for i in os.listdir(Path_to_db)]
        else:
            Path_to_db=Path_to_db.split()
         
        #parse output_path
        for input_path in Path_to_db:
           
            if is_fasta(input_path):
                
                suf=Path(input_path).suffix
                Output_path=input_path
                if Bacterial_only:   Output_path=Output_path.replace(suf,"_BacArch"+suf)
                if Remove_ambiguous: Output_path=Output_path.replace(suf,"_NoAmb"+suf)
                if No_Dump:          Output_path=Output_path.replace(suf,"_NoDump"+suf)
                if No_Fragments:     Output_path=Output_path.replace(suf,"_NoFrag"+suf)
                if Equate_IL:        Output_path=Output_path.replace(suf,"_IJeqL"+suf)
                if Add_decoy:        Output_path=Output_path.replace(suf,"_Decoy"+suf)
                if Add_taxid:        Output_path=Output_path.replace(suf,"_taxid"+suf)
        
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
                recs=SeqIO.parse(input_path,format="fasta")
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
            

prep_db(Path_to_db)