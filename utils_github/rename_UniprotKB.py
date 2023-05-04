# -*- coding: utf-8 -*-
"""
Created on Mon May  1 10:53:47 2023

@author: ZR48SA
"""


#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()



#%% Functions

def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g
        
#%%

#add something to add "CRAP"

#for Bacterial_only

Path_to_db="H:/Databases/UniprotKB/UniprotKB.fasta"
Path_to_taxonomy="C:/MP-CHEW/CHEW/parsed_ncbi_taxonomy.tsv" #unused
Taxid_delimiter="OX=" #unused

#for Remove_ambiguous
Ambiguous_AAs=["B","X","Z","[","("] #"O","U","*"

#Add decoy
decoy_delimiter="" #unused
decoy_method="" #unused


import Bio
from Bio import SeqIO
import pandas as pd
import itertools
import random
import subprocess

# Options
Bacterial_only=True    # retain only Bacteria (and Archaea(!)) in database  
Equate_IL=True         # change I and J into L 
Remove_ambiguous=True  # remove ambiguous amino acids "B","X","Z","[","(" , and J in case IL is not equated
No_Fragments=False     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
No_Dump=True           # remove dump taxa (unspecific filler names of NCBI with bloated annotations, like "uncultured" or "bacterium")
Add_decoy=False        # append decoy of reversed or scrambled peptides
Add_taxid=True         # append decoy of reversed or scrambled peptides


decoy_method="reverse" #or "scrambled

# Functions

def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g

#parse output_path
Output_path=Path_to_db
if Bacterial_only:   Output_path=Output_path.replace(".fasta","_BacArch.fasta")
if Remove_ambiguous: Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
if No_Dump:          Output_path=Output_path.replace(".fasta","_NoDump.fasta")
if No_Fragments:     Output_path=Output_path.replace(".fasta","_NoFrag.fasta")
if Equate_IL:        Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
if Add_decoy:        Output_path=Output_path.replace(".fasta","_Decoy.fasta")
if Add_taxid:         Output_path=Output_path.replace(".fasta","_taxid.fasta")


# tax database and files

if Bacterial_only or No_Dump:
    ranks=["superkingdom","phylum","class","order","family","genus","species"] 
    taxdf=pd.read_csv(Path_to_taxonomy,sep="\t")

if Bacterial_only:
    taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)

if No_Dump:
    taxdf=taxdf[taxdf["Dump_taxid"].astype(str)=="False"] 


#%%

if not Equate_IL: Ambiguous_AAs+=["J"]

#read    
recs=SeqIO.parse(Path_to_db,format="fasta")
chunks=chunk_gen(recs)
#write IL datbase
print("writing "+Path(Output_path).stem)
with open(Output_path,"w+") as f:

    for ic,c in enumerate(chunks):
        print("chunk "+ str(ic))

        #chunk_df=pd.DataFrame([[str(r.seq),r.description] for r in c],columns=["seq","description"])
        chunk_df=pd.DataFrame([[r.id,str(r.seq),r.description] for r in c],columns=["id","seq","description"])
        
        
        
        if Bacterial_only or No_Dump: chunk_df=chunk_df[chunk_df.description.str.split(Taxid_delimiter).apply(lambda x:x[-1]).str.split(" ").apply(lambda x: x[0]).isin(taxdf["OX"])]

        if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
        if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Ambiguous_AAs],axis=1).any(axis=1)]
        if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)]
        if No_Dump:          chunk_df=chunk_df[chunk_df.description.str.split(Taxid_delimiter).apply(lambda x:x[-1]).str.split(" ").apply(lambda x: x[0]).isin(taxdf["OX"])]
        if Add_taxid:        chunk_df["id"]=chunk_df["id"]+"|"+chunk_df.description.str.split(Taxid_delimiter).apply(lambda x: x[-1]).str.split(" ").apply(lambda x:x[0])
    


        if Add_decoy:
            decoy=chunk_df.copy()
            if decoy_method=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
            if decoy_method=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
            decoy["id"]=decoy_delimiter+decoy["id"]
            chunk_df=pd.concat([chunk_df,decoy])
            
        f.write("\n"+"\n".join(">"+chunk_df["id"]+"\n"+chunk_df["seq"]))

#%% construct diamond database

diamond_path="C:/MP-CHEW/CHEW/"

command="cd" +' "'+diamond_path+'" && '
#./diamond makedb --in reference.fasta -d reference
command+="diamond makedb --in "+str(Path(Output_path)) + " -d "+Path(Output_path).stem
print(command)
stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

print(stderr)