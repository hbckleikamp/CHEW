# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:24:45 2022

@author: ZR48SA
"""



#%% set base path
from pathlib import Path
import os
from inspect import getsourcefile
# change directory to script directory (should work on windows and mac)
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())
basedir=os.getcwd()


#Modules


#from SMSNet_final_database_search_HK import *
import Bio
from Bio import SeqIO

import subprocess
import pandas as pd
import psutil
import shutil
import numpy as np
import re
import math
import random
import time
from collections import Counter
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

import matplotlib.pyplot as plt
import itertools
#Config


############ Function Config #############

output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq","full_sseq"]
ranks=np.array(["superkingdom","phylum","class","order","family","genus","species"]) 
ncbi_taxdf=pd.read_csv("parsed_ncbi_taxonomy.tsv",sep="\t")

#Update metadata paths!
metadata_filepaths=["C:/MultiNovo/GTDB/Metadata/bac120_taxonomy.tsv",
"C:/MultiNovo/GTDB/Metadata/ar53_taxonomy.tsv"]
gtdb_taxdf=pd.concat([pd.read_csv(f,sep="\t",header=None,names=["OX","taxonomy"]) for f in metadata_filepaths]).set_index("OX")
gtdb_taxdf[ranks]=gtdb_taxdf["taxonomy"].str.rsplit(";",expand=True)
gtdb_taxdf=gtdb_taxdf[ranks]



### Paths ###
Temporary_directory  =basedir  #Writing directory for temporary indices of MSFragger and Diamond, make sure there is enough space here
Output_directory     =basedir  #Directory where all generated outputs are written to





############ Functions #############



def raw2mzML(raw_file): #

    output_folder=str(Path(Output_directory,"HB_mzML"))
    if not os.path.exists(output_folder): os.mkdir(output_folder)

    command="cd" +' "'+output_folder+'" && msconvert '
    command+='"'+raw_file+'"' 
    command+=' --mzML --filter "peakPicking vendor" --filter "zeroSamples removeExtra" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>"'
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    output_file=str(Path(output_folder,Path(raw_file).stem+".mzML"))
    return output_file

def read_pin(pinfile):
    with open(pinfile, "r") as f:
        lines=f.readlines()    
    header=lines[0].replace("\n","").split("\t")
    return pd.DataFrame([i.replace("\n","").split("\t",len(header)-1) for i in lines[1:]],columns=header)


def MSFragger_annotation(input_files,   #full_path, .mzML
                         database_path, #full_path, .fasta
                         output_folder,
                         params_path=str(Path(basedir,"closed_fragger.params")),
                         max_no_hits=5, #total number of top hits retained after mergeing database splits
                         no_splits=None, 
                         no_batches=None
                         ):


    MSFragger_jar_path=str(Path(basedir,"MSFragger-3.5.jar"))  
    pep_split_path=str(Path(basedir,"msfragger_pep_split_HK.py"))
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))

    
    #rewrite closed_fragger.params according to database path
    with open(params_path,"r+") as f:
        lines=f.readlines()
        lines=["database_name = "+database_path+" #database name here\n" if line.startswith("database_name =") else line for line in lines]
        f.seek(0)
        f.writelines(lines)
        
        
    # to kill or not to kill (avoid a bug and RAM overload in file mergeing)
    if sum(["pepxml_pin" in i.split("#")[0]  for i in lines]):
       pep_split_path=str(Path(basedir,"msfragger_pep_split_HK_nokill.py"))     
   
    os.chdir(Temporary_directory) #change to tempdir for writing indices
    stderr=""
    retry=0
    while True: #rerun with different settings untill settings are found that fit RAM
        
        #remove old peptide split indices
        if os.path.exists(str(Path(Temporary_directory,"split_peptide_index_tempdir"))): shutil.rmtree(str(Path(Temporary_directory,"split_peptide_index_tempdir"))) 
        JavaMem=int(psutil.virtual_memory().available/10**9*0.8) #check available RAM
        
        if no_splits==None:
            no_splits=math.ceil(os.path.getsize(database_path)/10**9*150/JavaMem) #starting number of splits for db splitting, database size*150 is an estimated scaling factor, which will be effected by the number of allowed modifications and missed cleavages
        
        if no_batches==None:
            no_batches=int(psutil.disk_usage(Path(pep_split_path).anchor).free/10**9/os.path.getsize(database_path)) #number of files annotated per time
        batches=np.array_split(input_files, math.ceil(no_batches))

        print("running MSFragger, total number of splits: "+str(no_splits)+" , total number of batches: "+str(len(batches))+", Java Heap: "+str(JavaMem))
        
        
        retry+=1
        if retry>1:
            print("retrying")
            time.sleep(30) #wait 30 seconds for ram comsuption to reset before allocating new java heap
        
        if retry>8:
            print("retried 8 times, unknown error") #>8 ives you 2**8 splits (256) maximum when starting from 2 splits
            "i"+1 #hard exit
            
        command="cd" +' "'+basedir+'" && '
        command+=" && ".join(["".join([' python ',
                                  ' "'+pep_split_path+'" ',
                                  str(no_splits),
                                  ' "'+"java -jar -Xmx"+str(JavaMem)+"G"+'" ',
                                  ' "'+MSFragger_jar_path+'" ',
                                  ' "'+params_path+'" ',
                                
                                  #" MSFragger-3.5.jar closed_fragger.params ", #can also be written as full paths without the "cd" 
                                  '"'+'" "'.join(batch)+'"']) for batch in batches])
        
        print(command)
        stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        
        
        if "No space left on device" in str(stderr): #retry with smaller file batch
            print("out of disk space, retrying with smaller batch size")
            no_batches+=1
            if no_batches>len(input_files):
                print("unsolvable space error, exiting")
                "i"+1 #hard exit
            continue
                
        if "java.lang.OutOfMemoryError" in str(stderr) or "Not enough memory allocated to MSFragger." in str(stderr): #retry with more splits
            print("out of Java memory, retrying with more splits")
            no_splits=no_splits*2
            continue
        
        if "returned non-zero exit status" not in str(stderr):
            break

    os.chdir(basedir) #change back to basedir
    

    #this step can take up a lot of RAM, depending on the number of Database splits
    folders=[i for i in Path(basedir,Temporary_directory,"split_peptide_index_tempdir").glob("*") if i.is_dir()]
    outpaths=[]
    for mzML_file in input_files:
        
        outpath=str(Path(basedir,Output_directory,output_folder,Path(mzML_file).stem+".pin"))
        outpaths.append(outpath)
        
        if len(folders):
            #! Doesnt work for final search, since pin files are cleared up automatically 
            pinfile=Path(mzML_file).stem+".pin"
            pepdf=pd.concat([read_pin(str(Path(folder,pinfile))) for folder in folders])
            pepdf=pepdf.sort_values(by=["ScanNr","ExpMass","log10_evalue"]).groupby(["ScanNr","ExpMass"],sort=False).head(max_no_hits)
            pepdf["Proteins"]=pepdf.Proteins.str.replace("\t"," ")
            pepdf.to_csv(outpath,sep="\t")

        else:
            pin=str(Path(Path(mzML_file).parents[0],Path(mzML_file).stem))+".pin" #files are written to the location of the mzml files
            if os.path.exists(pin):
                shutil.move(pin,str(Path(Output_directory,output_folder,Path(mzML_file).stem+".pin")))
        
        #if pepXML files exist, move them
        pepxml=str(Path(Path(mzML_file).parents[0],Path(mzML_file).stem))+".pepXML" #files are written to the location of the mzml files
        if os.path.exists(pepxml):
            shutil.move(pepxml,str(Path(Output_directory,output_folder,Path(mzML_file).stem+".pepXML")))

    return outpaths




#Write MSfragger and PepNet files to fasta files fo Diamond alignment
def write_to_Diamond_fasta(MSFragger_files, #filepath or list of filepaths
                           output_folder,
                           max_evalue=10,
                           Top_score_fraction=0.9,
                           header_info=["alignment_Target_Decoy",
                                        "hyperscore"],
                           output_file=None, #if all are to be written to the same output file, put filename here, otherwise they are written to separate files
                           ):
       
    if not os.path.exists(str(Path(Output_directory ,output_folder))): os.mkdir(str(Path(Output_directory ,output_folder)))
    
    #remove output file to prevent overappending to existing file
    if output_file!=None:
        out_path=str(Path(Output_directory ,output_folder,Path(output_file).stem+".fa"))
        if os.path.exists(out_path): os.remove(out_path)
    
    output_paths=[]
    if type(MSFragger_files)!=type(list()):
        MSFragger_files=[MSFragger_files]
    MSFragger_files.sort()

    for im,MSFragger_file in enumerate(MSFragger_files):
        print("writing "+MSFragger_file+" to fasta") #debug
    
        #write outputs
        if output_file==None:
            out_path=str(Path(Output_directory ,output_folder,Path(MSFragger_file).stem+".fa"))
            
        if MSFragger_file.endswith(".pin"):

            pepdf=pd.read_csv(MSFragger_file,sep="\t")
            pepdf[["hyperscore","log10_evalue"]]=pepdf[["hyperscore","log10_evalue"]].astype(float)
            pepdf["evalue"]=10**pepdf.log10_evalue
            if max_evalue:
                pepdf=pepdf[pepdf["evalue"]<=max_evalue]
            if Top_score_fraction: #pick best scoring candidate
                pepdf=pepdf[(pepdf["evalue"]/pepdf.groupby(["ScanNr","ExpMass"])["evalue"].transform('min'))<=(1/Top_score_fraction)]
        
            pepdf["peptide_neighbours"]=pepdf["Peptide"].str.replace("c","").str.replace("n","").str.replace("-","").str.replace(".","").apply(lambda x: re.sub("[\[\[].*?[\]\]]", "", x).replace(",","")) #remove ptms in peptides
            pepdf.loc[:,"alignment_Target_Decoy"]="Target"
            pepdf.loc[pepdf["Proteins"].str.contains("decoy"),"alignment_Target_Decoy"]="Decoy"
        
    
            pepdf=pepdf[["peptide_neighbours"]+[i for i in header_info if i in pepdf.columns]].drop_duplicates().reset_index(drop=True)
            hdict=pepdf[[i for i in header_info if i in pepdf.columns]].fillna("").astype(str).T.to_dict()
            shdict=[str(hdict[ix]).replace(" ","") for ix in range(len(pepdf))]
            heads=">"+pepdf["peptide_neighbours"]+";"+shdict
            
            with open(out_path,"a") as f: #create target
                f.write("\n".join([heads[ix]+"\n"+peptide for ix,peptide in enumerate(pepdf["peptide_neighbours"])])+"\n")
            output_paths.append(out_path)       
        
    return list(set(output_paths))


def Diamond_alignment(input_file, #.fa
                      database_path,
                      output_folder,
                      output_columns=["qseqid","sseqid","stitle","bitscore"],
                      select=" -k25 ", #-top or  -k + integer (see diamond docs)
                      block_size=5,
                      index_chunks=1,
                      minimum_pident=80,
                      minimum_coverage=80,
                      minimum_bitscore=20,
                      tempdir=Temporary_directory,
                      Tricks=True,
                      evalue=None,
                      mode=None,
                      other_args=""
                      ):
    
    
    if type(input_file)==list:
        input_file=input_file[0]

    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    output_file=str(Path(Output_directory ,output_folder,Path(input_file).stem+".tsv"))

    command="cd "+'"'+basedir +'"'+ " && " + \
            "".join(['"'+str(Path(basedir,"diamond"))+'"',   
            " blastp -q " +'"'+input_file+'"',
            " -d "+'"'+database_path+'"',
            " -o "+'"'+output_file+'"',
            " -c" + str(index_chunks), 
            " --log ",
            " -b "+ str(block_size),  #Main parameter that determines ram consumption and performance
            " "+select+" ",
            " --id "+          str(minimum_pident),
            " --min-score "+   str(minimum_bitscore),
            " --query-cover  "+str(minimum_coverage),
            " -f 6 qseqid "+" ".join(output_columns)+" ",
            " -t "+'"'+tempdir+'"'+other_args])
 
    if Tricks:
        command+=" --algo ctg --dbsize 1 "  #little hacks 
        
    if evalue!=None:
        command+=" -e "+str(evalue)+" "
    
    if mode!=None:
        command+=mode
    
    print(command)
    
    #do this via a "bat" file because of diamond bug that does not want to work with custom matrices
    batfile=str(Path(basedir,"alignment.bat"))
    with open(batfile,"w") as bat:
        bat.write("#!/bin/bash"+"\n"+command)

    stdout, stderr =subprocess.Popen(batfile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    shutil.move(str(Path(basedir,"diamond.log")), str(Path(Output_directory ,output_folder,Path(input_file).stem+".log")))

    return output_file

            

#read alignment file to generator and filter on top x% scoring
def Diamond_alignment_Reader(input_file,
                              output_columns=["qseqid","sseqid","stitle","bitscore"],
                              score_cutoff=0.9,
                              Read_batch=1000000): 
    
    cdf=pd.read_csv(input_file, sep='\t', chunksize=Read_batch,names=output_columns) #read to generator
    sc=[] #dummy
    for ix,c in enumerate(cdf):
        print(ix)
        
        _,index=np.unique(c.qseqid,return_index=True)
        d=c.iloc[0:index.max()]
        if ix>0:
            d=pd.concat([sc,d])

        d.loc[:,"Sequence"]=d.qseqid.str.rsplit(";",expand=True).iloc[:,0]
        d=d[(d["bitscore"]/d.groupby("qseqid")["bitscore"].max()>=score_cutoff).tolist()]
        yield d

        sc=c.iloc[index.max():]
        
    #last one
    d=sc
    d.loc[:,"Sequence"]=d.qseqid.str.rsplit(";",expand=True).iloc[:,0]
    d=d[(d["bitscore"]/d.groupby("qseqid")["bitscore"].max()>=score_cutoff).tolist()]
    
    yield d
    



import more_itertools as mit
def fill_g(x):
    cons=[np.array(list(g)) for g in mit.consecutive_groups(np.argwhere(x==""))]
    if cons:
        for c in cons:   
            c=c.flatten().tolist()
            if c[-1]!=len(x)-1: 
                for i in c:
                    x[i]="gap_"+str(x[c[-1]+1])+"_"+ranks[i]
    return list(x)


def Process_diamond_alignment(Alignment,
                              output_folder,
                              output_columns=["qseqid","sseqid","stitle","bitscore"],
                              taxonomy="NCBI",
                              Precision_prefilter=0.4, #Target Decoy precision based denoising pre LCA filter
                              Frequency_prefilter=2,
                              weight_rank="species",
                              weight_cutoff=0.6):
          
     
    iterable=Diamond_alignment_Reader(Alignment,output_columns=output_columns)
    al=pd.concat([i for i in iterable])
     
    if taxonomy=="NCBI":
        #add NCBI taxonomy based on UniprotKB alignment headers
        al["OX"]=al["stitle"].str.split("OX=").apply(lambda x: x[-1]).str.split(" ").apply(lambda x: x[0]).astype(int)
        al["gene"]=al["stitle"].str.split("OX=").apply(lambda x: x[-1]).str.split(" ").apply(lambda x: x[0])
        al=al.merge(ncbi_taxdf,how="left",on="OX").fillna("") #use index instead of merge
    
        #fill gaps (only necessary for NCBI taxonomy)
        gaps=al[(al[ranks]=="").any(axis=1)]
        if len(gaps):
            u,ix,inv=np.unique(gaps[ranks].apply(";".join,axis=1),return_inverse=True,return_index=True)
            u=gaps.iloc[ix][ranks].values
            gaps[ranks]=np.array(list((map(fill_g,u))))[inv]
            al.loc[gaps.index,ranks]=gaps
    
    if taxonomy=="GTDB":
        al["OX"]=al.loc[:,"sseqid"].apply(lambda x: "_".join(x.split("_")[-3:])) 
        al=al.merge(gtdb_taxdf.reset_index(),how="left",on="OX").fillna("") #GTDB
    
    al=pd.concat([al, #add information stored in the header of the fasta file to the read dataframe
              pd.DataFrame.from_dict(al.qseqid.str.split(";").apply(lambda x: eval(x[1])).values.tolist(),orient="columns")],axis=1)
    
    ### Filtering ###
     
    
    al["hyperscore"]=al["hyperscore"].astype(float) #make this work for different input headers
    al["decoy"]=al.sseqid.str.startswith("decoy_")    
    target=al[~al["decoy"]]
    decoy= al[ al["decoy"]]
    
    # precision pre-filter
    mw=pd.concat([d.groupby(weight_rank).size() for d in [target,decoy]],axis=1).fillna(0).reset_index() 
    mw.columns=[weight_rank,"target_count","decoy_count"]
    mw["precision"]=mw["target_count"]/(mw["decoy_count"]+mw["target_count"])
    mw.to_csv(str(Path(Output_directory,output_folder,"precision.tsv")),sep="\t")
    
    taxa=mw.loc[(mw["precision"]>=Precision_prefilter) & (mw["target_count"]>=Frequency_prefilter),weight_rank]
    taxids=al.loc[al[weight_rank].isin(taxa),"OX"].unique().tolist()
    target=target[target["OX"].isin(taxids)]
    
    target=target.merge(mw[[weight_rank,"precision"]],on=weight_rank,how="left")
    target["corrected_score"]=target["bitscore"]*target["precision"]
    
    tlca=weighted_lca(target,weight_column="corrected_score",weight_cutoff=weight_cutoff)
    tlca_file=str(Path(Output_directory,output_folder,"lca.tsv"))
    tlca.to_csv(tlca_file,sep="\t")
     
    
    #ps=tlca.loc[tlca["genus"]!="","proteins"].reset_index(drop=True).str.split(", ") #quick and dirty, genus only
    ps=tlca["proteins"].reset_index(drop=True).str.split(", ")
    
    ps=ps.explode().unique().tolist()
    

    
    return ps,tlca_file
            



def weighted_lca(df, #dataframe with at least a column called Peptide, and rank
        group_on="Sequence",
        weight_column="bitscore",  
        protein_column="sseqid", #name of column containing proteins
        weight_cutoff=0.8   # minimum fraction of total weights 
        ):

    
    
    if weight_column not in df.columns:
         print("weight column not detected, weighing on frequency")
         df[weight_column]=1
        
    df[group_on]=df[group_on].astype(str)


    lin=[]
    for rank in ranks:
        wc=df.groupby([group_on,rank]).agg({weight_column:['sum']})/df.groupby(group_on).agg({weight_column:['sum']})
        wc.columns=wc.columns.droplevel(1)
        wc=wc.reset_index(rank)
        lin.append(wc[wc[weight_column]>=weight_cutoff][rank])
        
    lin=pd.concat(lin,axis=1)
    lcas=pd.DataFrame(df[group_on]).drop_duplicates().merge(lin,on=group_on).set_index(group_on)
    last=lcas.fillna(method="ffill",axis=1).iloc[:,-1]
    lcas["proteins"]=df[df[ranks].add(df[group_on],axis=0).isin(last.tolist()+last.index).any(axis=1)].groupby(group_on)[protein_column].apply(lambda x: ", ".join(list(set(x))))
    
    #add back proteins with no common ancestor concensus
    no_lca=df[~df[group_on].isin(lcas.index)]
    no_lca=pd.DataFrame(no_lca.groupby(group_on)[protein_column].apply(lambda x: ", ".join(x)))
    no_lca.columns=["proteins"]
    no_lca[ranks.tolist()]=[""]*len(ranks)
    lcas=pd.concat([lcas,no_lca],axis=0)
    
    
    return lcas.fillna("")





def make_diamond_database(Database, #.fa
                          output_folder,
                          output_file):
    
    
    

    output_path=str(Path(Output_directory,output_folder,Path(output_file).stem))
    command="cd "+'"'+basedir +'"'+ " && "
    command+="diamond makedb --in "+'"'+Database+'"' + " -d "+'"'+output_path+'"'
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


    return output_path+".dmnd"




from Bio import SeqIO
def Unique_fasta(input_file,chunk_size=10**5):
#source: https://stackoverflow.com/questions/66462611/remove-duplicated-sequences-in-fasta-with-python   

    output_file=input_file.replace(".fa","_unique.fa")
    chunks=chunk_gen(SeqIO.parse(input_file, "fasta"),size=chunk_size)
    seen = set()
    
    
    with open(output_file,"w") as f: #clears file if exists
        pass
    with open(output_file, "a") as f:
 
        for chunk in chunks:
        
            records = []
            for record in chunk:  
    
                if record.name not in seen:
    
                    seen.add(record.name)
                    records.append(record)
                    
            df=pd.DataFrame([[rec.description,str(rec.seq)] for rec in records],columns=["description","seq"])
            f.write("\n".join(">"+df.description+"\n"+df.seq)+"\n")

    return output_file, len(seen) 



import itertools
def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g
        







def Write_alignment_to_database(Alignment,output_folder):
    
    iterable=Diamond_alignment_Reader(Alignment,output_columns=["qseqid","sseqid","stitle","bitscore","full_sseq"],Read_batch=10**6)
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    file=str(Path(Output_directory,output_folder,"target.fa"))

    
    with open(file,"w") as f:
    
        for ix,batch in enumerate(iterable):
            batch.full_sseq=batch.full_sseq.str.replace("*","",regex=False) #remove protein ends
            f.write("\n".join((">"+batch.stitle+"\n"+batch.full_sseq))+"\n")

    return file

def slice_gen(gen,indices): #indices should be sorted
    ixx=0
    for ix,i in enumerate(gen):
        if ix==indices[ixx]:
            ixx+=1
            yield i
        if ixx==len(indices):
            break

def pseudo_randomize(x,steps=10):
    r=random.sample(range(steps),steps)
    return "".join(sum([list(x)[i::len(r)] for i in r],[]))


def write_decoy(input_file,output_folder,method="pseudo_random" #pseudo_random, reverse or random
                  ,Write_batch=10**6
                  ,decoy_prefix="decoy_",
                 
                  #this section is for adding a scaling decoy form the initial database
                  Sample=False,
                  Initial_entries=False):
 
    
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    decoy=str(Path(Output_directory,output_folder,"decoy.fa"))
    
    with open(input_file) as f: #in case of massive files there could also be a batched reader added here
        lines=f.readlines()
    
    #here a scaled decoy is added
    if Sample and Initial_entries:
        if Sample>Initial_entries:
            indices=list(set(random.choices(range(Initial_entries),k=Sample))) 
        else:
            indices=list(set(random.sample(range(Initial_entries),Sample))) 
        indices.sort()
        chunks=chunk_gen(slice_gen(SeqIO.parse(input_file, "fasta"),indices),size=Write_batch)
    
    else:
        chunks=chunk_gen(SeqIO.parse(input_file, "fasta"),size=Write_batch)
    
    with open(decoy,"w") as f: #clears file if exists
        pass
    
    with open(decoy,"a") as f:
    
        for chunk in chunks:
    
            df=pd.DataFrame([[rec.description,str(rec.seq)] for rec in chunk],columns=["description","seq"])
            df.description=decoy_prefix+df.description
            
            if method=="reverse":
                df.seq=df.seq.str[::-1]
                
            if method=="pseudo_random": #random staggered join (faster) 
                df.seq=df.seq.apply(pseudo_randomize)
    
            if method=="random":
                df.seq=df.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
        
            f.write("\n".join(">"+df.description+"\n"+df.seq)+"\n")
    
    return decoy
    
def merge_target_decoy(target,decoy):
    
    
    target_decoy=str(Path(Path(target).parents[0],"target_decoy.fa"))
    with open(target_decoy,'wb') as o:
        shutil.copyfileobj(open(target,'rb'), o)
        shutil.copyfileobj(open(decoy ,'rb'), o)

    return target_decoy



def write_database_composition(Database,
                               output_folder,
                               taxonomy="NCBI"):
    

    recs=SeqIO.parse(Database,format="fasta")
    chunks=chunk_gen(recs,size=10**9)
    

    entries=0
    taxa=[]
    
    for chunk in chunks:
        df=pd.DataFrame([[r.id,r.description,str(r.seq)] for r in chunk],columns=["id","description","seq"])

        if taxonomy=="NCBI":                
            taxa.extend(df.description.str.split("OX=").apply(lambda x: x[-1]).str.split(" ").apply(lambda x: x[0]).astype(int))
        if taxonomy=="GTDB":
            taxa.extend(df["id"].str.split("_").apply(lambda x: "_".join(x[-3:])).tolist())
        
        entries+=len(df)

    tax_counts=pd.DataFrame.from_dict(Counter(taxa),orient="index",columns=["Count"])
    
    if taxonomy=="NCBI":
        tax_counts=tax_counts.merge(ncbi_taxdf.set_index("OX"),left_index=True,right_index=True,how="left").dropna().sort_values(by="Count",ascending=False)
    if taxonomy=="GTDB":
        tax_counts=tax_counts.merge(gtdb_taxdf,how="left",left_index=True,right_index=True).dropna().sort_values(by="Count",ascending=False)
    
    tax_counts.to_csv(str(Path(Output_directory,output_folder,"database_composition.tsv")),sep="\t")
           
    
    return tax_counts,len(tax_counts),entries #composition,"richness",total entries
    




### Final section ###
from scipy.optimize import curve_fit
def monod(x,a,b):
    return a*x/(x+b)

def extrapolate_database(folders,
                         output_folder):


    folders=output_folders
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))

    names=[Path(folder).name for folder in folders]
    dfs=[]
    
   
    for folder in folders:

    
        df=pd.read_csv(str(Path(folder,"database_composition.tsv")),sep="\t")
        dfs.append(df.groupby("species")["Count"].sum())
        

    dfs=pd.concat(dfs,axis=1).dropna()#.fillna(0)#.dropna()
    dfs.columns=names
    

    dfs=1/dfs #make inverse for fitting
    dfs=dfs.sort_values(by=names[-1],ascending=False)#.reset_index()

    #Monod curve fit
    x=np.arange(1,len(folders)+1)
    opts=[]
    for r,i in dfs.iterrows():
        y=i.values
        try:
            popt, pcov = curve_fit(monod, x, y,maxfev=1000,p0=[y[-1],0] , bounds=(np.array([y[-1], -0.1]),
                                                                                  np.array([1,100])))
        except:
            popt=[1,1]
        opts.append(popt)
    
    dfs[["a","b"]]=opts
    dfs[names]=(1/dfs[names]) 
    dfs["a_inv"]=1/dfs["a"]
    dfs["above_cutoff"]=dfs["a_inv"]>=2

    dfs.to_csv(str(Path(basedir,Output_directory,output_folder,"extrapolated_species.tsv")),sep="\t")
    taxa=dfs[dfs["above_cutoff"]].index.tolist()

    return taxa


def filter_Database_taxonomy(Database,
                    taxids,
                    output_folder,
                    taxonomy="NCBI",
                    output_file=None
                        ):
        
        if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
        
        if output_file==None:
            output_file="target.fa" #refactor later
        
        output_path=str(Path(Output_directory,output_folder,output_file))
                         
        recs=SeqIO.parse(Database,format="fasta")
        chunks=chunk_gen(recs)
        
        with open(output_path,"w") as f: #clears file if exists
            pass
        with open(output_path,"a") as f:
            
            for ic,c in enumerate(chunks):
                print("chunk "+ str(ic))
        
                chunk_df=pd.DataFrame([[str(r.seq),r.description,r.id] for r in c],columns=["seq","description","id"])
                
                if taxonomy=="NCBI":
                    chunk_df=chunk_df[chunk_df["description"].str.split("OX=").apply(lambda x: x[-1].split(" ")[0]).astype(int).isin(taxids)]
                if taxonomy=="GTDB":
                    chunk_df=chunk_df[chunk_df["id"].str.split("_").apply(lambda x: "_".join(x[-3:])).isin(taxids)]
    
    
                f.write("\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"])+"\n")
    
        return output_path
    
def load_full_db(Database,taxonomy=""): #load in memory (works only for small databases)
    recs=SeqIO.parse(Database,format="fasta")
    rdf=pd.DataFrame([[str(r.seq),r.description,r.id] for r in recs],columns=["seq","description","id"])
    

    if taxonomy=="NCBI":                
        rdf["OX"]=rdf.description.str.split("OX=").apply(lambda x: x[-1]).str.split(" ").apply(lambda x: x[0]).astype(int)
        rdf=rdf.merge(ncbi_taxdf.set_index("OX"),left_on="OX",right_index=True,how="left").dropna()
    if taxonomy=="GTDB":
        rdf["OX"]=rdf["id"].str.split("_").apply(lambda x: "_".join(x[-3:])).tolist()
        rdf=rdf.merge(gtdb_taxdf,how="left",left_on="OX",right_index=True).dropna()
    
    return rdf.set_index("id")

def filter_database_proteins(in_memory_database,proteins,
                             output_folder,
                            ):

  

    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    output_path=str(Path(Output_directory,output_folder,"target.fa"))
    
    in_memory_database=in_memory_database.loc[proteins,:]    
    
    
    cols=['OX', 'superkingdom', 'phylum', 'class', 'order','family', 'genus', 'species', 'Dump_taxid', 'OS']
    cols=[c for c in cols if c in in_memory_database.columns]
    comp=in_memory_database.groupby(cols).size().rename("Count").sort_values(ascending=False).reset_index()
    comp.to_csv(str(Path(basedir,Output_directory,output_folder,"database_composition.tsv")),sep="\t")
        
    
    with open(output_path,"w") as f: #clears file if exists
        pass
    with open(output_path,"a") as f:
        f.write("\n".join(">"+in_memory_database["description"]+"\n"+in_memory_database["seq"])+"\n")
    
    

    return in_memory_database,output_path,len(comp),len(in_memory_database) #mem_db, written_db, richess, entries
    
    


    

#%%

#Script

#This example workflow is written for 3 cycles of MSFragger, +PEPNET in cycle 1, +SMSNET in cycle 3
#Decoy sequences are only in MSFragger annotations during the last cycle.


############ Input Raw Files #############

# raw_files=["C:/MP-CHEW/Datasets/PXD023217 CAMPI/Raw/S01.raw",
# "C:/MP-CHEW/Datasets/PXD023217 CAMPI/Raw/S02.raw"]

############ Prep Raw Files #############


#convert raw to mzML
#mzML_files=[raw2mzML(file) for file in raw_files]

mzML_files=["C:/MP-CHEW/HB_mzML/S01.mzML",
"C:/MP-CHEW/HB_mzML/S02.mzML"]






############ parameters ############


taxonomy="GTDB"

#databases
clustered_database="H:/Databases/Uniref50/uniref50_BacArch_NoAmb_IJeqL.fasta" 



unclustered_diamond_database="H:/Databases/UniprotKB/UniprotKB_BacArch_NoAmb_NoDump_IJeqL.dmnd"
unclustered_database_fasta="H:/Databases/UniprotKB/UniprotKB_BacArch_NoAmb_NoDump_IJeqL.fa"

unclustered_diamond_database="H:/Databases/GTDB/GTDB_merged_NoAmb_IJeqL.dmnd"
unclustered_diamond_fasta="H:/Databases/GTDB/GTDB_merged_NoAmb_IJeqL.fasta"

#MSfragger params files
params_fast          =str(Path(basedir,"closed_fragger_fast.params"))
params_mid           =str(Path(basedir,"closed_fragger_mid.params"))
params_final         =str(Path(basedir,"closed_fragger_final.params"))

#%%


# ############ Initialize #############


output_folder="cycle_0" 
output_folders=[]
output_folders.append(str(Path(Output_directory,output_folder)))

# #MSFragger annotation
MSFragger_files=MSFragger_annotation(input_files=mzML_files,   #full_path, .mzML
                                      database_path=clustered_database, #full_path, .fasta
                                      output_folder=output_folder,
                                      params_path=params_fast,
                                      no_splits=20, #50 #(1GB+0.1GB per file)*splits = ~temporary writing space (also depends on database size...)
                                      no_batches=1 #depends on number of samples
                                      )


#%%

output_folder="cycle_0"
MSFragger_files=["C:/MP-CHEW/CHEW/cycle_0/S02.pin",
"C:/MP-CHEW/CHEW/cycle_0/S01.pin"]

#Write to fasta for Diamond alignment
fasta_file=write_to_Diamond_fasta(MSFragger_files=MSFragger_files,
                                  Top_score_fraction=False, 
                                  output_folder=output_folder,
                                  output_file="Diamond_fasta.fa")




#Diamond alignment
Alignment=Diamond_alignment(fasta_file,
                            database_path=unclustered_diamond_database, #Database,#"UniprotKB_NoAmb_IJeqL.dmnd",
                            output_columns=["qseqid","sseqid","stitle","bitscore","full_sseq"], #add full_sseq for writing first database
                            output_folder=output_folder,
                            select=" -k25",  
                            block_size=3,  #decrease for lower RAM consumption
                            index_chunks=1) #increase for lower RAM consumption


target=Write_alignment_to_database(Alignment,output_folder=output_folder)

#%%

output_folder="cycle_0"
target="C:/MP-CHEW/CHEW/cycle_0/target_unique.fa"

utarget,entries=Unique_fasta(target) 
composition,richness,entries=write_database_composition(Database=utarget,
                                                    taxonomy=taxonomy,
                                                    output_folder=output_folder)





#frequency filtering
taxids=composition[composition["Count"]>20].index.tolist() #this step is kind of tricky, make it dynamic? 


target=filter_Database_taxonomy(utarget,taxids,
                        taxonomy=taxonomy,
                        output_file="filtered.fa",
                        output_folder=output_folder)


DB_in_mem=load_full_db(target,taxonomy=taxonomy)
composition,richness,entries=write_database_composition(Database=target,
                                    taxonomy=taxonomy,
                                    output_folder=output_folder)

Initial_Database,Initial_entries=target,entries #backup since cycle overwrites utarget,entries
decoy=write_decoy(Initial_Database,output_folder=output_folder)
Database=merge_target_decoy(Initial_Database,
                            decoy)



#%%
output_folders=[]


cycle=0

while True:

    cycle+=1
    print("Starting cycle: "+str(cycle))
    output_folder="cycle_"+str(cycle) 
    output_folders.append(str(Path(Output_directory,output_folder)))
    old_richness=richness



    MSFragger_files=MSFragger_annotation(mzML_files,
                                        Database, #unique, non-ambiguous
                                        output_folder=output_folder,
                                        params_path=params_mid,
                                        no_splits=2,
                                        no_batches=1)
    

    weight_rank="species"
    max_evalue=10
    Precision_prefilter=0.7 #Target Decoy precision based denoising pre LCA filter
    Frequency_prefilter=2 # static prefiler cutoff
    Top_score_fraction=0.9
    taxonomy="GTDB"
    
    

    prots=[]
    for ix,f in enumerate(MSFragger_files):
    

    
        pepdf=pd.read_csv(f,sep="\t")
        
        pepdf["u_ix"]=np.unique(pepdf["ScanNr"].astype(str)+"_"+pepdf["ExpMass"].astype(str),return_inverse=True)[1]
        pepdf["u_ix"]=pepdf["u_ix"].astype(str)+"_"+str(ix)
        #pepdf["Peptide"]=pepdf["Peptide"].str[1:-1].str.strip(".").str.replace("[","]").str.split("]").apply(lambda x: "".join(x[::2]))

        
      
        pepdf[["hyperscore","log10_evalue"]]=pepdf[["hyperscore","log10_evalue"]].astype(float)
        pepdf["evalue"]=10**pepdf.log10_evalue
        if max_evalue:
            pepdf=pepdf[pepdf["evalue"]<=max_evalue]
        if Top_score_fraction:
            pepdf=pepdf[(pepdf["evalue"]/pepdf.groupby("u_ix")["evalue"].transform('min'))<=(1/Top_score_fraction)]
    
    
        prots.append(pepdf[["u_ix","Proteins",'hyperscore']])

    prots=pd.concat(prots)
    prots["Proteins"]=prots["Proteins"].str.strip().str.split(" ")
    prots=prots.explode("Proteins")
    prots["Decoy"]=prots["Proteins"].str.startswith("decoy_")
    
 
    
    # s=p[p.index.str.contains("Streptomyces")]
    # b=p[p.index.str.contains("Bifido")]
    # l=p[p.index.str.contains("Lactiplanti")]
    
    # s["precision"].plot.hist(bins=10)
    # b["precision"].plot.hist(bins=10)
    # l["precision"].plot.hist(bins=10)
    


    
    #for UniprotKB this needs you to put the OX= directly in the header. implement this in future versions!
    
    if taxonomy=="GTDB":
        prots["OX"]=prots["Proteins"].str.replace("decoy_","").apply(lambda x: "_".join(x.split("_")[-3:]))
        prots[ranks]=gtdb_taxdf.loc[prots["OX"].tolist(),ranks].values
    
    g=prots.groupby([weight_rank,"Decoy"]).size()
    p=g.rename("Count").reset_index().pivot(index="species",columns="Decoy",values="Count").fillna(0)
    p.columns=["decoy" if i else "target" for i in p.columns]
    
    
    ### in fhte final pipeline you need to add a file specific index and then merge
    #because this way a taxa and protein only gets removed when it is not detected in all files. 
    
    
    p["precision"]=p["target"]/(p["target"]+p["decoy"])
    p.to_csv(str(Path(Output_directory,output_folder,"precision.tsv")),sep="\t")
    p=p[p["target"]>=Frequency_prefilter]
    p=p[p["precision"]>=Precision_prefilter]
    
    target=prots[~prots["Decoy"]]
    target=target.merge(p,how="inner",left_on=weight_rank,right_index=True)
    target["score"]=target['hyperscore']*target["precision"]
    
    tlca=weighted_lca(target,group_on="u_ix",weight_column="score",protein_column="Proteins",weight_cutoff=0.6)
    tlca.to_csv(str(Path(Output_directory,output_folder,"lca.tsv")),sep="\t")
    
    
    ps=tlca["proteins"].reset_index(drop=True).str.split(", ")
    proteins=ps.explode().unique().tolist()
        
 

    DB_in_mem, target, richness, entries=filter_database_proteins(DB_in_mem,proteins,output_folder=output_folder)
    decoy=write_decoy(Initial_Database,output_folder=output_folder,
                      
                      #here a scaled decoy is added
                      Sample=entries*cycle+1,
                      Initial_entries=Initial_entries)
    

    Database=merge_target_decoy(target,
                                decoy)




    # #Cleanup
    # for i in [decoy,fasta_file,Alignments]:
    #     try:
    #         os.remove(i)
    #     except:
    #         pass

    decrease_ratio=richness/old_richness

    print("Database size decrease fraction : "+str(decrease_ratio))
    if decrease_ratio>0.95: #0.95
        break

#%% pick leafs



df=pd.read_csv("C:/MP-CHEW/CHEW/cycle_12_lca.tsv",sep="\t")


df["Sequence"]=df["u_ix"] #bandaid


df["proteins"]=df["proteins"].str.split(", ")
edf=df.explode("proteins")
edf["taxa"]=edf["proteins"].str.split("_").apply(lambda x: "_".join(x[-3:]))

q=edf.groupby("taxa").size().sort_values(ascending=False)
edf=edf[edf["taxa"].isin(q[(q.cumsum()/q.sum())<0.9].index)] #remove tail
df=df.merge(edf.groupby("Sequence",sort=False)["taxa"].apply(list),left_on="Sequence",right_index=True,how="inner")
#df["taxa"]=edf.groupby("Sequence",sort=False)["taxa"].apply(list).values

#node strategy
ranks= [ 'phylum', 'class', 'order', 'family',
       'genus']#, 'species']
c=0
tax=[]
for r in ranks[::-1]:
    print(r)
    
    u=df[r].dropna().unique()
    
    
    
    for t in u:
        print(t)
        tu=df[df[r]==t]
        ul=len(tu)
        tue=tu.explode("taxa").reset_index()
        
        
        if len(tue["taxa"].unique())>5:
   
            
            g=pd.DataFrame(tue.groupby("taxa")["index"].apply(set))
            g["l"]=g["index"].apply(len)
            g=g.sort_values(by="l",ascending=False)
            
            t=set()
            ls=[0]
            d1=0
            for ix,i in enumerate(g["index"]):
                t.update(i)
                l=len(t)
             
                
                ls.append(l)
                d=ls[-1]-ls[-2]
                
                if d<2:   #too small increase
                    break
                if d==d1: #linear regime, should make this: "linearish"
                    break
                
                if (l/ul)>0.95:
                    break
                
                d1=d
                

            
            #tax.extend(g.index[(ix-1):].tolist()) #taxa to remove
            
            tax.extend(g.index[:ix].tolist()) #taxa to keep?
    
    #edf=edf[~edf.taxa.isin(tax)]
    edf=edf[edf.taxa.isin(tax)]
    df.pop("taxa")
    df=df.merge(edf.groupby("Sequence",sort=False)["taxa"].apply(list),left_on="Sequence",right_index=True,how="inner")
    

ja=edf[["genus","taxa"]].drop_duplicates().groupby("genus").size()

#%% Final search

# output_folder="final_search"
# unclustered_db="H:/Databases/GTDB/GTDB_merged_NoAmb_IJeqL.fasta"
# db=filter_Database_taxonomy(Database=unclustered_db,taxids=tax, taxonomy=taxonomy, output_file="FINAL.fa",
#             output_folder=output_folder)

write_decoy(db,method="reverse",output_folder=output_folder)

#%%
target="C:/MP-CHEW/CHEW/final_search/FINAL.fa"
decoy="C:/MP-CHEW/CHEW/final_search/decoy.fa"
d=merge_target_decoy(target, decoy)

#%%



MSFragger_files=MSFragger_annotation(mzML_files,
                                    database_path="C:/MP-CHEW/CHEW/final_search/target_decoy.fa",
                                    output_folder=output_folder,
                                    params_path=params_mid,
                                    no_splits=2,
                                    no_batches=1)

#%%
write_database_composition(target,output_folder=output_folder,taxonomy=taxonomy)

#%%
# #%%

# ############ Final round #############
# print("Starting Final search")
# output_folder="Final_full_db"


# taxonomy="GTDB"
# #%% bandaid

# #def extrapolate_database(folders,
# #                         output_folder):



# output_folders=[


# "C:/MP-CHEW/CHEW/cycle_1",
# "C:/MP-CHEW/CHEW/cycle_2",
# "C:/MP-CHEW/CHEW/cycle_3",
# "C:/MP-CHEW/CHEW/cycle_4",
# "C:/MP-CHEW/CHEW/cycle_5",
# "C:/MP-CHEW/CHEW/cycle_6",
# "C:/MP-CHEW/CHEW/cycle_7",
# "C:/MP-CHEW/CHEW/cycle_8",
# "C:/MP-CHEW/CHEW/cycle_9",
# "C:/MP-CHEW/CHEW/cycle_10",
# "C:/MP-CHEW/CHEW/cycle_11",
# "C:/MP-CHEW/CHEW/cycle_12"
#                 ]
# #output_folders.sort()

# output_folder="C:/MultiNovo/Final_full_db"

# rank="genus"
# folders=output_folders
# if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))

# names=[Path(folder).name for folder in folders]
# dfs=[]

   
# for folder in folders:


#     df=pd.read_csv(str(Path(folder,"database_composition.tsv")),sep="\t")
#     dfs.append(df.groupby(rank)["Count"].sum())
    

# dfs=pd.concat(dfs,axis=1).dropna()#.fillna(0)#.dropna()
# dfs.columns=names


# dfs=1/dfs #make inverse for fitting
# dfs=dfs.sort_values(by=names[-1],ascending=False)#.reset_index()

# #Monod curve fit
# x=np.arange(1,len(folders)+1)
# opts=[]
# for r,i in dfs.iterrows():
#     y=i.values
#     try:
#         popt, pcov = curve_fit(monod, x, y,maxfev=1000,p0=[y[-1],1000] , bounds=(np.array([y[-1], -0.1]),
#                                                                               np.array([1,np.inf])))
#     except:
#         popt=[1,1]
#     opts.append(popt)

# dfs[["a","b"]]=opts
# dfs[names]=(1/dfs[names]) 
# dfs["a_inv"]=1/dfs["a"]
# dfs["above_cutoff"]=dfs["a_inv"]>=2

# dfs.to_csv(str(Path(basedir,Output_directory,output_folder,"extrapolated_species.tsv")),sep="\t")
# taxa=dfs[dfs["above_cutoff"]].index.tolist()


#%%





#%%
dfs=dfs.sort_values(by=)


#%%

#%%

a=pd.read_csv("C:/MP-CHEW/CHEW/pretty good/cycle_8/Diamond_fasta_target_lca.tsv",sep="\t")

#%%

# species=extrapolate_database(folders=output_folders,
#                             output_folder=output_folder)


target=construct_final_db(Database=unclustered_database_fasta,
                           species=species,
                           taxonomy=taxonomy,
                           output_folder=output_folder)


decoy=write_decoy(target,output_folder=output_folder,method="reverse")


Database=merge_target_decoy(target,
                            decoy)


MSFragger_files=MSFragger_annotation(mzML_files,
                                    Database, #unique, non-ambiguous
                                    output_folder=output_folder,
                                    params_path=params_final_selenium,
                                    no_splits=32,#8
                                    no_batches=1
                                    )



#%%
MSFragger_files=["C:/MultiNovo/Final_full_db/T0_1_HK_JLN_UPLIFT.pin",

"C:/MultiNovo/Final_full_db/T0_A1_HK_JLN_UPLIFT.pin",
"C:/MultiNovo/Final_full_db/T1_2_HK_JLN_UPLIFT.pin",
"C:/MultiNovo/Final_full_db/T1_1_HK_JLN_UPLIFT.pin",
"C:/MultiNovo/Final_full_db/T0_A2_HK_JLN_UPLIFT.pin",
"C:/MultiNovo/Final_full_db/T1_A1_HK_JLN_UPLIFT.pin",
"C:/MultiNovo/Final_full_db/T0_2_HK_JLN_UPLIFT.pin",
"C:/MultiNovo/Final_full_db/T1_A2_HK_JLN_UPLIFT.pin"]

target="C:/MultiNovo/Final_full_db/final_database.fasta"

#% Post processing
import string
import csv
mfile="C:/MultiNovo/unimod.txt"
with open(mfile,"r") as f:
    m=f.read()

ms=m.split("[Term]")[2:]
unimod_names=[i.split("name: ")[1].split("\n")[0] for i in ms]
unimod_masses=[i.split('xref: delta_mono_mass "')[1].split('"')[0] for i in ms]
unimod_df=pd.DataFrame(list(zip(unimod_names,unimod_masses)),columns=["unimod_name","unimod_mass"])
unimod_df["unimod_mass"]=unimod_df["unimod_mass"].astype(float)

std_aa_mass = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841,
               'T': 101.04768,'C': 103.00919,'L': 113.08406,'I': 113.08406,'J': 113.08406,
               'N': 114.04293,'D': 115.02694,'Q': 128.05858,'K': 128.09496,'E': 129.04259,
               'M': 131.04049,'H': 137.05891,'F': 147.06841,'U': 150.95364,'R': 156.10111,
               'Y': 163.06333,'W': 186.07931,'O': 237.14773}  



MSFragger_files.sort()
pepXML_files=[i.replace(".pin",".pepXML") for i in MSFragger_files]

proteins=set()
pepdfs=[]
for ix,file in enumerate(MSFragger_files):
    print(ix)
    

    pepdf=read_pin(file)
    pepdf['Proteins']=pepdf['Proteins'].str.replace("\t"," ")  
    pepdf[["hyperscore","log10_evalue"]]=pepdf[["hyperscore","log10_evalue"]].astype(float)

    pepdf.pop("delta_hyperscore")
    pepdf=pepdf.groupby(pepdf.columns[1:-1].tolist())["Proteins"].apply(" ".join).reset_index()
    
    pepdf["evalue"]=10**pepdf.log10_evalue
    if max_evalue:
        pepdf=pepdf[pepdf["evalue"]<=max_evalue]
    if Top_score_fraction:
        pepdf=pepdf[(pepdf["evalue"]/pepdf.groupby(["ScanNr","ExpMass"])["evalue"].transform('min'))<=(1/Top_score_fraction)]
        
        
    pepdf=pepdf.groupby([i for i in pepdf.columns if i!="rank"],sort=False).nth(0).reset_index() #remove duplicate hits
    
    #Filter at 0.05 FDR
    pepdf["Decoy"]=pepdf["Proteins"].str.contains("decoy")
  

    pepdf=pepdf.sort_values(by="hyperscore",ascending=False).reset_index(drop=True)
    pepdf["FDR"]=pepdf["Decoy"].cumsum()/np.arange(1,len(pepdf)+1)
    if (pepdf["FDR"]<0.05).sum():
        minscore=pepdf["hyperscore"].iloc[(~(pepdf["FDR"]<0.05)).idxmax()]
        pepdf["passed_0.05_FDR_threshold"]=pepdf["hyperscore"]>=minscore
    
    pepdfs.append(pepdf)
    proteins.update(pepdf["Proteins"].str.split().explode().tolist())
    
#proteins.to_csv(str(Path(Output_directory,output_folder,"found_proteins.tsv")),sep="\t")

#retrieve fasta headers and full protein sequences from final database

recs=SeqIO.parse(Database,format="fasta")
chunks=chunk_gen(recs,size=10**9)
dfs=[]
for chunk in chunks:
    df=pd.DataFrame([[r.id,r.description,str(r.seq)] for r in chunk],columns=["Proteins","description","seq"])
    dfs.append(df[df["Proteins"].isin(proteins)])
dfs=pd.concat(dfs)

if taxonomy=="NCBI":
    dfs["OX"]=dfs.description.str.split("OX=").apply(lambda x: x[-1]).str.split(" ").apply(lambda x: x[0]).astype(int)
    dfs=dfs.merge(ncbi_taxdf,on="OX",how="left") 

#write found proteins
with open(str(Path(Output_directory,output_folder,"found_proteins.fa")),"w") as f:
    f.write("\n"+"\n".join(">"+dfs["description"]+"\n"+dfs["seq"]))


for ix,pepdf in enumerate(pepdfs):
    

    ### taxonomic annotation ###
    if taxonomy=="NCBI":
        lins=pd.DataFrame(pepdf["Proteins"].str.split().explode()).reset_index().merge(dfs,how="left",on="Proteins").set_index("index")
    if taxonomy=="GTDB":
        lins=pd.DataFrame(pepdf["Proteins"].str.split().explode().apply(lambda x: "_".join(x.split("_")[-3:]))).reset_index().merge(gtdb_taxdf,how="left",left_on="Proteins",right_on="OS").merge(dfs,how="left",on="Proteins").set_index("index")

    pepdf[ranks]=lins.groupby("index")[ranks].nth(0)[(lins.groupby("index")[ranks].nunique()==1)] #vectorized lca
   
    ### Selenium annotation

    # variable_mod_07 = 57.021464 UC 1           #Carbamidomethylation HK (now coded as optional modification)
    # variable_mod_08 = 104.965914 C 1           #Carbamidomethylation HK + S to Se
    # variable_mod_09 = 63.939365  M 1           #M oxidation HK + S to Se 
    pepdf["seleno_peptide"] = pepdf["Peptide"].str.contains("U")
    pepdf["cysteine_S_Se"]  =(pepdf['Peptide'].str.contains("C[57",regex=False)) | (pepdf['Peptide'].str.contains("C[104",regex=False))
    pepdf["methionine_S_Se"]=(pepdf['Peptide'].str.contains("M[57",regex=False)) | (pepdf['Peptide'].str.contains("M[63",regex=False))
    pepdf["selenoprotein"]  =lins.seq.str.contains("U").groupby("index").any()
    
    
    #write output
    pepdf.to_csv(MSFragger_files[ix].replace(".pin","_processed.pin"))
    
    ### CalisP formatting (convert pip format to sequest PSMs format)
    #Annotated Sequence	Confidence	PSM Ambiguity	Modifications	First Scan	Spectrum File	Charge	m/z [Da]	MH+ [Da]	RT [min]	    XCorr	Protein Accessions	Master Protein Accessions	# Proteins
    # PSMs=pd.DataFrame(pepdf[['ScanNr', "SpecId","ExpMass", "retentiontime", "hyperscore", "Proteins"]])
    # PSMs.columns=["First Scan","Spectrum File", "MH+ [Da]","RT [min]","XCorr","Protein Accessions"] 
    PSMs=pd.DataFrame(pepdf[['ScanNr',"ExpMass", "retentiontime", "hyperscore", "Proteins"]])
    PSMs.columns=["First Scan", "MH+ [Da]","RT [min]","XCorr","Protein Accessions"] 
    
    #reformat columns
    PSMs["Protein Accessions"]=PSMs["Protein Accessions"].str.rstrip().str.replace(" ", ", ") 
    #PSMs["Spectrum File"]=PSMs["Spectrum File"].str.split(".").apply(lambda x: x[0]+".raw")
    PSMs["Spectrum File"]=Path(MSFragger_files[ix]).stem+".raw"
    
    PSMs[["Confidence","PSM Ambiguity"]]="" #filler
    

    
    #add charge from pepXML file
    with open (pepXML_files[ix],"r") as f:
        
 
    # test="C:/MultiNovo/cycle_0/T0_1_HK_JLN_UPLIFT.pepXML"
    # with open (test,"r") as f:
        xml_string=f.read()
        ps=[i.split("</spectrum_query>")[0] for i in xml_string.split("<spectrum_query ")]
        parse_targets=["start_scan","assumed_charge"]
        s=[[i.split(p+'=')[1].split('"')[1] for i in ps[1:]] for p in parse_targets]
        specdf=pd.DataFrame(s,index=parse_targets).T.astype(int) 
    
    PSMs["Charge"]=specdf.merge(PSMs["First Scan"].astype(int),left_on="start_scan",right_on="First Scan",how="right")["assumed_charge"].fillna(1).astype(int)
    PSMs["m/z [Da]"]= (PSMs["MH+ [Da]"].astype(float)+PSMs["Charge"]*1.007825319)/PSMs["Charge"]
    PSMs["MH+ [Da]"]=PSMs["MH+ [Da]"].astype(float)+1.007825319 #add H+


    #%missing columns: Annotated Sequence	Confidence	PSM Ambiguity	Modifications Charge	m/z [Da] Master Protein Accessions # Proteins
    PSMs["Master Protein Accessions"]=PSMs["Protein Accessions"].copy()
    PSMs["# Proteins"]=PSMs["Protein Accessions"].str.count(",")+1
    
    #add modifications
    peptides=pepdf["Peptide"].apply(lambda x: x[2:-2])
    PSMs["Modifications"]=""
    

    #get known modifications
    has_mod=peptides[peptides.str.contains("[",regex=False)]
    unimods=has_mod.str.replace("[","]").str.split("]").apply(lambda x: x[1::2]).explode().drop_duplicates().astype(float).reset_index(drop=True).values
    #unimods=np.array([float(i.strip(string.ascii_letters)) for i in ";".join(pepdf.columns).split(";nmc;")[1].split(";Peptide;")[0].split(";")])
    um_dict=unimod_df.iloc[np.argmin(abs(unimod_df["unimod_mass"].values-unimods.reshape(-1,1)),axis=1)].set_index(unimods.astype(str))["unimod_name"].to_dict()

    ms=[]
    for i in has_mod.str.replace("[","]",regex=False).str.split("]"):
        m=[]
        aa_ix=0
        for oi,x in enumerate(i):
            if "." not in x:
                aa_ix+=len(x)
            else:
                if aa_ix:
                    m+=[i[oi-1][-1]+str(aa_ix)+"("+um_dict.get(x)+")"]
                
        ms.append("; ".join(m))
    PSMs.loc[has_mod.index,"Modifications"]=ms
    PSMs["Modifications"]=PSMs["Modifications"].str.replace("n1(Acetyl)","N-Term(Prot)(Acetyl)",regex=False)
    PSMs["Modifications"]=PSMs["Modifications"].str.replace("Cys->CamSec","Carbamidomethyl",regex=False) #for CalisP
 
    
    #lowercasing
    l_peptides=peptides.str.lstrip("n")
    for i in std_aa_mass.keys():
        l_peptides=l_peptides.str.replace(i+"[",i.lower()+"[",regex=False)
    l_peptides=l_peptides.apply(lambda x: re.sub("[\[\[].*?[\]\]]", "", x))
    l_peptides.update(l_peptides.loc[peptides.str.startswith("n")].apply(lambda x: x[1].lower()+x[2:]))
    PSMs["Annotated Sequence"]=l_peptides
    

    #write output
    PSMs=PSMs[[
          "Annotated Sequence" ,
          "Confidence"    ,
          "PSM Ambiguity" ,
          "Modifications" ,
          "First Scan",
          "Spectrum File"  ,
          "Charge" ,
          "m/z [Da]" ,
          "MH+ [Da]" ,
          "RT [min]" ,
          "XCorr" ,
          "Protein Accessions" ,
          "Master Protein Accessions" ,
          "# Proteins" ]].fillna("").drop_duplicates()
    
    PSMs.columns=['"'+c+'"' for c in PSMs.columns]
    PSMs='"'+PSMs.astype(str)+'"'
    PSMs.to_csv(MSFragger_files[ix].replace(".pin","_PSMs.txt"),sep="\t",index=False,quoting=csv.QUOTE_NONE,escapechar='\\')

    #still missing charge Ambigouus, mz and some others
    
