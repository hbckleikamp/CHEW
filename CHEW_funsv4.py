# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 13:26:43 2023

@author: ZR48SA
"""

#Modules



import Bio
from Bio import SeqIO

import subprocess
import pandas as pd
import psutil
import shutil
import sys
import numpy as np
import re
import math
import random
import time
from collections import Counter
import itertools
import more_itertools as mit
import inspect
import datetime
from pathlib import Path
import os

### only needed for post processing:
import seaborn as sns
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
### 

#local module
from load_vars import *
from config import *
#### function variables ###

from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0])) # set base path
basedir=os.getcwd()



#%% basic params


#SMSNet filepaths
simple_masks=pd.read_csv(str(Path(basedir,"simple_unmasks.tsv")),sep="\t")

#Base taxonomy filepaths
ranks=np.array(["superkingdom","phylum","class","order","family","genus","species"]) 
taxdf_path=str(Path(basedir,"Setup","parsed_taxonomy.tsv"))


IO_batch=10**6 #how much lines should be written or read at once from a fasta file 
diamond_output_columns=["qseqid","sseqid","stitle","bitscore","full_sseq"]



#%%
import glob
list_of_files = glob.glob(str(Path(basedir,"*.CHEW_params"))) 
latest_file   = max(list_of_files, key=os.path.getctime)
kws=load_variables(latest_file)

#%%
####



#### Overcomplicated stuff for handling arg passing into functions ###

class function_variables(): #function variables
    def __init__ (self):
        self.output_folder=basedir
        self.tmp_folder=basedir
        self.Output_directory=basedir
        self.Temporary_directory=basedir
        self.basedir=basedir


def dargs(func):
    signature = inspect.signature(func)
    return {
        k: v.default
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    }    


def passed_kwargs():

  

    def decorator(function,
                  kws=kws
                  ):
        def inner(*args, **kwargs):
            
            class_instance=function_variables()

            #update variables (hierachy: kws dict < defaults < passed)
            
            #can also try to update kws each time a function is run?
            class_instance.__dict__.update(dargs(function))
            class_instance.__dict__.update(kws)
            class_instance.__dict__.update(kwargs)
            
            ### some other general functions ###

            #parse output_folders
            output_folder,tmp_folder,Output_directory,Temporary_directory,basedir=class_instance.output_folder,class_instance.tmp_folder,class_instance.Output_directory,class_instance.Temporary_directory,class_instance.basedir
            
            if not os.path.isabs(output_folder): output_folder=str(Path(basedir,Output_directory,output_folder))
            if not os.path.isabs(tmp_folder):       tmp_folder=str(Path(basedir,Temporary_directory,tmp_folder))
            
            
            class_instance.output_folder=output_folder
            class_instance.tmp_folder=tmp_folder 
                
            #make missing directories
            for i in output_folder,tmp_folder:
                if not os.path.exists(i): os.makedirs(i, exist_ok=True) 


            #arg logger (add "write_args":True to kws)
            if hasattr(class_instance,"write_args"):
                if class_instance.write_args:
                    filename=datetime.now.strftime("%y_%m_%d_%H_%M_%S")+function.__name__+".txt"
                    df=pd.DataFrame.from_dict(class_instance.__dict__)
                    df.to_csv(filename,sep="\t")
                    

            inner.vars=class_instance

            return function(*args, **kwargs)
        return inner
    return decorator




def isiterable(e):
    try:
        iter(e)
        return 1
    except:
        return 0

def check_required_args(v,l): #list of strs
    for i in l:    
        try: 
            e=getattr(v,i)
            if isiterable(e):
                if len(e)==0:
                    raise Exception ("the variable '"+i+"' is a required argument, and is not defined properly! (variable is empty)" )
        except AttributeError as err:
            raise Exception ("the variable '"+i+"' is a required argument, and is not defined properly! (variable name does not exist)" )


   

    

#%% Decorated funs

@passed_kwargs()
def test(**kwargs):
    v=test.vars
    ja=v.ja
    
    check_required_args(v,["ja"])
    
    print("ja is "+str(ja))

@passed_kwargs()
def raw2mzML(*,output_folder="mzML",msconvert_filepath=msconvert_filepath,**kwargs): 
    
    v=raw2mzML.vars #parse function arguments
    output_folder=v.output_folder
    input_files=v.input_files
    check_required_args(v,["input_files"])

    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".raw")]
        else: input_files=input_files.split()

    
    

    output_files=[]
    for input_file in input_files:
        output_file=str(Path(output_folder,Path(input_file).stem+".mzML"))
        output_dir=str(Path(output_file).parents[0])
        
        if not input_file.endswith(".mzML"):
            output_files.append(output_file)
            
            command="cd" +' "'+output_folder+'" && '+msconvert_filepath
            command+='"'+input_file+'"' 
            command+=' --mzML --filter "peakPicking vendor" --filter "zeroSamples removeExtra" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>" -o '+'"'+output_dir+'"'
            #command+=' --mzML --filter "peakPicking vendor" --filter "zeroSamples removeExtra" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>"'
           
            print(command)
            stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            
        else:    
            output_files.append(input_file)
            
    return output_files

@passed_kwargs()
def raw2mgf(*,output_folder="mgf",msconvert_filepath=msconvert_filepath,**kwargs):
    
    v=raw2mgf.vars #parse function arguments
    output_folder=v.output_folder
    input_files=v.input_files
    check_required_args(v,["input_files"])
    
    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".raw")]
        else: input_files=input_files.split()


    #standard filter
    #filter_conditions=' --filter "peakPicking vendor" --filter "zeroSamples removeExtra" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>" '
    
    #more stringent filter to increase speed of SMSNet
    #filter_conditions=' --filter "peakPicking vendor" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>" --filter "threshold bpi-relative .005 most-intense" --filter "threshold count 50 most-intense" '
    
    #more stringent filter to increase speed of SMSNet
    filter_conditions=' --filter "peakPicking vendor" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>" --filter "threshold count 100 most-intense" --filter "zeroSamples removeExtra" '
    
    
    output_files=[]
    for input_file in input_files:
        output_file=str(Path(output_folder,Path(input_file).stem+".mgf"))
        output_dir=str(Path(output_file).parents[0])
        
        if not input_file.endswith(".mgf"):
            output_files.append(output_file)
        
            command="cd" +' "'+output_folder+'" && '+msconvert_filepath
            command+='"'+input_file+'"' 
            command+=' --mgf '+filter_conditions+' -o '+'"'+output_dir+'"'
            print(command)
            stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            
        
        else:    
            output_files.append(input_file)
        
        reformat_mfg(output_files[-1])
        
    return output_files



    

@passed_kwargs()
def MSFragger_annotation(*, 
                         params_path=params_fast,
                         max_no_hits=5, #total number of top hits retained after mergeing database splits
                         no_splits=None, 
                         no_batches=None,
                         output_folder="inital_annotation",
                         **kwargs):
    
   
    v=MSFragger_annotation.vars #parse function arguments
    
    #required arg
    input_files=v.input_files   
    database_path=v.database_path
    check_required_args(v,["input_files","database_path"])
    
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    params_path=v.params_path
    no_splits,no_batches=v.no_splits,v.no_batches
    max_no_hits=v.max_no_hits

    if type(input_files)==str:
        if os.path.isdir(input_files):
            x=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".mzML")]
            if not len(x):
                x=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".raw")]
            if len(x):
                input_files=x
    if type(input_files)==str:
        input_files=input_files.split()
    if type(input_files)==str:
        input_files=[input_files]
    
    #if raw, convert to mzML
    for ix,input_file in enumerate(input_files):
        if input_file.endswith(".raw"):
            input_files[ix]=raw2mzML(input_files=input_file)[0]
    input_files=[i for i in list(set(input_files)) if i.endswith(".mzML")]
    
    
    
    MSFragger_jar_path=str(Path(basedir,"MSFragger-3.5.jar"))  
    pep_split_path=str(Path(basedir,"msfragger_pep_split_HK.py"))


    
    #rewrite closed_fragger.params according to database path
    with open(params_path,"r+") as f:
        lines=f.readlines()
        lines=["database_name = "+database_path+" #database name here\n" if line.startswith("database_name =") else line for line in lines]
        f.seek(0)
        f.writelines(lines)
        
        
    # to kill or not to kill (avoid a bug and RAM overload in file mergeing)
    if sum(["pepxml_pin" in i.split("#")[0]  for i in lines]):
       pep_split_path=str(Path(basedir,"msfragger_pep_split_HK_nokill.py"))     
   
    os.chdir(tmp_folder) #change to tempdir for writing indices
    stderr=""
    retry=0
    while True: #rerun with different settings untill settings are found that fit RAM
        
        #remove old peptide split indices
        if os.path.exists(str(Path(tmp_folder,"split_peptide_index_tempdir"))): shutil.rmtree(str(Path(tmp_folder,"split_peptide_index_tempdir"))) 
        JavaMem=int(psutil.virtual_memory().available/10**9*0.8) #check available RAM
        
        if not no_splits:
            no_splits=math.ceil(os.path.getsize(database_path)/10**9*150/JavaMem) #starting number of splits for db splitting, database size*150 is an estimated scaling factor, which will be effected by the number of allowed modifications and missed cleavages
        
        if not no_batches:
            no_batches=int(psutil.disk_usage(Path(pep_split_path).anchor).free/10**9/os.path.getsize(database_path)) #number of files annotated per time
        if not no_batches>0: no_batches=1
        batches=np.array_split(input_files, math.ceil(no_batches))
        

        print("running MSFragger, total number of splits: "+str(no_splits)+" , total number of batches: "+str(len(batches))+", Java Heap: "+str(JavaMem))
        
        
        retry+=1
        if retry>1:
            print("error!, retrying")
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
    folders=[i for i in Path(basedir,tmp_folder,"split_peptide_index_tempdir").glob("*") if i.is_dir()]
    outpaths=[]
    for mzML_file in input_files:
        
        outpath=str(Path(output_folder,Path(mzML_file).stem+".pin"))
        outpaths.append(outpath)
        
        if len(folders):
            #! Doesnt work for final search, since pin files are cleared up automatically 
            pinfile=Path(mzML_file).stem+".pin"
            pepdf=pd.concat([read_pin(str(Path(folder,pinfile))) for folder in folders])
            if max_no_hits: pepdf=pepdf.sort_values(by=["ScanNr","ExpMass","log10_evalue"]).groupby(["ScanNr","ExpMass"],sort=False).head(max_no_hits)
                
            pepdf["Proteins"]=pepdf.Proteins.str.replace("\t"," ")
            pepdf.to_csv(outpath,sep="\t")

        else:
            pin=str(Path(Path(mzML_file).parents[0],Path(mzML_file).stem))+".pin" #files are written to the location of the mzml files
            if os.path.exists(pin):
                shutil.move(pin,str(Path(output_folder,Path(mzML_file).stem+".pin")))
        
        #if pepXML files exist, move them
        pepxml=str(Path(Path(mzML_file).parents[0],Path(mzML_file).stem))+".pepXML" #files are written to the location of the mzml files
        if os.path.exists(pepxml):
            shutil.move(pepxml,str(Path(output_folder,Path(mzML_file).stem+".pepXML")))

    return outpaths


@passed_kwargs()
def SMSNet_annotation(*,output_folder="inital_annotation",**kwargs):
    
    v=SMSNet_annotation.vars  #parse function arguments
    
    #required arguments
    input_files=v.input_files 
    check_required_args(v,["input_files"])
    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    if type(input_files)==str:
        if os.path.isdir(input_files):
            x=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".mgf")]
            if not len(x):
                x=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".mzML")]
            if not len(x):
                x=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".raw")]
            if len(x):
                input_files=x
    if type(input_files)==str:
        input_files=input_files.split()
    if type(input_files)==str:
        input_files=[input_files]
    
    for ix,input_file in enumerate(input_files):
        if input_file.endswith(".raw") or input_file.endswith(".mzML"):
            input_files[ix]=raw2mgf(input_files=input_file)[0]
    input_files=[i for i in list(set(input_files)) if i.endswith(".mgf")]
    
    
    output_files=[]
    for input_file in input_files:
        
        command="conda activate smsnet  && cd "+str(Path(basedir,"SMSNet-master"))+" && "
        command+="python run_HK.py --model_dir smsnet_phospho --inference_input_file "+input_file+" --rescore"
        print(command)
        stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        
        output_file=str(Path(output_folder,Path(input_file).stem+"_SMSNET.tsv"))
        
        SMSNet_output_file=str(Path(str(Path(input_file).parents[0])+"_output",'p-mod_fdr10.tsv'))

        
        shutil.move(SMSNet_output_file, output_file)

        output_files.append(output_file)
        #cleanup    
        #shutil.rmtree(str(Path(input_file).parents[0])+"_output")
    return output_files



from more_itertools import sliced
import more_itertools as mit
  

@passed_kwargs()
def write_to_Diamond_fasta(*,
                           
                           #MSFragger score filters
                           max_evalue=10,
                           Top_score_fraction=0.9, #in case of multiple top candidates retain the top scoring fraction
                           
                           #SMSNet score filters
                           simple_unmask=True,     #attempts to solve low complexity masked SMSNet regions
                           X_padding=False,        #replace gaps with X
                           SMSNet_ppm=False,       #max ppm tolerance
                           SMSNet_minscore=False,  #minum mean peptide score
                           
                           #fasta writing parameters
                           header_info=["Target_Decoy"],         #information columns that should be retained in the fasta header
                           add_decoy=True,         #decoy alignment for database QC
                           unique_peptides=True,   #only write unique combinations of header and peptide
                           min_length=5,           #minimum tag length
                           output_file="peplist.fa",  #output_file name
                           **kwargs
                           ):
       

    
    v=write_to_Diamond_fasta.vars #parse function arguments
    
    #required arguments
    input_files=v.input_files 
    check_required_args(v,["input_files"])

    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    max_evalue=v.max_evalue
    Top_score_fraction=v.Top_score_fraction
    simple_unmask=v.simple_unmask
    X_padding=v.X_padding
    SMSNet_ppm=v.SMSNet_ppm
    SMSNet_minscore=v.SMSNet_minscore
    header_info=v.header_info
    unique_peptides=v.unique_peptides
    min_length=v.min_length
    output_file=v.output_file

 
    #remove output file to prevent overappending to existing file
    out_path=str(Path(output_folder,Path(output_file).stem+".fa"))
    if os.path.exists(out_path): os.remove(out_path)
    
    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith(".pin") or i.endswith("SMSNET.tsv")]
        else:
            input_files=input_files.split()



    pepdfs=[]

    for im,input_file in enumerate(input_files):
        print("writing "+input_file+" to fasta") 

        pepdf=pd.read_csv(input_file,sep="\t")
        
        if input_file.endswith(".pin"):  #parse MSFragger outputs
            pepdf[["hyperscore","log10_evalue"]]=pepdf[["hyperscore","log10_evalue"]].astype(float)
            pepdf["evalue"]=10**pepdf.log10_evalue
            if max_evalue:
                pepdf=pepdf[pepdf["evalue"]<=max_evalue]
            if Top_score_fraction: #pick best scoring candidate
                pepdf=pepdf[(pepdf["evalue"]/pepdf.groupby(["ScanNr","ExpMass"])["evalue"].transform('min'))<=(1/Top_score_fraction)]
            
            pepdf["peptide_neighbours"]=pepdf["Peptide"].str.replace("c","").str.replace("n","").str.replace("-","").str.replace(".","").apply(lambda x: re.sub("[\[\[].*?[\]\]]", "", x).replace(",","")) #remove ptms in peptides
            pepdf["tag"]=pepdf["peptide_neighbours"]
        
        
        if input_file.endswith("SMSNET.tsv"): #parse SMSNet outputs
            pepdf=parse_SMSNet_output(input_file=input_file,
                                      SMSNet_ppm=SMSNet_ppm,
                                      SMSNet_minscore=SMSNet_minscore,
                                      simple_unmask=simple_unmask,
                                      X_padding=X_padding)


        if "tag" not in pepdf.columns: pepdf["tag"]=pepdf["Peptide"]

        pepdf=pepdf[pepdf["tag"].fillna("").apply(len)>=min_length]
        pepdf["tag"]=pepdf["tag"].str.replace("I","L").str.replace("J","L")
        
        if add_decoy: #in this pipeline this is only used for database QC
            pepdf["Alignment_Decoy"]=False
            decoy=pepdf.copy()
            decoy["tag"]=decoy.str[::-1] #reverse
            decoy["Alignment_Decoy"]=True
            pepdf=pd.concat([pepdf,decoy])
            
        
        #parse fasta header
        for i in header_info: #round numeric header info to save space
            try:
                pepdf[i]=pd.to_numeric(pepdf[i]).round(2)
            except:
                pass
        hdict=pepdf[[i for i in header_info if i in pepdf.columns ]].fillna("").astype(str) 
        pepdf["dict"]="{"+("'"+hdict.columns+"'"+":"+hdict).apply(",".join,axis=1)+"}" 
   
      
        #Writing output
        if unique_peptides:
            pepdfs.append(pepdf[["tag","dict"]])
        else:
            with open(out_path,"a") as f: #create target
                f.write("\n".join(">"+pepdf["tag"]+";"+pepdf["dict"]+"\n"+pepdf["tag"])+"\n")

    if unique_peptides:
        pepdf=pd.concat(pepdfs).drop_duplicates()
        with open(out_path,"a") as f: #create target
            f.write("\n".join(">"+pepdf["tag"]+";"+pepdf["dict"]+"\n"+pepdf["tag"])+"\n")

    return out_path

#function to parse tabular Peptide files
@passed_kwargs()
def parse_peplist(*,
                  header_info=[],         #information columns that should be retained in the fasta header
                  unique_peptides=True,   #only write unique combinations of header and peptide
                  min_length=4,           #minimum tag length
                  output_file=None,       #if all are to be written to the same output file, put filename here, otherwise they are written to separate files
                  add_decoy=False,
                  
                  **kwargs
                  ):

    v=parse_peplist.vars #parse function arguments
    input_file=v.input_file #required argument
    check_required_args(v,["input_file"])
    
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    header_info=v.header_info
    unique_peptides=v.unique_peptides
    min_length=v.min_length
    output_file=v.output_file
    
    
    if output_file==None: output_file=Path(input_file).stem+".fa"
    out_path=str(Path(output_folder,output_file))


    if isinstance(input_file, pd.DataFrame):
        pepdf=input_file
    elif isinstance(input_fle,pd.Series):
        pepdf=input_file.reset_index()

    elif  is_fasta(input_file): #do nothing
        if not add_decoy:
            return input_file
        else:
            pepdf=pd.DataFrame([[r.id,str(r.seq)] for r in Bio.SeqIO.parse(SMSNet_peplist,format="fasta")],columns=["id","Peptide"])
 
    elif os.path.exists(input_file):
        pepdf=read_table(input_file,Keyword="Peptide")
        pepdf=pepdf[pepdf["Peptide"].fillna("").apply(len)>=min_length]
        pepdf["id"]=pepdf["Peptide"]
   
    
    
    if add_decoy:
        header_info+=["Alignment_Decoy"]
        
        #remove palindromic peptides and peptides that exist in both forward and reverse within the same dataset
        rf=pepdf.merge(pepdf["Peptide"].str[::-1].rename("rPeptide"),how="inner",left_on="Peptide",right_on="rPeptide")
        if len(rf):
            rf=rf[["Peptide","rPeptide"]]
            ps=set(rf["Peptide"].tolist())
            rps=set(rf["rPeptide"].tolist())
            tp=rf["Peptide"].str[-1].isin(["K","R"]).sum()
            trp=rf["rPeptide"].str[-1].isin(["K","R"]).sum()
            if tp==trp: rm=list(ps)
            if trp>tp:  rm=list(tp-trp)
            if trp<tp:  rm=list(trp-tp) 
            pepdf=pepdf[~pepdf["Peptide"].isin(rm)]
        
        #add decoy
        d=pepdf.copy()
        d["Peptide"]=d["Peptide"].str[::-1]
        pepdf["Alignment_Decoy"],d["Alignment_Decoy"]=False,True
        pepdf=pd.concat([pepdf,d])

        
        
    #parse fasta header
    for i in header_info: #round numeric header info to save space
        try:
            pepdf[i]=pd.to_numeric(pepdf[i]).round(2)
        except:
            pass
    hdict=pepdf[[i for i in header_info if i in pepdf.columns ]].fillna("").astype(str) 
    pepdf["dict"]="{"+("'"+hdict.columns+"'"+":"+hdict).apply(",".join,axis=1)+"}" 
    
    pepdf=pepdf[["Peptide","dict"]]
    if unique_peptides: pepdf=pepdf.drop_duplicates()
    
    with open(out_path,"w") as f:
        pass
    
    with open(out_path,"a") as f: #create target
        f.write("\n".join(">"+pepdf["Peptide"]+";"+pepdf["dict"]+"\n"+pepdf["Peptide"])+"\n")

    return out_path
    
@passed_kwargs()    
def BlastP_alignment(*,
                     blast_folderpath=blast_folderpath,
                     word_size=5, #mimum word size
                     evalue=100, #maximum evalue
                     max_target_seqs=10, #max number of sequences returned
                     **kwargs):
    
    v=BlastP_alignment.vars #parse function arguments
    #required args
    input_file=v.input_file
    database_path=v.database_path
    check_required_args(v,["input_file","database_path"])
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    output_path=str(Path(output_folder,"blastp_alignment.txt"))
    command=str(Path(blast_folderpath,"blastp"))
    command+="".join([" -matrix PAM30 -gapopen 5 -gapextend 2 ", 
                      " -evalue "+str(evalue),
                      " -word_size "+str(word_size),
                      " -max_target_seqs "+str(max_target_seqs),
                      " -query "+'"'+input_file+'"',
                      " -db "+'"'+database_path+'"',
                      " -out "+'"'+output_path+'"',
                      ' -outfmt "6 qseqid sseqid pident bitscore sseq"'])
    
    print(command)
    stddout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    return output_path

@passed_kwargs()
def Diamond_alignment(*,
                      select=" -k25 ", #-top or  -k + integer (see diamond docs)
                      block_size=5,
                      index_chunks=1,
                      minimum_pident=80,
                      minimum_coverage=80,
                      minimum_bitscore=20,
                      other_args=" --algo ctg --dbsize 1 ",
                      diamond_output_columns=diamond_output_columns,
                      diamond_filepath=diamond_filepath,
                      **kwargs
                      ):
    
    v=Diamond_alignment.vars #parse function arguments
    
    #required args
    input_file=v.input_file
    database_path=v.database_path
    check_required_args(v,["input_file","database_path"])

    #default args
    select=v.select
    block_size=v.block_size
    index_chunks=v.index_chunks
    minimum_pident=v.minimum_pident
    minimum_coverage=v.minimum_coverage
    minimum_bitscore=v.minimum_bitscore
    other_args=v.other_args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    if type(input_file)==list:
        input_file=input_file[0]


    output_file=str(Path(output_folder,Path(input_file).stem+".tsv"))
    

    command="cd "+'"'+basedir +'"'+ " && " + \
            "".join([diamond_filepath,   
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
            " -f 6  qseqid "+" ".join(diamond_output_columns)+" ", #add qseqid?
            " -t "+'"'+tmp_folder+'"'+other_args])
    
    print(command)
    
    #do this via a "bat" file because of diamond bug that does not want to work with custom matrices
    batfile=str(Path(basedir,"alignment.bat"))
    with open(batfile,"w") as bat:
        bat.write("#!/bin/bash"+"\n"+command)

    stdout, stderr =subprocess.Popen(batfile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    shutil.move(str(Path(basedir,"diamond.log")), str(Path(output_folder,Path(input_file).stem+".log")))

    return output_file

   

@passed_kwargs()
def Write_alignment_to_database(*,output_file="target.fa",unique_accs=True,**kwargs):
    
    v=Write_alignment_to_database.vars #parse function arguments
    
    #required args
    input_file=v.input_file
    check_required_args(v,["input_file"])
    
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    output_file=v.output_file
    unique_accs=v.unique_accs
    
    iterable=Diamond_alignment_Reader(input_file)

    output_path=str(Path(output_folder,output_file))
    print(output_file)
    print(output_path)
    seen=set()  
   
    with open(output_path,"w") as f:
    
        for ix,batch in enumerate(iterable):
            
            if unique_accs:
                batch=batch[~batch.sseqid.isin(seen)]
                seen.update(batch.sseqid.tolist())
            
            batch.full_sseq=batch.full_sseq.str.replace("*","",regex=False) #remove protein ends
            f.write("\n".join((">"+batch.stitle+"\n"+batch.full_sseq))+"\n")

    return output_path

    
@passed_kwargs()
def add_proteins_SMSNet(*,
                    
                    #SMSnet params
                    simple_unmask=True,     #attempts to solve low complexity masked SMSNet regions
                    X_padding=False,        #filla remaining masks with X
                    SMSNet_ppm=False,       #max ppm tolerance
                    SMSNet_minscore=False,  #mininum mean peptide score
                    
                    Diamond=True,
                    max_entries=500000,     #if final db is larger than this, skip Exact_tag and BlastP because of performance loss
                    Exact_tag=True,
                    BlastP=True,
                    min_tag_length=5,
                    max_tag_length=10,
                    score_cutoff=0.95,
                    minimum_pident=80,
                    minimum_coverage=80,
                    minimum_bitscore=20,
                    max_no_proteins=10, #dont report tags with more than this nr of proteins, this reduces output file size and number of decoy matches.
                    
                    **kwargs):
 

    
    v=add_proteins_SMSNet.vars #parse function arguments
    
            
    #required arguments
    input_files=v.input_files 
    database_path=v.database_path
    check_required_args(v,["input_files","database_path"])

    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    simple_unmask=v.simple_unmask
    X_padding=v.X_padding
    SMSNet_ppm=v.SMSNet_ppm
    SMSNet_minscore=v.SMSNet_minscore

    Diamond=v.Diamond
    Exact_tag=v.Exact_tag
    BlastP=v.BlastP
    min_tag_length=v.min_tag_length
    max_tag_length=v.max_tag_length
    score_cutoff=v.score_cutoff
    minimum_pident=v.minimum_pident
    minimum_coverage=v.minimum_coverage
    max_no_proteins=v.max_no_proteins #dont report tags with more than this nr of proteins, this reduces output file size and number of decoy matches.

    ### Add proteins ###

    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith("SMSNET.tsv")]
        else:
            input_files=input_files.split()

    SMSNet_peplist=write_to_Diamond_fasta(input_files=input_files,X_padding=False)
    tags=pd.DataFrame([[r.id,str(r.seq)] for r in Bio.SeqIO.parse(SMSNet_peplist,format="fasta")],columns=["id","tag"])
    matched=[]
    
    #check entries
    with open(database_path, "rbU") as f:
        entries = sum(1 for _ in f)/2
    if entries>max_entries:
        print("final database too large!!, skipping Exact tag matching and BlastP")
        Exact_tag=False
        BlastP=False
    
    if Diamond:
    
        print("Diamond alignment")

        final_database_dmnd=make_diamond_database(input_file=database_path) #this is a diamond databse with decoy proteins
        Alignment=Diamond_alignment(input_file=SMSNet_peplist,
                                    database_path=final_database_dmnd,
                                    minimum_pident=0, 
                                    minimum_coverage=0,
                                    minimum_bitscore=0,
                                    other_args=" --algo ctg --masking 0 --dbsize 1 ",
                                    #other_args=" --gapopen 5 --gapextend 30 --algo ctg --masking 0 --dbsize 1 --custom-matrix "+'"'+custom_matrix_path+'"',
                                    diamond_output_columns=["qseqid","sseqid","pident","bitscore","sseq"])
        
        
        al=pd.concat([i for i in Diamond_alignment_Reader(Alignment,
                                                          diamond_output_columns=["qseqid","sseqid","pident","bitscore","sseq"],
                                                          score_cutoff=score_cutoff)])
        tags=tags[~tags.tag.isin(al["tag"])]
        
        
        
        
        al["Proteins"]=al["sseqid"]
        al["coverage"]=al["sseq"].apply(len).divide(al["tag"].apply(len),axis=0)*100
        al["pident"]=al["pident"]
        al=al[al["pident"]>=minimum_pident]
        al=al[al["coverage"]>=minimum_coverage]
        al=al[al["bitscore"]>=minimum_bitscore]
        
        #top score filter
        al["score"]=al["pident"]*al["coverage"]*al["bitscore"]
        b=al.set_index("tag")["score"]
        al=al[b.index.isin(((b/al.groupby("tag")["score"].max())>=score_cutoff).index)]
        
        m=al.groupby("tag").agg({"pident":max,"coverage":max,"Proteins":" ".join}).reset_index()
        matched.append(m)
        

    if Exact_tag: #Exact tag should only be effective for short-mid length peptides (5-10)
    
        print("Exact substring matching")
        tags["len"]=tags["tag"].apply(len)
        unitags=tags[(tags.len<=max_tag_length) & (tags.len>min_tag_length)][["tag","len"]].drop_duplicates()
        other_tags=tags[~tags.tag.isin(unitags.index)]

        wins=np.sort(unitags["len"].unique())
        wunitags=[unitags.loc[unitags["len"]==w,"tag"].tolist() for w in wins]
        unitags=unitags.set_index("tag")
        unitags["No_Proteins"]=0
        unitags["Proteins"]=""
        unitags.pop("len")
        
        recs=SeqIO.parse(database_path,format="fasta")
        chunks=chunk_gen(recs)
        
        for ic,chunk in enumerate(chunks):
            df=pd.DataFrame([[r.id,str(r.seq)] for r in chunk],columns=["id","tag"])
            df.tag+="*"

            for counter,index_slice in enumerate(sliced(range(len(df)), 10000)):
                cdf = df.iloc[index_slice]
                ids=cdf.id
        
                print(counter)
                
                mw=cdf.tag.apply(lambda x: ["".join(w) for w in  mit.windowed(x,int(wins[0]))]).explode()
        
                prev_win=wins[0]
                cum_win=0
                m=mw[mw.isin(wunitags[0])].reset_index()
                
                for iw,win in enumerate(wins[1:]):

                    m["Proteins"]=ids.loc[m["index"].tolist()].values
                    r=m.groupby("tag").agg({"Proteins":[len,lambda x: " ".join(x)+" "]})
                    unitags.loc[r.index]+=r.values
            
                    delta_win=win-prev_win
                    cum_win+=delta_win

                    mw+=mw.str[-delta_win:].iloc[delta_win:].tolist()+[""]*delta_win
                    prev_win=win
                    m=mw[mw.isin(wunitags[iw])].reset_index()
            
        m=unitags[unitags["No_Proteins"]>0].reset_index()[["tag","Proteins"]]
        m["pident"]=100
        m["coverage"]=100
        matched.append(m)
        tags=pd.concat([other_tags,unitags[unitags["No_Proteins"]==0].reset_index()])
    
      
   
    if BlastP:
        
        print("BlastP alignment")
        utags=tags.tag.drop_duplicates()
        blastp_tags=str(Path(output_folder,"blastP_tags.fa"))
        with open(blastp_tags,"w") as f:
            f.write("\n".join(">"+utags+"\n"+utags)+"\n")
      
        blastp_db=make_blast_database(input_file=database_path)
        
        Alignment=BlastP_alignment(input_file=blastp_tags,database_path=blastp_db)
     

        al=pd.read_csv(Alignment,sep="\t",header=None,names=["tag","Proteins","pident","bitscore","sseq"])

        if len(al):
            al["coverage"]=(al["sseq"].apply(len)-al["sseq"].str.count("-").values).divide(al["tag"].apply(len),axis=0)*100
            al["pident"]=al["pident"]
            al=al[al["pident"]>=minimum_pident]
            al=al[al["coverage"]>=minimum_coverage]
            al=al[al["bitscore"]>=minimum_bitscore]
            
            #score filter
            al["score"]=al["pident"]*al["coverage"]*al["bitscore"]
            b=al.set_index("tag")["score"]
            al=al[b.index.isin(((b/al.groupby("tag")["score"].max())>=score_cutoff).index)]
            
            m=al.groupby("tag").agg({"pident":max,"coverage":max,"Proteins":" ".join}).reset_index()
            
            matched.append(m)
            

    matched_tags=pd.concat(matched)
    matched_tags["No_Proteins"]=1+matched_tags["Proteins"].str.strip().str.count(" ")
    matched_tags.loc[matched_tags["No_Proteins"]>max_no_proteins,"Proteins"]="unspecific"


    ### add back proteins
    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith("SMSNET.tsv")]
        else:
            input_files=input_files.split()
        
    output_paths=[]
    for input_file in input_files:
       
        
        pepdf=parse_SMSNet_output(input_file=input_file,simple_unmask=simple_unmask,X_padding=X_padding,SMSNet_ppm=SMSNet_ppm,SMSNet_minscore=SMSNet_minscore)
        pepdf=pepdf.merge(matched_tags,how="left",on="tag")
        output_path=str(Path(output_folder,Path(input_file).name))
        output_paths.append(output_path)
        pepdf.to_csv(output_path,sep="\t")
        
    return output_paths

@passed_kwargs()
def write_decoy(*,
                method="pseudo_random", #pseudo_random, reverse or random
                decoy_prefix="decoy_",
                output_file="decoy.fa",
                 
                #this section is for adding a scaling decoy form the initial database
                Sample=False,
                Initial_entries=False,**kwargs):
 
    v=write_decoy.vars #parse function arguments
    
    #required args
    input_file=v.input_file
    
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    output_file=v.output_file
    method=v.method
    decoy_prefix=v.decoy_prefix
    
    req_args=["input_file"]
    if v.Sample: req_args+=["Sample","Initial_entries"]
    check_required_args(v,req_args)
    
    output_path=str(Path(output_folder,output_file))
    
    #here a scaled decoy is added
    if Sample and Initial_entries:
        if Sample>Initial_entries:
            indices=list(set(random.choices(range(Initial_entries),k=Sample))) 
        else:
            indices=list(set(random.sample(range(Initial_entries),Sample))) 
        indices.sort()
        chunks=chunk_gen(slice_gen(SeqIO.parse(input_file, "fasta"),indices))
    
    else:
        chunks=chunk_gen(SeqIO.parse(input_file, "fasta"))
    
    with open(output_path,"w") as f: #clears file if exists
        pass
    with open(output_path,"a") as f:
    
        for chunk in chunks:
    
            df=pd.DataFrame([[rec.description,str(rec.seq)] for rec in chunk],columns=["description","seq"])
            df.description=decoy_prefix+df.description
            
            if method=="reverse":       df.seq=df.seq.str[::-1]
            if method=="pseudo_random": df.seq=df.seq.apply(pseudo_randomize)  #random staggered join (faster)    
            if method=="random":        df.seq=df.seq.apply(lambda x: ''.join(random.sample(x,len(x)))) #slow
        
            f.write("\n".join(">"+df.description+"\n"+df.seq)+"\n")
    
    return output_path
    


@passed_kwargs()
def write_database_composition(*,output_file="database_composition.tsv", **kwargs):
    
    v=write_database_composition.vars #parse function arguments
    
    #required args
    input_file=v.input_file
    check_required_args(v,["input_file"])
    taxdf=v.taxdf
      
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder


    recs=SeqIO.parse(input_file,format="fasta")
    chunks=chunk_gen(recs)
    
    output_path=str(Path(output_folder,output_file))

    entries=0
    taxa=[]
    
    for chunk in chunks:
        df=pd.DataFrame([[r.id,str(r.seq)] for r in chunk],columns=["id","seq"])
        taxa.extend(df.id.str.split("|").apply(lambda x: x[-1]))
        entries+=len(df)

    tax_counts=pd.DataFrame.from_dict(Counter(taxa),orient="index",columns=["Count"])
    tax_counts=tax_counts.merge(taxdf,how="left",left_index=True,right_index=True).dropna().sort_values(by="Count",ascending=False)
    
 
    tax_counts.to_csv(output_path,sep="\t")
           
    return tax_counts,len(tax_counts),entries #composition,"richness",total entries
    



@passed_kwargs()
def filter_Database_taxonomy(*,output_file="ft_target.fa",**kwargs):

    v=filter_Database_taxonomy.vars #parse function arguments
    
    #required args
    input_file=v.input_file
    taxids=v.taxids
    check_required_args(v,["input_file","taxids"]) 
    
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    output_path=str(Path(output_folder,output_file))
    

    #specifically designed for GTDB separate folder
    if os.path.isdir(input_file):

        fafiles=[]
        for root, dirs, files in os.walk(input_file, topdown=False):
            for file in files:
                if file.endswith(".fa") or file.endswith(".fasta"):
                    if "_".join(file.split("_")[:3]) in taxids:
                        fafiles.append(str(Path(root,file)))
            
        output_path=merge_files(fafiles,output_path=output_path)
    

    else:

        #standard for single database
        recs=SeqIO.parse(input_file,format="fasta")
        chunks=chunk_gen(recs)
        
        with open(output_path,"w") as f: #clears file if exists
            pass
        with open(output_path,"a") as f:
            
            for ic,c in enumerate(chunks):
                print("chunk "+ str(ic))
        
                chunk_df=pd.DataFrame([[str(r.seq),r.description,r.id] for r in c],columns=["seq","description","id"])
                chunk_df=chunk_df[chunk_df["id"].str.split("|").apply(lambda x: x[-1]).isin(taxids)]
    
                f.write("\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"])+"\n")

    return output_path
    
@passed_kwargs()
def filter_Database_proteins(*,output_file="fp_target.fa",**kwargs): #filters based on header

    v=filter_Database_proteins.vars #parse function arguments
    
    #required args
    input_file=v.input_file
    proteins=v.proteins 
    check_required_args(v,["input_file","proteins"]) 
    
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    output_path=str(Path(output_folder,output_file))
    

    #standard for single database
    recs=SeqIO.parse(input_file,format="fasta")
    chunks=chunk_gen(recs)
    
    with open(output_path,"w") as f: #clears file if exists
        pass
    with open(output_path,"a") as f:
        
        for ic,c in enumerate(chunks):
            print("chunk "+ str(ic))
    
            chunk_df=pd.DataFrame([[str(r.seq),r.description,r.id] for r in c],columns=["seq","description","id"])
            chunk_df=chunk_df[chunk_df["id"].isin(proteins)]

            f.write("\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"])+"\n")

    return output_path
    


@passed_kwargs()
def filter_Database_proteins_in_mem(*,output_file="fp_target.fa", #uses simple dataframe indexing (also writes db composition)
                            **kwargs):

    v=filter_Database_proteins_in_mem.vars #parse function arguments
    
    #required args
    input_file=v.input_file #this is a dataframe thats loaded in memory so technically not a file
    proteins=v.proteins
    check_required_args(v,["input_file","proteins"]) #taxids can be empty
    
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    output_path=str(Path(output_folder,output_file))
    
    input_file=input_file.loc[proteins,:]     
    
    
    cols=['OX', 'superkingdom', 'phylum', 'class', 'order','family', 'genus', 'species', 'Dump_taxid', 'OS']
    cols=[c for c in cols if c in input_file.columns]
    comp=input_file.groupby(cols).size().rename("Count").sort_values(ascending=False).reset_index()
    comp.to_csv(str(Path(output_folder,"database_composition.tsv")),sep="\t")
        
    
    with open(output_path,"w") as f: #clears file if exists
        pass
    with open(output_path,"a") as f:
        f.write("\n".join(">"+input_file["description"]+"\n"+input_file["seq"])+"\n")
    
    return input_file,output_path,len(comp),len(input_file) #mem_db, written_db, richess, entries
             

@passed_kwargs()
def filter_Database_taxonomy_in_mem(*,output_file="ft_target.fa", #uses simple dataframe indexing (also writes db composition)
                            **kwargs):

    v=filter_Database_taxonomy_in_mem.vars #parse function arguments
    
    #required args
    input_file=v.input_file #this is a dataframe thats loaded in memory so technically not a file
    proteins=v.proteins
    check_required_args(v,["input_file","proteins"]) #taxids can be empty
    
    #default args
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    output_path=str(Path(output_folder,output_file))
    
    input_file=input_file.loc[proteins,:]     
    
    
    cols=['OX', 'superkingdom', 'phylum', 'class', 'order','family', 'genus', 'species', 'Dump_taxid', 'OS']
    cols=[c for c in cols if c in input_file.columns]
    comp=input_file.groupby(cols).size().rename("Count").sort_values(ascending=False).reset_index()
    comp.to_csv(str(Path(output_folder,"database_composition.tsv")),sep="\t")
        
    
    with open(output_path,"w") as f: #clears file if exists
        pass
    with open(output_path,"a") as f:
        f.write("\n".join(">"+input_file["description"]+"\n"+input_file["seq"])+"\n")
    
    return input_file,output_path,len(comp),len(input_file) #mem_db, written_db, richess, entries
             


@passed_kwargs()
def load_full_db(Database,**kwargs): #load in memory (works only for small databases)

    v=load_full_db.vars #parse function arguments
    taxdf=v.taxdf

    recs=SeqIO.parse(Database,format="fasta")
    rdf=pd.DataFrame([[str(r.seq),r.description,r.id] for r in recs],columns=["seq","description","id"])
    
    rdf["OX"]=rdf.id.str.split("|").apply(lambda x: x[-1])
    rdf=rdf.merge(taxdf,how="left",left_on="OX",right_index=True).dropna()
    
    return rdf.set_index("id")

@passed_kwargs()
def make_blast_database(*,
       
                        blast_folderpath=blast_folderpath,
                        **kwargs):

     v=make_blast_database.vars #parse function arguments
    
     #required args
     input_file=v.input_file #this is a dataframe thats loaded in memory so technically not a file
     check_required_args(v,["input_file"]) 
     output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    
     output_path=str(Path(output_folder,"BlastPdb_"+Path(input_file).name))


     command=str(Path(blast_folderpath,"makeblastdb"))
     command+=" -in "+'"'+input_file+'"'+" -dbtype prot -out "+'"'+output_path+'"'
     print(command)
     stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


     return output_path



@passed_kwargs()
def make_diamond_database(*, #.fa
                          diamond_filepath=diamond_filepath,
                          output_file=False,
                          **kwargs):
    
    
    v=make_diamond_database.vars #parse function arguments
    
    #required args
    input_file=v.input_file #this is a dataframe thats loaded in memory so technically not a file
    check_required_args(v,["input_file"]) 
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    
    if output_file: output_path=str(Path(output_folder,Path(output_file).stem))
    else:           output_path=str(Path(output_folder,Path(input_file).stem))


    command="cd "+'"'+basedir +'"'+ " && "
    command+="diamond makedb --in "+'"'+input_file+'"' + " -d "+'"'+output_path+'"'
    print(command)
    stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()


    return output_path+".dmnd"

@passed_kwargs()
def refine_database(*,
                    
                    #MSFragger params
                    params_path=params_mid,
                    no_splits=2,
                    no_batches=1,    

                    #MSFragger score filters
                    max_evalue=10,
                    Top_score_fraction=0.9, #in case of multiple top candidates retain the top scoring fraction
                    
                    #Pre LCA filter
                    Frequency_prefilter=0, # 2   Static prefiler cutoff
                    Precision_prefilter=0, # 0.7 Target Decoy precision based denoising pre LCA filter

                    #focusing LCA parameters
                    weight_rank="species",
                    weight_cutoff=0.6,   
                    
                    #Post LCA representative picking
                    min_count=2,
                    min_ratio=0.99, 
                    prefilter_remove=False,
                    denoise_remove=True,
                    denoise_ranks=ranks.tolist(),
                    
                    #general parameters
                    min_rate=0.95,
                    output_folder="refine",
                    
                    **kwargs):

    #required files mzML files, database

    v=refine_database.vars #parse function arguments
    
    input_files=v.input_files     #mzML files
    database_path=v.database_path
    check_required_args(v,["input_files","database_path"]) #taxids can be empty

    #parse non-default arguments
    params_path=v.params_path
    no_splits=v.no_splits
    no_batches=v.no_batches
    max_evalue=v.max_evalue
    Top_score_fraction=v.Top_score_fraction
    Frequency_prefilter=v.Frequency_prefilter
    Precision_prefilter=v.Precision_prefilter
    weight_rank=v.weight_rank
    min_count=v.min_count
    min_ratio=v.min_ratio
    prefilter_remove=v.prefilter_remove
    denoise_remove=v.denoise_remove
    denoise_ranks=v.denoise_ranks
    min_rate=v.min_rate
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    taxdf=v.taxdf


    input_files=raw2mzML(input_files=input_files)     #make raw to mzml if supplied files are raw

    DB_in_mem=load_full_db(database_path)
    composition,richness,entries=write_database_composition(input_file=database_path)
    
    #variable name backup since cycle overwrites utarget,entries
    Initial_Database,Initial_entries=database_path,entries
    
    decoy=write_decoy(input_file=Initial_Database,output_folder=output_folder)
    database_path=merge_files([Initial_Database,decoy])

    cycle=0    
    while True:
    
        cycle+=1
        print("Starting cycle: "+str(cycle))
        folder_name=output_folder+"_"+str(cycle) 
        old_richness=richness
    
    
    
        MSFragger_files=MSFragger_annotation(input_files=input_files,
                                            database_path=database_path, 
                                            output_folder=folder_name,
                                            params_path=params_path, 
                                            no_splits=no_splits,
                                            no_batches=no_batches)
        
        prots=[]
        for ix,f in enumerate(MSFragger_files):
        
    
            pepdf=pd.read_csv(f,sep="\t")
            
            pepdf["u_ix"]=np.unique(pepdf["ScanNr"].astype(str)+"_"+pepdf["ExpMass"].astype(str),return_inverse=True)[1]
            pepdf["u_ix"]=pepdf["u_ix"].astype(str)+"_"+str(ix)
            pepdf[["hyperscore","log10_evalue"]]=pepdf[["hyperscore","log10_evalue"]].astype(float)
            pepdf["evalue"]=10**pepdf.log10_evalue
            if max_evalue:         pepdf=pepdf[pepdf["evalue"]<=max_evalue]
            if Top_score_fraction: pepdf=pepdf[(pepdf["evalue"]/pepdf.groupby("u_ix")["evalue"].transform('min'))<=(1/Top_score_fraction)]
        
        
            prots.append(pepdf[["u_ix","Proteins",'hyperscore']])
    
        prots=pd.concat(prots)
        prots["Proteins"]=prots["Proteins"].str.strip().str.split(" ")
        prots=prots.explode("Proteins")
        prots["Decoy"]=prots["Proteins"].str.startswith("decoy_")
        prots["OX"]=prots["Proteins"].str.replace("decoy_","").str.split("|").apply(lambda x: x[-1]) 
        
        #prots[ranks]=taxdf.loc[prots["OX"].tolist(),ranks].values 
        prots=prots.merge(taxdf,left_on="OX",right_index=True,how="inner") #deal with erroneous taxids  
 
        
    
        #filter proteins based on frequency and precision of their taxonomy
        g=prots.groupby([weight_rank,"Decoy"]).size()
        p=g.rename("Count").reset_index().pivot(index="species",columns="Decoy",values="Count").fillna(0)
        p.columns=["decoy" if i else "target" for i in p.columns]
        p["precision"]=p["target"]/(p["target"]+p["decoy"])
        p.to_csv(str(Path(folder_name,"precision.tsv")),sep="\t")
        rt=p[(p.target>=Frequency_prefilter) & (p.precision>=Precision_prefilter)].index #kept taxa
        rm=prots[prots[weight_rank].isin(rt)]
        if prefilter_remove: prots=rm
        else: prots=pd.concat([rm, prots[prots.u_ix.isin(list(set(prots.u_ix)-set(rm.u_ix)))]])

        target=prots[~prots["Decoy"]]
        decoy=prots[prots["Decoy"]]
        target=target.merge(p,how="inner",left_on=weight_rank,right_index=True)
        target["score"]=target['hyperscore']*target["precision"]
        tlca=weighted_lca(target,group_on="u_ix",weight_column="score",protein_column="Proteins",weight_cutoff=weight_cutoff)
        tlca.to_csv(str(Path(folder_name,"lca.tsv")),sep="\t")
        
        
        proteins,taxids=denoise_nodes(tlca,min_count=min_count,min_ratio=min_ratio,remove=denoise_remove,denoise_ranks=denoise_ranks)
        print("remaining proteins: "+str(len(proteins)))
        

        #DB_in_mem, target_db, richness, entries=filter_Database_proteins_in_mem(input_file=DB_in_mem,proteins=proteins,output_folder=folder_name)
        DB_in_mem, target_db, richness, entries=filter_Database_proteins_in_mem(input_file=DB_in_mem,proteins=proteins,output_folder=folder_name)
        
        
        decoy_db=write_decoy(input_file=Initial_Database,
                             output_folder=folder_name,
                             Sample=entries*cycle+1,
                             Initial_entries=Initial_entries)
        
        database_path=merge_files([target_db,decoy_db])

    
        decrease_ratio=richness/old_richness
        print("remaining taxa: "+str(richness))
        print("Database richness decrease fraction : "+str(decrease_ratio))
        if decrease_ratio>min_rate: 
            break

    
    return tlca,database_path,DB_in_mem

@passed_kwargs()
def Post_processing(*,
                    max_evalue=10,
                    Top_score_fraction=0.9,
                   
                    SMSNet_ppm=20,       #max ppm tolerance
                    SMSNet_minscore=False,  #minum mean peptide score
                    
                    FDR=0.05,
                    min_peptide_count=1,
                    remove_unannotated=False,

                    max_no_proteins=10,     #maximum different number of proteins linked to a single peptide
                    
                    mgf_files="",           #optionally add information from list of mgf files (adds intensity, charge and m/z)
                    database="",            #optionally add found protein output fasta file for further functional annotation  
                    
                    
                    **kwargs):
                    
    
    v=Post_processing.vars #parse function arguments
    
    #required arguments
    input_files=v.input_files 
    check_required_args(v,["input_files"])
    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    max_evalue=v.max_evalue
    Top_score_fraction=v.Top_score_fraction
    SMSNet_ppm=v.SMSNet_ppm
    SMSNet_minscore=v.SMSNet_minscore
    FDR=v.FDR
    min_peptide_count=v.min_peptide_count
    remove_unannotated=v.remove_unannotated
    max_no_proteins=v.max_no_proteins
    
    mgf_files=v.mgf_files
    database=v.database #database used for writing output protein fasta
    
    taxdf=v.taxdf
    
#     #%% Test
    
#     input_files=["C:/MP_CHEW/CHEW/SwissProt_Mix24dn/final/Q20518_Mix24X_SMSNET.tsv",
# "C:/MP_CHEW/CHEW/SwissProt_Mix24dn/final/Q20516_Mix24X_SMSNET.tsv",
# "C:/MP_CHEW/CHEW/SwissProt_Mix24dn/final/Q20517_Mix24X_160629055637_SMSNET.tsv"]
    
#     mgf_files="C:/MP_CHEW/CHEW/mgf"
    
    #### parsing input files 
    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith("SMSNET.tsv") or i.endswith(".pin")]
        else: input_files=input_files.split()

    input_files.sort()
    
    if len(mgf_files):
        if type(mgf_files)==str:
            if os.path.isdir(mgf_files):
                mgf_files=[str(Path(mgf_files,i)) for i in os.listdir(mgf_files) if i.endswith(".mgf")]
            else: mgf_files=mgf_files.split()
    
    
    pin_files=[i for i in input_files if i.endswith(".pin")]
    SMSNet_files=[i for i in input_files if i.endswith("SMSNET.tsv")]
    mgf_files=[i for i in mgf_files if i.endswith(".mgf")]
    
    input_files=pin_files+SMSNet_files+mgf_files
    
    if len(SMSNet_files):
        if len(pin_files)!=len(SMSNet_files):
            print("warning: unequal amount of SMSNet_files and MSFragger files detected")
            
    
    cols=['Peptide',"Proteins",
          'ScanNr','ExpMass','hyperscore',
          'log10_evalue',"Prediction","mean_score",
          "SMSNet_score", "evalue",'MassError(ppm)',
          "No_Proteins","pident","coverage"] #retained columns
          
    #outputs
    found_proteins=set()
    taxcounts=[]
    protcounts=[]
    unique_protcounts=[]
    taxints=[]
    protints=[]
    unique_protints=[]
    
    #make file pairs
    pairs=[]
    paired_up=set()
    for file in input_files:
        
        pair=[]
        s=Path(file).stem
        for file in input_files:
            if s in file:
                if file not in paired_up:
                
                    pair.append(file)
                    paired_up.update([file])
        if len(pair):
            pairs.append(pair)
        
    
    
    for pair in pairs:
        print(pair[0])
        pepdfs=[]
        pfile=""
        mgf_info=""
        for file in pair:

            if file.endswith(".mgf"):
                mgf_info=get_mgf_info(file)
            else:

                
                if file.endswith(".pin"):
                    pepdf=read_pin(file)
                    pfile=file
 
                else:
                    pepdf=pd.read_csv(file,sep="\t")
                    
        
                c=pepdf.columns
                if "Peptide" in c: pepdf["Peptide"]=pepdf["Peptide"].apply(lambda x: re.sub("[\[\[].*?[\]\]]", "", x)).str.split(".").apply(lambda x: x[1]).str.replace("I","L").str.replace("J","L") #remove ptms in peptides
                if "Proteins" in c: pepdf["Proteins"]=pepdf["Proteins"].str.replace("\t"," ",regex=True) #shouldnt Proteins alsways be a column?
                if "No_Proteins" not in c and "Proteins" in c : pepdf["No_Proteins"]=1+pepdf["Proteins"].str.strip().str.count(" ") 

                if "tag" in c: pepdf["Peptide"]=pepdf["tag"].str.replace("I","L").str.replace("J","L") #remove ptms in peptides
                if "Scores" in c: pepdf["SMSNet_score"]=pepdf["Scores"].str.rsplit(";",expand=True).astype(float).sum(axis=1)
                if "ScanNum" in c: pepdf["ScanNr"]=pepdf["ScanNum"]
                if 'ObservedM+H' in c: pepdf["ExpMass"]=pepdf['ObservedM+H']-1.007276466621 #-H+
                if  "log10_evalue" in c: pepdf["evalue"]=10**pepdf.log10_evalue.astype(float)
                
                c=pepdf.columns
                c=[i for i in c if i!="index" and "Unnamed:" not in i]
                
                
                pepdfs.append(pepdf[[i for i in cols if i in c]].reset_index(drop=True).drop_duplicates())
                
        if not(len(pepdfs)):
            continue
        
        pepdfs=pd.concat(pepdfs,axis=0,ignore_index=True).reset_index(drop=True)
        

        c=pepdfs.columns
        pepdfs["ScanNr"]=pepdfs["ScanNr"].astype(int)
        if "ExpMass" in c: pepdfs["ExpMass"]=pepdfs["ExpMass"].fillna(0).astype(float).round(2)
        if "No_Proteins" in c: 
            pepdfs["No_Proteins"]=pepdfs["No_Proteins"].fillna(0)
            pepdfs.loc[pepdfs["No_Proteins"].fillna(0)>max_no_proteins,"Proteins"]=""
            
            pepdfs["Decoy"]=pepdfs["Proteins"].str.count("decoy_")/pepdfs["No_Proteins"].fillna(0) #decoy fraction
            
        else:
            pepdfs["Decoy"]=pepdfs["Proteins"].str.contains("decoy_").fillna(False)
            
        
        if len(mgf_info): pepdfs=pepdfs.merge(mgf_info,on="ScanNr",how="left")
        
        for i in ["ExpMass","hyperscore","evalue"]:
            if i in c:
                pepdfs[i]=pepdfs[i].astype(float) 


    
        ##### De novo and non-denovo comparison
        if not pfile:
            pfile=file  
      
        basename=Path(pfile).stem
        #basepath=str(Path(Path(pfile).parents[0],basename))  
        basepath=str(Path(output_folder,basename))  
    
        if ("hyperscore" in c) and ("SMSNet_score" in c):
        
            mf_pepdf=pepdfs.loc[pepdfs.hyperscore.notnull(),["Peptide","ScanNr","hyperscore"]]
            sm_pepdf=pepdfs.loc[pepdfs.SMSNet_score.notnull(),["Peptide","ScanNr","SMSNet_score"]]
            shared_peptides=mf_pepdf.merge(sm_pepdf,on=["Peptide","ScanNr"],how="inner").drop_duplicates()
            shared_peptides=shared_peptides.groupby('ScanNr').apply(lambda x: x.sort_values(by="hyperscore",ascending=False).iloc[0]).reset_index(drop=True)
    
            x=shared_peptides["SMSNet_score"].values.reshape(-1,1)
            y=shared_peptides["hyperscore"].astype(float).values
            
            #lasso regression
            X_train, X_test, y_train, y_test = train_test_split(x,y,test_size=0.10,random_state=42)
            lasso_model = Lasso().fit(X_train,y_train)
            y_pred = lasso_model.predict(X_test)
            r2=round(r2_score(y_test, y_pred),2)
            
    
            pepdfs.loc[pepdfs["hyperscore"].isnull(),"hyperscore"]=lasso_model.predict(pepdfs["SMSNet_score"].dropna().values.reshape(-1,1))
    
            #plotting 
            fig,ax=plt.subplots()
            plt.xlabel("SMSNet_score")
            plt.ylabel("hyperscore")
            plt.scatter(x,y,s=0.1)
            plt.title(Path(pfile).stem+" "+"r2: "+str(r2))
            plt.savefig(basepath+"_dn_corr.png",dpi=400)
    
            #### plot high scoring ####
              
    
            
            #assign categories
            pepdfs["cat"]="MSFragger_unique"
            q=(pepdfs["Peptide"]+pepdfs["ScanNr"].astype(str)).isin(shared_peptides["Peptide"]+shared_peptides["ScanNr"].astype(str))
            pepdfs.loc[q,"cat"]="Shared"
            pepdfs.loc[~q & (pepdfs["SMSNet_score"].notnull()) & (pepdfs["No_Proteins"]!=0),"cat"]="SMSNet_unique_matched"
            pepdfs.loc[~q & (pepdfs["SMSNet_score"].notnull()) & (pepdfs["No_Proteins"]==0),"cat"]="SMSNet_unique_unmatched"
            pepdfs.loc[pepdfs["Proteins"].fillna("").str.contains("decoy_"),"cat"]="Decoy"
            ucats=pepdfs["cat"].unique()
            cats=["MSFragger_unique","Decoy","Shared","SMSNet_unique_unmatched","SMSNet_unique_matched"]
            cats=[i for i in cats if i in ucats]
            
            pepdfs["mass_corr_hyperscore"]=pepdfs["hyperscore"]/pepdfs["ExpMass"]
            colors=[sns.color_palette()[0]]+["#9FA5A9"]+sns.color_palette()[1:4]
            metric="mass_corr_hyperscore"
            x=round(pepdfs[metric].max()/20,3)
            bins=np.arange(0,pepdfs[metric].max()+x,x) #maybe instead bin this into quantiles?
            bw=0.001
            d=pepdfs[[metric,"cat"]]
            d["bin"]=np.digitize(d[metric],bins)
            
            g=d.groupby(["bin","cat"]).size()
            pv=pd.pivot_table(g.reset_index(),values=[0], index=['bin'],columns=['cat'])
            pv.columns=pv.columns.droplevel(0)
            pv=pv[cats[::-1]]
            pv.index=bins[pv.index]
    
            matched=pv[[i for i in cats if i in ["SMSNet_unique_matched","Shared"]]].sum(axis=1)
            total=pv[[i for i in cats if i in ["SMSNet_unique_matched","Shared","SMSNet_unique_unmatched"]]].sum(axis=1)
            q09=pv.index>d[metric].quantile(0.90)
            tqr=matched[q09].sum()/total[q09].sum()
            d=pepdfs[["cat",metric]].reset_index()
       
            fig,ax=plt.subplots()
            g=sns.histplot(data=d,x=metric,hue="cat",stat="count",multiple="stack",kde=False,binwidth=bw,
                          hue_order=cats,
                          palette=colors,
                          legend=True)
            
            #add scatter
            ax2 = ax.twinx()
            axd=matched/total
            axd=axd.fillna(method="ffill").fillna(0)
            line,=ax2.plot(axd.index, axd.values, 'b--')
            line.set_color("black")
            plt.text(ax.get_xticks()[-1],0.1,"SMSNet_09.score_match_ratio:"+str(round(tqr,3)))
            ax2.set_ylim([0,1])
            ax2.set_ylabel("SMSNet_match ratio",rotation=270,labelpad=-30)
            plt.title(Path(pfile).stem)
            plt.legend( loc='upper right',bbox_to_anchor=(1.6, 0.9), labels=['SMSNet_match ratio'])
            sns.move_legend(g, "upper left", bbox_to_anchor=(1.1, .7), title='category')
    
            plt.savefig(basepath+"dn_QC.png",dpi=300,bbox_inches="tight")
            pv.to_csv(basepath+"_dn_QC.tsv",sep="\t")
            
        #Filtering
        if remove_unannotated:                       pepdfs=pepdfs[pepdfs["No_Proteins"]!=0]
        if "evalue" in c and max_evalue:             pepdfs=pepdfs[(pepdfs["evalue"].isnull()) | (pepdfs["evalue"]<=max_evalue)]
        if "hyperscore" in c and Top_score_fraction: pepdfs=pepdfs[(pepdfs["hyperscore"]/pepdfs.groupby(["ScanNr","ExpMass"])["hyperscore"].transform('max')>=Top_score_fraction)]
        if 'MassError(ppm)' in c and SMSNet_ppm:     pepdfs=pepdfs[pepdfs["MassError(ppm)"].fillna(0).abs()<=SMSNet_ppm]
        if 'SMSNet_score' in c and SMSNet_minscore:  pepdfs=pepdfs[pepdfs['SMSNet_score'].fillna(0).abs()>=SMSNet_minscore]
        if 'SMSNet_score' in c and "ExpMass" in c:   pepdfs["mass_cor_denovo_score"]=pepdfs["SMSNet_score"]/pepdfs["ExpMass"].values
        
        
        
        if min_peptide_count: 
            s=pepdfs.groupby("Peptide").size()
            pepdfs=pepdfs[pepdfs["Peptide"].isin(s[s>=min_peptide_count].index)]
        
        c=pepdfs.columns
        if "hyperscore" in c: 
            FDR_col="hyperscore"
        elif "mass_cor_denovo_score" in c: 
            FDR_col="mass_cor_denovo_score"
        else: FDR_col="SMSNet_score"
       
      
        
        if FDR:
            pepdfs=pepdfs.sort_values(by=FDR_col,ascending=False)
            pepdfs["FDR"]=pepdfs["Decoy"].cumsum()/(pepdfs["Proteins"]!="").cumsum()
            
            #write hits
            dec=pepdfs["Decoy"].sum()
            t=(1- pepdfs["Decoy"]).sum()#-len(pepdfs["Proteins"]=="")
      
            pepdfs["passed FDR"]=pepdfs["FDR"]<=FDR
            
            
            du=pepdfs[pepdfs["passed FDR"]]["Decoy"].sum()
            tu=(1- pepdfs[pepdfs["passed FDR"]]["Decoy"]).sum()#-len(pepdfs[ (pepdfs["passed FDR"]) & (pepdfs["Proteins"]=="")])
            ddf=pd.DataFrame([[dec,t],[du,tu]],columns=["Decoy","Target"],index=["pre","post"])
            ddf.to_csv(basepath+"_performance.tsv",sep="\t")
    
        #vectorized lca     
        edf=pepdfs.copy()
        edf["Proteins"]=edf["Proteins"].str.split()
        edf=edf.explode("Proteins")
        edf["taxids"]=edf["Proteins"].fillna("|").str.split("|").apply(lambda x: x[-1])
        
        edf=edf.merge(taxdf,left_on="taxids",right_index=True,how="inner") #deal with erroneous taxids  
        
        #edf[taxdf.columns]=taxdf.loc[edf["taxids"].tolist()].values 
        lins=edf[ranks].reset_index().fillna("")
        
       
        pepdfs["No_Proteins"]=pepdfs["No_Proteins"].fillna(0)
        pepdfs[ranks]=lins.groupby("index")[ranks].nth(0)[(lins.groupby("index")[ranks].nunique()==1)].fillna("")   
        pepdfs=pepdfs.fillna("")
        pepdfs.to_csv(basepath+"_PSMs.tsv",sep="\t")
        
        found_proteins.update(pepdfs.Proteins.str.split().explode().tolist())
        
        #make sure only unique counts are stored (in case of hybrid annotation shared PSMs would be counted twice)
        pepdfs=pepdfs.groupby(["ScanNr","ExpMass"]).nth(0)  
        edf=pepdfs.copy()
        edf["Proteins"]=edf["Proteins"].str.split()
        edf=edf.explode("Proteins")
        
        #store here only those above FDR
        protcounts.append(edf[edf["passed FDR"]].groupby("Proteins").size().rename(basename))
        unique_protcounts.append(edf[(edf["No_Proteins"]==1) & (edf["passed FDR"])].groupby("Proteins").size().rename(basename))
        taxcounts.append(pepdfs[pepdfs["passed FDR"]].groupby(ranks.tolist()).size().rename(basename))
        
        if "intensity" in c:
            protints.append(edf[edf["passed FDR"]].groupby("Proteins")["intensity"].sum().rename(basename))
            unique_protints.append(edf[(edf["No_Proteins"]==1) & (edf["passed FDR"])].groupby("Proteins")["intensity"].sum().rename(basename))
            taxints.append(pepdfs[pepdfs["passed FDR"]].groupby(ranks.tolist())["intensity"].sum().rename(basename))
    
    #write combined outputs
    cnames=["tax_counts","prot_counts","unique_prot_counts","tax_ints","prot_ints","unique_prot_ints"]
    for ix_c,i in enumerate([taxcounts,protcounts,unique_protcounts,taxints,protints,unique_protints]):
        if len(i):
            pd.concat(i,axis=1).to_csv(str(Path(output_folder,cnames[ix_c]+".tsv")),sep="\t")
    
    
    if len(database):
        print("writing found proteins")
        filter_Database_proteins(input_file=database,proteins=found_proteins,output_file="found_proteins.fa")



def compare_dbs(input_file, #aligned fasta
                db1,  #diamond database 1
                db2): #diamond database 2
    
    
    
    dbs=[db1,db2]
    results=[]
    
    for db_ix,db in enumerate(dbs):
        if not db.endswith(".dmnd"): db=make_diamond_database(input_file=db)
        Alignment=Diamond_alignment(input_file=input_file,
                                    database_path=db)
        iterable=Diamond_alignment_Reader(Alignment,output_columns=diamond_output_columns)
        al=pd.concat([i for i in iterable])
        results.append(al)

        
    rs=[]
    for r in results:
        i=r.copy()
        i["Alignment_Decoy"]=i.qseqid.str.contains(":'Decoy'")
        i=i[["Peptide","sseqid","bitscore","Alignment_Decoy"]].drop_duplicates()
        i.loc[i["Alignment_Decoy"],"Peptide"]=i.loc[i["Alignment_Decoy"],"Peptide"].str[::-1]
        i=i.sort_values(by=["Peptide","Alignment_Decoy","bitscore"],ascending=False)
        x=i.groupby(["Peptide","Alignment_Decoy"],sort=False)["bitscore"].nth(0).reset_index()
        x["count"]=i.groupby(["Peptide","Alignment_Decoy"],sort=False).size().values
        rs.append(x.set_index(["Peptide","Alignment_Decoy"]))
    
    rs=pd.concat(rs,axis=1)
    rs.columns=["db1_score","db1_count","db2_score","db2_count"]
    rs=rs.reset_index()

    return rs 
    


       
def database_qc_per_score(pepdf,scoring_metric,bins=10):
    
    pepdf["quantile"]=pd.cut(pepdf[scoring_metric],bins)

    qdf=[]
    for n,g in pepdf.groupby("quantile"):
        
        c=pepdf[pepdf[scoring_metric]>n.left]
    
        cf=rs[rs.Peptide.isin(c.Peptide)]
        dc=cf[cf["Alignment_Decoy"]]
        tc=cf[~cf["Alignment_Decoy"]]
        completeness=tc.db2_score.sum()/tc.db1_score.sum()
        precision=(tc.db2_score.sum()/(tc.db2_score.sum()+dc.db2_score.sum()))/(tc.db1_score.sum()/(tc.db1_score.sum()+dc.db1_score.sum()))
        redundancy=tc.db1_count.mean()/tc.db2_count.mean()
        
        qdf.append([n.left,len(c),len(tc),completeness,precision,redundancy])
        qdf=pd.DataFrame(qdf,columns=["minimum "+scoring_metric,"peptide_count","matched_peptide_count","completeness ratio","precision ratio","redundancy ratio"]).fillna(0)
    
        return qdf

@passed_kwargs()
def database_QC_from_peplist(*,
                             placeholder="",
                             **kwargs):
    
    v=database_QC_from_peplist.vars #parse function arguments
    
    #required arguments
    input_file=v.input_file 
    db1=v.database_1
    db2=v.database_2
    check_required_args(v,["input_file","database_1","database_2"])


    output_folder,tmp_folder=v.output_folder,v.tmp_folder

    peplist=parse_peplist(input_file)
    rs=compare_dbs(peplist,db1,db2)
    
    rs.columns=["Peptide","Alignment_Decoy",Path(db1).stem+"_bitscore",Path(db1).stem+"_count",Path(db2).stem+"_score",Path(db2).stem+"_count"]
    
    rs.to_csv(str(Path(output_folder,"database_comparison_peptides.tsv")),sep="\t")
    
@passed_kwargs()
def database_QC_from_CHEW_PSMs(
                                *,
                               scoring_metric='mass_corr_hyperscore',
                               peplist="",
                               **kwargs):
    
    
    v=database_QC_from_CHEW_PSMs.vars #parse function arguments
    
    
    
    #required arguments
    input_files=v.input_files 
    check_required_args(v,["input_files","database_1","database_2"])
    
    input_files=v.input_files
    db1=v.database_1
    db2=v.database_2
    scoring_metric=v.scoring_metric
    peplist=v.peplist
    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    
    if len(peplist) and os.path.exists(peplist):
        peplist=parse_peplist(peplist)        
    else:            
        peplist=write_to_Diamond_fasta(input_files)
    rs=compare_dbs(peplist,db1,db2)
    
    #CHEW PSMs specific files 
    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith("_PSMs.tsv")]
        else: input_files=input_files.split()


    for input_file in input_files:
        
        pepdf=pd.read_csv(input_file,sep="\t")
        if "Peptide" not in pepdf.columns and "tag" in pepdf.columns: pepdf["Peptide"]=pepdf["tag"]
        qdf=database_qc_per_score(pepdf,scoring_metric=scoring_metric)
    
    rs.columns=["Peptide","Alignment_Decoy",Path(db1).stem+"_bitscore",Path(db1).stem+"_count",Path(db2).stem+"_score",Path(db2).stem+"_count"]
    
    rs.to_csv(str(Path(output_folder,"database_comparison_peptides.tsv")),sep="\t")
    qdf.to_csv(str(Path(output_folder,"database_comparison_per_score.tsv")),sep="\t")
        
    
    
    
@passed_kwargs()
def database_QC_from_SMSNet(
                                *,
                               peplist="",
                               **kwargs):
    
    v=database_QC_from_SMSNet.vars #parse function arguments
    
    
    
    #required arguments
    input_files=v.input_files 
    check_required_args(v,["input_files","database_1","database_2"])
    
    input_files=v.input_files
    db1=v.database_1
    db2=v.database_2

    peplist=v.peplist
    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    
    
    if len(peplist) and os.path.exists(peplist):
        peplist=parse_peplist(peplist)        
    else:            
        peplist=write_to_Diamond_fasta(input_files)
    rs=compare_dbs(peplist,db1,db2)
    
    #CHEW PSMs specific files 
    if type(input_files)==str:
        if os.path.isdir(input_files):
            input_files=[str(Path(input_files,i)) for i in os.listdir(input_files) if i.endswith("_SMSNet.tsv")]
        else: input_files=input_files.split()


    for input_file in input_files:
        
        pepdf=parse_SMSNet_output(t,simple_unmask=True,X_padding=False,SMSNet_ppm=False,SMSNet_minscore=False)
        pepdf["mean_score"]=pepdf.Scores.str.rsplit(";",expand=True).fillna(0).astype(float).mean(axis=1)
        if "Peptide" not in pepdf.columns and "tag" in pepdf.columns: pepdf["Peptide"]=pepdf["tag"]
        
        qdf=database_qc_per_score(pepdf,scoring_metric="mean_score")
    
    rs.columns=["Peptide","Alignment_Decoy",Path(db1).stem+"_bitscore",Path(db1).stem+"_count",Path(db2).stem+"_score",Path(db2).stem+"_count"]

    rs.to_csv(str(Path(output_folder,"database_comparison_peptides.tsv")),sep="\t")
    qdf.to_csv(str(Path(output_folder,"database_comparison_per_score.tsv")),sep="\t")
    
    
    #%%
    
    t="C:/MP_CHEW/CHEW/S56dn/peplist/S05_SMSNET.tsv"
    r=parse_SMSNet_output(t,simple_unmask=True,X_padding=False,SMSNet_ppm=False,SMSNet_minscore=False)
    
    #%%
    
#%%
@passed_kwargs()
def database_QC(initial_alignment="", #optional
                final_alignment="",   #optional
                
                initial_peplist="",   #required
                final_peplsit="" ,    #required
                
                initial_database="",  #required
                final_database="",    #required
                
                SMSnet_files=""       #required 
                ):
    
    
    v=database_QC.vars #parse function arguments
    
    #required arguments
    input_files=v.input_files 
    check_required_args(v,["initial_peplist",
                           "final_peplist",
                           "initial_database",
                           "final_database",
                           "SMSNet_files"])
    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    


#%%  Test

    initial_alignment="" #C:/MP_CHEW/CHEW/SwissProt_Mix24dn/peplist/peplist.tsv"
    final_alignment=""  #"C:/MP_CHEW/CHEW/SwissProt_Mix24dn/final/peplist.tsv"
    
    output_folder="C:/MP_CHEW/CHEW/db_comp/dbqc_test"
    kws.update(output_folder=output_folder)
    peplist="C:/MP_CHEW/CHEW/Mix24dn/peplist/peplist.fa"
    initial_database="H:/Databases/Swiss-Prot/Swiss-Prot/uniprot_sprot_BacArch_NoAmb_NoDump_NoDesc_IJeqL_taxid.dmnd"
    final_database="C:/MP_CHEW/CHEW/SwissProt_Mix24dn/final/ft_target.fa"
    
    initial_database="H:/Databases/UniprotKB/UniprotKB_BacArch_NoAmb_NoDump_IJeqL_taxid.dmnd"
    final_database="C:/MP_CHEW/CHEW/UniprotKB_Mix24dn/final/ft_target.fa"
    
    SMSNet_files=["C:/MP_CHEW/CHEW/Mix24/Initial_annotation/Q20518_Mix24X_SMSNET.tsv",
"C:/MP_CHEW/CHEW/Mix24/Initial_annotation/Q20516_Mix24X_SMSNET.tsv",
"C:/MP_CHEW/CHEW/Mix24/Initial_annotation/Q20517_Mix24X_160629055637_SMSNET.tsv"]

    initial_database,final_database=[  make_diamond_database(input_file=i )  if not i.endswith(".dmnd") else i.replace(".dmnd","") for i in [initial_database, final_database]    ] #check if databases are diamond databases
    # if len(peplist):
    # elif:
    # else:
    #     "I"+1
        
    
    
    ### check if initial  alignment has been done with or without decoy (or not)
    
    initial_al=pd.DataFrame() #placeholder
    if len(initial_alignment):

        #read initial alignment
        al=Diamond_alignment_Reader(initial_alignment,diamond_output_columns=diamond_output_columns)
        al["dict"]=al.qseqid.str.split(";").apply(lambda x: x[1]).str.strip("{}") #strip does not have regex
        al[al["dict"].iloc[0,:].str.split(":").apply(lambda x: x[0::2]).tolist()]=al["dict"].str.split(":").apply(lambda x: x[1::2]) #clumsy dictionary eval
        initial_al=al.copy()


    #read inital peplist (required argument)
    if not is_fasta(peplist): peplist=parse_peplist(peplist)  
    pepdf=pd.DataFrame([[r.description,str(r.seq)] for r in SeqIO.parse(peplist,"fasta")],columns=["description","seq"])
    
    if "Target_Decoy" not in initial_al.columns:
        
        #write decoy sequences and align 
        if not pepdf.description.str.contains(":Decoy}").sum():
           
            d=pepdf.copy()
            if pepdf.description.str.endswith("}").sum(): #they have a dict   
               pepdf.description=pepdf.description.str.replace("}",",Target_Decoy:Target}",regex=False).str.replace("{,","{",regex=False)
               d.description=d.seq.str[::-1]+";"+d.description.str.split(";").apply(lambda x: x[-1]).str.replace("}",",Target_Decoy:Decoy}",regex=False).str.replace("{,","{",regex=False)
            else:
               pepdf.description=pepdf.description+";{Target_Decoy:Target}"    
               d.description=d.seq.str[::-1]+";{Target_Decoy:Decoy}"  
            d.seq=d.seq.str[::-1]
        
            if len(initial_al):
                peplist=str(Path(output_folder,"decoy_peplist.fa"))
                pepdf=d.copy()
            else:
                peplist=str(Path(output_folder,"target_decoy_peplist.fa"))
                pepdf=pd.concat([pepdf,d])
            
            with open(peplist,"w") as f:
                f.write("\n".join(">"+pepdf.description+"\n"+pepdf.seq)+"\n")
            
                
        Alignment=Diamond_alignment(input_file=peplist,database_path=initial_database, other_args=" --algo ctg --dbsize 1 ")
        al=pd.concat([i for i in Diamond_alignment_Reader(Alignment)])
        initial_al=pd.concat([initial_al,al])
        
    
    ### check if final alignment has been done with or without decoy (or not)
    
    final_al=pd.DataFrame() #placeholder
    if len(final_alignment):

        #read final alignment
        al=Diamond_alignment_Reader(final_alignment,diamond_output_columns=diamond_output_columns)
        al["dict"]=al.qseqid.str.split(";").apply(lambda x: x[1]).str.strip("{}") #strip does not have regex
        al[al["dict"].iloc[0,:].str.split(":").apply(lambda x: x[0::2]).tolist()]=al["dict"].str.split(":").apply(lambda x: x[1::2]) #clumsy dictionary eval
        final_al=al.copy()
    
    if "Target_Decoy" not in final_al.columns:
    
        Alignment=Diamond_alignment(input_file=peplist,database_path=final_database)
        al=pd.concat([i for i in Diamond_alignment_Reader(Alignment)])
        final_al=pd.concat([final_al,al])
    
    
    
    #%%
    rs=[]
    for i in initial_al, final_al:
        i["Alignment_Decoy"]=i.qseqid.str.contains(":Decoy")
        i=i[["tag","sseqid","bitscore","Alignment_Decoy"]].drop_duplicates()
        i.loc[i["Alignment_Decoy"],"tag"]=i.loc[i["Alignment_Decoy"],"tag"].str[::-1]
        i=i.sort_values(by=["tag","Alignment_Decoy","bitscore"],ascending=False)
        r=i.groupby(["tag","Alignment_Decoy"],sort=False)["bitscore"].nth(0).reset_index()
        r["count"]=i.groupby(["tag","Alignment_Decoy"],sort=False).size().values
        rs.append(r.set_index(["tag","Alignment_Decoy"]))
    
    rs=pd.concat(rs,axis=1)
    rs.columns=["initial_score","initial_count","final_score","final_count"]
    rs=rs.reset_index()
    #%%
    # fig,ax=plt.subplots()
    # rs["initial_score"].plot.hist(bins=20)
    # rs["final_score"].plot.hist(bins=20)
    
    fig,ax=plt.subplots()
    rs["initial_count"].plot.hist(bins=20)
    fig,ax=plt.subplots()
    rs["final_count"].plot.hist(bins=20)
    
    fig,ax=plt.subplots()
    d=rs[rs["Alignment_Decoy"]]
    d["initial_score"].plot.hist(bins=20)
    d["final_score"].plot.hist(bins=20)
    
    fig,ax=plt.subplots()
    t=rs[~rs["Alignment_Decoy"]]
    t["initial_score"].plot.hist(bins=20)
    t["final_score"].plot.hist(bins=20)
    
    #something is wrong with the counts!
    

    
    #overall
    
    print(t.final_score.sum()/t.initial_score.sum())
    #print((d.initial_score.sum()/t.initial_score.sum())/(d.final_score.sum()/t.final_score.sum())) #decoy fraction instead of decoy count
    #print(d.final_score.sum()/d.initial_score.sum()) #decoy fraction instead of decoy count
    print((t.final_score.sum()/(t.final_score.sum()+d.final_score.sum()))/(t.initial_score.sum()/(t.initial_score.sum()+d.initial_score.sum())))
    print(rs["initial_count"].mean()/rs["final_count"].mean())
    
    
    #edit this back in CHEW funs later!
    
    #read each SMSNet file
    #for each scoring bin
    
    
    #for top scoring....
    
    #%% Compare to SMSNet score
    input_files=["C:/MP_CHEW/CHEW/Mix24/Initial_annotation/Q20518_Mix24X_SMSNET.tsv",
    "C:/MP_CHEW/CHEW/Mix24/Initial_annotation/Q20516_Mix24X_SMSNET.tsv",
    "C:/MP_CHEW/CHEW/Mix24/Initial_annotation/Q20517_Mix24X_160629055637_SMSNET.tsv"]
    for input_file in input_files: 
        if input_file.endswith("SMSNET.tsv"): #parse SMSNet outputs
            pepdf=parse_SMSNet_output(input_file=input_file,
                                      SMSNet_ppm=False,
                                      SMSNet_minscore=False,
                                      simple_unmask=True,
                                      X_padding=False)



        pepdf=pepdf[pepdf["tag"].fillna("").apply(len)>=min_length]
        pepdf["tag"]=pepdf["tag"].str.replace("I","L").str.replace("J","L")
        
        
        break
    
    m=pepdf.merge(rs,on="tag",how="left")
    m["total_score"]=m["score"]
    
    #%% Do it with PSMs files instead
    
    input_files=["C:/MP_CHEW/CHEW/UniprotKB_Mix24dn/output/Q20516_Mix24X_PSMs.tsv",
"C:/MP_CHEW/CHEW/UniprotKB_Mix24dn/output/Q20517_Mix24X_160629055637_PSMs.tsv",
"C:/MP_CHEW/CHEW/UniprotKB_Mix24dn/output/Q20518_Mix24X_PSMs.tsv"]
    
    for input_file in input_files:
        break
    pepdf=pd.read_csv(input_file,sep="\t")
    # pepdf=simple_unmask_fun(pepdf)
    # m=pepdf.merge(rs,on="tag",how="left")
    
    #%%

    
    #DBQC peplist (more flexible input)
    
    
    #DBQC CHEW PSMs (more detailed output)
    
    
    
    #%% 
    
    
    
    #%% compute bitscore ratios between final and initial
    
    
    
    #for each seq count top bitscore
    
    
    #for each seq, count length in target and length in decoy
    #or for each seq count top bitscore in target and top bitscore in decoy
    #top bitscore decoy/ top bitscore target
    # if missing ->1
    
    
    
        
    #%%
    #read SMSNet files
    
    
    
    #scoring metric
    
    #if EXPMass


#%% Helper funs



def get_mgf_info(file):
    
    
    with open(file,"r") as f:
        t=f.readlines()



    scans=[]
    for i in t:
        i=i.split("\n")[0]
        
        
        if i.startswith("TITLE="):
            scan=i.split("scan=")[1]
        
        elif i.startswith("PEPMASS"):
       
            s=i.split("PEPMASS=")[1].split()
            if len(s)==1:
                s=s+[0]
                
            mass,intensity=s
        
        elif i.startswith("CHARGE="):
            charge=i.split("=")[1].split("+")[0]
            scans.append([scan,mass,intensity,charge])
            
    scandf=pd.DataFrame(scans,columns=["ScanNr","m/z","intensity","charge"]).astype(float)
    scandf["ScanNr"]=scandf["ScanNr"].astype(int)
        
    return scandf

def parse_SMSNet_output(input_file,simple_unmask,X_padding,SMSNet_ppm,SMSNet_minscore):
    
    pepdf=""
    if type(input_file)==str:
        pepdf=pd.read_csv(input_file,sep="\t")
    if isinstance(input_file, pd.DataFrame):
        pepdf=input_file
        
    pepdf["Prediction"]=pepdf["Prediction"].str.replace("I","L").str.replace("J","L") #equate I and J to L
    if SMSNet_ppm: pepdf=pepdf[pepdf['MassError(ppm)'].abs()<=SMSNet_ppm]
    if SMSNet_minscore: pepdf=pepdf[pepdf.Scores.str.rsplit(";",expand=True).astype(float).mean(axis=1)>=SMSNet_minscore]
    if simple_unmask: #unmask simple combinations(1 or 2)
        preds=pepdf["Prediction"].drop_duplicates() 
        
        masks=preds.str.replace(")","(",regex=False).str.split("(",regex=False).explode() 
        masses=list(set(masks[masks.str.contains(".",regex=False)].astype(str).tolist()))
        masses.sort()
        
        for m in masses:

            match=simple_masks[(simple_masks["upper"]>float(m)) & (simple_masks["lower"]<float(m))]
            combs=sum([["".join(p) for p in itertools.permutations(c)] for c in match["compositions"]],[])
            if len(combs):
                
                preds=preds.str.replace("("+m+")",'"],'+str(combs)+',["',regex=False)
            else:
                if X_padding:
                    preds=preds.str.replace("("+m+")",'"],'+str(["X"*math.ceil(float(m)//110)])+',["',regex=False)
                    
        preds='[["'+preds.str.lstrip('"],').str.rstrip('"[,')+'"]]'
        preds=preds.str.replace('[["[','[[',regex=False).str.replace(']"]]',']]',regex=False).apply(eval)
        preds=preds.apply(lambda x:list(itertools.product(*x))).explode().apply(lambda x: "".join(x)).rename("Prediction")
      
        pepdf=pepdf.rename(columns={"Prediction":"Original_Prediction"})
        pepdf=pepdf.merge(preds,left_index=True,right_index=True,how="right").reset_index()

    #select longest tag for alignment
    pepdf["tag"]=pd.Series([max(pep[::2],key=len) for pep in pepdf["Prediction"].str.upper().str.replace("(",")",regex=False).str.split(")",regex=False)],name="longest_tag").str.strip("X")

    return pepdf

def is_fasta(input_file):
    fasta=SeqIO.parse(input_file,"fasta")
    return any(fasta)


#make sure that mgf files are properly formatted and filter mgf scans to process only HQ scans with SMSNet
def reformat_mfg(mgf_file,
                 minimum_peaks=20,
                 max_charge=6,
                 minimum_ion_pairs=2):


    proton_mass=1.007276
    water_mass=18.01056
    with open(mgf_file,"r+") as f:
        lines=f.readlines()
        lines=[line.replace("Scan: ","scan=") if line.startswith("TITLE=") else line for line in lines] #add Scan and SEQ=
        

        batches=[]
        
        t,s,rt,m,c="TITLE=Run: run, Index: 1, scan=1\n","SEQ=\n","RTINSECONDS=0\n","PEPMASS=0 0\n","CHARGE=1+\n"
        batch=[]
        peaks=[]
        
        for line in lines:
            if    line.startswith("TITLE="):        t=line
            elif  line.startswith("SEQ="):          s=line
            elif  line.startswith("RTINSECONDS="):  rt=line   
            elif  line.startswith("PEPMASS="):      m=line    
            elif  line.startswith("CHARGE="):       c=line
            elif  line.startswith("BEGIN IONS"):    pass
            elif  line.startswith("END IONS"):      pass
            
            else: peaks.append(line)    
                
            if line.startswith("END IONS"):

               if len(peaks)>minimum_peaks:
                   charge=int(c.split("CHARGE=")[-1].strip("+\n"))
                   if not charge>max_charge:
               
                       #check for min fragment pairs 
                       mz=(float(m.split()[0][8:])-proton_mass)*charge 
                       ms=pd.Series(peaks).str.split().apply(lambda x: x[0]).astype(float).values
                       md=abs(ms.reshape(-1,1)-(mz-ms+proton_mass))
                       mdw=abs(md-water_mass)
                       ion_pairs=(sum(mdw<0.1).sum()+sum(mdw<0.1).sum())/2
                       
                       if ion_pairs>=minimum_ion_pairs:
                
                           batch.extend(["BEGIN IONS\n"]) 
                           for i in [t,s,rt,m,c]: batch.extend([i])
                           batch.extend(peaks)
                           batch.extend([line])
                           batches.extend(batch)
               
                    

               t,s,rt,m,c="TITLE=Run: run, Index: 1, scan=1\n","SEQ=\n","RTINSECONDS=0\n","PEPMASS=0 0\n","CHARGE=1+\n"
               peaks=[]
               batch=[] 
               
       
        f.truncate(0)
        f.seek(0)
        f.writelines(batches)




def read_pin(pinfile):
    with open(pinfile, "r") as f:
        lines=f.readlines()    
    header=lines[0].replace("\n","").split("\t")
    df=pd.DataFrame([i.replace("\n","").split("\t",len(header)-1) for i in lines[1:]],columns=header)

    return df[[i for i in header if i!=""]].drop_duplicates()


#read alignment file to generator and filter on top x% scoring
def Diamond_alignment_Reader(input_file,*,
                             diamond_output_columns=diamond_output_columns,
    
                              score_cutoff=0.9): 
    
    cdf=pd.read_csv(input_file, sep='\t', chunksize=IO_batch,names=diamond_output_columns) #read to generator
    sc=[] #dummy
    for ix,c in enumerate(cdf):
        

        
        print(ix)
        

        
        _,index=np.unique(c.qseqid,return_index=True)
        d=c.iloc[0:index.max()]
        if ix>0:
            d=pd.concat([sc,d])

        d.loc[:,"tag"]=d.qseqid.str.rsplit(";",expand=True).iloc[:,0]
        d=d[(d["bitscore"]/d.groupby("qseqid")["bitscore"].max()>=score_cutoff).tolist()]
        yield d

        sc=c.iloc[index.max():]
        
    #last one
    d=sc
    d.loc[:,"tag"]=d.qseqid.str.rsplit(";",expand=True).iloc[:,0]
    d=d[(d["bitscore"]/d.groupby("qseqid")["bitscore"].max()>=score_cutoff).tolist()]
    
    yield d
    

def fill_g(x):
    cons=[np.array(list(g)) for g in mit.consecutive_groups(np.argwhere(x==""))]
    if cons:
        for c in cons:   
            c=c.flatten().tolist()
            if c[-1]!=len(x)-1: 
                for i in c:
                    x[i]="gap_"+str(x[c[-1]+1])+"_"+ranks[i]
    return list(x)


def weighted_lca(df,*, #dataframe with at least a column called Peptide, and rank
        group_on="tag",
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


def denoise_nodes(df, #lca df 
                  min_count=2,
                  min_ratio=0.99, 
                  remove=False,
                  denoise_ranks=["phylum","class","order","family","genus","species"]):

    if type(min_ratio)!=list: min_ratio=[min_ratio]
    if len(min_ratio)<len(denoise_ranks):
        min_ratio+=[max(min_ratio)]*(len(denoise_ranks)-len(min_ratio)) #pad max
    
    df=df.copy()
    df["proteins"]=df["proteins"].str.split(", ")
    edf=df.explode("proteins").reset_index()
    edf["OX"]=edf["proteins"].str.split("|").apply(lambda x: x[-1])
    dfs=[df]
    
    for ir,r in enumerate(denoise_ranks[::-1]): #reverse order from specific to unspecific
        # print("denoising: "+r)
        taxids=[]
        
        for n,tu in edf.groupby(r):
            if len(tu)<min_count or n=="":
                continue
            
            if tu["OX"].nunique()==1:
                taxids.extend([tu["OX"].iloc[0]]) 
            else:
                g=pd.DataFrame(tu.groupby("OX")["u_ix"].apply(set))
                g["l"]=g["u_ix"].apply(len)
                g=g.sort_values(by="l",ascending=False)
                
                u=set()
                ls=[0]
                s=0
                for ix,i in enumerate(g["u_ix"]): 
                    u.update(i)
                    l=len(u)
                    ls.append(l)
                    d=ls[-1]-ls[-2]
                    
                    if d<min_count: #absolute  
                        break
                    s+=l
                    if 1-(d/s)>min_ratio[ir]: #relative
                        break
    
                taxids.extend(g.index[:ix].tolist()) 
        edf=edf[edf["OX"].isin(taxids)]
        dfs.append(edf.groupby("u_ix")["proteins"].apply(list).rename("proteins"+str(ir)))
        #print(edf["u_ix"].nunique())
    
    if remove:
        proteins=edf["proteins"]
        
    else:
        cdf=pd.concat(dfs,axis=1)
        protcols=[i for i in cdf.columns if i.startswith("proteins")]
        df["proteins"]=cdf[protcols].ffill(axis=1)[protcols[-1]] #if keep
        proteins=df.explode("proteins")["proteins"]
    
    
    taxids=proteins.str.split("|").apply(lambda x: x[-1])

    return list(set(proteins)), list(set(taxids))


def chunk_gen(it,*kwargs):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//IO_batch):
        yield g
        

def slice_gen(gen,indices): #indices should be sorted
    ixx=0
    for ix,i in enumerate(gen):
        if ix==indices[ixx]:
            ixx+=1
            yield i
        if ixx==len(indices):
            break

def pseudo_randomize(x,steps=10): #not true random but stepped slicing (is faster)
    r=random.sample(range(steps),steps)
    return "".join(sum([list(x)[i::len(r)] for i in r],[]))

def merge_files(files,output_path=False):
    
    if not output_path:
        output_path=str(Path(Path(files[0]).parents[0],"target_decoy.fa"))
        
    with open(output_path,'wb') as o:
        for file in files:
            shutil.copyfileobj(open(file,'rb'), o)
            
    return output_path

#%%




