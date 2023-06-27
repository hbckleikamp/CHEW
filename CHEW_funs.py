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

#local module
from load_vars import *

#### function variables ###

from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0])) # set base path
basedir=os.getcwd()


#SMSnet filepaths
simple_masks=pd.read_csv(str(Path(basedir,"simple_unmasks.tsv")),sep="\t")

#MSFragger filepaths
MSFragger_jar_path=str(Path(basedir,"MSFragger-3.5.jar"))  
pep_split_path=str(Path(basedir,"msfragger_pep_split_HK.py"))
params_fast =str(Path(basedir,"closed_fragger_fast.params"))  #faster initial search
params_mid  =str(Path(basedir,"closed_fragger_mid.params"))   #medium search during refinement
params_final=str(Path(basedir,"closed_fragger_final.params")) #detailed search for smaller db



#Base taxonomy filepaths
ranks=np.array(["superkingdom","phylum","class","order","family","genus","species"]) 
taxdf_path=str(Path(basedir,"parsed_ncbi_taxonomy.tsv"))


IO_batch=10**6 #how much lines should be written or read at once from a fasta file 
diamond_output_columns=["qseqid","sseqid","stitle","bitscore","full_sseq"]





#kws="python scope is stupid"

#from raw2fasta_test import kws  
#kws=""


import glob
list_of_files = glob.glob(str(Path(basedir,"*.CHEW_params"))) 
latest_file   = max(list_of_files, key=os.path.getctime)
kws=load_variables(latest_file)


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
def raw2mzML(*,output_folder="mzML",**kwargs): 
    
    v=raw2mzML.vars #parse function arguments
    output_folder=v.output_folder
    input_files=v.input_files
    check_required_args(v,["input_files"])

    if type(input_files)==str:
        input_files=input_files.split()
    if type(input_files)==str:
        input_files=[input_files]
    
    output_files=[]
    for input_file in input_files:
        
        output_file=str(Path(output_folder,Path(input_file).stem+".mzML"))
        
        
        if not input_file.endswith(".mzML"):
            output_files.append(output_file)
            
            command="cd" +' "'+output_folder+'" && msconvert '
            command+='"'+input_file+'"' 
            command+=' --mzML --filter "peakPicking vendor" --filter "zeroSamples removeExtra" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>"'
            print(command)
            stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            
        else:    
            output_files.append(input_file)
            
    return output_files

@passed_kwargs()
def raw2mgf(*,output_folder="mgf",**kwargs):
    
    v=raw2mgf.vars #parse function arguments
    output_folder=v.output_folder
    input_files=v.input_files
    check_required_args(v,["input_files"])
    
    if type(input_files)==str:
        input_files=input_files.split()
    if type(input_files)==str:
        input_files=[input_files]
    
    output_files=[]
    for input_file in input_files:
        output_file=str(Path(output_folder,Path(input_file).stem+".mgf"))
        output_files.append(output_file)
        
        if not input_file.endswith(".mgf"):
            output_files.append(output_file)
        
            command="cd" +' "'+output_folder+'" && msconvert '
            command+='"'+input_file+'"' 
            command+=' --mgf --filter "peakPicking vendor" --filter "zeroSamples removeExtra" --filter "titleMaker Run: <RunId>, Index: <Index>, Scan: <ScanNumber>"'
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
        input_files=input_files.split()
    if type(input_files)==str:
        input_files=[input_files]
    
    input_files=[raw2mzML(input_file) if input_file.endswith(".raw") else input_file for input_file in input_files] #if raw, convert to mzML

    
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
def SMSnet_annotation(**kwargs):
    
    v=SMSnet_annotation.vars  #parse function arguments
    
    #required arguments
    input_files=v.input_files 
    check_required_args(v,["intput_file"])
    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder

    if type(input_files)==str:
        input_files=input_files.split()
    if type(input_files)==str:
        input_files=[input_files]
    
    input_files=[raw2mgf(input_file) if input_file.endswith(".raw") or input_file.endswith(".mzML")  else input_file for input_file in input_files] #if raw or mzML, convert to mzML
    
    output_files=[]
    for input_file in input_files:
        
        command="conda activate smsnet  && cd "+str(Path(basedir,"SMSNet-master"))+" && "
        command+="python run_HK.py --model_dir smsnet_phospho --inference_input_file "+input_file+" --rescore"
        print(command)
        stdout, stderr =subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        
        output_file=str(Path(output_folder,Path(input_file).stem+"_SMSNET.tsv"))
        shutil.move(str(Path(str(Path(input_file).parents[0])+"_output",'p-mod_fdr10.tsv')), output_file)

        output_files.append(output_file)
        #cleanup    
        shutil.rmtree(str(Path(input_file).parents[0])+"_output")
    return output_files

#Write MSfragger and PepNet files to fasta files fo Diamond alignment




def parse_SMSnet_output(input_file,simple_unmask,SMSnet_ppm,SMSnet_minscore):
    
    pepdf=""
    if type(input_file)==str:
        pepdf=pd.read_csv(input_file,sep="\t")
    if isinstance(input_file, pd.DataFrame):
        pepdf=input_file
        
    
    if SMSnet_ppm: pepdf=pepdf[pepdf['MassError(ppm)'].abs()<=SMSnet_ppm]
    if SMSnet_minscore: pepdf=pepdf[pepdf.Scores.str.rsplit(";",expand=True).astype(float).mean(axis=1)>=SMSnet_minscore]
    if simple_unmask: #unmask simple combinations(1 or 2)
        preds=pepdf["Prediction"].str.replace("I","L").str.replace("J","L").drop_duplicates() #equate I and J to L
        
        masks=preds.str.replace(")","(",regex=False).str.split("(",regex=False).explode() 
        masses=list(set(masks[masks.str.contains(".",regex=False)].astype(float).tolist()))
        masses.sort()
        
        for m in masses:
            match=simple_masks[(simple_masks["upper"]>m) & (simple_masks["lower"]<m)]
            combs=sum([["".join(p) for p in itertools.permutations(c)] for c in match["compositions"]],[])
            if len(combs): preds=preds.str.replace("("+str(m)+")",'"],'+str(combs)+',["',regex=False)
                
        preds='[["'+preds.str.lstrip('"],').str.rstrip('"[,')+'"]]'
        preds=preds.str.replace('[["[','[[',regex=False).str.replace(']"]]',']]',regex=False).apply(eval)
        preds=preds.apply(lambda x:list(itertools.product(*x))).explode().apply(lambda x: "".join(x)).rename("Prediction")
       
        pepdf.pop("Prediction")
        pepdf=pepdf.merge(preds,left_index=True,right_index=True,how="right").reset_index()

    #select longest tag for alignment
    pepdf["tag"]=pd.Series([max(pep[::2],key=len) for pep in pepdf["Prediction"].str.upper().str.replace("(",")",regex=False).str.split(")",regex=False)],name="longest_tag") 

    return pepdf

@passed_kwargs()
def add_proteins_SMSnet(*,
                    
                    #SMSnet score filters
                    simple_unmask=True,     #attempts to solve low complexity masked SMSnet regions
                    SMSnet_ppm=False,       #max ppm tolerance
                    SMSnet_minscore=False,  #minum mean peptide score
                    
                    **kwargs):
 
    
    v=add_proteins_SMSnet.vars #parse function arguments
    
    #default arguments
    output_folder,tmp_folder=v.output_folder,v.tmp_folder
    simple_unmask=v.simple_unmask
    SMSnet_ppm=v.SMSnet_ppm
    SMSnet_minscore=v.SMSnet_minscore
    
    #required arguments
    input_files=v.input_files 
    Alignment=v.Alignment
    check_required_args(v,["input_files","Alignment"])
    
    al=pd.concat([i for i in Diamond_alignment_Reader(Alignment)])
    al["Proteins"]=al["sseqid"]
    al=al[["tag","Proteins"]].groupby("tag")["Proteins"].apply(lambda x: " ".join(x))

    output_paths=[]
    if type(input_files)!=type(list()):
        input_files=[input_files]
    input_files.sort()
    
    for input_file in input_files:
        
        output_path=str(Path(output_folder,Path(input_file).name.replace("_SMSNET.tsv","_processed_SMSNET.tsv")))
        output_paths.append(output_path)
        pepdf=parse_SMSnet_output(input_file,simple_unmask,SMSnet_ppm,SMSnet_minscore)
        pepdf.merge(al,on="tag",how="left").to_csv(output_path,sep="\t")
        
    return output_paths

@passed_kwargs()
def write_to_Diamond_fasta(*,
                           
                           #MSFragger score filters
                           max_evalue=10,
                           Top_score_fraction=0.9, #in case of multiple top candidates retain the top scoring fraction
                           
                           #SMSnet score filters
                           simple_unmask=True,     #attempts to solve low complexity masked SMSnet regions
                           SMSnet_ppm=False,       #max ppm tolerance
                           SMSnet_minscore=False,  #minum mean peptide score
                           
                           #fasta writing parameters
                           header_info=[],         #information columns that should be retained in the fasta header
                           unique_peptides=True,   #only write unique combinations of header and peptide
                           min_length=4,           #minimum tag length
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
    SMSnet_ppm=v.SMSnet_ppm
    SMSnet_minscore=v.SMSnet_minscore
    header_info=v.header_info
    unique_peptides=v.unique_peptides
    min_length=v.min_length
    output_file=v.output_file
 
 
    #remove output file to prevent overappending to existing file
    out_path=str(Path(output_folder,Path(output_file).stem+".fa"))
    if os.path.exists(out_path): os.remove(out_path)
    
    output_paths=[]
    if type(input_files)!=type(list()):
        input_files=[input_files]
    input_files.sort()

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
            pepdf["tag"]=pepdf["peptide_neighbours"].str.replace("I","L").str.replace("J","L")
        
        
        if input_file.endswith("SMSNET.tsv"): #parse SMSNet outputs
            pepdf=parse_SMSnet_output(input_file=input_file,
                                      SMSnet_ppm=SMSnet_ppm,SMSnet_minscore=SMSnet_minscore,simple_unmask=simple_unmask)

        pepdf=pepdf[pepdf["tag"].fillna("").apply(len)>=min_length]
        
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
                  **kwargs
                  ):

    v=parse_peplist.vars #parse function arguments
    input_file=v.input_file #required argument
    check_required_args(v,["input_file"])
    
    
    #rewrite this to parse multiple files
    
    if  is_fasta(input_file): #do nothing
        return input_file
    
    else: #parse tabular input
        
        output_folder,tmp_folder=v.output_folder,v.tmp_folder
        header_info=v.header_info
        unique_peptides=v.unique_peptides
        min_length=v.min_length
        output_file=v.output_file
        
        if output_file==None: output_file=Path(input_file).stem+".fa"
        out_path=str(Path(output_folder,output_file))
        
        pepdf=read_table(input_file,Keyword="Peptide")
        pepdf=pepdf[pepdf["Peptide"].fillna("").apply(len)>=min_length]
        
        
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
def Diamond_alignment(*,
                      select=" -k25 ", #-top or  -k + integer (see diamond docs)
                      block_size=5,
                      index_chunks=1,
                      minimum_pident=80,
                      minimum_coverage=80,
                      minimum_bitscore=20,
                      other_args=" --algo ctg --dbsize 1 ",
                      diamond_output_columns=diamond_output_columns,
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
            " -f 6  qseqid "+" ".join(diamond_output_columns)+" ", #add qseqid?
            " -t "+'"'+tmp_folder+'"'+other_args])
 

    command+=other_args
        

    
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

    v=filter_Database_taxonomy.vars #parse function arguments
    
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
def load_full_db(Database,**kwargs): #load in memory (works only for small databases)

    v=load_full_db.vars #parse function arguments
    taxdf=v.taxdf

    recs=SeqIO.parse(Database,format="fasta")
    rdf=pd.DataFrame([[str(r.seq),r.description,r.id] for r in recs],columns=["seq","description","id"])
    
    rdf["OX"]=rdf.id.str.split("|").apply(lambda x: x[-1])
    rdf=rdf.merge(taxdf,how="left",left_on="OX",right_index=True).dropna()
    
    return rdf.set_index("id")


@passed_kwargs()
def make_diamond_database(*, #.fa
          
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
        prots[ranks]=taxdf.loc[prots["OX"].tolist(),ranks].values
    

    
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

#%% Helper funs

def is_fasta(input_file):
    fasta=SeqIO.parse(input_file,"fasta")
    return any(fasta)

def reformat_mfg(mgf_file):
    with open(mgf_file,"r+") as f:
        lines=f.readlines()
        lines=[line.replace("Scan: ","scan=")+"SEQ=\n" if line.startswith("TITLE=") else line for line in lines] #add Scan and SEQ=
        f.seek(0)
        f.writelines(lines)

def read_pin(pinfile):
    with open(pinfile, "r") as f:
        lines=f.readlines()    
    header=lines[0].replace("\n","").split("\t")
    return pd.DataFrame([i.replace("\n","").split("\t",len(header)-1) for i in lines[1:]],columns=header)


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
    edf["OX"]=edf["proteins"].str.split("_").apply(lambda x: "_".join(x[-3:]))
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


def chunk_gen(it):
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




