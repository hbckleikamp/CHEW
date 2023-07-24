
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



#%% parameters

svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])


GTDB=True
NCBI=True
rm=False #remove taxonomy downloads after parsing

### update kws from parsed arguments
parser = argparse.ArgumentParser(description="Create taxdf Input arguments",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()

### Define keyword dictionary 
cvars=set([str(k)+":#|%"+str(v)  for k,v in locals().copy().items() if k!="svars"])
kws={i.split(":#|%")[0]:i.split(":#|%")[1]   for i in list(cvars-svars)}
for k in kws.keys(): parser.add_argument("-"+k)
args = {k:v for k,v in vars(parser.parse_args()).items() if v is not None}
kws.update(args)
locals().update(kws)


#%%

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
        

   


def merge(outpath,files):
    with open(outpath,'wb') as wfd:
        for f in files:
            print(f)
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return

def download_extract(urls,path):
    
    download(urls,path)
    extract_subfolders(path)
    



#%%

gtdb_taxdf=pd.DataFrame()
ncbi_taxdf=pd.DataFrame()
paths=[]
ranks=["superkingdom","phylum","class","order","family","genus","species"] #ranks to be included
Dump_keywords=[" bacterium"," archaeon", "unclassified","uncultured","unidentified","environmental samples"] #flagging of dump taxa

#%% Download GTDB taxonomy

if GTDB:

    base_url="https://data.gtdb.ecogenomic.org/releases/latest/"
    r=requests.get(base_url).text
    
    files=[i.split('">')[0] for i in r.split('href="')[1:]]
    taxfiles=[i for i in files if i.endswith("_taxonomy.tsv")]
    
    path=str(Path(basedir,"gtdb_taxonomy"))
    paths.append(path)
    
    output_paths=[]
    for file in taxfiles:
        
        output_path=str(Path(path,file))
        output_paths.append(output_path)
        download_extract(base_url+file,path)
        
    gtdb_taxdf=pd.concat([pd.read_csv(path,sep="\t",header=None) for path in output_paths])
    gtdb_taxdf.columns=["OX","lineage"]
    gtdb_taxdf[ranks]=gtdb_taxdf["lineage"].str.rsplit(";",expand=True)
    gtdb_taxdf=gtdb_taxdf[["OX"]+ranks]
    gtdb_taxdf["Dump_taxid"]=False
    gtdb_taxdf["OS"]=gtdb_taxdf["OX"]

#%%

if NCBI:
        
    url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"    
        
    path=str(Path(basedir,"taxdump"))
    paths.append(path)
    
    download_extract(url,path)
    
    
    nodes=str(Path(path,"nodes.dmp")) #full path location to nodes.dmp
    names=str(Path(path,"names.dmp")) #full path location to names.dmp
    
    with open(names,"r") as f: lines=pd.DataFrame([l.split("\t")[0:7] for l in f.readlines()])
    namesdf=lines.iloc[:,[0,2,4,6]]
    namesdf.columns=["taxid","name","x","type"]
    namesdf=namesdf.loc[namesdf["type"]=="scientific name",["taxid","name"]]
    namesdf=namesdf.set_index("taxid")
    
                                                
    dump_taxids=list(set(sum([ namesdf[namesdf.name.str.contains(i)].index.tolist() for i in Dump_keywords],[]))) #taxids that contain a dump keyword
    
    
    #iteratively construct taxonomies from nodes.dmp
    
    with open(nodes,"r") as f: lines=pd.DataFrame([l.split("\t")[0:5] for l in f.readlines()])
    nodedf=lines.iloc[:,[0,2,4]]
    nodedf.columns=["taxid","parent_taxid","rank"]
    
    
    has_rank=nodedf

    inc=0               #used for renaming columns
    count=len(has_rank) #used for breaking option 1 (no more remaining candidates)
    break_counter=10     #used for breaking option 2 (3 loops same output)
    
    
    parent=nodedf.copy()
    
    completed=[]
    while count!=0:
        
        
    
        #renaming
        parent.columns=[str(inc)+"_taxid",str(inc)+"_parent_taxid",str(inc)+"_rank",]
        
        if inc:
            has_rank=has_rank.merge(parent,left_on=str(inc-1)+"_parent_taxid",
                                    right_on=str(inc)+"_taxid",how="left")
            
        else:
            has_rank=has_rank.merge(parent,left_on="parent_taxid",
                                    right_on=str(inc)+"_taxid",how="left")
            
    
    
    
        d=(has_rank[str(inc)+"_rank"]==ranks[0]).sum() #to check progress, substract those that have the first rank
        
        ix=has_rank[has_rank[str(inc)+"_rank"]=="superkingdom"].index
        
        
        completed.append(has_rank.loc[ix.tolist(),:])
        has_rank=has_rank.drop(ix)
        #this is to break out of the loop if there are no new hits within x iterations
        if d==0:
            break_counter+=1
        else:
            break_counter=0
        if break_counter>=10:
            break
        
        count-=d
        inc+=1
        print(count)
        

    rsp=[]
    for ix,has_rank in enumerate(completed):
        print("parsing batch "+str(ix)+"/"+str(len(completed)-1))
        rs=[]
        for r in ranks:
            print(r)
            v=np.argwhere((has_rank==r).values)
            has_rank.values[v[:,0],v[:,1]]
            rdf=pd.DataFrame([has_rank.values[v[:,0],0], has_rank.values[v[:,0],v[:,1]-2]]).T
            rdf.columns=["idx",r]
            rdf=rdf.set_index("idx")
            rs.append(rdf)
        
        
        p=rs[0]
        for r in rs[1:]:
            p=p.merge(r,how="left",on="idx")
            
        tcols=[i for i in has_rank.columns if "taxid" in i] #flag dump taxa
        
        if len(tcols)>1:
            p["Dump_taxid"]=pd.concat([has_rank[t].isin(dump_taxids) for t in tcols],axis=1).any(axis=1).values
        else:
            p["Dump_taxid"]=has_rank[tcols[0]].isin(dump_taxids).values 
        
        
        
        
        rsp.append(p)
    
    p=pd.concat(rsp)
    
    
    #root is used as placeholder for nans
    ncbi_taxdf=pd.DataFrame(namesdf.loc[p[ranks].fillna("1").values.flatten()].values.reshape(-1,len(ranks)),columns=ranks) 
    ncbi_taxdf[ncbi_taxdf=="root"]=""
    ncbi_taxdf.index=p.index
    ncbi_taxdf["Dump_taxid"]=p["Dump_taxid"]
    
    
    
    ncbi_taxdf.index.name="OX"
    ncbi_taxdf.index=ncbi_taxdf.index.astype(float).astype(int).astype(str)
    
    namesdf.index.name="OX"
    namesdf.columns=["OS"]
    ncbi_taxdf=ncbi_taxdf.merge(namesdf,on="OX").reset_index()

#write outputs
taxdf=pd.concat([gtdb_taxdf,ncbi_taxdf]).set_index("OX")
taxdf.index=taxdf.index.astype(str)
taxdf.loc["",:]=[""]*len(taxdf.columns)
taxdf.to_csv( "parsed_taxonomy.tsv",sep="\t")

#%% Cleanup
if rm:
    for path in paths:
        shutil.rmtree(path)

