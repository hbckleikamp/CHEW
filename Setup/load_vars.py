# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:19:24 2023

@author: ZR48SA
"""


import argparse
from datetime import datetime
import os
from pathlib import Path
import pandas as pd
from collections import Counter


def parse_kws(cvars,svars):
    kws=dict()
    new_vars=list(cvars-svars)
    for i in new_vars:
        s=i.split(":#|%")
        k=s[0]
        v=s[1].strip()
        try:
            v=eval(s[1])
        except:
            pass
        kws.update({k:v})
    
    
    parser = argparse.ArgumentParser(description="input arguments",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser._action_groups.pop()
    for k in kws.keys(): parser.add_argument("-"+k)
    args=vars(parser.parse_args())
    
    for k,v in args.items():
        if v!=None:
            v=v.strip()
            try:
                v=eval(v)     
            except:
                pass
            kws.update({k:v})
    
    return kws



def read_table(tabfile,*,
               Keyword=False, # "Value" or "Peptide" (Keyword is used in dynamic delimiter deetection)
               ):
    
    
    
    try:
        tab=pd.read_excel(tabfile,engine='openpyxl')
        if Keyword in tab.columns:
            return tab
    except:
        pass
        
    try: 
        tab=pd.read_csv(tabfile,sep=",")
        if Keyword in tab.columns:
            return tab
    except:
        pass
        
        
    try: 
        tab=pd.read_csv(tabfile,sep="\t")
        if Keyword in tab.columns:
            return tab
    except:
        pass    


    with open(tabfile,"r") as f: 
        tab=pd.DataFrame(f.readlines())
    
    
   
        
    #dynamic delimiter detection
    #if file delimiter is different, split using different delimiters until the desired column name is found
    if Keyword:
        if Keyword not in tab.columns: 
            delims=[i[0] for i in Counter([i for i in str(tab.iloc[0]) if not i.isalnum()]).most_common()]
            for delim in delims:
                if delim==" ": delim="\s"
                try:
                    tab=pd.read_csv(tabfile,sep=delim)
                    if Keyword in tab.columns:
                        break
                except:
                    pass
                
        
        #2nd try but with f.readlines
        if Keyword not in tab.columns: 
            
            with open(tabfile,"r") as f: 
                tab=pd.DataFrame(f.readlines())
                
                
            for delim in delims:
                if delim==" ": delim="\s"
                try:
                    tab=pd.read_csv(tabfile,sep=delim)
                    if Keyword in tab.columns:
                        break
                except:
                    pass
                  
    return tab
    
def load_variables(file):
    
    tab=read_table(file,Keyword="Value").fillna("")
    
    vs=[]
    
    for i in tab.Value:
        if i=="None":
            vs.append(None)
        else:
            try: 
                vs.append(eval(i))
            except:
                vs.append(i)
            
    return dict(zip(tab.Key, vs))
    
