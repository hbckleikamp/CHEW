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

def read_table(tabfile,*,
               Keyword=False, # "Value" or "Peptide" (Keyword is used in dynamic delimiter deetection)
               ):
    
    
    if tabfile.endswith('.xlsx') or tabfile.endswith('.xls'): tab=pd.read_excel(tabfile,engine='openpyxl')   
    elif tabfile.endswith('.csv'): tab=pd.read_csv(tabfile,sep=",")
    elif tabfile.endswith('.tsv'):                               tab=pd.read_csv(tabfile,sep="\t")
    else:

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
    
