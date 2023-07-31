# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:36:32 2023

@author: ZR48SA
"""

import pandas as pd
import numpy as np
import itertools

#%% set base path
from pathlib import Path
import os
from inspect import getsourcefile
# change directory to script directory (should work on windows and mac)
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())
basedir=os.getcwd()

## mass calculation      
         
std_aa_mass = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841,
               'T': 101.04768,'C': 103.00919,'L': 113.08406,'I': 113.08406,'J': 113.08406,
               'N': 114.04293,'D': 115.02694,'Q': 128.05858,'K': 128.09496,'E': 129.04259,
               'M': 131.04049,'H': 137.05891,'F': 147.06841,'U': 150.95364,'R': 156.10111,
               'Y': 163.06333,'W': 186.07931,'O': 237.14773}

#add fixed mods
fixed_mods={'C' : 57.021464,'U' : 57.021464} #Cysteine is alkylated
[std_aa_mass.update({k:(std_aa_mass.get(k)+v)})     for k,v in fixed_mods.items()]

std_aa_df=pd.DataFrame.from_dict(std_aa_mass,orient="index")
std_aa_df.index=std_aa_df.index.str.upper() #add modifications or not?

std_aa_df.loc["",:]=np.nan



#add variable mods (core, mass_shift, new name)
#here you can do some complicated stuff and add acetylation
#%%
#Standard SMSnet modifications
variable_mods=pd.DataFrame([["M",15.994915,"m"],
                            ["Y",79.96633,"y"], #phosphorylation
                            ["S",79.96633,"s"], #phosphorylation
                            ["T",79.96633,"t"]],#phosphorylation 
                           

                           columns=["Target","Mass","name"])
#Add acetylations?

[std_aa_mass.update({i["name"]:(std_aa_mass.get(i["Target"])+i["Mass"])}) for n,i in variable_mods.iterrows()]




if not os.path.exists(str(Path(basedir,"simple_unmasks.tsv"))):
# #precompute all combinations of up to x amino acids (for unmasking SMSnet)

    p_aas=np.array([i for i in std_aa_mass.keys() if i not in ["I","J"]])#,"O"]]) #no I J O
    p_mass=np.array([std_aa_mass.get(i) for i in p_aas])
    
    permutated_aa=[]
    for i in range(8): #what range is reasonable to unmask? increasing leads to slower unmasking
        
        m=np.array(list(itertools.combinations_with_replacement(p_mass, i+1))).sum(axis=1)
        aa=list(itertools.combinations_with_replacement(p_aas, i+1))
        df=pd.DataFrame(m,columns=["mass"])
        df["composition"]=aa
        permutated_aa.append(df)
    permutated_aa=pd.concat(permutated_aa)
    permutated_aa["composition"]=permutated_aa["composition"].apply(lambda x: "".join(list(x)))
    
    permutated_aa["mass"]=permutated_aa["mass"].round(4)
    permutated_aa=permutated_aa.drop_duplicates()
    permutated_aa=permutated_aa.set_index("mass")
    
    #remove unlikely stuff:
        #more than one of O,U,y,s,t
    permutated_aa=permutated_aa[~pd.concat([permutated_aa["composition"].str.contains(i) for i in ["OO","UU","yy","ss","tt"]],axis=1).any(axis=1)]
    permutated_aa=permutated_aa[pd.concat([permutated_aa["composition"].str.contains(i) for i in ["O","U","y","s","t"]],axis=1).sum(axis=1)<2]
        
    permutated_aa=permutated_aa.groupby("mass")["composition"].apply(list)
    permutated_aa.to_csv(str(Path(basedir,"permutated_aa.tsv")),sep="\t")

permutated_aa=pd.read_csv(str(Path(basedir,"permutated_aa.tsv")),sep="\t")

# #find permutated_aas with only one or two possibilities
if not os.path.exists(str(Path(basedir,"simple_unmasks.tsv"))):

    spermuta=permutated_aa.explode("composition").sort_index()
    
    
    spermuta["combinations"]=spermuta["composition"].apply(lambda x:1+ (len(list(set(x)))-1)*(len(x)-1)) #fast approximation of the amount of permutations 
    g=spermuta.groupby(spermuta.index)["combinations"].sum()<3
    xspermuta=spermuta[spermuta.index.isin(g[g].index.tolist())]
    
    ix_counter=0
    ppm50_merged=[]
    ls=spermuta.index*(1-50/1000000)
    hs=spermuta.index*(1+50/1000000)
    
    lowest=0
    for ix,i in enumerate(xspermuta.index):
        df=spermuta[(spermuta.index>=ls[ix]) & (spermuta.index<=hs[ix])]
    
        df["combinations"]=df["composition"].apply(lambda x:1+ (len(list(set(x)))-1)*(len(x)-1)) #fast approximation of the amount of permutations for 2 amino acid combinations     
        s=df["combinations"].sum()
        if s<3:
            ppm50_merged.append([ls[ix],hs[ix],", ".join(df["composition"]),s])
    
    ppm50_merged=pd.DataFrame(ppm50_merged,columns=["lower","upper","compositions","combinations"]).drop_duplicates()
    ppm50_merged["compositions"]=ppm50_merged["compositions"].str.split(", ")
    ppm50_merged=ppm50_merged.explode("compositions")
    ppm50_merged["unique_aas"]=ppm50_merged["compositions"].apply(lambda x: len(list(set(x))))
    ppm50_merged["length"]=ppm50_merged["compositions"].apply(len)
    ppm50_merged=ppm50_merged[ppm50_merged["unique_aas"]<3]
    ppm50_merged=ppm50_merged[(ppm50_merged["unique_aas"]==1) | (ppm50_merged["length"]<3)] #either unique_aas=1 or len >3
    ppm50_merged.to_csv(str(Path(basedir,"simple_unmasks.tsv")),sep="\t")

simple_masks=pd.read_csv(str(Path(basedir,"simple_unmasks.tsv")),sep="\t",usecols=["lower","upper","compositions"])