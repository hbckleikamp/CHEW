# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:42:24 2023

@author: ZR48SA
"""
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0])) # set base path
basedir=os.getcwd()





#%% Filepaths

msconvert_filepath="msconvert"  #msconvert is assumed to be added to environment variables as path
blast_folderpath=""             #the folder of local blast+ installation  is assumed to be added to environment variables as path
diamond_filepath=str(Path(basedir,"diamond"))

#MSFragger filepaths
MSFragger_jar_path=str(Path(basedir,"MSFragger-3.5.jar"))  
pep_split_path=str(Path(basedir,"msfragger_pep_split_HK.py"))
params_fast =str(Path(basedir,"closed_fragger_fast.params"))  #faster initial search
params_mid  =str(Path(basedir,"closed_fragger_mid.params"))   #medium search during refinement
params_final=str(Path(basedir,"closed_fragger_final.params")) #detailed search for smaller 


#%% prep paths
msconvert_filepath=' "'+str(Path(Path(msconvert_filepath).parents[0],Path(msconvert_filepath).name))+'" '
diamond_filepath=  ' "'+str(Path(Path(  diamond_filepath).parents[0],Path(  diamond_filepath).name))+'" '
#%%

