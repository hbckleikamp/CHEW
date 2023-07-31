# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 13:28:44 2023

@author: ZR48SA
"""

#%% clear variables and console, stor current variables

try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

from load_vars import *
from config import *

import warnings
warnings.filterwarnings("ignore") #remove when debugging!


from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0])) # set base path
basedir=os.getcwd()
print(basedir)
svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])



#%% Define parameter dict (kws)

#required arguments
input_files=""     # list or space delimited string of filepaths, or folder
mgf_files=""
SMSNet_files=""    #
database_path=""   #fasta file 

#Optional arguments
Temporary_directory=""
Output_directory=""
output_folder="final"  
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value

#default arguments
params_path=params_mid                                 # path to params file with detailed MSFragger parameters
max_no_hits=5                                               # max number of hits retained from each database split
no_splits=None                                              # number of database splits, determines performance and temporary index size
no_batches=None                                             # number of file splits,     determines performance and temporary index size

#Example syntax
#Annotate_final.py -SMSNet_fies " your smsnet files " -input_files "  your input folder  " -database_path " your database " -Output_directory "  your output folder   "



#%% Update parameter dict (kws)
           
### Define keyword dictionary 
cvars=set([str(k)+":#|%"+str(v)  for k,v in locals().copy().items() if k!="svars"])
kws={i.split(":#|%")[0]:i.split(":#|%")[1]   for i in list(cvars-svars)}
kws=parse_kws(cvars,svars) #use this everywhere

#log keyword dictionary
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_Annotate_final.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)

taxdf=read_table(taxdf_path,Keyword="OX").set_index("OX")
taxdf.index=taxdf.index.astype(str)
kws.update({"taxdf":taxdf}) 

#%% Annotate MSFragger

final_annotations=MSFragger_annotation(input_files=input_files,database_path=database_path)
    
if len(SMSNet_files): #identify de novo tags

    final_annotations+=add_proteins_SMSNet(input_files=SMSNet_files,database_path=database_path)    
    


Post_processing(input_files=final_annotations,database=database_path,no_splits=4,mgf_files=mgf_files)
