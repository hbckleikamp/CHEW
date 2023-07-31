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


import warnings
warnings.filterwarnings("ignore") #remove when debugging!


from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0])) # set base path
basedir=os.getcwd()
print(basedir)
svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])

#%% Define parameter dict (kws)

#required arguments


databases="" #required

peplist=""        #Optional: tabular or fasta format file that will be aligned
input_files=""    #Optional: instead of supplying a peplist, a peplist is constructed from CHEW PSMs or SMSNet files 



#Optional arguments
output_folder=""
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value



#Example syntax
#raw2mgf.py -input_files "  your input folder  " -Output_directory "  your output folder   "


#%% Update parameter dict (kws)
           
### Define keyword dictionary 
cvars=set([str(k)+":#|%"+str(v)  for k,v in locals().copy().items() if k!="svars"])
kws={i.split(":#|%")[0]:i.split(":#|%")[1]   for i in list(cvars-svars)}
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()

#log keyword dictionary
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_raw2mgf.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)



#%%


# dynamic reading of databases and input files
#add option for mutliple database comparison
#add function: database qc from alignment
#input would either be databases or alignments




if input_files="":
    
    
@passed_kwargs()
def database_QC_from_peplist(*,
                             placeholder="",
                             **kwargs):

    check_required_args(v,["input_file","database_1","database_2"])


def database_QC_from_CHEW_PSMs(
                                *,
                               scoring_metric='mass_corr_hyperscore',
                               peplist="",
                               **kwargs):
    

    check_required_args(v,["input_files","database_1","database_2"])

    
    
    
@passed_kwargs()
def database_QC_from_SMSNet(
                              

#make 

mzML_files=raw2mgf() 






