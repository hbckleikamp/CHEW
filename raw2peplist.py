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
input_files=[]    # list or space delimited string of filepaths (.mgf for SMSNet .mzML for MSFragger, .raw for both)
clustered_database_fasta=""   #fasta file 

#Which annotation should be used? (one or more should be True)
MSFragger=False
SMSNet=False

#Example syntax:
#raw2peplist.py -input_files "  your input folder  "  -database_path " your database " -MSFragger 1 -SMSNet 1




#output folders
Temporary_directory=""
Output_directory=""
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value


### annotate_MSFragger
params_path=str(Path(basedir,"closed_fragger_fast.params")) # path to params file with detailed MSFragger parameters
max_no_hits=5                                       # max number of hits retained from each database split
no_splits=None                                      # number of database splits, determines performance and temporary index size
no_batches=None                                     # number of file splits,     determines performance and temporary index size

### write_to_Diamond_fasta

#MSFragger score filters
max_evalue=10                                       # maximum allowed evalue score of peptides 
Top_score_fraction=0.9                              # in case of multiple top candidates retain the top scoring fraction 

#SMSNet score filters
simple_unmask=True     #attempts to solve low complexity masked SMSNet regions
SMSNet_ppm=False       #max ppm tolerance
SMSNet_minscore=False  #minimum mean peptide score

#Database construction parameters
header_info=[]          #information columns that should be retained in the fasta header
unique_peptides=True    #only write unique combinations of header and peptide
min_length=4            #minimum tag length



#%% Update parameter dict (kws)
           
### Define keyword dictionary 
cvars=set([str(k)+":#|%"+str(v)  for k,v in locals().copy().items() if k!="svars"])
kws={i.split(":#|%")[0]:i.split(":#|%")[1]   for i in list(cvars-svars)}


### update kws from parsed arguments
parser = argparse.ArgumentParser(description="CHEW Input arguments",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
for k in kws.keys(): parser.add_argument("-"+k)
args = {k:v for k,v in vars(parser.parse_args()).items() if v is not None}
kws.update(args)

### update kws from variable_tab 
if variable_tab: kws.update(load_variables(variable_tab)) #uses columns: Key, Value
if "variable_tab" in kws.keys():
    if kws.get("variable_tab"): 
        kws.update(load_variables(kws.get("variable_tab")))

#log keyword dictionary
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_raw2peplist.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)



#%% Raw2peplist

kws.update({"output_folder":"annotations"}) #update output_subfolder

annotations=[]

mzML_files=raw2mzML() 


#annotate de novo
if SMSNet: 
    mgf_files=raw2mgf()
    SMSNet_files=SMSNet_annotation(input_files=mgf_files)
    annotations+=SMSNet_files
    
#annotate with clustered database
if MSFragger: 
    MSFragger_files=MSFragger_annotation(input_files=mzML_files,database_path=clustered_database_fasta)
    annotations+=MSFragger_files

kws.update({"output_folder":"peplist"}) #update output_subfolder

Diamond_fasta=write_to_Diamond_fasta(input_files=annotations)


    

