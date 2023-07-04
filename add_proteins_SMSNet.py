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
input_files=[]     # list or space delimited string of filepaths SMSNet files
database_path=""   #fasta file: database used to cunstrcut diamond database and for exact string searching
peplist=""         #fasta file: aligned against diamond database

#output folders
Temporary_directory=""
Output_directory=""
output_file="final"
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value


#SMSNet score filters
simple_unmask=True      #attempts to solve low complexity masked SMSNet regions
SMSNet_ppm=False        #max ppm tolerance
SMSNet_minscore=False   #minimum mean peptide score
Exact_tag=database_path #exact srting searching of SMSnet against a database

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
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_add_proteins_SMSNet.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)



#%% Add proteins SMSnet

if type(input_files)==str:
    if os.path.isdir(input_files):
        x=[str(Path(input_files,i)) for i in os.listdir(input_files)]
        if len(x):
            input_files=x
if type(input_files)==str:
    input_files=input_files.split()
if type(input_files)==str:
    input_files=[input_files]

SMSNet_files=[i for i in input_files if i.endswith("SMSNET.tsv")]


#Aligning of SMSNet tags against diamond database
final_database_dmnd=make_diamond_database(input_file=database_path) 
if not len(peplist): peplist=write_to_Diamond_fasta(SMSNet_files)
Alignment=Diamond_alignment(input_file=peplist,database_path=final_database_dmnd)  

#adding proteins including Exact tag matches
add_proteins_SMSNet(input_files=SMSNet_files,Alignment=Alignment)    

