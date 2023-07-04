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
input_file=""    #fasta file
unclustered_database_dmnd=""  #diamond database

#Example syntax:
#peplist2Initial_DB.py -input_file "  your input file  " -unclustered_database_dmnd " your database "

#Optional arguments
Temporary_directory=""
Output_directory=""
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value

#Default arguments

#Diamond specific arguments
select=" -k25 " #-top or  -k + integer (see diamond docs)
block_size=5
index_chunks=1
minimum_pident=80
minimum_coverage=80
minimum_bitscore=20
other_args=" --algo ctg --dbsize 1 "

initial_minimum_taxid_frequency=20 #minimum taxid frequency that Diamond database should have




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
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_peplist2Initial_DB.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)

#load taxdf
taxdf=read_table(taxdf_path,Keyword="OX").set_index("OX")
taxdf.index=taxdf.index.astype(str)
kws.update({"taxdf":taxdf}) 


#%%  Peplist2Initial_DB

#align peptide list with DIAMOND
Alignment=Diamond_alignment(input_file=input_file,
                            database_path=unclustered_database_dmnd)


#write matched sequences to database
initial_target=Write_alignment_to_database(input_file=Alignment)

#filter database based on frequency
if initial_minimum_taxid_frequency:
    composition,richness,entries=write_database_composition(input_file=initial_target)
    taxids=composition[composition["Count"]>=initial_minimum_taxid_frequency].index.tolist() 
    initial_target=filter_Database_taxonomy(input_file=initial_target,taxids=taxids)

