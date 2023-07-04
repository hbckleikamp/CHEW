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
input_files=""    # list or space delimited string of filepaths, or folder (mzML files)
database_path=""   #fasta file 

#Example syntax:
#Refine_DB.py -input_files "  your input files  " -database_path " your database "

#Optional arguments
Temporary_directory=""
Output_directory=""
output_folder="refine"
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value

#Default arguments

#MSFragger params
params_path=params_mid # path to params file with detailed MSFragger parameters
no_splits=2            # number of database splits, determines performance and temporary index size
no_batches=1           # number of file splits,     determines performance and temporary index size

#MSFragger score filters
max_evalue=10
Top_score_fraction=0.9 #in case of multiple top candidates retain the top scoring fraction

#Pre LCA filter
Frequency_prefilter=2,  # 2   Static prefiler cutoff, taxa should have more PSMs 
Precision_prefilter=0.7 # 0.7 Target Decoy precision based denoising pre LCA filter
prefilter_remove=False  # if prefilter completly removes a scan, retain taxa?

#focusing LCA parameters
weight_rank="species"  # rank used for calculating weights during focusing lca
weight_cutoff=0.6      # weight cutoff for focusing lca 

#Post LCA representative picking
min_count=2
denoise_ranks=["genus","species"]
min_ratio=[0.99,0.95]
denoise_remove=True

min_rate=0.95 #break out of refinement when database richness decreases less than min_rate

### post refinement filtering ###
final_denoise_ranks=["class","order","family","genus","species"]
final_min_ratio=0.99
final_minimum_taxid_frequency=5




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
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_RefineDB.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)

taxdf=read_table(taxdf_path,Keyword="OX").set_index("OX")
taxdf.index=taxdf.index.astype(str)
kws.update({"taxdf":taxdf}) 




#%% Refine_db

tlca,database_path,DB_in_mem=refine_database(input_files=input_files,
                                              database_path=database_path)

kws.update({"output_folder":"refine_final"})

proteins,taxids=denoise_nodes(tlca,min_ratio=final_min_ratio,denoise_ranks=final_denoise_ranks,remove=True)

if final_minimum_taxid_frequency:
    DB_in_mem,database_path,richness,entries=filter_Database_proteins_in_mem(input_file=DB_in_mem,proteins=proteins)
    composition,richness,entries=write_database_composition(input_file=database_path)
    taxids=composition[composition["Count"]>=final_minimum_taxid_frequency].index.tolist() 

#write results
pd.DatFrame(taxids,columns=["Taxids"]).to_csv(database_path.replace(".fa","_taxids.tsv"),sep="\t")
pd.DatFrame(proteins,columns=["Proteins"]).to_csv(database_path.replace(".fa","_proteins.tsv"),sep="\t")

print("final number of taxa in db: "+str(len(taxids)))


