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

from CHEW_funs import *
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0])) # set base path
basedir=os.getcwd()
print(basedir)
svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])


#%% ### Description

#This script serves as the 1st part of the CHEW workflow.
#By annotation with de novo sequencing (SMSNet) or a clusterd database (MSFragger),
#a peptide list is generated for alignment with Diamond


#construction of parameter dict (kws) is as follows:
    #lowest priority: params defined inside script
    #next: params from argparse (only during command line usage)
    #highest priority: params passed from file (variable tab)

#Agrument handling is as follows:
    #lowest priority: default function arguments
    #next: arguments inside parameter dict (kws)
    #highest priority: explicitly passed arguments in function call (will orverwrite arguments in kws)


#%% Define parameter dict (kws)

#required arguments
input_files=[]    # list or space delimited string of filepaths (.mgf for SMSnet .mzML for MSFragger, .raw for both)
clustered_database_fasta=""   #fasta file 
unclustered_database_fasta="" #fasta file (or folder with fasta files in case of GTDB)
unclustered_database_dmnd=""  #diamond database

#variable handling and logging
variable_tab=""   # Optional: supply parameters from a file
write_vars=False  # Optional: log parameter dict (kws) used within script to file

#output folders
Temporary_directory=""
Output_directory=""
output_folder=""
output_file=""

#### Section 1: raw2peplist

#Which annotation should be used? (one or more should be True)
MSFragger=False
SMSNet=False

### annotate_MSFragger
params_path=str(Path(basedir,"params_fast.params")) # path to params file with detailed MSFragger parameters
max_no_hits=5                                       # max number of hits retained from each database split
no_splits=None                                      # number of database splits, determines performance and temporary index size
no_batches=None                                     # number of file splits,     determines performance and temporary index size

### write_to_Diamond_fasta

#MSFragger score filters
max_evalue=10                                       # maximum allowed evalue score of peptides 
Top_score_fraction=0.9                              # in case of multiple top candidates retain the top scoring fraction 

#SMSnet score filters
simple_unmask=True     #attempts to solve low complexity masked SMSnet regions
SMSnet_ppm=False       #max ppm tolerance
SMSnet_minscore=False  #minimum mean peptide score

#Database construction parameters
header_info=[]          #information columns that should be retained in the fasta header
unique_peptides=True    #only write unique combinations of header and peptide
min_length=4            #minimum tag length

#### Section 2: peplist2db

### Diamond_alignment

minimum_taxid_frequency=20

### Write_alignment_to_database



#### Section 3: refine_db


#### Section 4: construct_db


#### Section 5: final annotation


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
if variable_tab: kws.update(load_variables(variable_tab))
if "variable_tab" in kws.keys():
    if kws.get("variable_tab"): 
        kws.update(load_variables(kws.get("variable_tab")))

#log keyword dictionary
if write_vars:
    filename=datetime.now.strftime("%y_%m_%d_%H_%M_%S")+"raw2fasta.txt" #the basename needs to be changed manually per script
    pd.DataFrame.from_dict(kws).to_csv(filename,sep="\t")


#%% Raw2peplist

annotations=[]

#annotate de novo
if SMSNet: annotations+=annotate_SMSnet()
    
#annotate with clustered database
if MSFragger: annotations+=annotate_MSFragger(database_path=clustered_dabase_fasta)

Diamond_fasta=write_to_Diamond_fasta(input_files=annotations)

#%% peplist2db

#parse peplist (if file is not fasta)
Diamond_fasta=parse_peplist(Diamond_fasta) 

#align peptide list with DIAMOND
Alignment=Diamond_alignment(input_file=Diamond_fasta)

#write matched sequences to database
target=Write_alignment_to_database(input_file=Alignment)

#filter database based on frequency
if minimum_taxid_frequency:
    composition,richness,entries=write_database_composition(input_file=target)
    taxids=composition[composition["Count"]>=minimum_taxid_frequency].index.tolist() 
    target=filter_Database_taxonomy(input_file=target,taxids=taxids)


#%% Refine_db

#load in mem

#

#output_folders=

#%% Construct_db


#%%

#annotate_MSFragger(database_path=clustered_dabase_fasta)
#post processing (FDR, lca)



#%% Cleanup



