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
input_files=[]    # list or space delimited string of filepaths (.mgf for SMSNet .mzML for MSFragger, .raw for both)
clustered_database_fasta=""   #fasta file 
unclustered_database_fasta="" #fasta file (or folder with fasta files in case of GTDB)
unclustered_database_dmnd=""  #diamond database

#variable handling and logging
variable_tab="C:/MP-CHEW/CHEW/test_params.xlsx"   # Optional: supply parameters from a file, uses columns: Key, Value
#write_vars=True   # Optional: log parameter dict (kws) used within script to file

#output folders
Temporary_directory=""
Output_directory=""
#output_folder=""
#output_file=""

#### Section 1: raw2peplist

#Which annotation should be used? (one or more should be True)
MSFragger=False
SMSNet=False

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

#### Section 2: peplist2db

### Diamond_alignment

initial_minimum_taxid_frequency=20

### Write_alignment_to_database



#### Section 3: refine_db
final_minimum_taxid_frequency=5

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
if variable_tab: kws.update(load_variables(variable_tab)) #uses columns: Key, Value
if "variable_tab" in kws.keys():
    if kws.get("variable_tab"): 
        kws.update(load_variables(kws.get("variable_tab")))

#log keyword dictionary
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_raw2fasta.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)

taxdf=read_table(taxdf_path,Keyword="OX").set_index("OX")
taxdf.index=taxdf.index.astype(str)

#important for taxdf, add this to setup! 
taxdf.loc["",:]=""

kws.update({"taxdf":taxdf}) 

#%%



#%% Raw2peplist

kws.update({"output_folder":"peplist"}) #update output_subfolder

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

Diamond_fasta=write_to_Diamond_fasta(input_files=annotations)



#%% Initial_db

#parse peplist (if file is not fasta)
Diamond_fasta=parse_peplist(input_file=Diamond_fasta) 

#align peptide list with DIAMOND
Alignment=Diamond_alignment(input_file=Diamond_fasta,
                            database_path=unclustered_database_dmnd)


#write matched sequences to database
initial_target=Write_alignment_to_database(input_file=Alignment)

#filter database based on frequency
if initial_minimum_taxid_frequency:
    composition,richness,entries=write_database_composition(input_file=initial_target)
    taxids=composition[composition["Count"]>=initial_minimum_taxid_frequency].index.tolist() 
    initial_target=filter_Database_taxonomy(input_file=initial_target,taxids=taxids)


#%% Refine_db

kws.update({"output_folder":"refine"})



tlca,database_path,DB_in_mem=refine_database(input_files=mzML_files,
                                              database_path=initial_target,
                                              denoise_ranks=["genus","species"],
                                              min_ratio=[0.99,0.95],
                                              Frequency_prefilter=2,
                                              Precision_prefilter=0.7,
                                              prefilter_remove=False,
                                              denoise_remove=True,
                                              )

kws.update({"output_folder":"refine_final"})

proteins,taxids=denoise_nodes(tlca,min_ratio=0.99,denoise_ranks=["class","order","family","genus","species"],remove=True)

if final_minimum_taxid_frequency:
    DB_in_mem,database_path,richness,entries=filter_Database_proteins_in_mem(input_file=DB_in_mem,proteins=proteins)
    composition,richness,entries=write_database_composition(input_file=database_path)
    taxids=composition[composition["Count"]>=final_minimum_taxid_frequency].index.tolist() 

print("final number of taxa in db: "+str(len(taxids)))

#%% Final_db

kws.update({"output_folder":"final"})

final_target=filter_Database_taxonomy(input_file=unclustered_database_fasta,taxids=taxids)
final_decoy=write_decoy(input_file=final_target,method="reverse")
final_database=merge_files([final_target,final_decoy])

if SMSNet:  #construct diamond DB for tag alignment
    final_database_dmnd=make_diamond_database(input_file=final_database) #this is a diamond databse with decoy proteins
    



#%% Final_annotation



final_annotations=[]

if SMSNet: #align de novo tags
    Alignment=Diamond_alignment(input_file=Diamond_fasta,database_path=final_database_dmnd)  #if diamond_fasta exists..
    final_annotations+=add_proteins_SMSNet(input_files=SMSNet_files,Alignment=Alignment,
                                           Exact_tag=final_target)    

if MSFragger: 
    final_annotations+=MSFragger_annotation(input_files=mzML_files,database_path=final_database,params_path=params_final)




#%% Post processing

kws.update({"output_folder":"output"})

Post_processing(input_files=final_annotations)
    

