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

import time
start_time=time.time()

import warnings
warnings.filterwarnings("ignore") #remove when debugging!


from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0])) # set base path
basedir=os.getcwd()
print(basedir)
svars=set([str(k)+":#|%"+str(v) for k,v in locals().copy().items()])







#%% Define parameter dict (kws)

#required arguments
#This subroutine has the most rquired arguments out of all, therefore it is recommended to work with paramter files

#should either supply a list of mzML (and or) SMSNet file
mzML_files=""                 # list or space delimited string or folder of filepaths mzML files, used for MSFragger
peplist=""                    # fasta file: aligned against diamond database
clustered_database_fasta=""   # fasta file 
unclustered_database_fasta="" # fasta file (or folder with fasta files in case of GTDB)
unclustered_database_dmnd=""  # diamond database


#Example syntax:
#peplist2Annotation.py -mzML_files "path" -peplist "path" -clustered_database_fasta "path" -unclustered_database_fasta "path" -unclustered_database_dmnd "path"

SMSNet_files=""    #optional: list or space delimited string or folder of filepaths SMSNet files, used for SMSNet

#output folders
Temporary_directory=""
Output_directory=""
output_folder="peplist"
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value, recommended for large scripts

#%% Default arguments

#### Section 1: raw2peplist

#Which annotation should be used? (one or more should be True)
MSFragger=False
SMSNet=False

### annotate_MSFragger


#MSFragger score filters
max_evalue=10
Top_score_fraction=0.9 #in case of multiple top candidates retain the top scoring fraction

#SMSNet score filters
simple_unmask=True     #attempts to solve low complexity masked SMSNet regions
SMSNet_ppm=20          #max ppm tolerance
SMSNet_minscore=False  #minum mean peptide score

#Database construction parameters
header_info=[]          #information columns that should be retained in the fasta header
unique_peptides=True    #only write unique combinations of header and peptide
min_length=4            #minimum tag length

#### Section 2: peplist2db

#Diamond specific arguments
select=" -k25 " #-top or  -k + integer (see diamond docs)
block_size=5
index_chunks=1
minimum_pident=80
minimum_coverage=80
minimum_bitscore=20
other_args=" --algo ctg --dbsize 1 "

initial_minimum_taxid_frequency=20 #minimum taxid frequency that Diamond database should have

#### Section 3: refine_db
#Refine MSFragger params
refine_params=str(Path(basedir,"closed_fragger_mid.params")) 
refine_no_splits=2            # number of database splits, determines performance and temporary index size
refine_no_batches=1           # number of file splits,     determines performance and temporary index size

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
final_params=str(Path(basedir,"closed_fragger_final.params")) #detailed search for smaller db
final_min_ratio=0.99
final_minimum_taxid_frequency=5
final_denoise_remove=True

#### Section 4: final annotation
final_MSFragger=True
mgf_files=""
database=""
FDR=0.05            #false discovery rate
min_peptide_count=1 #minimum occurrence for each unique peptide
remove_unannotated=False 

#%%



#%% Update parameter dict (kws)
           
### Define keyword dictionary 
cvars=set([str(k)+":#|%"+str(v)  for k,v in locals().copy().items() if k!="svars"])
kws={i.split(":#|%")[0]:i.split(":#|%")[1]   for i in list(cvars-svars)}
kws=parse_kws(cvars,svars) #use this everywhere

print(kws)
print()
#log keyword dictionary
kws_filename=datetime.now().strftime("%y_%m_%d_%H_%M_%S")+"_peplist2Annotation.CHEW_params" #the basename needs to be changed manually per script, since inspect only works form CLI
kws_df=pd.DataFrame.from_dict(kws,orient="index").fillna("").reset_index()
kws_df.columns=["Key","Value"]
kws_df.set_index("Key").to_csv(kws_filename,sep="\t")


from CHEW_funs import *

locals().update(kws)

taxdf=read_table(taxdf_path,Keyword="OX").set_index("OX")
taxdf.index=taxdf.index.astype(str)
kws.update({"taxdf":taxdf}) 


#%% parse input files

all_inp=[]


for i in [mzML_files,SMSNet_files]:

    if type(i)==str:
        if os.path.isdir(i):
            x=[str(Path(i,f)) for f in os.listdir(i)]
            if len(x):
                i=x
    if type(i)==str:
        i=i.split()
    if type(i)==str:
        i=[i]

    all_inp+=i

SMSNet_files=[i for i in all_inp if i.endswith("SMSNET.tsv")]
mzML_files=[i for i in all_inp if i.endswith(".mzML" or ".raw")]

SMSNet=len(SMSNet_files)
  


#%% Initial_db

#parse peplist (if file is not fasta)
Diamond_fasta=parse_peplist(input_file=peplist) 

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
                                              params_path=refine_params,
                                              no_splits=refine_no_splits,
                                              no_batches=refine_no_batches
                                              
                                              )

kws.update({"output_folder":"refine_final"})

proteins,taxids=denoise_nodes(tlca,min_ratio=final_min_ratio,denoise_ranks=final_denoise_ranks,remove=True)

if final_minimum_taxid_frequency:
    DB_in_mem,database_path,richness,entries=filter_Database_proteins_in_mem(input_file=DB_in_mem,proteins=proteins)
    composition,richness,entries=write_database_composition(input_file=database_path)
    taxids=composition[composition["Count"]>=final_minimum_taxid_frequency].index.tolist() 

print("final number of taxa in db: "+str(len(taxids)))

#%% Final_db

kws.update({"output_folder":"final"})
print("retrieving proteomes")
final_target=filter_Database_taxonomy(input_file=unclustered_database_fasta,taxids=taxids)
final_decoy=write_decoy(input_file=final_target,method="reverse")
final_database=merge_files([final_target,final_decoy])




#%% Final_annotation

MSFragger=final_MSFragger #Always do final annotation hybrid!

final_annotations=[]

if SMSNet: #identify de novo tags

    final_annotations+=add_proteins_SMSNet(input_files=SMSNet_files,database_path=final_database)    

if MSFragger: 
    final_annotations+=MSFragger_annotation(input_files=mzML_files,database_path=final_database,params_path=final_params)




#%% Post processing

kws.update({"output_folder":"output"})

Post_processing(input_files=final_annotations,database=final_database,no_splits=4)
    
#%%
import time
print("total elapsed time: "+str(time.time()-start_time))


