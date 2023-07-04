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

database_path=""   #fasta file: database used to cunstruct diamond database and for exact string searching

mzML_files=""      #optional: list or space delimited string or folder of filepaths mzML files, used for MSFragger
SMSNet_files=""    #optional: list or space delimited string or folder of filepaths SMSNet files, used for SMSNet
#should either supply a list of mzML (and or) SMSNet files
peplist=""         #ofasta file: aligned against diamond database


#Optional
Temporary_directory=""
Output_directory=""
output_file="final"
variable_tab=""   # Optional: supply parameters from a file, uses columns: Key, Value



#default arguments

#MSFragger score cutoffs
max_evalue=10
Top_score_fraction=0.9  
   
#SMSNet score cutoffs
simple_unmask=True      #attempts to solve low complexity masked SMSNet regions
SMSNet_ppm=20          #max ppm tolerance
SMSNet_minscore=False  #minum mean peptide score
Exact_tag=database_path #exact srting searching of SMSnet against a database

#post processing
FDR=0.05            #false discovery rate
min_peptide_count=1 #minimum occurrence for each unique peptide
remove_unannotated=False 




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

taxdf=read_table(taxdf_path,Keyword="OX").set_index("OX")
taxdf.index=taxdf.index.astype(str)
kws.update({"taxdf":taxdf}) 

#%% Annotate final

all_inp=[]


for i in [mzML_files,SMSNet_files]:

    if type(i)==str:
        if os.path.isdir(i):
            x=[str(Path(i,i)) for i in os.listdir(i)]
            if len(x):
                i=x
    if type(i)==str:
        i=i.split()
    if type(i)==str:
        i=[i]

    all_inp+=i

SMSNet_files=[i for i in all_inp if i.endswith("SMSNET.tsv")]
mzML_files=[i for i in all_inp if i.endswith(".mzML" or ".raw")]

annotations=[]

if len(SMSNet_files):
    
    #Aligning of SMSNet tags against diamond database
    final_database_dmnd=make_diamond_database(input_file=database_path) 
    if not len(peplist): peplist=write_to_Diamond_fasta(SMSNet_files)
    Alignment=Diamond_alignment(input_file=peplist,database_path=final_database_dmnd)  
    
    #adding proteins including Exact tag matches
    annotations+=add_proteins_SMSNet(input_files=SMSNet_files,Alignment=Alignment)    


if len(mzML_files): 
    annotations+=MSFragger_annotation(input_files=mzML_files,database_path=database_path,params_path=params_final)




#%% Post processing

kws.update({"output_folder":"output"})

Post_processing(input_files=final_annotations)
    
