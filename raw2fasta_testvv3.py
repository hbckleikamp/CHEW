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

#remove later
import warnings
warnings.filterwarnings("ignore")
#

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
kws.update({"taxdf":taxdf}) 

#%%



#%% Raw2peplist

kws.update({"output_folder":"peplist"}) #update output_subfolder

annotations=[]

mzML_files=raw2mzML() 


#annotate de novo
if SMSNet: 
    # mgf_files=raw2mgf()
    # SMSnet_files=SMSnet_annotation(input_files=mgf_files)
    SMSnet_files=["C:/MP-CHEW/CHEW/SMSnet/Q20516_Mix24X_SMSNET.tsv",
    "C:/MP-CHEW/CHEW/SMSnet/Q20517_Mix24X_160629055637_SMSNET.tsv",
    "C:/MP-CHEW/CHEW/SMSnet/Q20518_Mix24X_SMSNET.tsv"]
    annotations+=SMSnet_files
    
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
    final_annotations+=add_proteins_SMSnet(input_files=SMSnet_files,Alignment=Alignment)    

if MSFragger: 
    final_annotations+=MSFragger_annotation(input_files=mzML_files,database_path=final_database,params_path=params_final)




#%% Post processing


FDR=0.05
min_peptide_count=1
remove_unannotated=True
input_files=final_annotations

input_files.sort()



pin_files=[i for i in input_files if i.endswith(".pin")]
smsnet_files=[i for i in input_files if i.endswith("SMSNET.tsv")]

if len(smsnet_files):
    if len(pin_files)!=len(smsnet_files):
        print("warning: unequal amount of SMSnet_files and MSFragger files detected")
        

cols=['Peptide',"Proteins",'ScanNr','ExpMass','hyperscore','log10_evalue',"Prediction","mean_score","SMSNet_score", "evalue"] #retained columns
      

#make file pairs
pairs=[]
paired_up=set()
for file in input_files:
    
    pair=[]
    s=Path(file).stem
    for file in input_files:
        if s in file:
            if file not in paired_up:
            
                pair.append(file)
                paired_up.update([file])
    if len(pair):
        pairs.append(pair)
    


for pair in pairs:
    print(len(pair))


    pepdfs=[]
    pfile=""
    for file in pair:
        
        
        
        try:
            pepdf=pd.read_csv(file,sep="\t")
        except:
            pepdf=read_pin(file)
            
        if file.endswith(".pin"):
            pfile=file


        c=pepdf.columns
        if "Peptide" in c: pepdf["Peptide"]=pepdf["Peptide"].apply(lambda x: re.sub("[\[\[].*?[\]\]]", "", x)).str.split(".").apply(lambda x: x[1]).str.replace("I","L").str.replace("J","L") #remove ptms in peptides
        if "Proteins" in c: pepdf["Proteins"]=pepdf["Proteins"].str.replace("\t"," ",regex=True)
        if "tag" in c: pepdf["Peptide"]=pepdf["tag"].str.replace("I","L").str.replace("J","L") #remove ptms in peptides
        if "Scores" in c: pepdf["SMSNet_score"]=pepdf["Scores"].str.rsplit(";",expand=True).astype(float).sum(axis=1)
        if "ScanNum" in c: pepdf["ScanNr"]=pepdf["ScanNum"]
        if 'ObservedM+H' in c: pepdf["ExpMass"]=pepdf['ObservedM+H']-1.007276466621 #-H+
        if  "log10_evalue" in c: pepdf["evalue"]=10**pepdf.log10_evalue.astype(float)
        
        c=pepdf.columns
        
        pepdfs.append(pepdf[[i for i in cols if i in c]].reset_index())

    pepdfs=pd.concat(pepdfs,axis=0)
    pepdfs["Decoy"]=pepdfs["Proteins"].str.contains("decoy_").fillna(False)
    c=pepdfs.columns
    pepdfs["ScanNr"]=pepdfs["ScanNr"].astype(int)
    for i in ["ExpMass","hyperscore","evalue"]:
        if i in c:
            pepdfs[i]=pepdfs[i].astype(float) 
  

    

    from sklearn.linear_model import Lasso
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import r2_score
    import matplotlib.pyplot as plt

    #De novo and non-denovo score correlation

    if ("hyperscore" in c) and ("SMSNet_score" in c):
    
        mf_pepdf=pepdfs.loc[pepdfs.hyperscore.notnull(),["Peptide","ScanNr","hyperscore"]]
        sm_pepdf=pepdfs.loc[pepdfs.SMSNet_score.notnull(),["Peptide","ScanNr","SMSNet_score"]]
        shared_peptides=mf_pepdf.merge(sm_pepdf,on=["Peptide","ScanNr"],how="inner").drop_duplicates()
        shared_peptides=shared_peptides.groupby('ScanNr').apply(lambda x: x.sort_values(by="hyperscore",ascending=False).iloc[0]).reset_index(drop=True)

        x=shared_peptides["SMSNet_score"].values.reshape(-1,1)
        y=shared_peptides["hyperscore"].astype(float).values
        
        #lasso regression
        X_train, X_test, y_train, y_test = train_test_split(x,y,test_size=0.10,random_state=42)
        lasso_model = Lasso().fit(X_train,y_train)
        y_pred = lasso_model.predict(X_test)
        r2=round(r2_score(y_test, y_pred),2)
        

        pepdfs.loc[pepdfs["hyperscore"].isnull(),"hyperscore"]=lasso_model.predict(pepdfs["SMSNet_score"].dropna().values.reshape(-1,1))

        #plotting 
        fig,ax=plt.subplots()
        plt.xlabel("SMSnet_score")
        plt.ylabel("hyperscore")
        plt.scatter(x,y,s=0.1)
        plt.title(Path(pfile).stem+" "+"r2: "+str(r2))
        plt.savefig(basepath+"_dn_corr.png",dpi=400)

        #plot high scoring?
  
    if not pfile:
        pfile=file  
  
    basepath=str(Path(Path(pfile).parents[0],Path(pfile).stem))  
  
    #Filtering
    if remove_unannotated: pepdfs=pepdfs[pepdfs["Proteins"].notnull()]
    if "evalue" in c and max_evalue: pepdfs=pepdfs[pepdfs["evalue"]<=max_evalue]
    if "hyperscore" in c and Top_score_fraction: pepdfs=pepdfs[(pepdfs["hyperscore"]/pepdfs.groupby(["ScanNr","ExpMass"])["hyperscore"].transform('max'))<=(1/Top_score_fraction)]
    if min_peptide_count: 
    
        s=pepdfs.groupby("Peptide").size()
        pepdfs=pepdfs[pepdfs["Peptide"].isin(s[s>=min_peptide_count].index)]
    
    if "hyperscore" in c: FDR_col="hyperscore"
    else: FDR_col="SMSNet_score"
    
    
    if FDR:
  

        pepdfs=pepdfs.sort_values(by=FDR_col,ascending=False)
        pepdfs["FDR"]=pepdfs["Decoy"].cumsum()/np.arange(1,len(pepdfs)+1)
        
        #write hits
        d=pepdfs["Decoy"].sum()
        t=(~pepdfs["Decoy"]).sum()
        pepdfs=pepdfs[pepdfs["FDR"]<=FDR]
        du=pepdfs["Decoy"].sum()
        tu=(~pepdfs["Decoy"]).sum()
        ddf=pd.DataFrame([[d,t],[du,tu]],columns=["Decoy","Target"],index=["pre","post"])
        ddf.to_csv(basepath+"_performance.tsv",sep="\t")

    edf=pepdfs.copy()
    edf["Proteins"]=edf["Proteins"].str.split()
    edf=edf.explode("Proteins")
    edf["taxids"]=edf["Proteins"].str.split("|").apply(lambda x: x[-1])
    

 
    edf[taxdf.columns]=taxdf.loc[edf["taxids"].tolist()].values
    lins=edf[ranks].reset_index().fillna("")
    

    pepdfs[ranks]=lins.groupby("index")[ranks].nth(0)[(lins.groupby("index")[ranks].nunique()==1)].fillna("") #vectorized lca        
    print(pepdfs.groupby("genus").size())
    pepdfs.to_csv(basepath+"_PSMs",sep="\t")
