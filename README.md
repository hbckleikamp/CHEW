# CHEW


CHEW (cyclic hybrid/homology environmental workflow)
is a data analysis pipeline that facilitates accesible metaprotoemics database construction and annotation.



To improve accessibility CHEW is constructed as a collection of modular routines.
CHEW can be executed as a single script, or if only parts of the CHEW pipeline should be executed, subroutines can be executed instead.
For example, if the user only wishes to convert raw to mzML, then the subroutine raw2mzML can be executed.


To further improve accesibility, each subroutine can be executed in a variety of ways.
A subroutine can be run directly from a python interpreter (Such as Spyder), or a subroutine can be callied from a command line interface. This also means that users can decide which parts of the CHEW pipeline they use in their workflows.
For example when the user wants to construct a focused databse from their own custom peplist, if other softwares are used than SMSNet or MSFragger, they can simply run Peplist2finalDB or Peplist2Annotation with their own custom peplist, and only use the latter half of the CHEW pipeline.

Default arguments are listed within the subroutine itself.
These can be edited manually, or be supplied from a tabular parameter file.
Each of default arguments can also be supplied from the command line, using the same syntax.


CHEW consists of 3 main routines, which are subdivided into several subroutines.


1. Raw2peplist:
    - This routine constructs an initial peptide list, which is later used to construct a final database
    - The routine optionally uses SMSNet for de novo annotation, and MSFragger for database searching against a clustered database   
2. Peplist2finalDB:
    - This routine constructs a focused database from a user-supplied peplist, through an initial alignment step with DIAMOND, followed by database searching with MSFragger and filtering cycles.
    - The user wants to use different softwares to generate their own peptides, they can supply a peplist in fasta or tabular format
    - The focused database can be used to annotate other files, or can be used with different database searching tools.
3. Annotate_final
    - This routine uses a focused database to perform final annoatation with MSFragger
    - Resulting files can be filtered on false discovery rate, and will contain lowest common ancestor annotation
    - Additionally a fasta file of found proteins is constructed for later functional annotation by the user.
    - If de novo sequencing data is included with SMSNet, quality control plots are generated for the focused database.



![alt text](https://github.com/hbckleikamp/CHEW/blob/main/CHEW_workflow.PNG)

Fig. 1 CHEW workflow subroutines





### CHEW Setup

CHEW relies on external softwares for  raw file conversion (msconvert) spectrum annotation (MSFragger, SMSNet) and alignment (DIAMOND).
The folder Setup contains individual scripts to setup up CHEW on your own system.
The script `run_Setup.py` executes each setup script in their default mode.

Since MSFragger relies on an individual download liscence it should be manually downloaded:
download MSFragger 3.5 .jar file from http://msfragger-upgrader.nesvilab.org/upgrader/, move to CHEW directory and extract contents
CHEW is not compatible with mac/osx systems, as they require DIAMOND to be run externally.

msconvert must also be downloaded manually from https://proteowizard.sourceforge.io/download.html
and added to the default path of the operating system

To allow easy lowest common ancestor annotation and compatibility with GTDB and NCBI taxonomies,
the script `2_Create_taxdf.py` creates a rank normalized file containing all lineages fro GTDB and NCBI, and also flags NCBI taxids
based for Dump-taxa (unclassified or uncultured sequences)

Additionally, it needs a clustered database, an unclustered database, and a unclusterd diamond database.
To function properly, headers of the unclustered database need to be prepped to a CHEW compatible format,
which means the id of each fasta header should contain the "taxonomy id", as such >Accession|...|"taxonomy id".
This header renaming is automated during database prepping.

Downloading, prepping of databases and diamond database construction is done automatically with `4_Download_unclustered_database.py`
This script can download a variety of different databases UniProtKB (default) Swiss-Prot, TrEMBL, Uniref100,90,50 RefSeq, NCBI_NR and GTDB
The database prepping function is also included as standalone in `db_prepper.py` 
Additionally, there are several arguments supplied when prepping the CHEW database, which as default includes:
removal of Eukaryotic sequences, removal of Dump taxa, euqating I and J to L, removal of Ambiguous amino acids.
Eeach of these parameters can be changed manually, or be changed when running the setup scripts from the command line.

Disclaimer: <br>
While CHEW is designed to work with low memory requirements and high speed, it can write large temporary files to memory when using MSFragger and DIAMOND.
For example: MSFragger annotation of 10 files against a 20 GB database could lead to ~200 GB of temporary indices. 
It is therefore recommended to operate CHEW on a local desktop or laptop, instead of a computer cluster, which often have limited memory and can crash when full.
CHEW can be installed on an external HDD or SSD file to ensure sufficient storage space for temporary indices.
In each subroutine an alternative directory can be specified for writing temporary indices using the argument `-Temporary_directory "your directory"`

### Example CLI commands
Further documentation and syntax is supplied within each individual script.

Download RefSeq, prep database and make DIAMOND database: <br>
`4_Download_unclustered_database.py -DB RefSeq`

just Download GTDB, without preping database nor making DIAMOND database: <br>
`4_Download_unclustered_database.py -DB GTDB -make_dmnd 0 -prepdb 0` 

Prep existing database for CHEW, retain non-prokaryotic sequences and Dump taxids: <br>
`db_prepper.py "path to db" -BacArch_only 0 -No_Dump 0` 

Convert raw files to mzML: <br>
`raw2mzML -input_files "folder with raw files"`

Construct a peplist with only SMSNet: <br>
`raw2peplist.py -input_files "  your input folder  "  -database_path " your database "  -SMSNet 1`

Run database construction from peplist: <br>
`peplist2Final_DB.py -mzML_files "path" -peplist "path" -clustered_database_fasta "path" -unclustered_database_fasta "path" -unclustered_database_dmnd "path"`

Focus an existing database: <br>
`Refine_DB.py -input_files "  mzML folder  " -database_path " your database "`

Focus a database from a list of proteins or taxonomy ids: <br>
`Construct_DB.py -input_file "  tsv file with proteins or taxonomy ids  " -database_path " your database "` 

Run complete CHEW pipeline with MSFragger (db search) and SMSNet (de novo): <br>
`CHEW.py -input_files "folder with raw files"  -clustered_database_fasta "path" -unclustered_database_fasta "path" -unclustered_database_dmnd "path" -MSFragger 1 -SMSNet 1`

Run  CHEW pipeline with only MSFragger (db search) without final annotation (up to database construction): <br>
`CHEW.py -input_files "folder with raw files"  -clustered_database_fasta "path" -unclustered_database_fasta "path" -unclustered_database_dmnd "path" -MSFragger 1 -final_annotation 0`


#### Licensing

The pipeline is licensed with standard MIT-license. <br>
If you would like to use this pipeline in your research, please cite the following papers: 
      
- Kong, Andy T., et al. "MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics." Nature methods 14.5 (2017): 513-520.        

- Karunratanakul, Korrawe, et al. "Uncovering thousands of new peptides with sequence-mask-search hybrid de novo peptide sequencing framework." Molecular & Cellular Proteomics 18.12 (2019): 2478-2491.

- Buchfink, Benjamin, Chao Xie, and Daniel H. Huson. "Fast and sensitive protein alignment using DIAMOND." Nature methods 12.1 (2015): 59-60.



#### Contact:
-Hugo Kleikamp (Developer) @HugoKleikamp (Twitter)



#### Recommended links to other repositories:
https://github.com/bbuchfink/diamond<br>
https://github.com/Nesvilab/MSFragger<br>
https://github.com/cmb-chula/SMSNet
