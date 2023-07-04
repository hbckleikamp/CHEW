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
3. Annotate_final
    - This routine uses a focused database to perform final annoatation with MSFragger
    - Resulting files can be filtered on false discovery rate, and will contain lowest common ancestor annotation
    - Additionally a fasta file of found proteins is constructed for later functional annotation by the user.
    - If de novo sequencing data is included with SMSNet, quality control plots are generated for the focused database.



![alt text](https://github.com/hbckleikamp/CHEW/blob/main/CHEW_workflow.PNG)
Fig. 1 CHEW workflow subroutines





### Running Setup



### Example CLI commands
