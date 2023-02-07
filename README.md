# SCD_NGS_Analysis
 An extension of CRISPResso2 for SCD amplicon-seq
 
 The purpose of this package is to analyze next-generation sequencing data using Basespace and the Crispresso2 pipeline
 This pipeline will  genotype the beta-hemoglobin (HBB) gene and off-target alleles in genomic DNA collected from ficoll-separate blood, 
 bone marrow, or gene products from patients living with Sickle Cell Disease (SCD). 
 

 ## Installation
 git clone the SCD_NGS_repository into a designated folder using the following command
 
 ```
 git clone https://github.com/ACTG802/SCD_NGS_Analysis
 ```
 
 Create a conda environment and import all required packages in one line below:
 
 ```
 conda env create -n Python_3_9 -f environment.yaml
 ```
 
 These two lines of code eliminates the need to conda install individual packages and copy script files into a folder in biospace one-by-one.
 
 ## How to run an analysis
 
 Follow the steps in the SCD ipeline plan of making a directory, import fastq files from basespace and running using the commmand
 
 ```
 sbatch  /groups/clinical/projects/crispresso/crispresso_wrapper.sh /groups/clinical/projects/run_folder
 ```
 
 Make sure fastq file names are consistent. Control samples are detected by containg the words 'Water', 'H2O' or 'NTC'. The amplicons of the HBB gene must contains 'HBB' in the file name and likewise for 'OT1' for the OT1 samples.
 
 ## Output
 
 
 ## Flags and QC
 
 
 
 
 
