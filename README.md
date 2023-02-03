# SCD_NGS_Analysis
 An extension of CRISPResso2 for SCD amplicon-seq
 
 The purpose of this package is to analyze next-generation sequencing data using Basespace and the Crispresso2 pipeline
 This pipeline will  genotype the beta-hemoglobin (HBB) gene and off-target alleles in genomic DNA collected from ficoll-separate blood, 
 bone marrow, or gene products from patients living with Sickle Cell Disease (SCD). 
 
 
 ## Requirements
 CRISPResso2, fastp and basespace
 
 The following commands can be used for installation of fastp and basespace. It is assumed CRISPResso2 is redily available as shared software 
 ```
 conda install -c bioconda fastp -y
 conda install -c hcc basespace-cli -y
 ```
 
 Please see these links for more information:
 
 [https://github.com/OpenGene/fastp]
 
 [http://crispresso2.pinellolab.org/submission]
 
 
 
 ## Installation
 git clone the SCD_NGS_repository using the following command
 
 ```
 git clone https://github.com/ACTG802/SCD_NGS_Analysis
 conda create -n Python_3_9 python=3.9
 conda activate Python_3_9
 installation.sh
 ```
 
 ## Flags and QC
 
 
 
 
 
