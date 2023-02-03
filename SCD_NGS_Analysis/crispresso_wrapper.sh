#!/bin/bash
#SBATCH -p standard
#SBATCH -c 48
#SBATCH --job-name CRISPResso
#SBATCH -o %j.out
#SBATCH -e %j.err


##############################################################################################################
# EXECUTION NOTES:                                                                                           #
# - Run within data directory                                                                                #
# BEFORE RUNNING:                                                                                            #
# -Update everything under "Ref file Paths"                                                                  #
# -Check Fastq Nomeclature/Pattern under "GET SAMPLE NAMES" and "Setup file names"                           #
##############################################################################################################


## ACTIVATE CONDA
eval "$(conda shell.bash hook)"
conda activate Test_Py_3_9

pwd=$(pwd)
data_dir=$1

## PERMISSIONS
chmod -R 775 $pwd

## SETUP INPUTS, VARIABLES AND PATHS

# Variables
threads=8

# Directory paths
results_dir=${data_dir}Merged
qc_dir=${data_dir}QC_Metrics

mkdir -m777 $results_dir
mkdir -m777 $qc_dir

echo $data_dir
echo $results_dir

## GET SAMPLE NAMES
#samples=($(ls ${data_dir}/*_L001_R1_001.fastq.gz | sed -e 's/\_L001_R1_001.fastq.gz$//'))
samples=($(ls -R ${data_dir}*/*_L001_R1_001.fastq.gz | sed -e 's/\_L001_R1_001.fastq.gz$//'))
echo ${samples[*]}

## QC, TRIM, MERGE
for SAMPLE in "${samples[@]}"; do
        # Print Variables
        echo $SAMPLE
        base="$(basename -- $SAMPLE)"
        echo $base

        # Setup file names
        r1=${SAMPLE}_L001_R1_001.fastq.gz
        r2=${SAMPLE}_L001_R2_001.fastq.gz
        merged=${results_dir}/${base}_merged.fastq
        matches=${results_dir}/${base}_matches.txt

        # FASTQC
	echo "fastqc -t $threads -o $qc_dir ${r1} ${r2}"
        #fastqc -t $threads -o $qc_dir ${r1} ${r2}

        # Trim and QC:
        echo "fastp --in1 $r1 --in2 $r2 --detect_adapter_for_pe --merge --merged_out $merged -h $qc_dir/${base}_merged_fastp.html -j $qc_dir/${base}_merged_fastp.json --overrepresentation_analysis -A &> $qc_dir/${base}_merged_fastp.log"
	#fastp --in1 $r1 --in2 $r2 --detect_adapter_for_pe --merge --merged_out $merged -h $qc_dir/${base}_merged_fastp.html -j $qc_dir/${base}_merged_fastp.json --overrepresentation_analysis -A &> $qc_dir/${base}_merged_fastp.log


        # FASTQC
        echo "fastqc -t $threads -o $qc_dir ${merged}"
	#fastqc -t $threads -o $qc_dir ${merged}

done

#multiqc *

## CRISPRESSO
eval "$(conda shell.bash hook)"
conda activate crispresso2_env

python /groups/doudna/projects/CRISPResso_Wrapper/Scripts/run_agg_crispresso.py ${data_dir}
