#### IMPORT ####
from datetime import datetime
import pandas as pd
import numpy as np
import glob
import os
import subprocess
import sys


#### SETUP ####
# CLI INPUTS
data_dir = sys.argv[1]#"/groups/doudna/projects/CRISPResso_Wrapper/Data/"
amplicon_data = "/groups/doudna/projects/CRISPResso_Wrapper/Data/Test_Data/Sample_Data_Info.xlsx"

# GET AMPLICON DETAILS
amplicon_df = pd.read_excel(amplicon_data, header = None, names = ["Amplicon_ID", "Amplicon"])
amplicon_df["Amplicon_ID"] = amplicon_df["Amplicon_ID"].apply(lambda x: x.replace(" ", "_"))
amplicon_dict = dict(zip(amplicon_df.Amplicon_ID, amplicon_df.Amplicon))

# PATHS, FOLDERS AND FILES
pwd = os.getcwd()
fasta_subfolders = [x for x in os.listdir(data_dir) if ".xlsx" not in x]
merge_dir = data_dir + "Merged/" # pwd + "Merged/"
r1_suffix = "_L001_R1_001.fastq.gz"
fastas = glob.glob(data_dir + "**/*gz", recursive=True)
prefixes = [x.replace(data_dir, "").replace(r1_suffix, "").split("/")[-1] for x in fastas if r1_suffix in x]


#### RUN CRISPRESSO ####

for prefix in prefixes:
    print(prefix)
    merged = merge_dir + prefix + "_merged.fastq"
    
    if "HBB" in prefix:
        amplicon = amplicon_dict['HBB_Amplicon']
        hdr = amplicon_dict['HBB_HDR']
        guide = amplicon_dict['HBB_sgRNA_sequence']
        cli =  f"CRISPResso --fastq_r1 {merged} --amplicon_seq {amplicon} -n {prefix} " +\
                 f"-g {guide} -e {hdr}"
        print(cli)
        #subprocess.run(cli, shell = True)
    elif "OT1" in prefix:
        amplicon = amplicon_dict['OT1_Amplicon']
        guide = amplicon_dict['OT1_sgRNA_sequence']

        cli = f"CRISPResso --fastq_r1 {merged} --amplicon_seq {amplicon} -n {prefix} " +\
                 f"-g {guide}"
        print(cli)
        #subprocess.run(cli, shell = True)

#### AGGREGATE RESULTS ####
#list folders  of each sample type
hbb_dirs = [data_dir+x for x in os.listdir(data_dir) if "HBB" in x and "Water" not in x and "html" not in x]
ot1_dirs = [data_dir+x for x in os.listdir(data_dir) if "OT1" in x and "Water" not in x and "html" not in x]
hbb_water = [data_dir+x for x in os.listdir(data_dir) if "HBB" in x and "Water" in x and "html" not in x]
ot1_water = [data_dir+x for x in os.listdir(data_dir) if "OT1" in x and "Water" in x and "html" not in x]
failed_runs = []
hbb_rows = []
ot1_rows = []
hbb_ctrls = []
ot1_ctrls = []
ctrl_names = []

## HBB
for hbb in hbb_dirs:
    run = "_".join(hbb.split("/")[-1].replace("CRISPResso_on_", "").split("_")[:-1])
    crispresso_files = [x for x in os.listdir(hbb)]
    print(hbb, run)
    if "CRISPResso_quantification_of_editing_frequency.txt" in crispresso_files:
        summary_df = pd.read_csv(f"{hbb}/CRISPResso_quantification_of_editing_frequency.txt", sep ="\t")
        # HBB HDR
        hdr = summary_df.loc[summary_df.Amplicon == "HDR"]
        hdr_percent = float(hdr.Unmodified/hdr.Reads_aligned_all_amplicons*100)
        # HBB INDELS
        ref =  summary_df.loc[summary_df.Amplicon == "Reference"]
        indel_percent = float(ref.Modified/ref.Reads_aligned_all_amplicons*100)
        #unmod_percent = ref["Unmodified%"]
        unmod_percent = float(ref.Unmodified/ref.Reads_aligned_all_amplicons*100)
        hbb_row = [run, hdr_percent, indel_percent, unmod_percent,hdr.Unmodified,ref.Modified,ref.Unmodified]
        hbb_rows.append(hbb_row)
    else:
        print(f"!!!!! {run} RUN FAILED !!!!!")
        failed_runs.append(run)
        
        
#for c in hdr.columns:
#    print(c)
 #   print(hdr[c])
    
    
## OT1 INDELS
for ot1 in ot1_dirs:
    run = "_".join(ot1.split("/")[-1].replace("CRISPResso_on_", "").split("_")[:-1])
    crispresso_files = [x for x in os.listdir(ot1)]
    print(ot1, run)
    if "CRISPResso_quantification_of_editing_frequency.txt" in crispresso_files:
        summary_df = pd.read_csv(f"{ot1}/CRISPResso_quantification_of_editing_frequency.txt", sep ="\t")
        # OT1 INDELS
        ref =  summary_df.loc[summary_df.Amplicon == "Reference"]
        indel_percent = float(ref.Modified/ref.Reads_aligned_all_amplicons*100)
        #unmod_percent = ref["Unmodified%"]
        unmod_percent = float(ref.Unmodified/ref.Reads_aligned_all_amplicons*100)
        ot1_row = [run, hdr_percent, indel_percent, unmod_percent,hdr.Unmodified,ref.Modified,ref.Unmodified]
        ot1_rows.append(ot1_row)
    else:
        print(f"!!!!! {run} RUN FAILED !!!!!")
        failed_runs.append(run)
        
## HBB WATER CTRL
for hbb_w in hbb_water:
    run = "_".join(hbb_w.split("/")[-1].replace("CRISPResso_on_", "").split("_")[:-1])
    crispresso_files = [x for x in os.listdir(hbb_w)]
    print(hbb_w, run)
    ctrl_names.append(run)
    if "CRISPResso_quantification_of_editing_frequency.txt" in crispresso_files:
        summary_df = pd.read_csv(f"{hbb_w}/CRISPResso_quantification_of_editing_frequency.txt", sep ="\t")
        # INDELS
        ref =  summary_df.loc[summary_df.Amplicon == "Reference"]
        indel_percent = float(ref["Modified%"].values[0])
        #unmod_percent = ref["Unmodified%"]
        unmod_percent = float(ref.Unmodified/ref.Reads_aligned_all_amplicons*100)
        hbb_w_row = [run, hdr_percent, indel_percent, unmod_percent,hdr.Unmodified,ref.Modified,ref.Unmodified]
        hbb_ctrls.append(hbb_w_row)
    else:
        print(f"!!!!! {run} RUN FAILED !!!!!")
        failed_runs.append(run)
        
## OT1 WATER CTRL
for ot1_w in ot1_water:
    run = "_".join(ot1_w.split("/")[-1].replace("CRISPResso_on_", "").split("_")[:-1])
    crispresso_files = [x for x in os.listdir(ot1_w)]
    print(ot1_w, run)
    ctrl_names.append(run)
    if "CRISPResso_quantification_of_editing_frequency.txt" in crispresso_files:
        summary_df = pd.read_csv(f"{ot1_w}/CRISPResso_quantification_of_editing_frequency.txt", sep ="\t")
        # INDELS
        ref =  summary_df.loc[summary_df.Amplicon == "Reference"]
        indel_percent = float(ref["Modified%"].values[0])
        #unmod_percent = ref["Unmodified%"]
        unmod_percent = float(ref.Unmodified/ref.Reads_aligned_all_amplicons*100)
        ot1_w_row = [run, hdr_percent, indel_percent, unmod_percent,hdr.Unmodified,ref.Modified,ref.Unmodified]
        ot1_ctrls.append(ot1_w_row)
    else:
        print(f"!!!!! {run} RUN FAILED !!!!!")
        failed_runs.append(run)
        
hdr_df = pd.DataFrame(hbb_rows, columns = ["HBB_Sample", 
                                           "HBB_HDR_Percent", 
                                           "HBB_Indel_Percent", "HBB_Unmodified_Percent", 
                                           'HBB_HDR_N', 'HBB_Indel_N', 'HBB_Unmodified_N']).sort_values("HBB_Sample")
hdr_df["Condition"] = hdr_df.HBB_Sample.apply(lambda x: x[:-4])
try:
    hbb_ctrl = [x for x in ctrl_names if "HBB" in x][0]
    if hbb_ctrl in failed_runs:
        hdr_df[hbb_ctrl] = "Passed: No Amplicon Detected"
    else:
        hdr_df[hbb_ctrl] = "Failed: Correct Amplicon Detected"
    col_reorder = ['Condition', 'HBB_Sample', 'HBB_HDR_Percent', 'HBB_Indel_Percent', 'HBB_Unmodified_Percent',
                   'HBB_HDR_N', 'HBB_Indel_N', 'HBB_Unmodified_N',
               'WaterControl-HBB', 'OT1_Sample', 'OT1_Indel_Percent', 'OT1_Unmodified_Percent', 
               'OT1_Indel_N', 'OT1_Unmodified_N','WaterControl-OT1']
except IndexError:
    print("!!!!! NO HBB CONTROL DETECTED !!!!!")
    col_reorder = ['Condition', 'HBB_Sample', 'HBB_HDR_Percent', 'HBB_Indel_Percent', 'HBB_Unmodified_Percent',
                   'HBB_HDR_N', 'HBB_Indel_N', 'HBB_Unmodified_N',
            'OT1_Sample', 'OT1_Indel_Percent', 'OT1_Unmodified_Percent',
          'OT1_Indel_N', 'OT1_Unmodified_N']
    pass 

ot1_df = pd.DataFrame(ot1_rows, columns = ["OT1_Sample", "OT1_Indel_Percent", 
                                           "OT1_Unmodified_Percent",'OT1_Indel_N', 'OT1_Unmodified_N']).sort_values("OT1_Sample")
ot1_df["Condition"] = ot1_df.OT1_Sample.apply(lambda x: x[:-4])
try:
    ot1_ctrl = [x for x in ctrl_names if "OT1" in x][0]
    if ot1_ctrl in failed_runs:
        ot1_df[ot1_ctrl] = "Passed: No Amplicon Detected"
    else:
        ot1_df[ot1_ctrl] = "Failed: Correct Amplicon Detected"
        col_reorder = ['Condition', 'HBB_Sample', 'HBB_HDR_Percent', 'HBB_Indel_Percent', 'HBB_Unmodified_Percent',
                       'HBB_HDR_N', 'HBB_Indel_N', 'HBB_Unmodified_N',
                   'WaterControl-HBB', 'OT1_Sample', 'OT1_Indel_Percent', 'OT1_Unmodified_Percent', 
                   'OT1_Indel_N', 'OT1_Unmodified_N','WaterControl-OT1']
except IndexError:
    print("!!!!! NO OT1 CONTROL DETECTED !!!!!")
    col_reorder = ['Condition', 'HBB_Sample', 'HBB_HDR_Percent', 'HBB_Indel_Percent', 'HBB_Unmodified_Percent',
                   'HBB_HDR_N', 'HBB_Indel_N', 'HBB_Unmodified_N',
               'WaterControl-HBB', 'OT1_Sample', 'OT1_Indel_Percent', 'OT1_Unmodified_Percent', 
               'OT1_Indel_N', 'OT1_Unmodified_N','WaterControl-OT1']
    pass

##output of all data
all_data = pd.merge(hdr_df, ot1_df, on="Condition")
#all_data = pd.concat([hdr_df, ot1_df], axis = 0)
#all_data = all_data[col_reorder]
#all_data.drop("Condition",  axis=1, inplace = True)
t#ime = datetime.now().strftime("%m_%d_%Y_%Hhr_%Mmin_%Ssec")
all_data.to_excel(f"{data_dir}Results_Summary_{time}.xlsx", index = False)
