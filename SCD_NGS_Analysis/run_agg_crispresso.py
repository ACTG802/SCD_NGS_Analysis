'''

Created by M. Trinidad
Edits T. Hudson
Last updated: 2/06/23
'''


#### IMPORT ####
from datetime import datetime
import pandas as pd
import glob
import os
import sys
from fpdf import FPDF


#### SETUP ####
# CLI INPUTS
data_dir = sys.argv[1]#'C:\\Users\\thuds\\Dropbox\\SCD_NGS_Assay\\SCD_17_Results\\'

amplicon_data = 'Sample_Data_Info.xlsx'

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

def agg_results(amplicon_type, data_dir):
    '''

    Parameters
    ----------
    amplicon_type : 'HBB' or 'OT1'
        Amplicons types found in CRISPResso Directory
    data_dir : Directory of CRISPResso output


    Returns
    -------
    df : pandas dataframe
        dataframe of all samples for amplicon specified
        Includes total reads and percent edited
    ct_dict : dictionary
        Displays the water controls found and whether amplicons were detected
    failed_runs : list
        All runs that did not output crisspresso files
      .
    min_error : list
        List of samples that total aligned reads (or sequence depth) did not reach 1000reads.

    '''
    
    dirs = [data_dir+x for x in os.listdir(data_dir) if amplicon_type in x and "Water" not in x and "H2O" not in x and "html" not in x]
    water_control =  [data_dir+x for x in os.listdir(data_dir) if amplicon_type in x and "Water" in x  and "html" not in x]
    if len(water_control)==0:
        water_control =  [data_dir+x for x in os.listdir(data_dir) if amplicon_type in x and "H2O" in x  and "html" not in x]
        
    rows = []
    ct_dict = {}
    failed_runs = []
    min_error = []
    ctrl_names = []
    
    for i in dirs:
        run = i[i.index("CRISPResso_on_")+len("CRISPResso_on_"):]
        crispresso_files = [x for x in os.listdir(i)]
        if "CRISPResso_quantification_of_editing_frequency.txt" in crispresso_files:
            summary_df = pd.read_csv(f"{i}/CRISPResso_quantification_of_editing_frequency.txt", sep ="\t")
            
            ref =  summary_df.loc[summary_df.Amplicon == "Reference"]
            tot_reads_aligned = float(ref.Reads_aligned_all_amplicons)
            unmod_percent = float(ref.Unmodified/ref.Reads_aligned_all_amplicons*100)
            indel_percent = float(ref.Modified/ref.Reads_aligned_all_amplicons*100)
            
            # INDELS
            if amplicon_type == 'HBB':
                hdr = summary_df.loc[summary_df.Amplicon == "HDR"]
                hdr_percent = float(hdr.Unmodified/tot_reads_aligned*100)
                
                row = [run, round(hdr_percent,2), indel_percent, unmod_percent,hdr.Unmodified.iloc[0],
                       ref.Modified.iloc[0],ref.Unmodified.iloc[0], tot_reads_aligned]
                 
            if amplicon_type == 'OT1':
                
                row = [run, indel_percent, unmod_percent,ref.Modified.iloc[0],ref.Unmodified.iloc[0], tot_reads_aligned]
            if tot_reads_aligned < 1000:
                min_error.append(run)
            rows.append(row)
        else:
            print(f"!!!!! {run} RUN FAILED !!!!!")
            print('CRISPResso files not present in directory')
            failed_runs.append(run)
            
    ## WATER CTRL
    for w in water_control:
        run = w[w.index("CRISPResso_on_")+len("CRISPResso_on_"):]
        crispresso_files = [x for x in os.listdir(w)]
        ctrl_names.append(run)
        
        if "CRISPResso_quantification_of_editing_frequency.txt" in crispresso_files:
            summary_df = pd.read_csv(f"{w}/CRISPResso_quantification_of_editing_frequency.txt", sep ="\t")
            
            # INDELS
            ref =  summary_df.loc[summary_df.Amplicon == "Reference"]
            tot_reads_aligned = float(ref.Reads_aligned_all_amplicons)
            indel_percent = float(ref["Modified%"].values[0])
            
            #unmod_percent = ref["Unmodified%"]
            unmod_percent = float(ref.Unmodified/ref.Reads_aligned_all_amplicons*100)
            row = [run, 'nan', indel_percent, unmod_percent, 'nan' , ref.Modified.iloc[0],ref.Unmodified.iloc[0], tot_reads_aligned]
            rows.append(row)
        else:
            print(f'{run} - No Amplicon Detected')
            failed_runs.append(run)
            
    ##make dataframe
    if amplicon_type == 'HBB':
        
        df = pd.DataFrame(rows, columns = ["HBB_Sample", 
                                               "HBB_HDR_Percent", 
                                               "HBB_Indel_Percent", "HBB_Unmodified_Percent", 
                                               'HBB_HDR_N', 'HBB_Indel_N', 'HBB_Unmodified_N', 'HBB_Total_Reads']).sort_values("HBB_Sample")
    if amplicon_type == 'OT1':
        df = pd.DataFrame(rows, columns = ["OT1_Sample", "OT1_Indel_Percent", 
                                                   "OT1_Unmodified_Percent",'OT1_Indel_N', 'OT1_Unmodified_N', 'OT1_Total_Reads']).sort_values("OT1_Sample")
    
    try:
        
        for ctrl in ctrl_names:
            
            if ctrl in failed_runs:
                
                ct_dict[ctrl] = "Passed: No Amplicon Detected"
                
            else:
                ct_dict[ctrl]  = "Failed: Correct Amplicon Detected"

    except IndexError:
        print(f"!!!!! NO {amplicon_type} CONTROL DETECTED !!!!!")
        pass 


    return df, ct_dict, failed_runs, min_error

hdr_df, hbb_ct_dict, hbb_failed_runs, hbb_min_error = agg_results(amplicon_type = 'HBB', data_dir = data_dir)
ot1_df, ot1_ct_dict, ot1_failed_runs, ot1_min_error = agg_results(amplicon_type = 'OT1', data_dir = data_dir)
min_error = hbb_min_error + ot1_min_error

def sort_samples(df):
    try:
        
       names = [name[name.index('_S')+2:] for name in df.iloc[:,0] if name]
       df['sorted'] = names
       df = df.sort_values('sorted')
       df = df.drop(columns = 'sorted')
    except:
         None
    return df

hdr_df = sort_samples(hdr_df)
ot1_df = sort_samples(ot1_df)
all_data = pd.concat([hdr_df, ot1_df], axis = 0)
time = datetime.now().strftime("%m_%d_%Y_%Hhr_%Mmin_%Ssec")
all_data.to_excel(f"{data_dir}Results_Summary_{time}.xlsx", index = False)



##-----------------------------
#Generate PDF


hbbdf = pd.DataFrame()
ot1df = pd.DataFrame()

x,y = len(hdr_df), len(ot1_df)

hbbdf['Sample Name'] = hdr_df['HBB_Sample']
ot1df['Sample Name'] = ot1_df['OT1_Sample']
hbbdf['HDR % Aligned (# of reads)'] = [f'{hdr_df.HBB_HDR_Percent.iloc[row]} %({hdr_df.HBB_HDR_N.iloc[row]})' for row in range(x)]

hbbdf['NHEJ % Aligned (# of reads)'] = [f'{round(hdr_df.HBB_Indel_Percent.iloc[row],2)} %({hdr_df.HBB_Indel_N.iloc[row]})' for row in range(x)]
ot1df['NHEJ % Aligned (# of reads)'] = [f'{round(ot1_df.OT1_Indel_Percent.iloc[row],2)} %({ot1_df.OT1_Indel_N.iloc[row]})' for row in range(y)]

hbbdf['Unmodified % Aligned (# of reads)'] = [f'{round(hdr_df.HBB_Unmodified_Percent.iloc[row],2)} %({hdr_df.HBB_Unmodified_N.iloc[row]})' for row in range(x)]
ot1df['Unmodified % Aligned (# of reads)'] = [f'{round(ot1_df.OT1_Unmodified_Percent.iloc[row],2)} %({ot1_df.OT1_Unmodified_N.iloc[row]})' for row in range(y)]

        
    
    
def draw_table(df,pw,y):
    pdf.set_font('Times', style = 'B', size = 9)
    line_height = pdf.font_size * 2.5
    colw = pw/df.shape[1]
    
    for c in df.columns:
        pdf.cell(w =colw, h = line_height, txt=c, border=1, align = 'C')
    y+=line_height
    pdf.set_xy(ch,y)

    
    pdf.set_font('Times', style = '', size = 9)
    for i in df.index:
        for c in df.columns:
            pdf.cell(w =colw, h = line_height, txt=df[c].iloc[i], border=1, align = 'C')
        y+=line_height
        pdf.set_xy(ch,y)
        
    #y = og_y +(( df.shape[0]) * m)
    #if y > pdf.h:
     #   y = y-pdf.h
    return y
    
time = datetime.now().strftime("%m/%d/%Y")
# Margin
m = 10
# Page width: Width of A4 is 210mm
pw = 210 - 2*m
# Cell height
ch = 10

pdf = FPDF()
pdf.add_page()
x,y = m,m
pdf.set_xy(m, m)
pdf.set_font('Times', 'B', 12)

#title
pdf.cell(w=0, h=ch, txt="Sickele Cell Amplicon-seq Summary Report", align = 'C')
pdf.set_font('Times','', 12)

#date
y = y +m
pdf.set_xy(ch,y)
pdf.cell(w = 0,h = ch, txt=f'Report Date:  {time}' )

# HBB Table
pdf.set_font('Times', '', 12)
y = y + m
pdf.set_xy(ch,y)
pdf.cell(w = 0, h = ch, txt = 'HBB Summary of Aligned Reads')
y = y+m
pdf.set_xy(ch,y)
y = draw_table(df=hbbdf,pw=pw,y=y)


#OT1 report
pdf.set_font('Times', '', 12)
pdf.set_xy(ch,y)
pdf.cell(w = 0, h = ch, txt = 'OT1 Summary of Aligned Reads')
y = y+ch
pdf.set_xy(ch,y)
y =draw_table(df=ot1df,pw=pw,y=y)

#error and QC
pdf.set_font('Times', '', 12)
ch = pdf.font_size * 2.5
pdf.set_xy(10,pdf.y)

txt = 'Run Warnings and QC'
lines = 1
if len(min_error) == 0:
     txt = txt + '\n Passed: All Samples Met Minimum Read Requirments'
     lines +=.1
else:
    for r in min_error:
         txt= txt +  f'\n Warning: {r} did not meet a minimum of 1000 aligned reads'
         lines +=.1

for k,v in hbb_ct_dict.items():
    txt= txt +  f'\n{k}: {v}'
    lines +=.1

for k,v in ot1_ct_dict.items():
    txt= txt +  f'\n{k}: {v}'
    lines +=.1

h = pdf.font_size * lines
pdf.multi_cell(w = 0 , h = h, txt = txt)

pdf.output('test.pdf', 'F')