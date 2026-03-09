#!/usr/bin/env python

import os
import re
import csv
import sys
import glob
import argparse
import pandas as pd

def get_sample_list(directory):
    dir_names = [x for x in os.listdir(directory) if os.path.isdir(os.path.join(directory, x))]
    sample_names = []
    pattern = r'(^\S+)_[ATCG]{10}'
    for dir_name in dir_names:
        match = re.match(pattern, dir_name)
        if match:
            sample_names.append(match.group(1))
    return sample_names

parser = argparse.ArgumentParser(description='Make QC Excel Spreadsheet')
parser.add_argument('-d','--batchdir',required=True,help='workflow batch output dir')
parser.add_argument('-s','--samplesheet',required=False,help='input sample spreadsheet with 2 tabs')

args = parser.parse_args()

in_ss      = args.samplesheet
out_dir    = os.path.abspath(args.batchdir)
batch_name = os.path.basename(out_dir)

if not os.path.isdir(out_dir):
    sys.exit("Workflow output dir does not exist: " + out_dir)

if in_ss:
    in_df = pd.read_excel(in_ss, sheet_name='QC Metrics')
    sample_list = in_df['SAMPLE ID'].tolist()
else:
    sample_list = get_sample_list(out_dir)
    in_df = pd.DataFrame({'Samples' : sample_list})

hap_scores      = []
hap_sites       = []
total_bases     = []
total_reads     = []
pct_map_reads   = []
error_rate_1    = []
error_rate_2    = []
pct_q30_1       = []
pct_q30_2       = []
min_ins_size    = []
pct_umi_dup     = []
avg_align_cov   = []
pct_align_reads = []
pct_target_20   = []
pct_target_100  = []
pct_target_1500 = []
total_giga_bases= []

for sample_name in sample_list:
    sample_name = sample_name.replace(" ", "")
    search = os.path.join(out_dir, f"{sample_name}*")
    sample_dir = glob.glob(search)[0]
    if not os.path.isdir(sample_dir):
        sys.exit(sample_dir + " is not a valid sample directory")

    dragen_dir = os.path.join(sample_dir, "dragen")
    mapping_metrics = glob.glob(os.path.join(dragen_dir, "*.mapping_metrics.csv"))[0]
    target_metrics = glob.glob(os.path.join(dragen_dir, "*.target_bed_coverage_metrics.csv"))[0]
    umi_metrics = glob.glob(os.path.join(dragen_dir, "*.umi_metrics.csv"))[0]
    haplotect_out = os.path.join(sample_dir, f"{sample_name}.haplotect.txt")

    if not (os.path.isfile(mapping_metrics) and os.path.isfile(target_metrics) and os.path.isfile(umi_metrics) and os.path.isfile(haplotect_out)):
        sys.exit(f"No dragen mapping and/or target metrics and/or umi metrics and/or haplotect out for {sample_dir}")

    with open(haplotect_out, "r") as haplotect_file:
        lines = haplotect_file.readlines()
        last_line = lines[-1]
        cols = last_line.split()
        hapscore = cols[6]
        if hapscore == 'NaN':
            hapscore = 0
        hap_scores.append(float(hapscore))
        hap_sites.append(int(cols[2]))

    with open(mapping_metrics, "r") as map_file:
        map_reader = csv.reader(map_file)
        for line in map_reader:
            if "MAPPING/ALIGNING PER RG" in line:
                break
            if "Total input reads" in line:
                total_reads.append(int(line[3]))
            elif "Mapped reads" in line:
                pct_map_reads.append(float(line[4]))
            elif "Total bases" in line:
                total_bases.append(int(line[3]))
                total_giga_bases.append(float(round(int(line[3])/1000000000, 2)))
            elif "Mismatched bases R1" in line:
                error_rate_1.append(float(line[4]))
            elif "Mismatched bases R2" in line:
                error_rate_2.append(float(line[4]))
            elif "Q30 bases R1" in line:
                pct_q30_1.append(float(line[4]))
            elif "Q30 bases R2" in line:
                pct_q30_2.append(float(line[4]))
            elif "Insert length: mean" in line:
                min_ins_size.append(float(line[3]))

    with open(target_metrics, "r") as target_file:
        target_reader = csv.reader(target_file)
        for line in target_reader:
            if "Average alignment coverage over target region" in line:
                avg_align_cov.append(float(line[3]))
            elif "PCT of target region with coverage [1500x: inf)" in line:
                pct_target_1500.append(float(line[3]))
            elif "PCT of target region with coverage [ 100x: inf)" in line:
                pct_target_100.append(float(line[3]))
            elif "PCT of target region with coverage [  20x: inf)" in line:
                pct_target_20.append(float(line[3]))
            elif "Aligned reads in target region" in line:
                pct_align_reads.append(float(line[4]))

    with open(umi_metrics, "r") as umi_file:
        umi_reader = csv.reader(umi_file)
        num_of_reads = 0
        num_of_cons  = 0
        for line in umi_reader:
            if "Number of reads" in line:
                 num_of_reads = int(line[3])
            if "Consensus pairs emitted" in line:
                 num_of_cons  = int(line[3])
                 break
        pct_umi_dup.append(float(round(100 - (num_of_cons * 2)/num_of_reads * 100, 2)))

qc_df = pd.DataFrame({
    "HAPLOTECT_SCORE"           : hap_scores, 
    "HAPLOTECT_SITES"           : hap_sites,
    "TOTAL_GIGA_BASES"          : total_giga_bases,
    "TOTAL_READS"               : total_reads,
    "PCT_MAPPED_READS"          : pct_map_reads,
    "MISMATCH_RATE_1"           : error_rate_1,
    "MISMATCH_RATE_2"           : error_rate_2,
    "PCT_Q30_BASES_1"           : pct_q30_1,
    "PCT_Q30_BASES_2"           : pct_q30_2,
    "MEAN_INSERT_SIZE"          : min_ins_size,
    "PCT_UMI_DUPLICATE_READS"   : pct_umi_dup,
    "AVG_ALIGN_TARGET_COVERAGE" : avg_align_cov,
    "PCT_TARGET_ALIGNED_READS"  : pct_align_reads,
    "PCT_TARGET_20x"            : pct_target_20,
    "PCT_TARGET_100x"           : pct_target_100,
    "PCT_TARGET_1500x"          : pct_target_1500
})

out_file1 = os.path.join(out_dir, batch_name + "_QC.xlsx")
all_df = pd.concat([in_df, qc_df], axis=1)
all_df.to_excel(out_file1, sheet_name='All QC', index=False)

if in_ss:
    sss_df = pd.DataFrame({
        'Library'          : sample_list,
        'Total Bases'      : total_bases,
        'Percent Q30 (R1)' : pct_q30_1,
        'Percent Q30 (R2)' : pct_q30_2
    })

    hap_scores_pct = ['{:.2f}%'.format(x * 100) for x in hap_scores]

    fcs_df = pd.DataFrame({
        'Sample'                : sample_list,
        'contaminationestimate' : hap_scores_pct
    })

    out_file2 = os.path.join(out_dir, batch_name + "_Genoox.xlsx")
    writer    = pd.ExcelWriter(out_file2)
    #writer = pd.ExcelWriter(out_file1, engine='xlsxwriter')

    in_df.to_excel (writer, sheet_name='QC Metrics - qPCR', index=False)
    sss_df.to_excel(writer, sheet_name='Single Sample Stats', index=False, float_format="%.2f")
    fcs_df.to_excel(writer, sheet_name='Final Coverage Stats - TCP', index=False, float_format="%.2f")
    writer.save()
