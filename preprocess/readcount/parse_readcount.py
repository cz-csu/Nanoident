#!/usr/bin/env python

import os
import argparse
import logging
import pandas as pd
import numpy as np
from scipy.stats import chisquare
# Per-base/indel data fields
# IMPORTANT: this relies on Python 3.6+ to maintain insertion order
# Each field is a key with value a function to convert to the
# appropriate data type
base_fields = {
    'base': str,
    'count': int,
    'avg_mapping_quality': float,
    'avg_basequality': float,
    'avg_se_mapping_quality': float,
    'num_plus_strand': int,
    'num_minus_strand': int,
    'avg_pos_as_fraction': float,
    'avg_num_mismatches_as_fraction': float,
    'avg_sum_mismatch_qualities': float,
    'num_q2_containing_reads': int,
    'avg_distance_to_q2_start_in_q2_reads': float,
    'avg_clipped_length': float,
    'avg_distance_to_effective_3p_end': float
}

# Create a simple dataframe
nat_df = {
    'position': [],
    'avg_basequality': [],
    #'avg_num_mismatches_as_fraction': [], #如果有甲基化，错配数量会更多一点
    #'avg_sum_mismatch_qualities': [], #如果有甲基化，错配质量会更大一点，猜测
    'base_num_fraction': [],
    'ATCGN': []
}
diff_df = {
    'position': [],
    'avg_basequality': [],
    'base_num_fraction': [],
    'g_num': [],
    'g_value': []
}
wga_df ={
    'position': [],
    'avg_basequality': [],
    #'avg_num_mismatches_as_fraction': [], #如果有甲基化，错配数量会更多一点
    #'avg_sum_mismatch_qualities': [], #如果有甲基化，错配质量会更大一点，猜测
    'base_num_fraction': [],
    'ATCGN': []
}
nat_output="/home/chenzh/readcount/HP_NAT_fq.tsv"
wga_output="/home/chenzh/readcount/HP_WGA_fq.tsv"
# Open the bam-readcount output file and read it line by line
# Note that the output has no header, so we consume every line
#"""
with open(nat_output) as nat_fh:
  for line in nat_fh:
    # Strip newline from end of line
    line = line.strip()
    # Fields are tab-separated, so split into a list on \t
    fields = line.split('\t')
    # The first four fields contain overall information about the position
    chrom = fields[0]  # Chromosome/reference
    position = int(fields[1])  # Position (1-based)
    nat_df['position'].append(position)
    reference_base = fields[2]  # Reference base
    depth = int(fields[3])  # Depth of coverage
    # The remaining fields are data for each base or indel
    # Iterate over each base/indel
    ref_num=0
    base_total_num=0
    ATCGN=[]
    for base_data_string in fields[4:]:
      # We will store per-base/indel data in a dict
      base_data = {}
      # Split the base/indel data on ':'
      base_values = base_data_string.split(':')
      # Iterate over each field of the base/indel data
      for i, base_field in enumerate(base_fields.keys()):
        # Store each field of data, converting to the appropriate
        # data type
        base_data[base_field] = base_fields[base_field](base_values[i])
      if base_data['base'] == reference_base:
        ref_num=base_data['count']
        nat_df['avg_basequality'].append(base_data['avg_basequality'])
        #nat_df['avg_num_mismatches_as_fraction'].append(base_data['avg_num_mismatches_as_fraction'])
        #nat_df['avg_sum_mismatch_qualities'].append(base_data['avg_sum_mismatch_qualities'])
      if base_data['base'] == "A" or base_data['base'] == "C" or base_data['base'] == "G" or base_data['base'] == "T" or base_data['base'] == "N":
        ATCGN.append(base_data['count'])
        base_total_num+=base_data['count']
    nat_df['ATCGN'].append((np.array(ATCGN)+0.5)/(base_total_num+2.5))
    if base_total_num==0:
      nat_df['base_num_fraction'].append(0)
    else:
      nat_df['base_num_fraction'].append(ref_num/base_total_num)

with open(wga_output) as wga_fh:
  for line in wga_fh:
    # Strip newline from end of line
    line = line.strip()
    # Fields are tab-separated, so split into a list on \t
    fields = line.split('\t')
    # The first four fields contain overall information about the position
    chrom = fields[0]  # Chromosome/reference
    position = int(fields[1])  # Position (1-based)
    wga_df['position'].append(position)
    reference_base = fields[2]  # Reference base
    depth = int(fields[3])  # Depth of coverage
    # The remaining fields are data for each base or indel
    # Iterate over each base/indel
    ref_num=0
    base_total_num=0
    ATCGN=[]
    for base_data_string in fields[4:]:
      # We will store per-base/indel data in a dict
      base_data = {}
      # Split the base/indel data on ':'
      base_values = base_data_string.split(':')
      # Iterate over each field of the base/indel data
      for i, base_field in enumerate(base_fields.keys()):
        # Store each field of data, converting to the appropriate
        # data type
        base_data[base_field] = base_fields[base_field](base_values[i])
      if base_data['base'] == reference_base:
        ref_num=base_data['count']
        wga_df['avg_basequality'].append(base_data['avg_basequality'])
      if base_data['base'] == "A" or base_data['base'] == "C" or base_data['base'] == "G" or base_data['base'] == "T" or base_data['base'] == "N":
        ATCGN.append(base_data['count'])
        base_total_num+=base_data['count']
    wga_df['ATCGN'].append((np.array(ATCGN)+0.5)/(base_total_num+2.5))
    if base_total_num==0:
      wga_df['base_num_fraction'].append(0)
    else:
      wga_df['base_num_fraction'].append(ref_num/base_total_num)
#nat_df=pd.DataFrame(nat_df)
#"""
row_num=len(nat_df['position'])
#"""
for i in range(row_num):
    diff_df['position'].append(nat_df['position'][i])
    diff_df['avg_basequality'].append(nat_df['avg_basequality'][i]-wga_df['avg_basequality'][i])
    diff_df['base_num_fraction'].append(nat_df['base_num_fraction'][i]-wga_df['base_num_fraction'][i])
    diff_df['g_num'].append(2 * np.sum(nat_df['ATCGN'][i] * np.log(nat_df['ATCGN'][i] / wga_df['ATCGN'][i])))
    diff_df['g_value'].append(chisquare(nat_df['ATCGN'][i], f_exp= wga_df['ATCGN'][i]).pvalue)
#"""
diff_df=pd.DataFrame(diff_df)    
diff_df.to_csv('HP_diff2.csv')
#wga_df = pd.read_csv('BA_wga.csv')
#print(nat_df['position'][0],row_num)
#print(df_loaded.head(3))