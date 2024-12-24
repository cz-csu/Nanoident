#!/bin/bash

# 参数设置默认值
Bactrial_PATH=${1:-"NO"}
DIFF_PATH=${2:-"analysis/NO_subset_difference.RDS"}
OUTPUT_PATH=${3:-"analysis/NO/NO_motifs_my_train"}
Motif=${4:-"analysis/NO/NO_motifs_my_train"}
REF_PATH=${5:-"reference/NO/NO_sequence.fasta"}
nanodisco characterize -p 4 -b $Bactrial_PATH -d $DIFF_PATH -o $OUTPUT_PATH -m $Motif -t nn -r $REF_PATH
#Change result path of function classify.Detected.Motifs.Get_test_data.basic_group in code/analysis_functions.R to the path you want.
