#!/bin/bash

# 参数
DATASET_WGA_PATH=$1
DATASET_NAT_PATH=$2
SAMPLE_WGA_NAME=$3
SAMPLE_NAT_NAME=$4
REF_PATH=$5
PRE_OUTPUT_DIR=$6
DIFF_OUTPUT_DIR=$7
MERGE_OUTPUT_DIR=$8
MERGE_NAME=$9

# 数据预处理
nanodisco preprocess -p 20 -f $DATASET_WGA_PATH -s $SAMPLE_WGA_NAME -o $PRE_OUTPUT_DIR -r $REF_PATH

# 数据预处理 for NAT
nanodisco preprocess -p 20 -f $DATASET_NAT_PATH  -s $SAMPLE_NAT_NAME -o $PRE_OUTPUT_DIR -r $REF_PATH

# 计算特征
nanodisco difference -nj 10 -nc 1 -p 2 -i $PRE_OUTPUT_DIR -o $DIFF_OUTPUT_DIR -w $SAMPLE_WGA_NAME -n $SAMPLE_NAT_NAME -r $REF_PATH

# 合并结果
nanodisco merge -d $DIFF_OUTPUT_DIR -o $MERGE_OUTPUT_DIR -b $MERGE_NAME


# 安装 bam-readcount
# 请确保 bam-readcount 已安装并在 PATH 中

# 计算 basecall 特征
bam-readcount -f /home/chenzh/Project1_1/nd_example2/home/nanodisco/reference/BA/BA_sequence.fasta /home/chenzh/Project1_1/nd_example2/home/nanodisco/analysis/BA/preprocessed_subset/BA_NAT.fasta > /home/chenzh/readcount/BA_NAT_fq.tsv
