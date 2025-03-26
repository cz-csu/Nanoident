
#!/bin/bash

# 参数
Bactrial_name=$1
DIFF_PATH=$2
MOTIF_ENRICH_OUTPUT_PATH=$3
REF_PATH=$4
# Motif 富集模块
nanodisco motif -p 20 -b $Bactrial_name -d $DIFF_PATH -o $MOTIF_ENRICH_OUTPUT_PATH -r $REF_PATH -a
