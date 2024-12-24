
#!/bin/bash

# 参数
REF_PATH=$1
DATA_PATH=$2
OUTPUT_PATH=$3
#For the basecall feature, we first need to install bam-readcount (https://github.com/genome/bam-readcount)
bam-readcount -f $REF_PATH  $DATA_PATH > $OUTPUT_PATH

#Then we need to use MotifNet/preprocess/readcount/parse_readcount.py change the output from TSV to CSV format we can use, you may need to change the path in the file.

python readcount/parse_readcount.py