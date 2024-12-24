#!/bin/bash

fasta=$1
output_dir=$2
/home/nanodisco/meme_4.11.4/build/bin/dreme -p $fasta  -oc $output_dir
