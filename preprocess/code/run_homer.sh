#!/bin/bash
fasta=$1
output_dir=$2
export PERL5LIB=$PERL5LIB:/home/nanodisco/homer/bin/
#singularity exec --env PATH=/home/chenzh/Project1_1/nd_example2/home/nanodisco/homer/.//bin/:$PATH nanodisco.sif /home/chenzh/Project1_1/nd_example2/home/nanodisco/homer/bin/findMotifs.pl $fasta fasta $output_dir -len 8,10,12
/home/nanodisco/homer/bin/findMotifs.pl $fasta fasta $output_dir -len 10,12,14