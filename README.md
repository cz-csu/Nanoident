# MotifNet
A multi-column convolutional neural network model for typing and fine mapping of various bacterial methylation motifs using Nanopore Sequencing.

## Getting Started

Replace the folders in the original nanodisco repository with this repository's code folder.

### Prerequisites

There are two environments to build, an R environment for data preprocessing and motif enrichment, and a python environment for MotifNet model training.

**R:**

follow the nanodisco (https://github.com/fanglab/nanodisco/tree/master) to build the virtual environment using Singularity (v3.2.1 and above).

**python:**

Use MotifNet/env/freeze.yml to generate the environment directly under conda.

```bash
conda env create -f freeze.yml
```
### Data Preprocess

```bash
#data preprocess for WGA
nanodisco preprocess -p 20 -f dataset/unzip_data/MinION_NO_WGA -s NO_WGA -o analysis/NO/preprocessed_subset -r reference/NO/NO_sequence.fasta
#data preprocess for NAT
nanodisco preprocess -p 20 -f dataset/unzip_data/MinION_NO_NAT -s NO_NAT -o analysis/NO/preprocessed_subset -r reference/NO/NO_sequence.fasta
#calculate the features
nanodisco difference -nj 10 -nc 1 -p 2 -i analysis/NO/preprocessed_subset -o analysis/NO/difference_subset -w NO_WGA -n NO_NAT -r reference/NO/NO_sequence.fasta
# merge results
nanodisco merge -d analysis/NO/difference_subset -o analysis -b HP_
#motif enrichment module
nanodisco motif -p 20 -b NO -d analysis/NO_subset_difference.RDS -o analysis/NO -r reference/NO/NO_sequence.fasta -a 

#For the basecall feature, we first need to install bam-readcount (https://github.com/genome/bam-readcount)
bam-readcount -f /home/chenzh/Project1_1/nd_example2/home/nanodisco/reference/BA/BA_sequence.fasta  /home/chenzh/Project1_1/nd_example2/home/nanodisco/analysis/BA/preprocessed_subset/BA_NAT.fasta > /home/chenzh/readcount/BA_NAT_fq.tsv

#Then we need to use MotifNet/preprocess/readcount/parse_readcount.py change the output from TSV to CSV format we can use, you may need to change the path in the file.

python parse_readcount.py
```

### MotifNet train and test

**Get train and valid data:**
```bash
Rscript code/train.R
```
In MotifNet/code/train.R using function prepare.meta.classification.data in MotifNet/code/analysis_functions.R to get the BA_classification_data.basic_group_22.rds of each bacterium. Then through the function evaluate.classifiers.performance.basic_group.ALL.time to obtain the data needed for training, you may need to change my results path in the function.

**Get test data:**

```bash
nanodisco characterize -p 4 -b NO -d analysis/NO_subset_difference.RDS -o analysis/NO/NO_motifs_my_train -m CAANNNNNNNCTGG,CCAGNNNNNNNTTG,CTCGAG,GCGGCCGC -t nn -r reference/NO/NO_sequence.fasta
```

Change result path of function classify.Detected.Motifs.Get_test_data.basic_group in code/analysis_functions.R to the path you want.

**Change data form:**

Since the data we are working with is in an R environment, we need to convert it to a pkl file that can be used in a python environment using data_from_R_to_python.ipynb.

**Then we can train our MotifNet model:**

```bash
python train.py
```
This file also contains a section for test results.

**If you want to do a loocv evaluation and see the results, you can run:**

```bash
python train_loocv_model.py
python train_loocv_test.py
```
### Pretrained model

you can just use the model we provide in best_model/MotifNet.pt, you can run postprocess.py achieve the best results with  module.