# MotifNet
A weighted multi-feature convolutional neural network model for fine mapping of various bacterial methylation motifs

## Getting Started

Replace the folders in the original nanodisco repository with this repository's code folder.

### Prerequisites

There are two environments to build, an R environment for data preprocessing and motif enrichment, and a python environment for MotifNet model training.

R:

follow the nanodisco (https://github.com/fanglab/nanodisco/tree/master) to build the virtual environment using Singularity (v3.2.1 and above).

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

#Then we need to use MotifNet/the preprocess readcount/parse_readcount.py the output by TSV into CSV format we can use,you may need to change the path in the file.

python parse_readcount.py
```

### MotifNet train and test

Get train and valid data:
```bash
Rscript code/train.R
```
In code/train.R using function code/analysis_functions.R/prepare.meta.classification.data get BA_classification_data.basic_group_22.rds of each bacterium.And then through the evaluate.classifiers.performance.basic_group.ALL.time to obtain the data needed for training, you may need to change my results path in function.

Get test data:

```bash
nanodisco characterize -p 4 -b NO -d analysis/NO_subset_difference.RDS -o analysis/NO/NO_motifs_my_train -m CAANNNNNNNCTGG,CCAGNNNNNNNTTG,CTCGAG,GCGGCCGC -t nn -r reference/NO/NO_sequence.fasta
```

Change result path of function code/analysis_functions.R/classify.Detected.Motifs.Get_test_data.basic_group to the path you want.

Change data form:

Since the data we are working with is in an R environment, we need to convert it to a pkl file that can be used in a python environment using data_from_R_to_python.ipynb.

Then we can train our MotifNet model:

```bash
python train.py
```
This file also contains a section for test results

If you want to do a loocv evaluation and see the results, you can run it

```bash
python train_loocv_model.py
python train_loocv_test.py
```
### Pretrained model

you can just use the model we provide in best_model/MotifNet.pt, you can run postprocess.py achieve the best results with  module.