suppressMessages(library(optparse))
path_script <- "/home/nanodisco/code/"
source(paste0(path_script,"analysis_functions.R"))
suppressMessages(load.libraries.characterize())
#BA
BA_motifs <- c("GCWGC","GGATCC")
BA_explicit_motifs <- c("G5mCWGC", "GGAT4mCC")
BA_mod_types <- c("5mC", "4mC")
BA_mod_poss <- c(2,5) #数的在explicit_motifs的位置,distance 是 0-base,mod_pos是1-base(更具如下)
BA_mod_poss_predict <- c(3,3)
#motif_summary$motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
BA_col_motifs<- c("#99FEF", "#99EC66")
BA_base_name<-"BA"
BA_genome<-"reference/BA/BA_sequence.fasta"
BA_path_diff_data<-"analysis/BA_subset_difference.RDS"
#BF
BF_motifs <- c("GATC")
BF_explicit_motifs <- c("GAT5mC")
BF_mod_types <- c("5mC")
BF_mod_poss <- c(4)
BF_mod_poss_predict <- c(4)
BF_col_motifs<- c("#99AC66")
BF_base_name<-"BF"
BF_genome<-"reference/BF/BF_sequence.fasta"
BF_path_diff_data<-"analysis/BF_subset_difference.RDS"
#CP
CP_motifs <- c("CCGG", "CACNNNNNRTAAA","GATC","GGWCC","GTATAC","TTTAYNNNNNGTG","VGACAT")
CP_explicit_motifs <- c("5mCCGG", "C6mACNNNNNRTAAA","GAT5mC","GGW5mCC","GTAT6mAC","TTT6mAYNNNNNGTG","VGAC6mAT")
CP_mod_types <- c("5mC", "6mA","5mC","5mC","6mA","6mA","6mA")
CP_mod_poss <- c(1,2,4,4,5,4,5)
CP_mod_poss_predict <- c(1,2,4,5,5,5,6)
CP_col_motifs<- c("#99FFFF", "#99CC66", "#99AC66", "#99BC66", "#89CC66", "#79CC66", "#69CC66")
CP_base_name<-"CP"
CP_genome<-"reference/CP/CP_sequence.fasta"
CP_path_diff_data<-"analysis/CP_subset_difference.RDS"
#EC
EC_motifs <- c("AACNNNNNNGTGC","CCWGG","GATC","GCACNNNNNNGTT")
EC_explicit_motifs <- c("A6mACNNNNNNGTGC","C5mCWGG","G6mATC","GC6mACNNNNNNGTT")
EC_mod_types <- c("6mA","5mC","6mA","6mA")
EC_mod_poss <- c(2,2,2,3)
EC_mod_poss_predict <- c(2,3,2,3)
EC_col_motifs<- c("#91FAF","#91FBF","#91FCF","#91FDF")
EC_base_name<-"EC"
EC_genome<-"reference/EC/EC_sequence.fasta"
EC_path_diff_data<-"analysis/EC_subset_difference.RDS"
#HP
HP_motifs <- c("CCGG","ATTAAT","CATG","CRTANNNNNNNWC","CSAG","CTRYAG","CYANNNNNNTTC","GCGC","GAGG","GANNNNNNNTAYG","GAANNNNNNTRG","GAATTC","GGCC","GMRGA","GTAC","GTNNAC","TCTTC","TCGA","TCNNGA","TGCA")
HP_explicit_motifs <- c("4mCCGG","ATTA6mAT","C6mATG","CRT6mANNNNNNNWC","CS6mAG","CTRY6mAG","CY6mANNNNNNTTC","G5mCGC","G6mAGG","G6mANNNNNNNTAYG","GA6mANNNNNNTRG","GA6mATTC","GG5mCC","GMRG6mA","GT6mAC","GTNN6mAC","T4mCTTC","TCG6mA","TCNNG6mA","TGC6mA")
HP_mod_types <- c("4mC","6mA","6mA","6mA","6mA","6mA","6mA","5mC","6mA","6mA","6mA","6mA","5mC","6mA","6mA","6mA","4mC","6mA","6mA","6mA")
HP_mod_poss <- c(1,5,2,4,3,5,3,2,2,2,3,3,3,5,3,5,2,4,6,4)
HP_mod_poss_predict <- c(1,6,3,4,3,5,3,3,2,2,3,4,3,5,3,5,3,4,6,4)
HP_col_motifs<- c("#31FAF","#31FBF","#31FCF","#31FDF","#32FAF","#32FBF","#32FCF","#32FDF","#33FAF","#33FBF","#33FCF","#33FDF","#81EDF","#34FBF","#34FCF","#34FDF","#35FAF","#35FBF","#35FCF","#35FDF")
HP_base_name<-"HP"
HP_genome<-"reference/HP/HP_sequence.fasta"
HP_path_diff_data<-"analysis/HP_subset_difference.RDS"
#MH
MH_motifs <- c("CTNAG","AGCT","CCACGK","GATC","GCYYGAT","GTAC")
MH_explicit_motifs <- c("4mCTNAG","AG4mCT","CCA4mCGK","G6mATC","GCYYG6mAT","GTA4mC")
MH_mod_types <- c("4mC","4mC","4mC","6mA","6mA","4mC")
MH_mod_poss <- c(1,3,4,2,6,4)
MH_mod_poss_predict <- c(1,3,5,2,7,4)
MH_col_motifs<- c("#91AAF","#91BBF","#91CCF","#91FCF","#91EDF","#91DDA")
MH_base_name<-"MH"
MH_genome<-"reference/MH/MH_sequence.fasta"
MH_path_diff_data<-"analysis/MH_subset_difference.RDS"
#NG
NG_motifs <- c("CCGCGG","GCCGGC","GAGNNNNNTAC","GCANNNNNNNNTGC","GGCC","GGNNCC","GGTGA","GTANNNNNCTC","RGCGCY")
NG_explicit_motifs <- c("C5mCGCGG","G5mCCGGC","G6mAGNNNNNTAC","GC6mANNNNNNNNTGC","GG5mCC","GGNN5mCC","GGTG6mA","GT6mANNNNNCTC","RG5mCGCY")
NG_mod_types <- c("5mC","5mC","6mA","6mA","5mC","5mC","6mA","6mA","5mC")
NG_mod_poss <- c(2,2,2,3,3,5,5,3,3)
NG_mod_poss_predict <- c(3,2,3,3,3,5,5,3,3)
NG_col_motifs<- c("#81AAF","#81BBF","#81CCF","#81DDF","#81EDF","#81DDA","#71DDA","#61DDA","#51DDA")
NG_base_name<-"NG"
NG_genome<-"reference/NG/NG_sequence.fasta"
NG_path_diff_data<-"analysis/NG_subset_difference.RDS"

list_classifier<-c("neural.network.predict")
length_vector<- 12
selected_starts<- -9
nb_motifs<-1
balancing<-0
nb_threads<-46
BA_difference_data <- readRDS(BA_path_diff_data) 
BF_difference_data <- readRDS(BF_path_diff_data) 
CP_difference_data <- readRDS(CP_path_diff_data) 
EC_difference_data <- readRDS(EC_path_diff_data) 
HP_difference_data <- readRDS(HP_path_diff_data) 
MH_difference_data <- readRDS(MH_path_diff_data) 
NG_difference_data <- readRDS(NG_path_diff_data) 
BA_motif_summary <- data.frame(motif=BA_motifs,explicit_motif=BA_explicit_motifs,mod_pos=BA_mod_poss_predict,mod_type=BA_mod_types,col_motif=BA_col_motifs, stringsAsFactors=FALSE)
BF_motif_summary <- data.frame(motif=BF_motifs,explicit_motif=BF_explicit_motifs,mod_pos=BF_mod_poss_predict,mod_type=BF_mod_types,col_motif=BF_col_motifs, stringsAsFactors=FALSE)
CP_motif_summary <- data.frame(motif=CP_motifs,explicit_motif=CP_explicit_motifs,mod_pos=CP_mod_poss_predict,mod_type=CP_mod_types,col_motif=CP_col_motifs, stringsAsFactors=FALSE)
EC_motif_summary <- data.frame(motif=EC_motifs,explicit_motif=EC_explicit_motifs,mod_pos=EC_mod_poss_predict,mod_type=EC_mod_types,col_motif=EC_col_motifs, stringsAsFactors=FALSE)
HP_motif_summary <- data.frame(motif=HP_motifs,explicit_motif=HP_explicit_motifs,mod_pos=HP_mod_poss_predict,mod_type=HP_mod_types,col_motif=HP_col_motifs, stringsAsFactors=FALSE)
MH_motif_summary <- data.frame(motif=MH_motifs,explicit_motif=MH_explicit_motifs,mod_pos=MH_mod_poss_predict,mod_type=MH_mod_types,col_motif=MH_col_motifs, stringsAsFactors=FALSE)
NG_motif_summary <- data.frame(motif=NG_motifs,explicit_motif=NG_explicit_motifs,mod_pos=NG_mod_poss_predict,mod_type=NG_mod_types,col_motif=NG_col_motifs, stringsAsFactors=FALSE)
motif_summary <- rbind(BA_motif_summary,BF_motif_summary,CP_motif_summary,EC_motif_summary,HP_motif_summary,MH_motif_summary,NG_motif_summary)
motif_summary=motif_summary[!duplicated(motif_summary),]

BA_classification_data = prepare.meta.classification.data(BA_difference_data, BA_base_name, BA_motif_summary, BA_genome, 0, TRUE, iupac_nc, nb_threads)
BF_classification_data = prepare.meta.classification.data(BF_difference_data, BF_base_name, BF_motif_summary, BF_genome, 0, TRUE, iupac_nc, nb_threads)
CP_classification_data = prepare.meta.classification.data(CP_difference_data, CP_base_name, CP_motif_summary, CP_genome, 0, TRUE, iupac_nc, nb_threads)
EC_classification_data = prepare.meta.classification.data(EC_difference_data, EC_base_name, EC_motif_summary, EC_genome, 0, TRUE, iupac_nc, nb_threads)
HP_classification_data = prepare.meta.classification.data(HP_difference_data, HP_base_name, HP_motif_summary, HP_genome, 0, TRUE, iupac_nc, nb_threads)
MH_classification_data = prepare.meta.classification.data(MH_difference_data, MH_base_name, MH_motif_summary, MH_genome, 0, TRUE, iupac_nc, nb_threads)
NG_classification_data = prepare.meta.classification.data(NG_difference_data, NG_base_name, NG_motif_summary, NG_genome, 0, TRUE, iupac_nc, nb_threads)

motif_center_summary<-NA
library(Rtsne)
ecp=evaluate.classifiers.performance.LOOCV.draw(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, "annotation_mod",motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nb_threads,balancing,list_classifier)
