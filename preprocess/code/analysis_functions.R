
load.libraries.all <- function(){
	library(optparse)
	library(ggplot2)
	library(zoo)
	library(metap) #Sumlog
	library(Biostrings)
	library(grid)
	library(BSgenome)
	library(gtable)
	library(RColorBrewer)
	library(foreach)
	library(doMC)
	library(plyr)
	library(dplyr)
	library(reshape)
	library(tidyr)
	# library(rgl) # Not needed
	library(pROC)
	library(data.table)
	library(xml2)
	library(scales)
	library(ggdendro) # heatmap
	# library(fpc) # heatmap
	# library(dendextend) # heatmap
	library(ggseqlogo) # heatmap
	library(gtools) # Ordering classifier
	library(ggbeeswarm)
	library(egg)
	library(hues)
	library(stringr)
	library(Rtsne)
	library(VGAM)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(earth)
	library(kernlab)
	library(randomForest)
	library(e1071)
}

load.libraries.motif <- function(){
	library(optparse)
	library(stringr)
	library(doMC)
	library(plyr)
	library(dplyr)
	library(zoo)
	library(metap) #Sumlog
	library(GenomicRanges)
	library(Biostrings)
	library(BSgenome)
	library(xml2)
	library(tidyr)
	library(ggplot2)
	library(RColorBrewer)
	library(egg)
}

load.libraries.characterize <- function(){
	library(optparse)
	library(stringr)
	library(foreach)
	library(doMC)
	library(Biostrings)
	library(GenomicRanges)
	library(dplyr)
	library(zoo)
	library(ggplot2)
	library(grid)
	library(gtable)
	library(hues)
	library(tidyr)
	library(RColorBrewer)
	library(caret)
	library(nnet)
	library(randomForest)
}

load.libraries.merge <- function(){
	library(optparse)
	library(foreach)
	library(plyr)
}

load.libraries.refine <- function(){
	library(stringr)
	library(doMC)
	library(Biostrings)
	library(GenomicRanges)
	library(plyr)
	library(dplyr)
	library(ggplot2)
	library(RColorBrewer)
	library(egg)
}

print_message <- function(message){
	# Print message to terminal
	cat(paste0("[",Sys.time(),"] ",message,".\n"))
}

check.input.motif <- function(opt){
	# Check if current differences file exist.
	if(is.null(opt$path_diff_data)){
		cat("Parameter -d/--path_diff_data is missing. Please provide the path to a current differences file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_diff_data)){
			cat(paste0("Current differences file doesn't exist (",opt$path_diff_data,").\n"))
			cat("Please check -d/--path_diff_data parameter and/or run nanodisco difference.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check list of motif
	if(!is.null(opt$list_motif)){
		# Check list of motif formating
		if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$list_motif)){
			cat(paste0("Unknown character found in comma separated list of motifs (",opt$list_motif,").\n"))
			cat("Please check -m/--list_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
			quit(save="no", status=4)
		}else{
			opt$list_motif <- str_split(opt$list_motif, ",", simplify=TRUE)[1,]
		}
	}
	# Check if reference genome file exist.
	if(is.null(opt$genome)){
		cat("Parameter -r/--genome is missing. Please provide the path to a reference genome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$genome)){
			cat(paste0("Reference genome file doesn't exist (",opt$genome,").\n"))
			cat("Please check -r/--genome parameter. Path to reference genome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}
	# Check if a subset of contigs should be analyzed.
	if(is.null(opt$list_contig) && is.null(opt$contigs_file)){
		opt$list_contig <- NA # No contigs 
	}else{
		if(!is.null(opt$list_contig) && !is.null(opt$contigs_file)){
			cat("Parameter -c/--list_contig and --contigs_file cannot be both supplied. Please select one of the two options.\n")
			quit(save="no", status=4)
		}else{
			if(!is.null(opt$list_contig)){
				# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
			}else if(!is.null(opt$contigs_file)){
				if(!file.exists(opt$contigs_file)){
					cat(paste0("List of contigs file doesn't exist (",opt$contigs_file,").\n"))
					cat("Please check --contigs_file parameter. Path to file with of contigs (one per line).\n")
					quit(save="no", status=4)
				}else{
					contigs_from_file <- read.table(opt$contigs_file, stringsAsFactors=FALSE)
					# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
					opt$list_contig <- paste0(contigs_from_file$V1, collapse=",")
				}
			}
		}
	}

	return(opt)
}

check.input.characterize <- function(opt){
	# Check if current differences file exist.
	if(is.null(opt$path_diff_data)){
		cat("Parameter -d/--path_diff_data is missing. Please provide the path to a current differences file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_diff_data)){
			cat(paste0("Current differences file doesn't exist (",opt$path_diff_data,").\n"))
			cat("Please check -d/--path_diff_data parameter and/or run nanodisco difference.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check list of motif
	if(is.null(opt$list_motif)){
		cat("Parameter -m/--list_motif is missing. Please provide a comma separated list of motifs (e.g. GATC,CCWGG).\n")
		quit(save="no", status=3)
	}else{
		# Check list of motif formating
		if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$list_motif)){
			cat(paste0("Unknown character found in comma separated list of motifs (",opt$list_motif,").\n"))
			cat("Please check -m/--list_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
			quit(save="no", status=4)
		}
	}
	# Check list of model
	if(is.null(opt$type_model)){
		cat("Parameter -t/--type_model is missing. Please provide a comma separated list of model to apply (nn: neural network, rf: random forest, or knn: k-nearest neighbor; e.g. nn,rf).\n")
		quit(save="no", status=3)
	}else{
		# Check list of motif formating
		if(!all(str_split(opt$type_model, ",", simplify=TRUE)[1,] %in% c("nn","rf","knn"))){
			cat(paste0("Unknown type found in comma separated list of model (",opt$type_model,").\n"))
			cat("Please check -t/--type_model parameter. Only the following models are recognized: nn, rf, knn (e.g. nn,rf).\n")
			quit(save="no", status=4)
		}
	}
	# Check if reference genome file exist.
	if(is.null(opt$genome)){
		cat("Parameter -r/--genome is missing. Please provide the path to a reference genome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$genome)){
			cat(paste0("Reference genome file doesn't exist (",opt$genome,").\n"))
			cat("Please check -r/--genome parameter. Path to reference genome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}
	# Check if a subset of contigs should be analyzed.
	if(is.null(opt$list_contig) && is.null(opt$contigs_file)){
		opt$list_contig <- NA # No contigs 
	}else{
		if(!is.null(opt$list_contig) && !is.null(opt$contigs_file)){
			cat("Parameter -c/--list_contig and --contigs_file cannot be both supplied. Please select one of the two options.\n")
			quit(save="no", status=4)
		}else{
			if(!is.null(opt$list_contig)){
				# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
			}else if(!is.null(opt$contigs_file)){
				if(!file.exists(opt$contigs_file)){
					cat(paste0("List of contigs file doesn't exist (",opt$contigs_file,").\n"))
					cat("Please check --contigs_file parameter. Path to file with of contigs (one per line).\n")
					quit(save="no", status=4)
				}else{
					contigs_from_file <- read.table(opt$contigs_file, stringsAsFactors=FALSE)
					# TODO could check if contig exist in metagenome.fasta and/or in methylation profile file
					opt$list_contig <- paste0(contigs_from_file$V1, collapse=",")
				}
			}
		}
	}

	return(opt)
}

check.model.version <- function(diff_data){
	basecaller <- attr(diff_data, "basecaller")
	bc_version <- attr(diff_data, "version")
	if(!basecaller %in% c("Albacore","Guppy")){ # Base caller is not recognized
		cat(paste0("Base caller used for raw *.fast5 processing isn't recognized (",basecaller,").\n"))
		cat("Matching models are not available but motif characterization will proceed with the default model.\n")
		cat("You might not obtain optimal results. Additional information can be found in our GitHub repository.\n")
		basecaller <- "Albacore"
		bc_version <- "v2.3.4"
	}else{
		model_name_pattern <- paste0("final_model_",basecaller,"_",bc_version,"_.*.RDS")
		if(!any(grepl(model_name_pattern,list.files("/home/nanodisco/models/")))){ # Exact model not available
			if(basecaller=="Albacore"){ # Base caller is recognized
				if(bc_version=="v2.3.4"){ # Default version
					model_name_pattern <- paste0("final_model_",basecaller,"_v2.3.4_.*.RDS")
					if(!any(grepl(model_name_pattern,list.files("/home/nanodisco/models/")))){ # Default version, models were removed
						cat(paste0("Models for ",basecaller," version ",bc_version," exist but are not in your container.\n"))
						cat("Please reach out for help on GitHub.\n")
						# TODO create function to retrieve models e.g. get_model Albacore vx.x.x
						quit(save="no", status=4)
					}
				}else{ # Version not supported
					cat(paste0("Models for ",basecaller," version ",bc_version," is not available.\n"))
					cat("Motif characterization will still proceed with the default model but obtained results might not be optimal.\n")
					cat("Additional information can be found in our GitHub repository.\n")
					bc_version <- "v2.3.4"
				}
			}else if(basecaller=="Guppy"){ # Base caller is recognized but not supported
				cat(paste0("Models for ",basecaller," version ",bc_version," is not yet available but we are working on it.\n"))
				cat("Motif characterization will still proceed with the default model but obtained results might not be optimal.\n")
				cat("Additional information can be found in our GitHub repository.\n")
				basecaller <- "Albacore"
				bc_version <- "v2.3.4"
			}
		}
	}

	return(list(basecaller=basecaller, version=bc_version))
}

check.input.merge <- function(opt){
	# Check if current differences file exist.
	if(is.null(opt$path_diff_data)){
		cat("Parameter -d/--path_diff_data is missing. Please provide the path to a current differences file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_diff_data)){
			cat(paste0("Current differences file doesn't exist (",opt$path_diff_data,").\n"))
			cat("Please check -d/--path_diff_data parameter and/or run nanodisco difference.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}

	return(opt)
}

check.input.refine <- function(opt){
	# Check if current differences file exist.
	if(is.null(opt$path_diff_data)){
		cat("Parameter -d/--path_diff_data is missing. Please provide the path to a current differences file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$path_diff_data)){
			cat(paste0("Current differences file doesn't exist (",opt$path_diff_data,").\n"))
			cat("Please check -d/--path_diff_data parameter and/or run nanodisco difference.\n")
			quit(save="no", status=4)
		}
	}
	# Check if output directory exist and create if not.
	if(!dir.exists(opt$path_output)){
		dir.create(opt$path_output, recursive=TRUE)
	}
	# Check list of motif
	if(!is.null(opt$list_motif)){
		# Check list of motif formating
		if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$list_motif)){
			cat(paste0("Unknown character found in comma separated list of motifs (",opt$list_motif,").\n"))
			cat("Please check -m/--list_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
			quit(save="no", status=4)
		}else{
			opt$list_motif <- str_split(opt$list_motif, ",", simplify=TRUE)[1,]
		}
	}else{
		cat("Parameter -m/--list_motif is missing. Please provide a comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG).\n")
		quit(save="no", status=3)
	}
	# Check list of motif
	if(!is.null(opt$candidate_motif)){
		# Check list of motif formating
		if(opt$candidate_motif=="all"){
			opt$candidate_motif <- opt$list_motif
		}else if(!grepl(pattern="^[ACGTRYSWKMBDHVN,]+$", x=opt$candidate_motif)){
			cat(paste0("Unknown character found in comma separated list of motifs (",opt$candidate_motif,").\n"))
			cat("Please check -M/--candidate_motif parameter. Only the following characters are recognized: 'ACGTRYSWKMBDHVN,' (e.g. GATC,CCWGG).\n")
			quit(save="no", status=4)
		}else{
			opt$candidate_motif <- str_split(opt$candidate_motif, ",", simplify=TRUE)[1,]
		}
	}else{
		cat("Parameter -M/--candidate_motif is missing. Please provide a comma separated list of motifs following IUPAC nucleotide code or 'all' (e.g. GATC,CCWGG).\n")
		quit(save="no", status=3)
	}
	# Check if reference genome file exist.
	if(is.null(opt$genome)){
		cat("Parameter -r/--genome is missing. Please provide the path to a reference genome file.\n")
		quit(save="no", status=3)
	}else{
		if(!file.exists(opt$genome)){
			cat(paste0("Reference genome file doesn't exist (",opt$genome,").\n"))
			cat("Please check -r/--genome parameter. Path to reference genome (.fasta or .fa).\n")
			quit(save="no", status=4)
		}
	}

	return(opt)
}

#   _____                           _                                        _                
#  |  __ \                         | |                                      | |               
#  | |  \/ ___ _ __   ___ _ __ __ _| |  _ __   __ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___ 
#  | | __ / _ \ '_ \ / _ \ '__/ _` | | | '_ \ / _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
#  | |_\ \  __/ | | |  __/ | | (_| | | | |_) | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
#   \____/\___|_| |_|\___|_|  \__,_|_| | .__/ \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
#                                      | |                                                    
#                                      |_|                                                    

options(scipen=20) # 999

seq_params <- list(fwd_llen=5, fwd_rlen=16, rev_llen=10, rev_rlen=11)
iupac_nc <- data.frame(
	code=c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N"),
	pattern=c("A","C","G","T","[AG]","[CT]","[CG]","[AT]","[GT]","[AC]","[CGT]","[AGT]","[ACT]","[ACG]","[ACGT]"),
	choice=c("A","C","G","T","AG","CT","CG","AT","GT","AC","CGT","AGT","ACT","ACG","ACGT")
) # [^] do not rev.comp easily

#   _                     _ _                   _       _        
#  | |                   | (_)                 | |     | |       
#  | |     ___   __ _  __| |_ _ __   __ _    __| | __ _| |_ __ _ 
#  | |    / _ \ / _` |/ _` | | '_ \ / _` |  / _` |/ _` | __/ _` |
#  | |___| (_) | (_| | (_| | | | | | (_| | | (_| | (_| | || (_| |
#  \_____/\___/ \__,_|\__,_|_|_| |_|\__, |  \__,_|\__,_|\__\__,_|
#                                    __/ |                       
#                                   |___/                        

merge.corrected.chunks <- function(data_path){ #TODO check chunk are well build, maybe confusion start and or multi mapping
	list_corrected_files <- paste0(list.files(data_path, pattern="chunk.*.corrected_data.rds"))
	idx <- as.integer(do.call(rbind,strsplit(list_corrected_files,".", fixed=TRUE))[,2])
	corrected_files_info <- data.frame(file_path=paste0(data_path,"/",list_corrected_files), chunk_idx=idx, stringsAsFactors=FALSE)
	corrected_files_info <- corrected_files_info[order(corrected_files_info$chunk_idx),]

	corrected_data <- foreach(corrected_file=corrected_files_info$file_path[seq(60,64)], .combine=rbind) %do% { # TODO remove head
		tmp_corrected_data <- readRDS(corrected_file)

		return(tmp_corrected_data)
	}

	return(corrected_data)
}

merge.stat.chunks <- function(data_path, inSilico=FALSE){
	list_stat_files <- paste0(list.files(data_path, pattern="chunk.[0-9]+.difference.rds"))
	idx <- as.integer(do.call(rbind,strsplit(list_stat_files,".", fixed=TRUE))[,2])
	stat_files_info <- data.frame(file_path=paste0(data_path,"/",list_stat_files), chunk_idx=idx, stringsAsFactors=FALSE)
	stat_files_info <- stat_files_info[order(stat_files_info$chunk_idx),]
	stat_data <- foreach(stat_file=stat_files_info$file_path, .combine=rbind.fill, .multicombine=TRUE, .maxcombine=500) %do% { # TODO progress
		tmp_stat_data <- readRDS(stat_file)
		#Reorder column TODO fix scoring.position when not enough data; default res
		if(inSilico){
			if(any(grepl("contig",colnames(tmp_stat_data)))){
				column_order <- match(c("contig","position","dir","strand","data_type","N_val","mean_diff","t_test_pval","u_test_pval"),colnames(tmp_stat_data))
			}else{
				column_order <- match(c("position","dir","strand","data_type","N_val","mean_diff","t_test_pval","u_test_pval"),colnames(tmp_stat_data))
			}
			tmp_stat_data <- tmp_stat_data[,column_order] # TODO can remove NA here later
		}else{
			if(any(grepl("contig",colnames(tmp_stat_data)))){
				column_order <- match(c("contig","position","dir","strand","N_wga","N_nat","mean_diff","t_test_pval","u_test_pval"),colnames(tmp_stat_data))
			}else{
				column_order <- match(c("position","dir","strand","N_wga","N_nat","mean_diff","t_test_pval","u_test_pval"),colnames(tmp_stat_data))
			}
			tmp_stat_data <- tmp_stat_data[,column_order] # TODO can remove NA here later
		}

		return(tmp_stat_data)
	}
	if(!is.null(attr(readRDS(stat_files_info$file_path[1]), "basecaller"))){
		tmp_stat_data <- readRDS(stat_files_info$file_path[1])
		attr(stat_data, "basecaller") <- attr(tmp_stat_data, "basecaller")
		attr(stat_data, "version") <- attr(tmp_stat_data, "version")
	}
	return(stat_data)
}

merge.stat.chunks.more_attr <- function(data_path, inSilico=FALSE){
	list_stat_files <- paste0(list.files(data_path, pattern="chunk.[0-9]+.difference.rds"))
	idx <- as.integer(do.call(rbind,strsplit(list_stat_files,".", fixed=TRUE))[,2])
	stat_files_info <- data.frame(file_path=paste0(data_path,"/",list_stat_files), chunk_idx=idx, stringsAsFactors=FALSE)
	stat_files_info <- stat_files_info[order(stat_files_info$chunk_idx),]
	stat_data <- foreach(stat_file=stat_files_info$file_path, .combine=rbind.fill, .multicombine=TRUE, .maxcombine=500) %do% { # TODO progress
		tmp_stat_data <- readRDS(stat_file)
		#Reorder column TODO fix scoring.position when not enough data; default res
		if(inSilico){
			if(any(grepl("contig",colnames(tmp_stat_data)))){ #,"std_diff","time_diff","s_level_diff"
				column_order <- match(c("contig","position","dir","strand","data_type","N_val","mean_diff","t_test_pval","u_test_pval","mean_nat","mean_wga"),colnames(tmp_stat_data))
			}else{
				column_order <- match(c("position","dir","strand","data_type","N_val","mean_diff","t_test_pval","u_test_pval","mean_nat","mean_wga"),colnames(tmp_stat_data))
			}
			tmp_stat_data <- tmp_stat_data[,column_order] # TODO can remove NA here later
		}else{
			if(any(grepl("contig",colnames(tmp_stat_data)))){
				column_order <- match(c("contig","position","dir","strand","N_wga","N_nat","mean_diff","t_test_pval","u_test_pval","mean_nat","mean_wga"),colnames(tmp_stat_data))
			}else{
				column_order <- match(c("position","dir","strand","N_wga","N_nat","mean_diff","t_test_pval","u_test_pval","mean_nat","mean_wga"),colnames(tmp_stat_data))
			}
			tmp_stat_data <- tmp_stat_data[,column_order] # TODO can remove NA here later
		}
		return(tmp_stat_data)
	}
	if(!is.null(attr(readRDS(stat_files_info$file_path[1]), "basecaller"))){
		tmp_stat_data <- readRDS(stat_files_info$file_path[1])
		attr(stat_data, "basecaller") <- attr(tmp_stat_data, "basecaller")
		attr(stat_data, "version") <- attr(tmp_stat_data, "version")
	}
	return(stat_data)
}

save.stat.chunks <- function(data_path, name_file, is_inSilico=FALSE){
	if(dir.exists(data_path)){
		stat_data <- merge.stat.chunks(data_path, is_inSilico)
		saveRDS(stat_data, file=paste0(dirname(data_path),"/",name_file))

		return(stat_data)
	}else{
		print("Directory do not exists.")

		return(NULL)
	}
}

save.merged.stat.chunks <- function(data_path, path_file){
	if(dir.exists(data_path)){
		stat_data <- merge.stat.chunks.more_attr(data_path)
		saveRDS(stat_data, file=path_file, compress="bzip2")

		return(stat_data)
	}else{
		print("Directory do not exists.")

		return(NULL)
	}
}


#####
## Missing values
#####
# Check for missing values which might create problem with rollingapply
check.missing.values <- function(stat_data){
	results <- foreach(contig_name=unique(stat_data$contig), .combine=rbind) %do% {
		fwd_subset <- subset(stat_data, contig==contig_name & dir=="fwd" & !is.na(mean_diff))
		fwd_expected_len <- max(fwd_subset$position) - (min(fwd_subset$position) - 1)
		fwd_measured_len <- nrow(fwd_subset)
		rev_subset <- subset(stat_data, contig==contig_name & dir=="rev" & !is.na(mean_diff))
		rev_expected_len <- max(rev_subset$position) - (min(rev_subset$position) - 1)
		rev_measured_len <- nrow(rev_subset)

		res <- data.frame(
			fwd_expected_len=fwd_expected_len, fwd_measured_len=fwd_measured_len,
			rev_expected_len=rev_expected_len, rev_measured_len=rev_measured_len)

		return(res)
	}

	missing_fwd <- sum(results$fwd_expected_len - results$fwd_measured_len)
	missing_rev <- sum(results$rev_expected_len - results$rev_measured_len)
	print(paste0("There is ",missing_fwd," missing values on fwd."))
	print(paste0("There is ",missing_rev," missing values on rev."))
	print(paste0("There is ",missing_fwd + missing_rev," missing values total."))
}

#   _____ _                   _                  
#  /  ___(_)                 | |                 
#  \ `--. _  __ _ _ __   __ _| |_ _   _ _ __ ___ 
#   `--. \ |/ _` | '_ \ / _` | __| | | | '__/ _ \
#  /\__/ / | (_| | | | | (_| | |_| |_| | | |  __/
#  \____/|_|\__, |_| |_|\__,_|\__|\__,_|_|  \___|
#            __/ |                               
#           |___/                                

define.ambiguous.genomic.ranges <- function(g_seq){
	gr_ambiguous_position <- foreach(direction=c("fwd","rev"), .combine=c) %do% {
		tmp_ambiguous_position <- vmatchPattern("N", g_seq, fixed=TRUE)
		tmp_ambiguous_position <- as(tmp_ambiguous_position, "GRanges")
		strand(tmp_ambiguous_position) <- ifelse(direction=="fwd","+","-")

		tmp_ambiguous_position <- reduce(tmp_ambiguous_position)

		return(tmp_ambiguous_position)
	}

	return(gr_ambiguous_position)
}

hide.ambiguous.matches <- function(gr_ambiguous_position, tmp_motif_matches, direction){
	gr_motif_matches <- as(tmp_motif_matches, "GRanges")
	strand(gr_motif_matches) <- ifelse(direction=="fwd","+","-")

	overlaps <- findOverlaps(gr_ambiguous_position, gr_motif_matches, type="any", select="all")
	gr_motif_matches <- gr_motif_matches[setdiff(seq_along(gr_motif_matches), overlaps@to)] # Select non overalpping motif matches

	return(gr_motif_matches)
}

find.isolated.motifs.sub <- function(motif_summary, idx_motif, direction, g_seq, gr_ambiguous_position){
	original_motif <- motif_summary$motif[idx_motif]
	mod_pos <- motif_summary$mod_pos[idx_motif]
	len_motif <- nchar(original_motif)

	if(direction=="rev"){
		motif <- reverseComplement(DNAString(original_motif)) # Double checked
		mod_pos <- (len_motif - mod_pos) + 1
	}else{
		motif <- original_motif
	}

	tmp_motif_matches <- vmatchPattern(motif, g_seq, fixed=FALSE)
	tmp_motif_matches <- hide.ambiguous.matches(gr_ambiguous_position, tmp_motif_matches, direction)

	motif_matches <- data.frame(contig_name=seqnames(tmp_motif_matches), contig_pos_motif=start(tmp_motif_matches) + (mod_pos - 1)) # Mark mod_pos
	if(nrow(motif_matches)>0){
		motif_matches$motif <- as.factor(as.character(original_motif)) # Add motif
	}

	return(motif_matches)
}

find.isolated.motifs <- function(genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, verbose=TRUE){
	g_seq <- readDNAStringSet(genome)
	names(g_seq) <- sapply(strsplit(names(g_seq), " "), `[`, 1) # Keep only first word like bwa index

	gr_ambiguous_position <- define.ambiguous.genomic.ranges(g_seq)

	pos_isolated_motifs <- foreach(direction=c("fwd","rev"), .combine=rbind) %do% { # Process both strand
		if(verbose){
			print(paste0("Processing ",direction," strand."))
			print(paste0("    Researching motifs from motif_summary."))
		}
		# List all modified bases and keep duplicate
		if(nbCPU>1){
			registerDoMC(nbCPU)
			motifs_matches <- foreach(idx_motif=seq(1,nrow(motif_summary)), .combine=rbind) %dopar% {
				motif_matches <- find.isolated.motifs.sub(motif_summary, idx_motif, direction, g_seq, gr_ambiguous_position)
				return(motif_matches)
			}
			registerDoSEQ()
		}else{
			motifs_matches <- foreach(idx_motif=seq(1,nrow(motif_summary)), .combine=rbind) %do% {
				motif_matches <- find.isolated.motifs.sub(motif_summary, idx_motif, direction, g_seq, gr_ambiguous_position)
				return(motif_matches)
			}
		}
		left_free <- left_signal + right_signal + error_margin  # same as right_free TODO switch strand if unbalanced value
		right_free <- left_signal + right_signal + error_margin # same as left_free
		if(left_signal!=-1 | right_signal!=-1 | error_margin!=-1){
			if(verbose){
				print(paste0("    Removing overlapping motifs"))
			}
			dir_pos_isolated_motifs <- motifs_matches %>%
				arrange(contig_name, contig_pos_motif) %>%
				mutate(dist_next_contig_motif=lead(contig_pos_motif,1) - contig_pos_motif, dist_prev_contig_motif=lag(dist_next_contig_motif,1)) %>%
				mutate(is_contig_end=ifelse(lead(contig_name,1)!=contig_name | is.na(dist_next_contig_motif), TRUE, FALSE)) %>% # Find contigs ends
				mutate(dist_next_contig_motif=ifelse(is_contig_end, nchar(g_seq)[match(contig_name, names(g_seq))]- contig_pos_motif, dist_next_contig_motif)) %>%
				mutate(is_contig_start=ifelse(lag(contig_name,1)!=contig_name | is.na(dist_prev_contig_motif), TRUE, FALSE)) %>% # Find contigs starts
				mutate(dist_prev_contig_motif=ifelse(is_contig_start, contig_pos_motif - 1, dist_prev_contig_motif)) %>%
				dplyr::select(-c(is_contig_start,is_contig_end)) %>%
				filter(dist_next_contig_motif>right_free & dist_prev_contig_motif>left_free) %>%
				mutate(dir=as.factor(direction)) %>%
				dplyr::select(-c(dist_next_contig_motif,dist_prev_contig_motif))
		}else{
			dir_pos_isolated_motifs <- motifs_matches %>%
				mutate(dir=as.factor(direction))
		}

		return(dir_pos_isolated_motifs)
	}

	return(pos_isolated_motifs)
}
extract.motifs.signature.old <- function(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU){
	gr_methylation <- GRanges(
		seqnames=methylation_signal$contig,
		ranges=IRanges(methylation_signal$position, methylation_signal$position),
		strand=ifelse(methylation_signal$dir=="fwd","+","-")
	)

	print_message("  Search motif(s) occurrences")
	motifs <- find.isolated.motifs(genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE)
	expected_signal <- motifs %>%
		mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
		mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin)
	gr_signal_motifs <- GRanges(
		seqnames=expected_signal$contig_name,
		ranges=IRanges(expected_signal$left_side, expected_signal$right_side),
		strand=as.factor(ifelse(expected_signal$dir=="fwd","+","-"))
	)

	if(filter_iso){
		overlapping_motifs <- overlapsAny(gr_signal_motifs, type="any", drop.self=TRUE)
		gr_signal_motifs <- gr_signal_motifs[!overlapping_motifs]
		expected_signal <- expected_signal[!overlapping_motifs,]
	}

	overlaps_motifs_methylation <- findOverlaps(gr_signal_motifs, gr_methylation, type="any", select="all")
	methylation_at_motifs <- data.frame(
		contig=expected_signal$contig_name[overlaps_motifs_methylation@from],
		pos_motif=expected_signal$contig_pos_motif[overlaps_motifs_methylation@from],
		motif=expected_signal$motif[overlaps_motifs_methylation@from],
		pos_signal=methylation_signal$position[overlaps_motifs_methylation@to],
		dir=methylation_signal$dir[overlaps_motifs_methylation@to],
		strand=methylation_signal$strand[overlaps_motifs_methylation@to],
		N_wga=methylation_signal$N_wga[overlaps_motifs_methylation@to],
		N_nat=methylation_signal$N_nat[overlaps_motifs_methylation@to],
		mean_diff=methylation_signal$mean_diff[overlaps_motifs_methylation@to]
	)

	methylation_at_motifs <- methylation_at_motifs %>% 
		mutate(distance=ifelse(dir=="fwd",pos_signal-pos_motif,(-(pos_signal-pos_motif)) - 7)) # Relative distance to mod_pos with strand correction

	filtered_methylation_to_motifs_distribution <- methylation_at_motifs %>%
		filter(N_wga>=min_cov & N_nat>=min_cov)

	return(filtered_methylation_to_motifs_distribution)
}
extract.motifs.signature <- function(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU){
	gr_methylation <- GRanges(
		seqnames=methylation_signal$contig,
		ranges=IRanges(methylation_signal$position, methylation_signal$position),
		strand=ifelse(methylation_signal$dir=="fwd","+","-")
	)

	print_message("  Search motif(s) occurrences")
	motifs <- find.isolated.motifs(genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE)
	print("cccccccccccccc")
	#expected_signal <- motifs %>%
	#	mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
	#	mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin)

	expected_signal <- motifs %>%
		#mutate(left_side=contig_pos_motif-10) %>%
		#mutate(right_side=contig_pos_motif+11)
		mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
		mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin)
		#mutate(left_side=ifelse(motifs$dir=="fwd",contig_pos_motif+expected_signal_left-signal_margin,contig_pos_motif-(expected_signal_right+signal_margin))) %>%
		#mutate(right_side=ifelse(motifs$dir=="fwd",contig_pos_motif+expected_signal_right+signal_margin,contig_pos_motif-(expected_signal_left-signal_margin)))
		#mutate(left_side=ifelse(motifs$dir=="fwd",contig_pos_motif+expected_signal_left-signal_margin-5,contig_pos_motif-(expected_signal_right+signal_margin+5))) %>%
		#mutate(right_side=ifelse(motifs$dir=="fwd",contig_pos_motif+expected_signal_right+signal_margin+5,contig_pos_motif-(expected_signal_left-signal_margin-5)))
	gr_signal_motifs <- GRanges(
		seqnames=expected_signal$contig_name,
		ranges=IRanges(expected_signal$left_side, expected_signal$right_side),
		strand=as.factor(ifelse(expected_signal$dir=="fwd","+","-"))
	)
	#filter_iso=FALSE
	if(filter_iso){
		overlapping_motifs <- overlapsAny(gr_signal_motifs, type="any", drop.self=TRUE)
		gr_signal_motifs <- gr_signal_motifs[!overlapping_motifs]
		expected_signal <- expected_signal[!overlapping_motifs,]
	}

	overlaps_motifs_methylation <- findOverlaps(gr_signal_motifs, gr_methylation, type="any", select="all")
	methylation_at_motifs <- data.frame(
		contig=expected_signal$contig_name[overlaps_motifs_methylation@from],
		pos_motif=expected_signal$contig_pos_motif[overlaps_motifs_methylation@from],
		motif=expected_signal$motif[overlaps_motifs_methylation@from],
		pos_signal=methylation_signal$position[overlaps_motifs_methylation@to],
		dir=methylation_signal$dir[overlaps_motifs_methylation@to],
		strand=methylation_signal$strand[overlaps_motifs_methylation@to],
		N_wga=methylation_signal$N_wga[overlaps_motifs_methylation@to],
		N_nat=methylation_signal$N_nat[overlaps_motifs_methylation@to],
		mean_diff=methylation_signal$mean_diff[overlaps_motifs_methylation@to],
		std_diff=methylation_signal$std_diff[overlaps_motifs_methylation@to],
		t_test_pval=methylation_signal$t_test_pval[overlaps_motifs_methylation@to],
		u_test_pval=methylation_signal$u_test_pval[overlaps_motifs_methylation@to],
		time_diff=methylation_signal$time_diff[overlaps_motifs_methylation@to]
		#sd_wga=methylation_signal$sd_wga[overlaps_motifs_methylation@to]
	)
	print("aaaaaaaaaaaa")
	#methylation_at_motifs$nat=methylation_signal$nat[overlaps_motifs_methylation@to]
	#methylation_at_motifs$wga=methylation_signal$wga[overlaps_motifs_methylation@to]
	methylation_at_motifs <- methylation_at_motifs %>% 
		mutate(distance=ifelse(dir=="fwd",pos_signal-pos_motif,(-(pos_signal-pos_motif)-7))) # Relative distance to mod_pos with strand correction，本来要-7反链
	#本来要-7
	#filtered_methylation_to_motifs_distribution <- methylation_at_motifs %>%
	#	filter(N_wga>=min_cov & N_nat>=min_cov)
	return(list(methylation_at_motifs=methylation_at_motifs))
}

draw.motifs.signatures <- function(methylation_signal, motif_summary, genome, filter_iso, split_dir, base_name_output, min_cov, nbCPU, iupac_nc){
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 4

	motifs_signature <- extract.motifs.signature(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU)

	height_fct <- length(unique(methylation_signal$strand))
	if(split_dir){
		width_fct <- length(unique(methylation_signal$dir))
	}else{
		width_fct <- 1
	}

	if(filter_iso & split_dir){
		pdf(paste0("Signatures_isolated_motifs_per_strand_x_dataTypes_",base_name_output,"_v7.pdf"),width=5*width_fct,height=5*height_fct)
	}else if(filter_iso & !split_dir){
		pdf(paste0("Signatures_isolated_motifs_",base_name_output,"_v7.pdf"),width=5*width_fct,height=5*height_fct)
	}else if(!filter_iso & split_dir){
		pdf(paste0("Signatures_motifs_per_strand_x_dataTypes_",base_name_output,"_v7.pdf"),width=5*width_fct,height=5*height_fct)
	}else if(!filter_iso & !split_dir){
		pdf(paste0("Signatures_motifs_",base_name_output,"_v7.pdf"),width=5*width_fct,height=5*height_fct)
	}

	stifle <- foreach(idx_motif=seq(1,nrow(motif_summary))) %do% { # Cannot %dopar% because drawing
		current_motif <- motif_summary$motif[idx_motif]
		mod_pos <- motif_summary$mod_pos[idx_motif]
		mod_type <- motif_summary$mod_type[idx_motif]
		len_motif <- nchar(current_motif)
		clean_motif <- paste0(substr(current_motif,1,mod_pos-1),mod_type,substr(current_motif,mod_pos+1,len_motif))

		print(paste0("Processing ",current_motif," motif (",idx_motif,"/",nrow(motif_summary),")."))
		signature_data <- subset(motifs_signature, motif==current_motif & !is.na(mean_diff))

		if(nrow(signature_data)>15){ # If motif found and more than 15 data point left TODO arbitrary threshold
			gp <- ggplot(signature_data) +
				geom_jitter(aes(x=as.factor(distance), y=mean_diff), size=1, height=0) + # TODO jitter x only
				geom_violin(aes(x=as.factor(distance), y=mean_diff), alpha=0.6) + # Throw warning if < one motif per strand|dir
				labs(title=paste0("Signature of ",clean_motif," motif by type of data")) +
				labs(y="Difference in mean event current (pA)") +
				geom_hline(yintercept=0, col="red") +
				coord_cartesian(ylim=c(-7,15))

			if(filter_iso){
				gp <- gp + labs(x=paste0("Distance from modified base at isolated ",clean_motif," motifs"))
			}else{
				gp <- gp + labs(x=paste0("Distance from modified base at ",clean_motif," motifs"))
			}

			if(split_dir){ # If fwd and rev
				gp <- gp + facet_grid(strand~dir)
			}
			print(gp)
		}
	}
	dev.off()
}

#   _____ _                   _   _             _ _         
#  /  ___(_)                 | | (_)           (_) |        
#  \ `--. _  __ _ _ __   __ _| |  _ _ __    ___ _| |_ _   _ 
#   `--. \ |/ _` | '_ \ / _` | | | | '_ \  / __| | __| | | |
#  /\__/ / | (_| | | | | (_| | | | | | | | \__ \ | |_| |_| |
#  \____/|_|\__, |_| |_|\__,_|_| |_|_| |_| |___/_|\__|\__,_|
#            __/ |                                          
#           |___/                                           

keep.local.peaks <- function(peak_pos, values, half_win_size){
	left_bound <- peak_pos - half_win_size
	left_bound <- ifelse(left_bound>0, left_bound, 1) #Limit to bound
	right_bound <- peak_pos + half_win_size
	right_bound <- ifelse(right_bound<length(values), right_bound, length(values)) #Limit to bound

	if(all(values[c(left_bound:(peak_pos-1), (peak_pos+1):right_bound)] <= values[peak_pos])){
		return(peak_pos) 
	}else{
		return(numeric(0))
	}
}

find.peaks <- function (values, half_win_size=2){
	shapes <- diff(sign(diff(values, na.pad=FALSE)))
	peaks <- which(shapes < 0) + 1
	local_peaks <- sapply(peaks, keep.local.peaks, values=values, half_win_size=half_win_size)

	return(unlist(local_peaks))
}
mark.peaks.origin <- function(subsample_data, half_win_size=2){
	values <- -log10(subsample_data$p)
	peaks <- find.peaks(values, half_win_size)
	subsample_data$peak <- ifelse(seq_along(values) %in% peaks,TRUE,FALSE)

	return(subsample_data)
}
mark.peaks <- function(subsample_data, half_win_size=2){
	values <- -log10(subsample_data$p)
	win_values <- -log10(subsample_data$win_p)
	peaks <- find.peaks(values, half_win_size)
	win_peaks <- find.peaks(win_values, half_win_size)
	subsample_data$peak <- ifelse(seq_along(values) %in% peaks,TRUE,FALSE)
	subsample_data$win_peak <- ifelse(seq_along(win_values) %in% win_peaks,TRUE,FALSE)
	return(subsample_data)
}

mark.modified.bases <- function(methylation_signal, motif_summary, genome, nbCPU, iupac_nc){
	left_signal <- -1 # Remove nothing
	right_signal <- -1
	error_margin <- -1
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 4

	gr_methylation <- GRanges(
		seqnames=methylation_signal$contig,
		ranges=IRanges(methylation_signal$position, methylation_signal$position),
		strand=ifelse(methylation_signal$dir=="fwd","+","-")
	)

	print(paste0("Searching motifs."))
	motifs <- find.isolated.motifs(genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE)
	gr_modified_position <- GRanges(
		seqnames=motifs$contig_name,
		ranges=IRanges(motifs$contig_pos_motif, motifs$contig_pos_motif),
		strand=as.factor(ifelse(motifs$dir=="fwd","+","-"))
	)

	overlaps_motifs_methylation <- findOverlaps(gr_modified_position, gr_methylation, type="any", select="all")
	methylation_signal$col_motif <- 0
	methylation_signal$col_motif[overlaps_motifs_methylation@to] <- as.integer(motifs$motif[overlaps_motifs_methylation@from])

	overlapping_motifs <- overlapsAny(gr_modified_position, type="any", drop.self=TRUE)
	gr_nonunique_modified_position <- gr_modified_position[overlapping_motifs]

	overlaps_nonunique_motifs_methylation <- findOverlaps(gr_nonunique_modified_position, gr_methylation, type="any", select="all")
	methylation_signal$col_motif[overlaps_nonunique_motifs_methylation@to] <- -1

	return(methylation_signal)
}

mySumlog <- function(x){
	p_values <- x[!is.na(x)]
	nval <- length(p_values)
	if(nval==0){
		return(NA)
	}else if(nval==1){
		return(p_values) # Not perfect
	}else{
		return(sumlog(p_values)$p) # Not perfect because != n
	}
}

myLogitp <- function(x){
	p_values <- x[!is.na(x)]
	p_values[p_values==1] <- p_values[p_values==1] - 10^-7 # Don't handle 1

	nval <- length(p_values)
	if(nval==0){
		return(NA)
	}else if(nval==1){
		return(p_values) # Not perfect
	}else{
		return(logitp(p_values)$p) # Not perfect because != n
	}
}

mySump <- function(x){
	p_values <- x[!is.na(x)]
	nval <- length(p_values)
	if(nval==0){
		return(NA)
	}else if(nval==1){
		return(p_values) # Not perfect
	}else{
		return(sump(p_values)$p) # Not perfect because != n
	}
}

mySumz <- function(x){
	p_values <- x[!is.na(x)]
	p_values[p_values==1] <- p_values[p_values==1] - 10^-7 # Don't handle 1

	nval <- length(p_values)
	if(nval==0){
		return(NA)
	}else if(nval==1){
		return(p_values) # Not perfect
	}else{
		return(sumz(p_values)$p) # Not perfect because != n
	}
}
rollingFunction2 <- function(sample_data, win_size, stat_val, smooth_func){
	# Setting missing value ?
	if(win_size==0){ #Do nothing
		res <- subset(sample_data,select=c("position","ks_pval"))
		colnames(res) <- c("position","win_p")
	}else{
		genomic_position <- sample_data$position
		genomic_position_range <- seq(min(genomic_position),max(genomic_position))

		window_p_values <- sample_data[,match("ks_pval",names(sample_data))]
		window_p_values[window_p_values==0] <- 10^-300
		print(head(window_p_values))
		gapped_dataset <- data.frame(genomic_position=genomic_position,window_p_values=window_p_values)
		filled_dataset <- data.frame(genomic_position=genomic_position_range, window_p_values=with(gapped_dataset, window_p_values[match(genomic_position_range, genomic_position)]))
		win_stat_dataserie <- zoo(filled_dataset$window_p_values, order.by=filled_dataset$genomic_position)
		win_res <- rollapply(win_stat_dataserie, width=win_size, by=1, FUN=get(smooth_func), align="center") 
		print("pppppppp")
		win_res <- as.data.frame(win_res)
		print(head(win_res))
		win_res$position <- as.numeric(rownames(win_res))
		colnames(win_res) <- c("win_p", "position")
		#res=merge(res,win_res,by="position")
	}
	res$win_p[res$win_p==0] <- 10^-300
	# Add back original values
	#res <- merge(res, sample_data, by="position") # Merge needed because some position are lost

	return(win_res)
}
rollingFunction.origin <- function(sample_data, win_size, stat_val, smooth_func){
	# Setting missing value ?
	if(win_size==0){ #Do nothing
		res <- subset(sample_data,select=c("position",stat_val))
		colnames(res) <- c("position","p")
	}else{
		p_values <- sample_data[,match(stat_val,names(sample_data))]
		p_values[p_values==0] <- 10^-300 # Add pseudovalue for sumlog; avoid warning("Some studies omitted")
		genomic_position <- sample_data$position
		genomic_position_range <- seq(min(genomic_position),max(genomic_position))

		gapped_dataset <- data.frame(genomic_position=genomic_position, p_values=p_values)
		filled_dataset <- data.frame(genomic_position=genomic_position_range, p_values=with(gapped_dataset, p_values[match(genomic_position_range, genomic_position)]))

		# TODO could try na.approx() or na.spline()
		stat_dataserie <- zoo(filled_dataset$p_values, order.by=filled_dataset$genomic_position)
		res <- rollapply(stat_dataserie, width=win_size, by=1, FUN=get(smooth_func), align="center") # sumlog return 0 if out of range ~>9.881313e-324

		res <- as.data.frame(res)
		res$position <- as.numeric(rownames(res))
		colnames(res) <- c("p", "position")
	}
	res$p[res$p==0] <- 10^-300 # Cap sumlog p.values; avoid missing large peaks; TODO cap 200? TODO multiple test correction

	# Add back original values
	res <- merge(res, sample_data, by="position") # Merge needed because some position are lost

	return(res)
}

rollingFunction <- function(sample_data, win_size, stat_val, smooth_func){
	# Setting missing value ?
	if(win_size==0){ #Do nothing
		res <- subset(sample_data,select=c("position",stat_val,"ks_pval"))
		colnames(res) <- c("position","p","win_p")
	}else{
		p_values <- sample_data[,match(stat_val,names(sample_data))]
		p_values[p_values==0] <- 10^-300 # Add pseudovalue for sumlog; avoid warning("Some studies omitted")

		genomic_position <- sample_data$position
		genomic_position_range <- seq(min(genomic_position),max(genomic_position))

		window_p_values <- sample_data[,match("ks_pval",names(sample_data))]
		window_p_values[window_p_values==0] <- 10^-300
		
		gapped_dataset <- data.frame(genomic_position=genomic_position, p_values=p_values)
		filled_dataset <- data.frame(genomic_position=genomic_position_range, p_values=with(gapped_dataset, p_values[match(genomic_position_range, genomic_position)]))
		#sink("test/test.txt")
		#print(filled_dataset)
		#sink()
		# TODO could try na.approx() or na.spline()
		stat_dataserie <- zoo(filled_dataset$p_values, order.by=filled_dataset$genomic_position)
		res <- rollapply(stat_dataserie, width=win_size, by=1, FUN=get(smooth_func), align="center") # sumlog return 0 if out of range ~>9.881313e-324

		gapped_dataset <- data.frame(genomic_position=genomic_position,window_p_values=window_p_values)
		filled_dataset <- data.frame(genomic_position=genomic_position_range, window_p_values=with(gapped_dataset, window_p_values[match(genomic_position_range, genomic_position)]))
		win_stat_dataserie <- zoo(filled_dataset$window_p_values, order.by=filled_dataset$genomic_position)
		win_res <- rollapply(win_stat_dataserie, width=win_size, by=1, FUN=get(smooth_func), align="center") 

		res <- as.data.frame(res)
		win_res <- as.data.frame(win_res)

		res$position <- as.numeric(rownames(res))
		win_res$position <- as.numeric(rownames(win_res))
		colnames(res) <- c("p", "position")
		colnames(win_res) <- c("win_p", "position")
		res=merge(res,win_res,by="position")
	}
	res$p[res$p==0] <- 10^-300 # Cap sumlog p.values; avoid missing large peaks; TODO cap 200? TODO multiple test correction
	res$win_p[res$win_p==0] <- 10^-300
	# Add back original values
	res <- merge(res, sample_data, by="position") # Merge needed because some position are lost

	return(res)
}

draw.local.signal <- function(base_name, final_stat_data, motif_summary, param_local_signal, genome){
	start <- param_local_signal$start
	end <- param_local_signal$end
	smooth_win_size <- param_local_signal$smooth_win_size
	peak_win_size <- param_local_signal$peak_win_size
	stat_val <- param_local_signal$stat_val
	smooth_func <- param_local_signal$smooth_func

	#Prepare statistics
	sample_data <- subset(final_stat_data, position >= start & position <= end)
	if(length(unique(sample_data$position))<50){return(NULL)} # Tweak limits; Stop if less than 50 data points overall

	sample_data <- sample_data %>% # Remove dir strand combination with less than n data points
		group_by(dir, strand) %>%
		filter(length(position)>10) # Tweak limits

	sample_data$strand <- factor(sample_data$strand, levels=c("t", "c1", "c2")) # Order strand

	graph_data <- ddply(sample_data, .(dir,strand), rollingFunction, win_size=smooth_win_size, stat_val=stat_val, smooth_func=smooth_func)
	graph_data <- ddply(graph_data, .(dir,strand), mark.peaks, half_win_size=peak_win_size)
	
	graph_data$grp <- with(graph_data, interaction(dir, strand, sep="_"))
	list_groups <- list(grp1=c("fwd_t","rev_c1","rev_c2"), grp2=c("fwd_c1","fwd_c2","rev_t"))
	levels(graph_data$grp)[levels(graph_data$grp)%in%list_groups[[1]]] <- names(list_groups)[1]
	levels(graph_data$grp)[levels(graph_data$grp)%in%list_groups[[2]]] <- names(list_groups)[2]

	graph_data$p[-log10(graph_data$p)>200] <- 10^-200 # Cap p.values; needed in addition to rollingFunction capping

	min_breaks <- seq(min(graph_data$position),max(graph_data$position),2)
	limits_y <- c(0,max(15,-log10(graph_data$p)))

	gp_fwd <- ggplot(subset(graph_data, grp=="grp1")) +
		geom_line(aes(x=position, y=-log10(p), colour=strand)) +
		geom_point(aes(x=position, y=-log10(p), colour=strand, fill=peak), shape=21) +
		scale_x_continuous(minor_breaks=min_breaks, limits=c(start, end)) +
		scale_y_continuous(limits=limits_y) +
		scale_colour_manual(values=c("t"="red","c1"="green","c2"="blue")) +
		scale_fill_manual(name="Peaks", values=c("FALSE"="black","TRUE"="yellow")) +
		labs(y="Forward Strand, -log10(p.value)")
	gp_rev <- ggplot(subset(graph_data, grp=="grp2")) +
		geom_line(aes(x=position, y=-log10(p), colour=strand)) +
		geom_point(aes(x=position, y=-log10(p), colour=strand, fill=peak), shape=21) +
		scale_x_continuous(minor_breaks=min_breaks, limits=c(start, end)) +
		scale_colour_manual(name="Data\ntypes",values=c("t"="red","c1"="green","c2"="blue")) + 
		scale_fill_manual(name="Peaks", values=c("FALSE"="black","TRUE"="yellow")) +
		scale_y_reverse(limits=rev(limits_y)) +
		labs(x="Genomic position", y="Reverse Strand, -log10(p.value)")

	#Prepare sequences
	genome_seq <- readDNAStringSet(genome)
	gr <- GRanges(	seqnames=Rle(names(genome_seq)),
					ranges=IRanges(start=start, end=end),
					strand=strand(c("+"))
	)
	sequence <- getSeq(genome_seq, gr)
	comp_sequence <- complement(sequence)
	sequence <- as.character(sequence)
	comp_sequence <- as.character(comp_sequence)
	len_seq <- end-start+1

	label_fwd <- data.frame(seq(start,end),strsplit(sequence,split="")[[1]],rep("fwd",len_seq))
	colnames(label_fwd) <- c("position","base","dir")
	label_rev <- data.frame(seq(start,end),strsplit(comp_sequence,split="")[[1]],rep("rev",len_seq))
	colnames(label_rev) <- c("position","base","dir")

	labels <- rbind(label_fwd,label_rev)
	labels <- merge(sample_data,labels, by.x=c("position","dir"), by.y=c("position","dir")) # Can add strand if needed

 	values_color <- motif_summary$col_motif
 	names(values_color) <- seq_along(motif_summary$col_motif)
 	values_color <- c("0"="#000000", "-1"="#FF0000", values_color) # When not modified; Add multi-motif at -1; #32CD32

	gp_lab_fwd <- ggplot(subset(labels, dir=="fwd")) +
		geom_text(aes(label=base, x=as.numeric(position), y=0, colour=as.factor(col_motif))) +
		scale_colour_manual(values=values_color) +
		scale_x_continuous(minor_breaks=seq(min(graph_data$position),max(graph_data$position),2), limits=c(start, end)) +
		scale_y_continuous(breaks=NULL)
	gp_lab_rev <- ggplot(subset(labels, dir=="rev")) +
		geom_text(aes(label=base, x=as.numeric(position), y=0, colour=as.factor(col_motif))) +
		scale_colour_manual(values=values_color) +
		scale_x_continuous(minor_breaks=seq(min(graph_data$position),max(graph_data$position),2), limits=c(start, end)) +
		scale_y_continuous(breaks=NULL)

	#Assemble panels
	pdf(file=NULL) # Avoid opening empty window
	tmp_graph <- ggplotGrob(gp_rev)
	tmp_graph <- gtable_add_rows(tmp_graph, unit(1,"null"), 5)

	ylab_gp_fwd <- ggplotGrob(gp_fwd)$grobs[grepl("ylab-l", ggplotGrob(gp_fwd)$layout$name)]
	yaxis_gp_fwd <- ggplotGrob(gp_fwd)$grobs[grepl("axis-l", ggplotGrob(gp_fwd)$layout$name)]
	panel_gp_fwd <- ggplotGrob(gp_fwd)$grobs[grepl("panel", ggplotGrob(gp_fwd)$layout$name)]
	tmp_graph <- gtable_add_grob(tmp_graph, ylab_gp_fwd, name="fwd_ylab", clip="off", t=6, l=2)
	tmp_graph <- gtable_add_grob(tmp_graph, yaxis_gp_fwd, name="fwd_yaxis", t=6, l=3)
	tmp_graph <- gtable_add_grob(tmp_graph, panel_gp_fwd, name="fwd_panel", t=6, l=4)

	tmp_graph <- gtable_add_rows(tmp_graph, unit(0.05,"null"), 6)
	panel_gp_lab_fwd <- ggplotGrob(gp_lab_fwd)$grobs[grepl("panel", ggplotGrob(gp_lab_fwd)$layout$name)]
	tmp_graph <- gtable_add_grob(tmp_graph, panel_gp_lab_fwd, name="lab_fwd_panel", t=7, l=4)
	tmp_graph <- gtable_add_rows(tmp_graph, unit(0.05,"null"), 7)
	panel_gp_lab_rev <- ggplotGrob(gp_lab_rev)$grobs[grepl("panel", ggplotGrob(gp_lab_rev)$layout$name)]
	tmp_graph <- gtable_add_grob(tmp_graph, panel_gp_lab_rev, name="lab_rev_panel", t=8, l=4)

	tmp_graph$layout[grepl("guide-box", tmp_graph$layout$name),c(1,3)] <- c(7,7)
	dev.off() # Close ggplotGrob windows

	pdf(paste0("Local_signal_",base_name,".pdf"), width=18)
	grid.draw(tmp_graph)
	dev.off()
}

#  ___  ___      _   _  __       _      _            _   _             
#  |  \/  |     | | (_)/ _|     | |    | |          | | (_)            
#  | .  . | ___ | |_ _| |_    __| | ___| |_ ___  ___| |_ _  ___  _ __  
#  | |\/| |/ _ \| __| |  _|  / _` |/ _ \ __/ _ \/ __| __| |/ _ \| '_ \ 
#  | |  | | (_) | |_| | |   | (_| |  __/ ||  __/ (__| |_| | (_) | | | |
#  \_|  |_/\___/ \__|_|_|    \__,_|\___|\__\___|\___|\__|_|\___/|_| |_|
#                                                                      
#                                                                      

# Define optimum region extraction window
explore.winParams.motif.detection <- function(genome, final_stat_data, iupac_nc, motif_summary, winParams){
	# Process statistical signal
	smooth_win_size <- winParams$smooth_win_size
	stat_val <- winParams$stat_val
	peak_win_size <- winParams$peak_win_size
	smooth_func <- winParams$smooth_func
	nbCPU <- winParams$nbCPU

	registerDoMC(nbCPU)
	final_stat_data <- subset(final_stat_data, !is.na(mean_diff)) # if smoothing needed TODO handle circular genome
	final_stat_data <- ddply(final_stat_data, .(contig,dir,strand), rollingFunction, win_size=smooth_win_size, stat_val=stat_val, smooth_func=smooth_func, .parallel=TRUE) # Warning here
	final_stat_data <- subset(final_stat_data, !is.na(p))
	final_stat_data <- ddply(final_stat_data, .(contig,dir,strand), mark.peaks, half_win_size=peak_win_size, .parallel=TRUE)
	registerDoSEQ()

	# Mark modified bases
	marked_final_stat_data <- mark.modified.bases(final_stat_data, motif_summary, genome, nbCPU, iupac_nc)

	# Compute motifs range
	motif_summary <- motif_summary %>% # TODO write begining
		mutate(motif_len=nchar(motif)) %>%
		mutate(llen=mod_pos-1, rlen=motif_len-mod_pos) # Relative length to modified position

	# Generate possible combination
	seq_param <- expand.grid(list(size=winParams$sizes,offset=winParams$offsets))
	seq_param <- seq_param %>% 
		mutate(llen=offset, rlen=size-(offset+1)) %>%
		arrange(size) %>%
		filter(rlen>=0)

	# Find fully contained motifs for any ranges of sequences tested
	registerDoMC(winParams$nbCPU)
	dir_seq_isolated_motif <- foreach(direction=c("fwd","rev"), .combine=rbind) %do% {
		sample_data <- subset(marked_final_stat_data, col_motif!=0 & dir==direction)
		list_DNA_mod <- sample_data$col_motif
		list_DNA_mod_param <- data.frame(pos=sample_data$position)
		list_DNA_mod_param$llen <- 0
		list_DNA_mod_param$rlen <- 0
		if(direction=="fwd"){
			list_DNA_mod_param$llen[list_DNA_mod>0] <- motif_summary$llen[list_DNA_mod[list_DNA_mod>0]]
			list_DNA_mod_param$rlen[list_DNA_mod>0] <- motif_summary$rlen[list_DNA_mod[list_DNA_mod>0]]
		}else{
			list_DNA_mod_param$llen[list_DNA_mod>0] <- motif_summary$rlen[list_DNA_mod[list_DNA_mod>0]]
			list_DNA_mod_param$rlen[list_DNA_mod>0] <- motif_summary$llen[list_DNA_mod[list_DNA_mod>0]]
		}
		ranges_DNA_mod <- IRanges(start=list_DNA_mod_param$pos-list_DNA_mod_param$llen, end=list_DNA_mod_param$pos+list_DNA_mod_param$rlen)

		peaks_data <- subset(final_stat_data, peak==TRUE & -log10(p)>winParams$threshold & dir==direction)
		list_peak_pos <- peaks_data$position

		seq_isolated_motif <- foreach(idx_param=seq(1,nrow(seq_param)), .combine=rbind) %dopar%{
			llen <- seq_param$llen[idx_param]
			rlen <- seq_param$rlen[idx_param]
			ranges_peaks <- IRanges(start=list_peak_pos-llen, end=list_peak_pos+rlen)

			# Find peaks above threshold in motif range
			cmp_peaks_ann <- findOverlaps(ranges_DNA_mod,ranges_peaks,type="within") # FP peaks non-overlapped
			detected_DNA_mod <- cmp_peaks_ann@from
			true_peaks <- cmp_peaks_ann@to

			list_DNA_mod_pos <- ranges_DNA_mod[detected_DNA_mod]@start + list_DNA_mod_param$llen[detected_DNA_mod]
			tmp <- subset(marked_final_stat_data, select=c("position","col_motif"), position %in% list_DNA_mod_pos & dir==direction)
			tmp2 <- merge(data.frame(position_DNA_mod=list_DNA_mod_pos), tmp, by.x="position_DNA_mod",by.y="position")
			tmp <- subset(peaks_data, position %in% list_peak_pos[true_peaks])
			tmp3 <- merge(data.frame(position_peak=list_peak_pos[true_peaks]), tmp, by.x="position_peak",by.y="position")
			DNA_mod_stat_data <- cbind(tmp2, tmp3)

			DNA_mod_stat_data$motif[DNA_mod_stat_data$col_motif>0] <- motif_summary$motif[DNA_mod_stat_data$col_motif[DNA_mod_stat_data$col_motif>0]]
			DNA_mod_stat_data$motif[DNA_mod_stat_data$col_motif==-1] <- "MultipleMotifs"
			DNA_mod_stat_data <- DNA_mod_stat_data %>%
				group_by(position_peak) %>%
				summarize(N_motifs=length(motif), motifs=paste(motif, collapse=" "), .groups="drop_last")

			tmp_seq_isolated_motif <- as.data.frame(table(DNA_mod_stat_data$motifs[DNA_mod_stat_data$N_motifs==1]))
			colnames(tmp_seq_isolated_motif) <- c("motif","freq")
			tmp_seq_isolated_motif <- merge(motif_summary,tmp_seq_isolated_motif, all.x=TRUE)
			tmp_seq_isolated_motif$llen <- llen
			tmp_seq_isolated_motif$rlen <- rlen
			tmp_seq_isolated_motif$offset <- seq_param$offset[idx_param]
			tmp_seq_isolated_motif$size <- seq_param$size[idx_param]
			tmp_seq_isolated_motif$dir <- as.factor(direction)

			tmp_seq_isolated_motif <- tmp_seq_isolated_motif %>%
				dplyr::select(motif, freq, llen, rlen, offset, size, dir) %>%
				arrange(motif) %>%
				spread(motif,freq)

			return(tmp_seq_isolated_motif)
		}

		return(seq_isolated_motif)
	}
	registerDoSEQ()

	results <- melt(dir_seq_isolated_motif, id.vars=c("dir","size","offset","llen","rlen"))
	colnames(results) <- c("dir","size","offset","llen","rlen","motif","N_detected")
	results$dir <- as.factor(results$dir)

	results <- results %>%
		group_by(dir, motif) %>%
		mutate(max_detected=max(N_detected, na.rm=TRUE), perc_detected=(N_detected/max_detected)*100) # TODO change to expect max N

	return(results)
}

# TODO match motif color
draw.winParams.motif.detection <- function(results, winParams){
	gp <- ggplot(results) +
		geom_line(aes(x=offset, y=perc_detected, colour=motif), alpha=0.7) +
		facet_wrap(dir~size) +
		labs(title="Explore sequence extraction parameters for motif detection") +
		labs(subtitle=paste0("Threshold: ",winParams$threshold,"; Smoothing Win. Size: ",winParams$smooth_win_size,"; Stats: ",winParams$stat_val,"; Peaks Win. Size: ",winParams$peak_win_size)) +
		labs(x="Offset (relative position of peaks)") +
		labs(y="Percentage of isolated motif detected")
	pdf(paste0("Motif_detection_parameters_",winParams$base_name,"_v1.pdf"), width=35, height=15)
	print(gp)
	dev.off()
}

# Volcano plots
mark.modified.peaks <- function(peaks_data, annotated_data, seq_params){
	annotated_data <- subset(annotated_data, col_motif!=0)

	gr_peaks <- GRanges( # TODO maybe reduce peaks window
		seqnames=peaks_data$contig,
		ranges=IRanges(peaks_data$position - ifelse(peaks_data$dir=="fwd",seq_params$fwd_llen, seq_params$rev_llen), peaks_data$position + ifelse(peaks_data$dir=="fwd",seq_params$fwd_rlen,seq_params$rev_rlen)),
		strand=ifelse(peaks_data$dir=="fwd","+","-")
	)
	gr_annotation <- GRanges(
		seqnames=annotated_data$contig,
		ranges=IRanges(annotated_data$position, annotated_data$position),
		strand=ifelse(annotated_data$dir=="fwd","+","-")
	)
	overlap_peaks_annotation <- findOverlaps(gr_peaks, gr_annotation, type="any")
	peaks_data$status <- factor("Unmodified", levels=c("Unmodified","Modified"))
	peaks_data$status[overlap_peaks_annotation@from] <- "Modified"

	return(peaks_data)
}

scatter.color.density <- function(x, y, group=NA){
	data <- data.frame(x=x, y=y, group=group) %>%
		group_by(group) %>%
		mutate(cols=densCols(x, y, nbin=1024, colramp=colorRampPalette(rev(rainbow(10, end=4/6)))))

	return(data$cols)
}

peaks.stat.volcano <- function(peaks_data, final_stat_data, base_name, detection_params, min_cov, motif_summary, genome, iupac_nc){
	# Annotate modified peaks
	annotated_data <- mark.modified.bases(final_stat_data, motif_summary, genome, detection_params$nbCPU, iupac_nc)
	annotated_peaks_data <- mark.modified.peaks(subset(peaks_data, N_wga>=min_cov & N_nat>=min_cov), annotated_data, detection_params$seq_params)

	# Color points by density
	annotated_peaks_data$cols <- scatter.color.density(annotated_peaks_data$mean_diff, -log10(annotated_peaks_data$p), annotated_peaks_data$status)

	# TODO add density panel if needed?
	pdf(paste0("Volcano_plots_peaks_",base_name,"_v1.pdf"), width=15, height=7)
	tmp_gp <- ggplot(annotated_peaks_data) +
		geom_hline(yintercept=quantile(-log10(annotated_peaks_data$p), probs=seq(0.1,0.9,0.1)), size=0.1, lty=2) +
		geom_point(aes(x=mean_diff, y=-log10(p), col=cols), pch=46) +
		geom_hline(yintercept=0, size=0.4) +
		geom_vline(xintercept=0, col=1) +
		scale_color_identity() +
		facet_grid(.~status) +
		labs(title=paste0("Distribution of peaks p values")) +
		labs(subtitle=paste0("p values from ",detection_params$stat_val," peaks (w=",detection_params$peak_win_size,") smoothed over ",detection_params$smooth_win_size," position(s)")) +
		labs(x="Measured mean pA differences", y="-log10(p.value)") +
		theme_bw()
	print(tmp_gp)
	dev.off()
}

# Motif detection
filter.best.peaks <- function(peaks_data, n_top_hits){
	peaks_data <- peaks_data[order(-log10(peaks_data$p), decreasing=TRUE),]
	best_peaks <- head(peaks_data, n_top_hits)
	nb_selected_rows <- nrow(best_peaks)
	nb_peaks <- nrow(peaks_data)
	min_threshold <- round(min(-log10(best_peaks$p)),1)
	print(paste0("Sampling ",nb_selected_rows," (th: ",min_threshold,") peaks over ",nb_peaks," possible ones."))

	return(best_peaks)
}
filter.best.win.peaks <- function(peaks_data, n_top_hits){
	peaks_data <- peaks_data[order(-log10(peaks_data$win_p), decreasing=TRUE),]
	best_peaks <- head(peaks_data, n_top_hits)
	nb_selected_rows <- nrow(best_peaks)
	nb_peaks <- nrow(peaks_data)
	min_threshold <- round(min(-log10(best_peaks$win_p)),1)
	print(paste0("Sampling ",nb_selected_rows," (th: ",min_threshold,") peaks over ",nb_peaks," possible ones."))

	return(best_peaks)
}
sample.best.peaks <- function(peaks_data, sample_size, threshold){
	peaks_data <- subset(peaks_data, -log10(p)>threshold)
	if(nrow(peaks_data)<sample_size){
		warning("Sample size too large. No random sampling performed.")
		random_peaks <- peaks_data
	}else{
		print(paste0("Sampling ",sample_size," peaks over ",nrow(peaks_data)," possible ones."))
		random_peaks <- sample_n(peaks_data, sample_size, replace=FALSE)
	}

	return(random_peaks)
}
extract.peak.sequences.homer <- function(selected_peaks_data, genome, seq_params, output_name){
    # Load reference genome information
    g_seq <- readDNAStringSet(genome)
    list_contig_name <- str_split(names(g_seq)," ", simplify=TRUE)[,1] # Convert contig names to mapped contig names

    # Generate data.frame of peaks coordinates
    df_selected_peaks_data <- data.frame(
        contig=names(g_seq)[match(selected_peaks_data$contig, list_contig_name)],
        contig_length=width(g_seq)[match(selected_peaks_data$contig, list_contig_name)],
        position=selected_peaks_data$position,
        dir=selected_peaks_data$dir
    )
    # Define genomic range for fasta sequence extraction
    df_selected_peaks_data <- df_selected_peaks_data %>%
        mutate(strand=ifelse(dir=="fwd","+","-")) %>%
        mutate(left_bound_size=ifelse(dir=="fwd",seq_params$fwd_llen,seq_params$rev_llen)) %>% # TODO Remove if change signal alignment
        mutate(right_bound_size=ifelse(dir=="fwd",seq_params$fwd_rlen,seq_params$rev_rlen)) %>% # TODO Remove if change signal alignment
        mutate(gr_start=position-left_bound_size, gr_end=position+right_bound_size) %>%
        filter(gr_start>0) %>% # Filter out out-of-bound ranges
        filter(gr_end<=contig_length) # Filter out out-of-bound ranges

    # Generate genomic ranges of peaks
    gr_selected_peaks_data <- GRanges(
        seqnames=Rle(df_selected_peaks_data$contig),
        ranges=IRanges(start=df_selected_peaks_data$gr_start, end=df_selected_peaks_data$gr_end),
        strand=df_selected_peaks_data$strand # Is strand specific
    )

    # Create a data frame for BED format
    bed_df <- data.frame(
        chrom = seqnames(gr_selected_peaks_data),
        start = start(gr_selected_peaks_data) - 1,  # BED format is 0-based
        end = end(gr_selected_peaks_data),
        name = paste0(str_split(df_selected_peaks_data$contig, " ", simplify=TRUE)[,1],"_",df_selected_peaks_data$position,"_",df_selected_peaks_data$dir),
        strand = strand(gr_selected_peaks_data)
    )

    # Write BED file
    write.table(bed_df, file = output_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

extract.peak.sequences <- function(selected_peaks_data, genome, seq_params, output_name){
	# Load reference genome information
	g_seq <- readDNAStringSet(genome)
	list_contig_name <- str_split(names(g_seq)," ", simplify=TRUE)[,1] # Convert contig names to mapped contig names

	# Generate data.frame of peaks coordinates
	df_selected_peaks_data <- data.frame(
		contig=names(g_seq)[match(selected_peaks_data$contig, list_contig_name)],
		contig_length=width(g_seq)[match(selected_peaks_data$contig, list_contig_name)],
		position=selected_peaks_data$position,
		dir=selected_peaks_data$dir
	)
	# Define genomic range for fasta sequence extraction
	df_selected_peaks_data <- df_selected_peaks_data %>%
		mutate(strand=ifelse(dir=="fwd","+","-")) %>%
		mutate(left_bound_size=ifelse(dir=="fwd",seq_params$fwd_llen,seq_params$rev_llen)) %>% # TODO Remove if change signal alignment
		mutate(right_bound_size=ifelse(dir=="fwd",seq_params$fwd_rlen,seq_params$rev_rlen)) %>% # TODO Remove if change signal alignment
		mutate(gr_start=position-left_bound_size, gr_end=position+right_bound_size) %>%
		filter(gr_start>0) %>% # Filter out out-of-bound ranges
		filter(gr_end<=contig_length) # Filter out out-of-bound ranges

	# Generate genomic ranges of peaks
	gr_selected_peaks_data <- GRanges(
		seqnames=Rle(df_selected_peaks_data$contig),
		ranges=IRanges(start=df_selected_peaks_data$gr_start, end=df_selected_peaks_data$gr_end),
		strand=df_selected_peaks_data$strand # Is strand specific
	)

	sequences <- getSeq(g_seq, gr_selected_peaks_data)
	names(sequences) <- paste0(str_split(df_selected_peaks_data$contig, " ", simplify=TRUE)[,1],"_",df_selected_peaks_data$position,"_",df_selected_peaks_data$dir)
	writeXStringSet(sequences, output_name)
}

execute.meme <- function(input_fasta_name, meme_output, mod, nbmotif, min_motif_size, max_motif_size, nbCPU){
	cmd_meme <- paste(paste0(path_script,"run_meme.sh"), input_fasta_name, meme_output, mod, nbmotif, min_motif_size, max_motif_size, nbCPU)

	system(cmd_meme)
}

read.meme.output <- function(meme_output_xml){
	motif_data <- tryCatch({
			read_xml(meme_output_xml) # Throw warning, duplicate maxsites definition
		}, warning=function(w) {
			if(grepl("Redefinition of element maxsites", w)){
				motif_data <- suppressWarnings(read_xml(meme_output_xml))

				return(motif_data)
			}else{
				print(w)

				return(NA)
			}
		}, error=function(e) {
			print("Unexpected error.")
			print(e)

			return(NA)
		}
	)

	# General
	motifs_info <- motif_data %>% xml_find_all("//motif")
	motifs <- motifs_info %>% xml_attr("name")
	motifs_length <- motifs_info %>% xml_attr("width") %>% as.integer()
	motifs_nbhits <- motifs_info %>% xml_attr("sites") %>% as.integer()
	motifs_evalue <- motifs_info %>% xml_attr("e_value") %>% as.numeric()
	motifs_llr <- motifs_info %>% xml_attr("llr") %>% as.integer()
	motifs_info <- data.frame(motifs, motifs_length, motifs_nbhits, motifs_evalue, motifs_llr)

	# Probabilities
	motifs_prob_detail <- motif_data %>% xml_find_all("motifs//probabilities") # or motifs//scores
	motifs_probs_results <- foreach(idx=seq(1, length(motifs_prob_detail)), .combine=rbind) %do% { # Not needed but easier to read
		letter_id <- motifs_prob_detail[idx] %>% xml_find_all("descendant::value") %>% xml_attr("letter_id")
		prob <- motifs_prob_detail[idx] %>% xml_find_all("descendant::value") %>% xml_contents() %>% xml_text() %>% as.numeric() # as.integer() if scores

		return(data.frame(motif=motifs[idx], letter_id=letter_id, prob=prob, pos=rep(seq(1,motifs_length[idx]),each=4)))
	}
	motifs_probs_results <- motifs_probs_results[!duplicated(motifs_probs_results),] # Remove duplicated motifs, appended when one motif ~= all matches
	motifs_probs_results <- motifs_probs_results %>%
		spread(key=letter_id, value=prob) %>%
		arrange(motif,pos)

	# Scores
	motifs_scores_detail <- motif_data %>% xml_find_all("motifs//scores") # or motifs//probabilities
	motifs_scores_results <- foreach(idx=seq(1, length(motifs_scores_detail)), .combine=rbind) %do% { # Not needed but easier to read
		letter_id <- motifs_scores_detail[idx] %>% xml_find_all("descendant::value") %>% xml_attr("letter_id")
		prob <- motifs_scores_detail[idx] %>% xml_find_all("descendant::value") %>% xml_contents() %>% xml_text() %>% as.integer() # as.numeric() if scores

		return(data.frame(motif=motifs[idx], letter_id=letter_id, prob=prob, pos=rep(seq(1,motifs_length[idx]),each=4)))
	}
	motifs_scores_results <- motifs_scores_results[!duplicated(motifs_scores_results),] # Remove duplicated motifs, appended when one motif ~= all matches
	motifs_scores_results <- motifs_scores_results %>%
		spread(key=letter_id, value=prob) %>%
		arrange(motif,pos)

	motifs_results <- merge(motifs_probs_results, motifs_scores_results, suffixes=c(".probs",".scores"), by=c("motif","pos"))
	motifs_results <- motifs_results %>% arrange(motif,pos)
	motifs_results$motif <- factor(motifs_results$motif, levels=levels(motifs_info$motifs))

	return(list(general=motifs_info, detail=motifs_results))
}

refine.meme.motif <- function(sig_motif_detail, meme_results, iupac_nc){
	DNA_base <- 4
	min_pos_h <- 0.5
	min_base_h <- 0.13
	probs_cols <- c("A.probs","C.probs","G.probs","T.probs")
	#print(meme_results)
	#print("kkkkkkkkkkkkkkkkkkkkk")
	#print(sig_motif_detail)
	#print("kkkkkkkkkkkkkkkkkkkkk")
	#sink("/home/nanodisco/train/refine.meme.motif.txt")
	nb_sequences <- subset(meme_results$general, motifs==unique(sig_motif_detail$motif))$motifs_nbhits
	correction <- (1/log(2)) * ((DNA_base - 1)/(2 * nb_sequences))

	rel_frequencies <- sig_motif_detail[,c("pos", probs_cols)] %>%
		gather(base, rel_freq, all_of(probs_cols)) %>%
		mutate(base=gsub(".probs","",base))

	uncertainties <- rel_frequencies %>%
		mutate(tmp_uncertainty=rel_freq * log2(rel_freq)) %>%
		group_by(pos) %>%
		summarize(uncertainty=-sum(tmp_uncertainty, na.rm=TRUE), .groups="drop_last")

	information_contents <- uncertainties %>%
		mutate(information_content=log2(DNA_base) - (uncertainty + correction))

	sequence_logo <- merge(rel_frequencies, information_contents) %>%
		arrange(pos, base) %>%
		mutate(base_h=rel_freq * information_content) %>%
		group_by(pos) %>%
		mutate(pos_h=sum(base_h))

	#print(sequence_logo)
	#sink()
	tmp_motif <- sequence_logo %>%
		filter(pos_h>min_pos_h) %>%
		filter(base_h>min_base_h) %>%
		group_by(pos) %>%
		summarize(bases=paste0(base, collapse=""), .groups="drop_last") %>%
		mutate(code=as.character(iupac_nc$code[match(bases,iupac_nc$choice)]))

	list_pos <- seq(min(tmp_motif$pos), max(tmp_motif$pos))
	if(nrow(tmp_motif)<length(list_pos)){ # Missing position (N)
		tmp_motif <- rbind(tmp_motif, data.frame(pos=list_pos[!list_pos %in% tmp_motif$pos], bases="ACGT", code="N"))
	}
	#print(tmp_motif)
	tmp_motif <- tmp_motif %>%
		arrange(pos) %>%
		summarize(motif=paste0(code, collapse=""), .groups="drop_last")

	motif <- tmp_motif$motif

	return(motif)
}

process.meme.output <- function(meme_results, nb_peaks){
	min_motif_length <- 4
	min_motif_evalue <- 10^-30
	# N_motif_hits <- 0.05*nb_peaks
	N_motif_hits <- 0
	# min_motif_llr <- 900 # TODO maybe depend on motif length*nb_seqs, so filter later
	min_motif_llr <- 0

	sig_motifs <- subset(meme_results$general, motifs_nbhits>=N_motif_hits & motifs_evalue<=min_motif_evalue & motifs_llr>=min_motif_llr)

	potential_motifs <- foreach(sig_motif=sig_motifs$motifs, .combine=c) %do% {
		sig_motif_detail <- subset(meme_results$detail, motif==sig_motif)

		motif <- refine.meme.motif(sig_motif_detail, meme_results, iupac_nc)

		return(motif)
	}
	potential_motifs <- potential_motifs[nchar(potential_motifs)>=min_motif_length]

	return(potential_motifs)
}
process.statistical.signal.origin <- function(final_stat_data, param_signal_processing){
	smooth_win_size <- param_signal_processing$smooth_win_size
	stat_val <- param_signal_processing$stat_val
	peak_win_size <- param_signal_processing$peak_win_size
	smooth_func <- param_signal_processing$smooth_func
	nbCPU <- param_signal_processing$nbCPU

	registerDoMC(nbCPU)
	final_stat_data <- subset(final_stat_data, !is.na(mean_diff)) # if smoothing needed TODO handle circular genome
	final_stat_data <- ddply(final_stat_data, .(contig,dir,strand), rollingFunction.origin, win_size=smooth_win_size, stat_val=stat_val, smooth_func=smooth_func, .parallel=TRUE) # Warning here
	final_stat_data <- subset(final_stat_data, !is.na(p))
	final_stat_data <- ddply(final_stat_data, .(contig,dir,strand), mark.peaks.origin, half_win_size=peak_win_size, .parallel=TRUE)
	registerDoSEQ()

	peaks_data <- final_stat_data[final_stat_data$peak==TRUE,]

	return(peaks_data)
}
process.statistical.signal <- function(final_stat_data, param_signal_processing){
	smooth_win_size <- param_signal_processing$smooth_win_size
	stat_val <- param_signal_processing$stat_val
	peak_win_size <- param_signal_processing$peak_win_size
	smooth_func <- param_signal_processing$smooth_func
	nbCPU <- param_signal_processing$nbCPU

	registerDoMC(nbCPU)
	final_stat_data <- subset(final_stat_data, !is.na(mean_diff)) # if smoothing needed TODO handle circular genome
	final_stat_data <- ddply(final_stat_data, .(contig,dir,strand), rollingFunction, win_size=smooth_win_size, stat_val=stat_val, smooth_func=smooth_func, .parallel=TRUE) # Warning here
	final_stat_data <- subset(final_stat_data, !is.na(p))
	final_stat_data <- subset(final_stat_data, !is.na(win_p))
	final_stat_data <- ddply(final_stat_data, .(contig,dir,strand), mark.peaks, half_win_size=peak_win_size, .parallel=TRUE)
	registerDoSEQ()
	peaks_data <- final_stat_data[final_stat_data$peak==TRUE,]
	win_peaks_data <- final_stat_data[final_stat_data$win_peak==TRUE,]
	#peaks_data <- rbind(peaks_data, win_peaks_data)
	print(nrow(peaks_data))
	print(nrow(win_peaks_data))
	return(list(peaks_data=peaks_data, win_peaks_data=win_peaks_data))
}

find.motifs.ranges <- function(genome, discovered_motifs, iupac_nc){
	g_seq <- readDNAStringSet(genome) # Not design for multiple contig
	names(g_seq) <- sapply(strsplit(names(g_seq), " "), `[`, 1) # Keep only first word like bwa index

	gr_ambiguous_position <- define.ambiguous.genomic.ranges(g_seq)

	genome_annotation <- foreach(direction=c("fwd","rev"), .combine=rbind) %do% { # Process both strand
		# List all modified bases and keep duplicate
		sub_genome_annotation <- foreach(idx_motif=seq(1,length(discovered_motifs)), .combine=rbind) %do% {
			original_motif <- discovered_motifs[idx_motif]
			len_motif <- nchar(original_motif)

			if(direction=="rev"){
				motif <- reverseComplement(DNAString(original_motif)) # Double checked
			}else{
				motif <- original_motif
			}

			tmp_motif_matches <- vmatchPattern(motif, g_seq, fixed=FALSE)
			tmp_motif_matches <- hide.ambiguous.matches(gr_ambiguous_position, tmp_motif_matches, direction)

			if(length(tmp_motif_matches)==0){ # If no motif occurrences
				motif_matches <- data.frame(
					contig_name=NA,
					motif=as.character(original_motif),
					start=NA,
					end=NA,
					tmp_motif=as.character(motif)
				)
			}else{
				motif_matches <- data.frame(
					contig_name=seqnames(tmp_motif_matches),
					motif=as.character(original_motif),
					start=start(tmp_motif_matches),
					end=end(tmp_motif_matches),
					tmp_motif=as.character(motif)
				)
			}

			return(motif_matches)
		}
		sub_genome_annotation$dir <- direction

		return(sub_genome_annotation)
	}
	genome_annotation$dir <- as.factor(genome_annotation$dir)
	genome_annotation <- subset(genome_annotation, !is.na(start))

	return(genome_annotation)
}

remove.explained.peaks <- function(peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU){
	arbitrary_mod_pos <- 1
	left_signal <- -1 # Remove nothing
	right_signal <- -1
	error_margin <- -1

	known_motifs_summary <- data.frame(motif=discovered_motifs, mod_pos=arbitrary_mod_pos, stringsAsFactors=FALSE) # Don't know mod. position yet
	known_motifs <- find.isolated.motifs(genome, known_motifs_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE)
	known_motifs <- known_motifs %>%
		mutate(start=ifelse(dir=="fwd", contig_pos_motif, contig_pos_motif - nchar(as.character(motif)) + 1)) %>%
		mutate(end=ifelse(dir=="fwd", contig_pos_motif + nchar(as.character(motif)) - 1, contig_pos_motif))

	g_seq <- readDNAStringSet(genome)
	list_contig_name <- str_split(names(g_seq)," ", simplify=TRUE)[,1]
	ranges_known_motifs <- GRanges(
		seqnames=Rle(names(g_seq)[match(known_motifs$contig_name, list_contig_name)]),
		ranges=IRanges(start=known_motifs$start, end=known_motifs$end),
		strand=ifelse(known_motifs$dir=="fwd","+","-") # Is strand specific
	)

	list_peak_contig <- peaks_data$contig
	list_peaks_pos <- peaks_data$position
	list_left_bound_size <- ifelse(peaks_data$dir=="fwd",seq_params$fwd_llen,seq_params$rev_llen) # TODO Remove if change signal alignment
	list_right_bound_size <- ifelse(peaks_data$dir=="fwd",seq_params$fwd_rlen,seq_params$rev_rlen)
	list_strand <- ifelse(peaks_data$dir=="fwd","+","-")

	ranges_peaks_pos <- GRanges(
		seqnames=Rle(names(g_seq)[match(list_peak_contig, list_contig_name)]),
		ranges=IRanges(start=list_peaks_pos-list_left_bound_size, end=list_peaks_pos+list_right_bound_size),
		strand=list_strand # Is strand specific
	)

	idx <- unique(findOverlaps(ranges_known_motifs, ranges_peaks_pos, type="within")@to) # TODO within? #, type="within"
	list_explained_peaks_pos <- list_peaks_pos[idx]
	list_explained_peaks_strand <- ifelse(list_strand[idx]=="+","fwd","rev")
	peaks_data <- subset(peaks_data, ! (position %in% list_explained_peaks_pos & dir %in% list_explained_peaks_strand))

	return(peaks_data)
}

# Replaced by tag.mutated.motifs
define.relative.position <- function(mutated_motifs_range){
	position <- seq(mutated_motifs_range$start_signal, mutated_motifs_range$end_signal)
	distance <- seq(1, (mutated_motifs_range$end_signal + 1) - mutated_motifs_range$start_signal)
	
	if(mutated_motifs_range$dir=="rev"){
		distance <- rev(distance)
	}

	return(data.frame(position=position, dir=mutated_motifs_range$dir, motif=mutated_motifs_range$motif, distance=distance))
}

remove.homopolymer <- function(peaks_data, genome, len_homo, iupac_nc, seq_params, nbCPU){
	library(stringi)

	homopolymers <- stri_dup(c("A","C","G","T"),len_homo)
	peaks_data <- remove.explained.peaks(peaks_data, homopolymers, seq_params, genome, iupac_nc, nbCPU)

	return(peaks_data)
}

# Refine motifs
generate.mutated.motifs <- function(motif){
	len_motif <- nchar(motif)
	splited_motif <- strsplit(motif,"")[[1]]
	list_mutated_motif <- foreach(idx=seq(1,len_motif), .combine=rbind) %do% {
		tmp_motif <- splited_motif
		sublist_mutated_motif <- foreach(nucleotide=c("A","C","G","T"), .combine=rbind) %do% {
			tmp_motif[idx] <- nucleotide

			return(data.frame(mutated_motif=paste0(tmp_motif,collapse=""), mutation_type=nucleotide))
		}
		sublist_mutated_motif$mutated_motif <- as.character(sublist_mutated_motif$mutated_motif)
		sublist_mutated_motif$pos_mutation <- idx

		return(sublist_mutated_motif)
	}
	list_mutated_motif$original_motif <- motif

	return(list_mutated_motif)
}

remove.explained.ranges <- function(mutated_motifs_expected_signal, discovered_motifs, genome, iupac_nc, left_signal, right_signal, error_margin, nbCPU){
	arbitrary_mod_pos <- 1
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 0

	known_motifs_summary <- data.frame(motif=discovered_motifs, mod_pos=arbitrary_mod_pos, stringsAsFactors=FALSE) # Don't know mod. position yet
	known_motifs <- find.isolated.motifs(genome, known_motifs_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE)
	known_motifs_expected_signal <- known_motifs %>%
		mutate(motif_length=nchar(as.character(motif))) %>%
		mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
		mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin) %>%
		mutate(left_side=ifelse(dir=="fwd", left_side - (arbitrary_mod_pos - 1), left_side - (motif_length - arbitrary_mod_pos))) %>% # Add full motif window with strand correction
		mutate(right_side=ifelse(dir=="fwd", right_side + (motif_length - arbitrary_mod_pos), right_side + (arbitrary_mod_pos - 1)))

	gr_signal_known_motifs <- GRanges(
		seqnames=known_motifs_expected_signal$contig_name,
		ranges=IRanges(known_motifs_expected_signal$left_side, known_motifs_expected_signal$right_side),
		strand=ifelse(known_motifs_expected_signal$dir=="fwd","+","-")
	)
	gr_signal_mutated_motifs <- GRanges(
		seqnames=mutated_motifs_expected_signal$contig_name,
		ranges=IRanges(mutated_motifs_expected_signal$left_side, mutated_motifs_expected_signal$right_side),
		strand=as.factor(ifelse(mutated_motifs_expected_signal$dir=="fwd","+","-"))
	)
	# gr_unexplained_signal_motifs <- GenomicRanges::setdiff(gr_signal_mutated_motifs, gr_signal_known_motifs)

	idx <- unique(findOverlaps(gr_signal_mutated_motifs, gr_signal_known_motifs, type="any")@from) # Find partial overlaps
	mutated_motifs_expected_signal <- mutated_motifs_expected_signal[-idx,] # Remove partial overlaps

	return(mutated_motifs_expected_signal)
}

tag.mutated.motifs <- function(potential_motif, discovered_motifs, final_stat_data, genome, nbCPU, iupac_nc){
	options(error = traceback)
	arbitrary_mod_pos <- 1
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 0

	mutated_potential_motif <- generate.mutated.motifs(potential_motif)
	print("bbbbbbbbb")
	mutated_motifs_summary <- data.frame(motif=unique(mutated_potential_motif$mutated_motif), mod_pos=arbitrary_mod_pos, stringsAsFactors=FALSE) # Don't know mod. position yet
	mutated_motifs <- find.isolated.motifs(genome, mutated_motifs_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE)
	print("eeeeee")
	mutated_motifs_expected_signal <- mutated_motifs %>%
		mutate(motif_length=nchar(as.character(motif))) %>%
		mutate(left_side=contig_pos_motif+expected_signal_left-signal_margin) %>%
		mutate(right_side=contig_pos_motif+expected_signal_right+signal_margin) %>%
		mutate(left_side=ifelse(dir=="fwd", left_side - (arbitrary_mod_pos - 1), left_side - (motif_length - arbitrary_mod_pos))) %>% # Add full motif window with strand correction
		mutate(right_side=ifelse(dir=="fwd", right_side + (motif_length - arbitrary_mod_pos), right_side + (arbitrary_mod_pos - 1)))

	# TODO maybe remove overlapping motif ranges: not sure.
	# if(filter_iso){
	# 	overlapping_motifs <- overlapsAny(gr_signal_motifs, type="any", drop.self=TRUE)
	# 	gr_signal_motifs <- gr_signal_motifs[!overlapping_motifs]
	# 	mutated_motifs_expected_signal <- mutated_motifs_expected_signal[!overlapping_motifs,]
	# }

	# Remove ranges explained in this batch (peaks are filtered after each batch)
	if(!is.null(discovered_motifs)){
		mutated_motifs_expected_signal <- remove.explained.ranges(mutated_motifs_expected_signal, discovered_motifs, genome, iupac_nc, left_signal, right_signal, error_margin, nbCPU)
	}
	print("ccccccccc")
	if(nrow(mutated_motifs_expected_signal)!=0){ # If a mutated motif remain
		gr_signal_mutated_motifs <- GRanges(
			seqnames=mutated_motifs_expected_signal$contig_name,
			ranges=IRanges(mutated_motifs_expected_signal$left_side, mutated_motifs_expected_signal$right_side),
			strand=as.factor(ifelse(mutated_motifs_expected_signal$dir=="fwd","+","-"))
		)

		gr_methylation <- GRanges(
			seqnames=final_stat_data$contig,
			ranges=IRanges(final_stat_data$position, final_stat_data$position),
			strand=ifelse(final_stat_data$dir=="fwd","+","-")
		)

		overlaps_motifs_methylation <- findOverlaps(gr_signal_mutated_motifs, gr_methylation, type="any", select="all")
		methylation_at_motifs <- data.frame(
			contig=mutated_motifs_expected_signal$contig_name[overlaps_motifs_methylation@from],
			pos_motif=mutated_motifs_expected_signal$contig_pos_motif[overlaps_motifs_methylation@from],
			motif=mutated_motifs_expected_signal$motif[overlaps_motifs_methylation@from],
			dir=mutated_motifs_expected_signal$dir[overlaps_motifs_methylation@from],
			position=final_stat_data$position[overlaps_motifs_methylation@to],
			strand=final_stat_data$strand[overlaps_motifs_methylation@to],
			N_wga=final_stat_data$N_wga[overlaps_motifs_methylation@to],
			N_nat=final_stat_data$N_nat[overlaps_motifs_methylation@to],
			mean_diff=final_stat_data$mean_diff[overlaps_motifs_methylation@to]
		)
		print("dddddddddd")
		methylation_at_motifs <- methylation_at_motifs %>%
			filter(!is.na(mean_diff)) %>%
			mutate(distance=ifelse(dir=="fwd",position-pos_motif,(-(position-pos_motif)) - 7)) %>% # Relative distance to mod_pos with strand correction
			mutate(distance=distance+(arbitrary_mod_pos-expected_signal_left)) # Correct to start at 1 

		# Filter on coverage
		# filtered_methylation_to_motifs_distribution <- methylation_at_motifs %>%
		# 	filter(N_wga>=min_cov & N_nat>=min_cov)

		tagged_final_stat_data <- merge(methylation_at_motifs, mutated_potential_motif, by.x=c("motif"), by.y="mutated_motif") # Add mutated motifs informations

		return(tagged_final_stat_data)
	}else{
		print(paste0("Potential motif ",potential_motif," do not have enough occurrences to check."))

		return(NULL)
	}
}

# TODO add progression
# TODO avoid generating refine plots twice
score.mutated.motifs <- function(tagged_final_stat_data){
	motif_summary <- tagged_final_stat_data %>% # TODO try to improve summary function
		group_by(mutation_type, pos_mutation, distance) %>%
		summarize(score2=abs(mean(mean_diff, na.rm=TRUE)), .groups="drop_last") %>%
		group_by(mutation_type, pos_mutation) %>%
		summarize(score=sum(score2), .groups="drop_last")
	motif_summary$mutation_type <- ordered(motif_summary$mutation_type, levels=c("T","G","C","A")) # Same order as facet_grid

	return(motif_summary)
}

refine.motif <- function(potential_motif, discovered_motifs, final_stat_data, genome, nbCPU, seq_params, iupac_nc, path_output, graph_subtitle, automated=TRUE){
	print_message(paste0("Generating refining plot for ",potential_motif))
	len_motif <- nchar(potential_motif)

	modification_at_motifs <- tag.mutated.motifs(potential_motif, discovered_motifs, final_stat_data, genome, nbCPU, iupac_nc)
	print("aaaaaa")
	if(!is.null(modification_at_motifs)){
		if(automated){
			compound_name <- paste0(potential_motif)
		}else{
			compound_name <- paste0(potential_motif,"_",graph_subtitle)
		}

		mutated_motif_score <- score.mutated.motifs(modification_at_motifs)
		myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

		# Plot mutated motif signal
		xmin_value <- min(modification_at_motifs$distance)
		xmax_value <- max(modification_at_motifs$distance)
		ymin_value <- max(-7,min(modification_at_motifs$mean_diff))
		ymax_value <- min(7,max(modification_at_motifs$mean_diff))
		cleanup_label <- function(x){ sprintf("%.0f", x) }

		gp_detail <- ggplot(modification_at_motifs) +
			geom_rect(data=mutated_motif_score, aes(xmin=xmin_value-0.5, xmax=xmax_value+0.5, ymin=ymin_value, ymax=ymin_value + 2, fill=score)) +
			geom_jitter(aes(x=distance, y=mean_diff), pch=46, height=0) +
			geom_violin(aes(x=distance, y=mean_diff, group=distance), alpha=0.6) + # Throw warning if < one motif per strand|dir
			geom_hline(yintercept=0, col="red") +
			facet_grid(mutation_type~pos_mutation) +
			scale_fill_gradientn(colours=myPalette(100), guide=FALSE) +
			labs(title=paste0("Refinement plot for ",compound_name," motifs")) +
			labs(y="Mean current differences (pA)") +
			coord_cartesian(xlim=c(xmin_value-0.5, xmax_value+0.5), ylim=c(ymin_value, ymax_value), expand=FALSE) +
			scale_x_continuous(breaks=seq(xmin_value, xmax_value, 2)) +
			theme_bw()

		filter_iso <- FALSE # Not used now
		if(filter_iso){
			gp_detail <- gp_detail + labs(x=paste0("Relative position from isolated ",potential_motif," motif occurrences (-3:motif:+2)"))
		}else{
			gp_detail <- gp_detail + labs(x=paste0("Relative position from ",potential_motif," motif occurrences (-3:motif:+2)"))
		}

		# Plot mutated motif scores
		gp_score <- ggplot(mutated_motif_score) +
			geom_tile(aes(x=pos_mutation, y=mutation_type, fill=score)) +
			geom_text(aes(x=pos_mutation, y=mutation_type, label=mutation_type)) +
			scale_fill_gradientn(colours=myPalette(100)) +
			# labs(title=paste0(compound_name," motif scores")) +
			labs(x="Mutated position", y="Mutated base", fill="Score") +
			coord_cartesian(expand=FALSE) +
			theme_bw() +
			theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(legend.position="right", legend.direction="horizontal") +
			guides(fill=guide_colourbar(title.vjust=0.5, title.position="top"))
		gp_empty <- ggplot() +
			geom_blank() +
			theme_bw() +
			theme(panel.border=element_blank())

		pdf(file=NULL) # Avoid opening empty window
		gp_combine <- arrangeGrob(grobs=list(gp_detail, gp_score, gp_empty), widths=c(1, 2), heights=c(5, 1.4), layout_matrix=rbind(c(1, 1), c(2, 3)))
		dev.off()

		# Handle warnings from mutated motifs with few/no data points
		withCallingHandlers({
			ggsave(paste0(path_output,"/","Refine_Motifs_",compound_name,".pdf"), gp_combine, width=len_motif*2.8, height=8.5, device="pdf")
		}, warning=function(w){
			if(grepl("no non-missing arguments to max; returning", w$message)){
				invokeRestart("muffleWarning") # Expected and "ok" warning
			}else if(grepl("Computation failed in .*stat_ydensity().*", w$message)){
				invokeRestart("muffleWarning") # Expected and "ok" warning
			}else{
				print(w)
			}
		})

		return(list(tag_mutated_motifs=modification_at_motifs, motif_summary=mutated_motif_score))
	}else{
		print("   Not enough occurrences to plot.")
		
		return(NULL)
	}
}

readline_clean <- function(prompt){
	cat(prompt,"\n")
	response <- readLines("stdin", n=1)

	return(response)
}

auto.confirm.motif <- function(discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle, score_threshold=2){
	refinement_data <- refine.motif(potential_motif, discovered_motifs, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle)
	if(!is.null(refinement_data)){
		res <- refinement_data$motif_summary %>%
			filter(score>score_threshold) %>%
			arrange(pos_mutation, desc(mutation_type)) %>%
			group_by(pos_mutation) %>%
			summarize(bases=paste0(mutation_type, collapse=""), .groups="drop_last") %>%
			group_by(pos_mutation) %>%
			mutate(code=as.character(iupac_nc$code[grepl(paste0("^",bases,"$"), as.character(iupac_nc$choice))]))
		if(nrow(res)>3){
			# Automated motif refinment results

			return(paste0(res$code, collapse=""))
		}else{
			# Automated motif refinment is too short
			return(NA)
		}
	}else{
		return(NA)
	}
}

# TODO shiny interface
user.confirm.motif <- function(discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle){
	has_responded0 <- 0
	while(has_responded0==0){
		response <- readline_clean(paste0("Do you want to refine ",potential_motif," motif? (Y/N)"))
		if(response=="Y"){
			stifle <- refine.motif(potential_motif, discovered_motifs, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle)
			
			has_responded1 <- 0
			while(has_responded1==0){
				response <- readline_clean(paste0("Does the signature is consistent with ",potential_motif," motif? (Y/N)"))
				if(response=="Y"){
					discovered_motifs <- unique(c(discovered_motifs, potential_motif))
					has_responded1 <- 1
				}else if(response=="N"){
					has_responded2 <- 0
					while(has_responded2==0){
						response <- readline_clean("Do you want to try an alternative motif? (N if no alternative)")
						if(response=="N"){
							has_responded2 <- 1
						}else if(nchar(response)>2){
							potential_motif <- response
							discovered_motifs <- user.confirm.motif(discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle)
							has_responded2 <- 1
						}else{
							has_responded2 <- 0
						}
					}
					has_responded1 <- 1
				}else{
					has_responded1 <- 0
				}
			}

			has_responded0 <- 1
		}else if(response=="N"){
			# Do nothung and return original discovered_motifs
			has_responded0 <- 1
		}else{
			has_responded0 <- 0
		}
	}

	return(discovered_motifs)
}

rev.quantile <- function(dataset, value){
	percentage <- ecdf(dataset)(value)

	return(percentage)
}

# 定义一个函数来合并相似的模体
merge_similar_motifs <- function(df, threshold) {
  i <- 1
  while (i <= nrow(df)) {
    # 找到与当前模体相似的模体
    similar_indices <- agrep(df$motif[i], df$motif, max.distance = threshold)
    
    # 如果找到了相似的模体
    if (length(similar_indices) > 1) {
      # 合并相似的模体
      df$motif[i] <- paste(df$motif[similar_indices], collapse = "/")
      df$log_p_value[i] <- min(df$log_p_value[similar_indices])
      
      # 去除相似的模体的行，除了当前行
      similar_indices <- similar_indices[similar_indices != i]
      df <- df[-similar_indices, ]
    }
    
    i <- i + 1
  }
  
  return(df)
}
motif.detection.origin <- function(motifs_peaks_data,final_stat_data, genome, iupac_nc, detection_params, discovered_motifs=NULL){
	threshold <- detection_params$threshold
	nb_peaks <- detection_params$nb_peaks
	seq_params <- detection_params$seq_params
	output_fasta_name <- detection_params$output_fasta_name
	output_bed_name <- detection_params$output_bed_name
	meme_main_dir <- detection_params$meme_main_dir
	meme_analysis_dir <- detection_params$meme_analysis_dir
	nbCPU <- detection_params$nbCPU

	smooth_func <- detection_params$smooth_func
	smooth_win_size <- detection_params$smooth_win_size
	peak_win_size <- detection_params$peak_win_size
	stat_val <- detection_params$stat_val

	automated <- detection_params$automated
	score_threshold <- detection_params$score_threshold

	if(!is.null(discovered_motifs)){
		print_message("Removed explained peaks")
		motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
	}

	nbLoops <- 1
	continue_research <- TRUE
 	no_progress <- FALSE
	final_result <- data.frame()
	while(continue_research){ # TODO Stop when top hits quantile below x threshold? When no more motifs?
		if(is.na(threshold)){
			selected_peaks_data <- filter.best.peaks(motifs_peaks_data, 2000)
		}else{
			selected_peaks_data <- sample.best.peaks(motifs_peaks_data, nb_peaks, threshold)
		}

		# TODO check and create motif_detection
		if(!file.exists(meme_main_dir)){
			dir.create(meme_main_dir)
		}
		extract.peak.sequences(selected_peaks_data, genome, seq_params, output_fasta_name)
		#extract.peak.sequences.homer(selected_peaks_data, genome, seq_params, output_bed_name)
		if(!file.exists(paste0(meme_main_dir,meme_analysis_dir))){
			dir.create(paste0(meme_main_dir,meme_analysis_dir))
		}
		meme_output <- paste0(meme_main_dir,meme_analysis_dir,"meme_",nbLoops)
		dreme_output <- paste0(meme_main_dir,meme_analysis_dir,"dreme_",nbLoops)
		homer_output <- paste0(meme_main_dir,meme_analysis_dir,"homer_",nbLoops)
		cmd_homer <- paste(paste0(path_script,"run_homer.sh"), output_fasta_name, homer_output)
		system(cmd_homer)
		
		execute.meme(output_fasta_name, meme_output, "zoops", 5, 4, 14, nbCPU)

		cmd_dreme <- paste(paste0(path_script,"run_dreme.sh"), output_fasta_name, dreme_output)
		system(cmd_dreme)

		meme_results <- read.meme.output(paste0(meme_output,"/meme.xml"))

		nb_current_peaks <- nrow(selected_peaks_data)
		potential_motifs <- process.meme.output(meme_results, nb_current_peaks)

#--------------------------------------------------------------------------------------------------------------
		doc <- xmlParse(paste0(meme_main_dir,meme_analysis_dir,"dreme_",nbLoops,"/dreme.xml"))

		motif_nodes <- getNodeSet(doc, "//motif")
		motifs <- sapply(motif_nodes, function(node) {xmlGetAttr(node, "seq")})
		p_value <- sapply(motif_nodes, function(node) {xmlGetAttr(node, "pvalue")})

		dreme_result <- data.frame(
		motif = motifs,
		log_p_value = -log10(as.numeric(p_value))
		)


		data <- read.csv(paste0(meme_main_dir,meme_analysis_dir,"homer_",nbLoops,"/knownResults.txt"), header = TRUE, sep = "\t",fileEncoding="utf-8")
		motifs <- data$Consensus
		log_p_value <-data$Log.P.value
		homer_old_result <- data.frame(
		motif = motifs,
		log_p_value = -as.numeric(log_p_value)
		)
		homer_old_result <- homer_old_result[homer_old_result$log_p_value > 30, ]
		#print(homer_old_result)

		# 读取文件
		lines <- readLines(paste0(meme_main_dir,meme_analysis_dir,"homer_",nbLoops,"/homerMotifs.all.motifs"))
		motifs <- c() 
		log_p_values <- c()
		# 遍历每一行
		for (line in lines) {
		# 如果这一行是一个模体的描述信息
		if (startsWith(line, ">")) {
			# 提取模体的名称
			name <- strsplit(line, "\t")[[1]][1]
			motif <- strsplit(name, ">")[[1]][2]
			log_p_value<- strsplit(line, "\t")[[1]][4]
			motifs <- c(motifs, motif)
			log_p_values <- c(log_p_values, log_p_value)
		} 
		}
		homer_new_result <- data.frame(
		motif = motifs,
		log_p_value = -as.numeric(log_p_values)
		)
		homer_new_result <- homer_new_result[homer_new_result$log_p_value > 300, ]
		#print(homer_new_result)

		final_result <- rbind(final_result,dreme_result, homer_old_result, homer_new_result)

#------------------------------------------------------------------------------------------------------------------------
		print_message(paste0(c("Potentials motifs:", paste(potential_motifs)), collapse=" "))
		if(!is.null(potential_motifs)){
			if(automated){
				current_discovered_motifs <- discovered_motifs

				newly_discovered_motifs <- foreach(potential_motif=potential_motifs, .combine=c) %do% {
					graph_subtitle <- paste0("sf: ",smooth_func,", sw: ",smooth_win_size,", pw: ",peak_win_size,", stat: ",stat_val,", np: ",nb_peaks,", th: ",threshold)
					tmp_newly_discovered_motifs <- auto.confirm.motif(current_discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle, score_threshold)
					current_discovered_motifs <- unique(c(current_discovered_motifs, tmp_newly_discovered_motifs))
					current_discovered_motifs <- current_discovered_motifs[!is.na(current_discovered_motifs)] # Remove potential unrefined motif

					return(tmp_newly_discovered_motifs)
				}

 				if(all(is.na(newly_discovered_motifs))){
					# No more motif found
 					no_progress <- TRUE
 				}else{
					discovered_motifs <- unique(c(discovered_motifs, newly_discovered_motifs))	
					discovered_motifs <- discovered_motifs[!is.na(discovered_motifs)] # Remove potential unrefined motif
 				}
			}else{
				newly_discovered_motifs <- foreach(potential_motif=potential_motifs, .combine=c) %do% {
					graph_subtitle <- paste0("sf: ",smooth_func,", sw: ",smooth_win_size,", pw: ",peak_win_size,", stat: ",stat_val,", np: ",nb_peaks,", th: ",threshold)
					tmp_newly_discovered_motifs <- user.confirm.motif(discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle)

					return(tmp_newly_discovered_motifs)
				}
				discovered_motifs <- unique(c(discovered_motifs, newly_discovered_motifs))
			}
		}

		print_message(paste0(c("Discovered motifs:", paste(discovered_motifs)), collapse=" "))

		if(automated && is.null(potential_motifs) ){
			# No more motif found in automated processing
			continue_research <- FALSE
		}else if(automated && no_progress){
			# No more motif found in automated processing
			continue_research <- FALSE
		}else if(automated && !is.null(potential_motifs)){
			# More motifs potentially exist
			print_message("Removing explained peaks")
			motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
			
			nbLoops <- nbLoops + 1
		}else{
			has_responded1 <- 0
			while(has_responded1==0){
				response <- readline_clean("Do you want to continue the research for motifs? (Y/N)")
				if(response=="Y"){
					if(!is.null(discovered_motifs)){
						print_message("Removing explained peaks")
						motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
					}

					nbLoops <- nbLoops + 1
					has_responded1 <- 1
					if(!is.na(threshold)){
						has_responded2 <- 0
						while(has_responded2==0){
							response2 <- readline_clean(paste0("Do you want to continue with threshold fixed at ",threshold,"? (Y/N/th)"))
							if(response2=="Y"){
								# Keep current threshold

								has_responded2 <- 1
							}else if(response2=="N"){
								threshold <- NA

								has_responded2 <- 1
							}else if(!is.na(as.numeric(response2))){ # TODO maybe allow float
								threshold <- as.numeric(response2)

								has_responded2 <- 1
							}
						}
					}
				}else if(response=="N"){
					continue_research <- FALSE

					has_responded1 <- 1
				}
			}
		}
	}
	final_result$motif=as.character(final_result$motif)
	df <- merge_similar_motifs(final_result, threshold = 0.2)
	sink(paste0(meme_main_dir,meme_analysis_dir,"/motif_find_result.txt")) 
	print("DREME and HOMER result:")
	print(df)
	print("-----------------------------------")
	print("MEME result:")
	print(discovered_motifs)
	sink()
	return(discovered_motifs)
}
motif.detection <- function(motifs_peaks_data,flag, final_stat_data, genome, iupac_nc, detection_params, discovered_motifs=NULL){
	threshold <- detection_params$threshold
	nb_peaks <- detection_params$nb_peaks
	seq_params <- detection_params$seq_params
	output_fasta_name <- detection_params$output_fasta_name
	output_fasta_name_win <- detection_params$output_fasta_name_win
	output_bed_name <- detection_params$output_bed_name
	meme_main_dir <- detection_params$meme_main_dir
	meme_analysis_dir <- detection_params$meme_analysis_dir
	nbCPU <- detection_params$nbCPU

	smooth_func <- detection_params$smooth_func
	smooth_win_size <- detection_params$smooth_win_size
	peak_win_size <- detection_params$peak_win_size
	stat_val <- detection_params$stat_val

	automated <- detection_params$automated
	score_threshold <- detection_params$score_threshold

	if(!is.null(discovered_motifs)){
		print_message("Removed explained peaks")
		motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
		#motifs_win_peaks_data <- remove.explained.peaks(motifs_win_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
	}

	nbLoops <- 1
	continue_research <- TRUE
 	no_progress <- FALSE
	final_result <- data.frame()
	while(continue_research){ # TODO Stop when top hits quantile below x threshold? When no more motifs?
		if(is.na(threshold)){
			if(flag==0){
				selected_peaks_data <- filter.best.peaks(motifs_peaks_data, 2000)
			}
			else{
				selected_peaks_data <- filter.best.win.peaks(motifs_peaks_data, 2000)
			}
			#selected_peaks_data <- rbind(selected_peaks_data, selected_win_peaks_data)
		}else{
			selected_peaks_data <- sample.best.peaks(motifs_peaks_data, nb_peaks, threshold)
		}

		# TODO check and create motif_detection
		if(!file.exists(meme_main_dir)){
			dir.create(meme_main_dir)
		}
		extract.peak.sequences(selected_peaks_data, genome, seq_params, output_fasta_name)
		#extract.peak.sequences.homer(selected_peaks_data, genome, seq_params, output_bed_name)
		if(!file.exists(paste0(meme_main_dir,meme_analysis_dir))){
			dir.create(paste0(meme_main_dir,meme_analysis_dir))
		}
		meme_output <- paste0(meme_main_dir,meme_analysis_dir,"meme_",nbLoops)
		dreme_output <- paste0(meme_main_dir,meme_analysis_dir,"dreme_",nbLoops)
		homer_output <- paste0(meme_main_dir,meme_analysis_dir,"homer_",nbLoops)


		execute.meme(output_fasta_name, meme_output, "zoops", 5, 4, 14, nbCPU)
		meme_results <- read.meme.output(paste0(meme_output,"/meme.xml"))
		nb_current_peaks <- nrow(selected_peaks_data)
		potential_motifs <- process.meme.output(meme_results, nb_current_peaks)

		print_message(paste0(c("Potentials motifs:", paste(potential_motifs)), collapse=" "))
		if(!is.null(potential_motifs)){
			if(automated){
				current_discovered_motifs <- discovered_motifs

				newly_discovered_motifs <- foreach(potential_motif=potential_motifs, .combine=c) %do% {
					graph_subtitle <- paste0("sf: ",smooth_func,", sw: ",smooth_win_size,", pw: ",peak_win_size,", stat: ",stat_val,", np: ",nb_peaks,", th: ",threshold)
					tmp_newly_discovered_motifs <- auto.confirm.motif(current_discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle, score_threshold)
					current_discovered_motifs <- unique(c(current_discovered_motifs, tmp_newly_discovered_motifs))
					current_discovered_motifs <- current_discovered_motifs[!is.na(current_discovered_motifs)] # Remove potential unrefined motif

					return(tmp_newly_discovered_motifs)
				}

 				if(all(is.na(newly_discovered_motifs))){
					# No more motif found
 					no_progress <- TRUE
 				}else{
					discovered_motifs <- unique(c(discovered_motifs, newly_discovered_motifs))	
					discovered_motifs <- discovered_motifs[!is.na(discovered_motifs)] # Remove potential unrefined motif
 				}
			}else{
				newly_discovered_motifs <- foreach(potential_motif=potential_motifs, .combine=c) %do% {
					graph_subtitle <- paste0("sf: ",smooth_func,", sw: ",smooth_win_size,", pw: ",peak_win_size,", stat: ",stat_val,", np: ",nb_peaks,", th: ",threshold)
					tmp_newly_discovered_motifs <- user.confirm.motif(discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle)

					return(tmp_newly_discovered_motifs)
				}
				discovered_motifs <- unique(c(discovered_motifs, newly_discovered_motifs))
			}
		}

		print_message(paste0(c("Discovered motifs:", paste(discovered_motifs)), collapse=" "))

		if(automated && is.null(potential_motifs) ){
			# No more motif found in automated processing
			continue_research <- FALSE
		}else if(automated && no_progress){
			# No more motif found in automated processing
			continue_research <- FALSE
		}else if(automated && !is.null(potential_motifs)){
			# More motifs potentially exist
			print_message("Removing explained peaks")
			motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
			nbLoops <- nbLoops + 1
		}else{
			has_responded1 <- 0
			while(has_responded1==0){
				response <- readline_clean("Do you want to continue the research for motifs? (Y/N)")
				if(response=="Y"){
					if(!is.null(discovered_motifs)){
						print_message("Removing explained peaks")
						motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
					}

					nbLoops <- nbLoops + 1
					has_responded1 <- 1
					if(!is.na(threshold)){
						has_responded2 <- 0
						while(has_responded2==0){
							response2 <- readline_clean(paste0("Do you want to continue with threshold fixed at ",threshold,"? (Y/N/th)"))
							if(response2=="Y"){
								# Keep current threshold

								has_responded2 <- 1
							}else if(response2=="N"){
								threshold <- NA

								has_responded2 <- 1
							}else if(!is.na(as.numeric(response2))){ # TODO maybe allow float
								threshold <- as.numeric(response2)

								has_responded2 <- 1
							}
						}
					}
				}else if(response=="N"){
					continue_research <- FALSE

					has_responded1 <- 1
				}
			}
		}
	}


	return(discovered_motifs)
}
motif.detection.mix <- function(motifs_peaks_data,motifs_win_peaks_data, final_stat_data, genome, iupac_nc, detection_params, discovered_motifs=NULL){
	threshold <- detection_params$threshold
	nb_peaks <- detection_params$nb_peaks
	seq_params <- detection_params$seq_params
	output_fasta_name <- detection_params$output_fasta_name
	output_bed_name <- detection_params$output_bed_name
	meme_main_dir <- detection_params$meme_main_dir
	meme_analysis_dir <- detection_params$meme_analysis_dir
	nbCPU <- detection_params$nbCPU

	smooth_func <- detection_params$smooth_func
	smooth_win_size <- detection_params$smooth_win_size
	peak_win_size <- detection_params$peak_win_size
	stat_val <- detection_params$stat_val

	automated <- detection_params$automated
	score_threshold <- detection_params$score_threshold

	if(!is.null(discovered_motifs)){
		print_message("Removed explained peaks")
		motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
		motifs_win_peaks_data <- remove.explained.peaks(motifs_win_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
	}

	nbLoops <- 1
	continue_research <- TRUE
 	no_progress <- FALSE
	final_result <- data.frame()
	while(continue_research){ # TODO Stop when top hits quantile below x threshold? When no more motifs?
		if(is.na(threshold)){
			selected_peaks_data <- filter.best.peaks(motifs_peaks_data, 1500)
			selected_win_peaks_data <- filter.best.win.peaks(motifs_win_peaks_data, 500)
			selected_peaks_data <- rbind(selected_peaks_data, selected_win_peaks_data)
		}else{
			selected_peaks_data <- sample.best.peaks(motifs_peaks_data, nb_peaks, threshold)
		}

		# TODO check and create motif_detection
		if(!file.exists(meme_main_dir)){
			dir.create(meme_main_dir)
		}
		extract.peak.sequences(selected_peaks_data, genome, seq_params, output_fasta_name)
		#extract.peak.sequences.homer(selected_peaks_data, genome, seq_params, output_bed_name)
		if(!file.exists(paste0(meme_main_dir,meme_analysis_dir))){
			dir.create(paste0(meme_main_dir,meme_analysis_dir))
		}
		meme_output <- paste0(meme_main_dir,meme_analysis_dir,"meme_",nbLoops)
		dreme_output <- paste0(meme_main_dir,meme_analysis_dir,"dreme_",nbLoops)
		homer_output <- paste0(meme_main_dir,meme_analysis_dir,"homer_",nbLoops)


		execute.meme(output_fasta_name, meme_output, "zoops", 5, 4, 14, nbCPU)
		meme_results <- read.meme.output(paste0(meme_output,"/meme.xml"))
		nb_current_peaks <- nrow(selected_peaks_data)
		potential_motifs <- process.meme.output(meme_results, nb_current_peaks)

		print_message(paste0(c("Potentials motifs:", paste(potential_motifs)), collapse=" "))
		if(!is.null(potential_motifs)){
			if(automated){
				current_discovered_motifs <- discovered_motifs

				newly_discovered_motifs <- foreach(potential_motif=potential_motifs, .combine=c) %do% {
					graph_subtitle <- paste0("sf: ",smooth_func,", sw: ",smooth_win_size,", pw: ",peak_win_size,", stat: ",stat_val,", np: ",nb_peaks,", th: ",threshold)
					tmp_newly_discovered_motifs <- auto.confirm.motif(current_discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle, score_threshold)
					current_discovered_motifs <- unique(c(current_discovered_motifs, tmp_newly_discovered_motifs))
					current_discovered_motifs <- current_discovered_motifs[!is.na(current_discovered_motifs)] # Remove potential unrefined motif

					return(tmp_newly_discovered_motifs)
				}

 				if(all(is.na(newly_discovered_motifs))){
					# No more motif found
 					no_progress <- TRUE
 				}else{
					discovered_motifs <- unique(c(discovered_motifs, newly_discovered_motifs))	
					discovered_motifs <- discovered_motifs[!is.na(discovered_motifs)] # Remove potential unrefined motif
 				}
			}else{
				newly_discovered_motifs <- foreach(potential_motif=potential_motifs, .combine=c) %do% {
					graph_subtitle <- paste0("sf: ",smooth_func,", sw: ",smooth_win_size,", pw: ",peak_win_size,", stat: ",stat_val,", np: ",nb_peaks,", th: ",threshold)
					tmp_newly_discovered_motifs <- user.confirm.motif(discovered_motifs, potential_motif, final_stat_data, genome, nbCPU, seq_params, iupac_nc, meme_output, graph_subtitle)

					return(tmp_newly_discovered_motifs)
				}
				discovered_motifs <- unique(c(discovered_motifs, newly_discovered_motifs))
			}
		}

		print_message(paste0(c("Discovered motifs:", paste(discovered_motifs)), collapse=" "))

		if(automated && is.null(potential_motifs) ){
			# No more motif found in automated processing
			continue_research <- FALSE
		}else if(automated && no_progress){
			# No more motif found in automated processing
			continue_research <- FALSE
		}else if(automated && !is.null(potential_motifs)){
			# More motifs potentially exist
			print_message("Removing explained peaks")
			motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
			motifs_win_peaks_data <- remove.explained.peaks(motifs_win_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
			nbLoops <- nbLoops + 1
		}else{
			has_responded1 <- 0
			while(has_responded1==0){
				response <- readline_clean("Do you want to continue the research for motifs? (Y/N)")
				if(response=="Y"){
					if(!is.null(discovered_motifs)){
						print_message("Removing explained peaks")
						motifs_peaks_data <- remove.explained.peaks(motifs_peaks_data, discovered_motifs, seq_params, genome, iupac_nc, nbCPU)
					}

					nbLoops <- nbLoops + 1
					has_responded1 <- 1
					if(!is.na(threshold)){
						has_responded2 <- 0
						while(has_responded2==0){
							response2 <- readline_clean(paste0("Do you want to continue with threshold fixed at ",threshold,"? (Y/N/th)"))
							if(response2=="Y"){
								# Keep current threshold

								has_responded2 <- 1
							}else if(response2=="N"){
								threshold <- NA

								has_responded2 <- 1
							}else if(!is.na(as.numeric(response2))){ # TODO maybe allow float
								threshold <- as.numeric(response2)

								has_responded2 <- 1
							}
						}
					}
				}else if(response=="N"){
					continue_research <- FALSE

					has_responded1 <- 1
				}
			}
		}
	}


	return(discovered_motifs)
}
wrapper.motif.detection.origin <- function(final_stat_data, genome, path_output, detection_params, iupac_nc, known_motifs=NULL, rm_homopolymer=NA){
	library(XML)
	initial_dir <- getwd()
	setwd(path_output)

	print_message("Processing statistical signal")
	peaks_data <- process.statistical.signal.origin(final_stat_data, detection_params)

	if(!is.na(rm_homopolymer)){
		print_message(paste0("Remove peaks explained by homopolymer of length ",rm_homopolymer," or more."))
		peaks_data <- remove.homopolymer(peaks_data, genome, rm_homopolymer, iupac_nc, detection_params$seq_params, detection_params$nbCPU)
	}

	discovered_motifs <- tryCatch({
			print_message("Searching potential modified motifs")
			motif.detection.origin(peaks_data, final_stat_data, genome, iupac_nc, detection_params, known_motifs)
		}, error=function(e) {
			print("Unexpected error.")
			print(e)

			return(NA)
		}
	)

	setwd(initial_dir)

	return(discovered_motifs)
}
wrapper.motif.detection <- function(final_stat_data, genome, path_output, detection_params,detection_params_win, iupac_nc, known_motifs=NULL, rm_homopolymer=NA){
	library(XML)
	initial_dir <- getwd()
	setwd(path_output)
	#把位点保留下来
	print_message("Processing statistical signal")
	data <- process.statistical.signal(final_stat_data, detection_params)

	if(!is.na(rm_homopolymer)){
		print_message(paste0("Remove peaks explained by homopolymer of length ",rm_homopolymer," or more."))
		data$peaks_data <- remove.homopolymer(data$peaks_data, genome, rm_homopolymer, iupac_nc, detection_params$seq_params, detection_params$nbCPU)
		data$win_peaks_data<- remove.homopolymer(data$win_peaks_data, genome, rm_homopolymer, iupac_nc, detection_params$seq_params, detection_params$nbCPU)
	}
	#print(head(peaks_data))
	#print(head(win_peaks_data))
	discovered_motifs <- tryCatch({
			print_message("Searching potential modified motifs")
			motif.detection.origin(data$peaks_data, final_stat_data, genome, iupac_nc, detection_params, known_motifs)
			#motif.detection(data$win_peaks_data,1, final_stat_data, genome, iupac_nc, detection_params_win, known_motifs)
			motif.detection.mix(data$peaks_data,data$win_peaks_data,final_stat_data, genome,iupac_nc,detection_params_win,known_motifs)
		}, error=function(e) {
			print("Unexpected error.")
			print(e)

			return(NA)
		}
	)

	setwd(initial_dir)

	return(discovered_motifs)
}

#  ___  ___      _   _  __       _____                   _   _             
#  |  \/  |     | | (_)/ _|     /  __ \                 | | (_)            
#  | .  . | ___ | |_ _| |_ ___  | /  \/ ___  _   _ _ __ | |_ _ _ __   __ _ 
#  | |\/| |/ _ \| __| |  _/ __| | |    / _ \| | | | '_ \| __| | '_ \ / _` |
#  | |  | | (_) | |_| | | \__ \ | \__/\ (_) | |_| | | | | |_| | | | | (_| |
#  \_|  |_/\___/ \__|_|_| |___/  \____/\___/ \__,_|_| |_|\__|_|_| |_|\__, |
#                                                                     __/ |
#                                                                    |___/ 

count.motifs <- function(genome, motif_summary, iupac_nc, nbCPU){
	left_signal <- -1 # Remove nothing
	right_signal <- -1
	error_margin <- -1

	print(paste0("Searching motifs."))
	motifs <- find.isolated.motifs(genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, nbCPU, FALSE)
	motifs_count <- motifs %>%
		group_by(motif) %>%
		summarize(n=n(), .groups="drop_last")

	return(motifs_count)
}

#  ___  ___      _   _  __       _____                _             
#  |  \/  |     | | (_)/ _|     /  ___|              (_)            
#  | .  . | ___ | |_ _| |_ ___  \ `--.  ___ ___  _ __ _ _ __   __ _ 
#  | |\/| |/ _ \| __| |  _/ __|  `--. \/ __/ _ \| '__| | '_ \ / _` |
#  | |  | | (_) | |_| | | \__ \ /\__/ / (_| (_) | |  | | | | | (_| |
#  \_|  |_/\___/ \__|_|_| |___/ \____/ \___\___/|_|  |_|_| |_|\__, |
#                                                              __/ |
#                                                             |___/ 

# Give similar results than signature; large overlapping considered
tag.motifs <- function(motif_to_tag, discovered_motifs, methylation_signal, genome, iupac_nc, expected_signal_left, expected_signal_right, signal_margin, min_cov){
	if(!is.null(discovered_motifs)){
		motifs_ranges <- find.motifs.ranges(genome, c(motif_to_tag, discovered_motifs), iupac_nc)
	}else{
		motifs_ranges <- find.motifs.ranges(genome, motif_to_tag, iupac_nc)
	}
	# Methylation position unknown between motif start and end
	# Nanopolish do not produce synmetric signal alignments, it's aligned on 6-mer start position - 1 with fwd strand referential
	# e.g. 1000:CGATCT pA is 999 on fwd strand; 1000:CGATCT pA is 994 on rev strand
	# Strand don't mater for expected signal because looks for overlaps between consecutive motifs on same strand
	motifs_ranges <- motifs_ranges %>%
		mutate(left_side=start + expected_signal_left - signal_margin) %>% #-12
		mutate(right_side=end + expected_signal_right + signal_margin) #5

	gr_motifs_ranges <- GRanges(
		seqnames=motifs_ranges$contig_name,
		ranges=IRanges(start=motifs_ranges$left_side, end=motifs_ranges$right_side),
		strand=ifelse(motifs_ranges$dir=="fwd","+","-") # Is strand specific
	)

	# Remove overlapping motifs
	overlapping_motifs_ranges <- overlapsAny(gr_motifs_ranges, type="any", drop.self=TRUE)
	gr_motifs_ranges <- gr_motifs_ranges[!overlapping_motifs_ranges]
	motifs_ranges <- motifs_ranges[!overlapping_motifs_ranges,]

	isolated_motifs_ranges <- subset(motifs_ranges, motif==motif_to_tag)
	# Correct motif range according to strand, needed to compute distance and align signal
	isolated_motifs_ranges <- isolated_motifs_ranges %>%
		mutate(start=ifelse(dir=="fwd",start,end)) %>%
		mutate(end=ifelse(dir=="fwd",end,start))

	gr_isolated_motifs_ranges <- GRanges(
		seqnames=isolated_motifs_ranges$contig_name,
		ranges=IRanges(start=isolated_motifs_ranges$left_side, end=isolated_motifs_ranges$right_side),
		strand=ifelse(isolated_motifs_ranges$dir=="fwd","+","-") # Is strand specific
	)

	gr_methylation <- GRanges(
		seqnames=methylation_signal$contig,
		ranges=IRanges(methylation_signal$position, methylation_signal$position),
		strand=ifelse(methylation_signal$dir=="fwd","+","-")
	)
	#find origin signal
	overlaps_motifs_methylation <- findOverlaps(gr_isolated_motifs_ranges, gr_methylation, type="any", select="all")
	methylation_at_motifs_ranges <- data.frame(
		contig=isolated_motifs_ranges$contig_name[overlaps_motifs_methylation@from],
		pos_motif=isolated_motifs_ranges$start[overlaps_motifs_methylation@from],
		motif=isolated_motifs_ranges$motif[overlaps_motifs_methylation@from],
		pos_signal=methylation_signal$position[overlaps_motifs_methylation@to],
		dir=methylation_signal$dir[overlaps_motifs_methylation@to],
		strand=methylation_signal$strand[overlaps_motifs_methylation@to],
		N_wga=methylation_signal$N_wga[overlaps_motifs_methylation@to],
		N_nat=methylation_signal$N_nat[overlaps_motifs_methylation@to],
		mean_diff=methylation_signal$mean_diff[overlaps_motifs_methylation@to]
	)

	methylation_at_motifs_ranges <- methylation_at_motifs_ranges %>% 
		mutate(distance=ifelse(dir=="fwd",pos_signal-pos_motif,(-(pos_signal-pos_motif)) - 7)) # Relative distance to mod_pos with strand correction

	filtered_methylation_to_motifs_ranges <- methylation_at_motifs_ranges %>%
		filter(N_wga>=min_cov & N_nat>=min_cov)

	return(filtered_methylation_to_motifs_ranges)
}

gtable_delete_row <- function(gtable, start_row, end_row){
	toRemove <- gtable$layout$t %in% seq(start_row,end_row)
	original_names <- gtable$layout$name[!toRemove]
	gtable$layout$name[!toRemove] <- paste0("toKeep_",original_names)
	gtable <- gtable_filter(gtable, "toKeep_") # Remove grobs
	gtable$heights <- gtable$heights[-seq(start_row,end_row)] # Remove hights
	gtable$layout$t <- ifelse(gtable$layout$t>end_row, start_row + (gtable$layout$t - end_row - 1), gtable$layout$t) # Correct t
	gtable$layout$b <- ifelse(gtable$layout$b>end_row, start_row + (gtable$layout$b - end_row - 1), gtable$layout$b) # Correct b
	gtable$layout$name <- original_names

	return(gtable)
}

draw.signature.center.detection <- function(tagged_final_stat_data, motif_score, modified_base_prediction, motif, expected_signal_left, expected_signal_right, signal_margin, signal_len, mod_pos, plot_name){
	len_motif <- nchar(motif)
	nb_bases_pre_motif <- abs(expected_signal_left - signal_margin + 4) # 4 is correction to align signal to 6-mer peaks (-4 to 0)
	nb_bases_post_motif <- abs(expected_signal_right + signal_margin + 4)

	gp_score <- ggplot(motif_score) +
		geom_bar(aes(x=as.factor(distance), weight=score_window)) +
		labs(x=paste0("Relative position from motif (-",nb_bases_pre_motif,":motif:",(len_motif-1)+nb_bases_post_motif,")"), y="Signal strength") +
		theme_bw() +
		theme(legend.position="none") +
		theme(plot.title=element_blank())

	# Convert to include pre/post-motif window
	predicted_mod_relative_pos <- modified_base_prediction$predicted_mod_pos + nb_bases_pre_motif + 1 # Because predicted_mod_pos 0-based
	chosen_mod_relative_pos <- modified_base_prediction$possible_mod_pos + nb_bases_pre_motif + 1 # Because possible_mod_pos 0-based

	motif_labelling <- c(rep("n",nb_bases_pre_motif),strsplit(motif,"")[[1]],rep("n",nb_bases_post_motif))
	motif_labelling_position <- seq(1,nchar(motif) + signal_len + (signal_margin*2 - 1))
	motif_labelling_fontface <- rep("plain",length(motif_labelling_position))
	motif_labelling_fontface[predicted_mod_relative_pos] <- "italic" # Predicted signature center
	motif_labelling_fontface[chosen_mod_relative_pos] <- "bold" # Possible signature center
	motif_labelling_color <- rep("#000000",length(motif_labelling_position))
	if(!is.na(mod_pos)){
		motif_labelling_color[nb_bases_pre_motif + mod_pos] <- "#FF0000" # Expected signature center
	}
	motif_information <- data.frame(
		x=motif_labelling_position,
		y=max(-15,min(tagged_final_stat_data$mean_diff, na.rm=TRUE)),
		label=motif_labelling,
		fontface=motif_labelling_fontface,
		color=motif_labelling_color
	)
	gp_motif <- ggplot(motif_information) +
		geom_text(aes(x=as.factor(x), y=y, label=label, fontface=fontface, col=color)) +
		scale_colour_manual(values=c("#000000"="black", "#FF0000"="red")) +
		scale_y_continuous(breaks=motif_labelling_position) +
		theme_bw() +
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
		theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
		theme(plot.title=element_blank()) +
		theme(legend.position="none") +
		theme(panel.border=element_rect(colour="black", fill=NA, size=1))

	gp_signal <- ggplot(subset(tagged_final_stat_data,!is.na(mean_diff))) + # Remove remaining NA
		geom_jitter(aes(x=as.factor(distance), y=mean_diff), pch=46, height=0) +
		geom_violin(aes(x=as.factor(distance), y=mean_diff), alpha=0.5) +
		geom_hline(yintercept=0, col="red") +
		# facet_grid(dir~.) + # Check strand consitency
		coord_cartesian(ylim=c(max(-15,min(tagged_final_stat_data$mean_diff, na.rm=TRUE)),min(15,max(tagged_final_stat_data$mean_diff, na.rm=TRUE)))) +
		labs(title=paste0("Research signature center for ",motif," motif")) +
		labs(y="Mean current differences (pA)") +
		theme_bw() +
		theme(legend.position="none") +
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

	pdf(file=NULL) # Avoid opening empty window	
	grob_gp_signal <- ggplotGrob(gp_signal)
	grob_gp_motif <- ggplotGrob(gp_motif)
	grob_gp_score <- ggplotGrob(gp_score)
	dev.off() # Remove ggplotGrob window
	combine_graph <- rbind(grob_gp_signal, grob_gp_motif, grob_gp_score, size="first")
	combine_graph$widths <- unit.pmax(grob_gp_signal$widths, grob_gp_score$widths) # Match largest margin

	tmp <- gtable_delete_row(combine_graph, 7, 15)
	# gtable_show_layout(tmp)
	clean_combine_graph <- gtable_delete_row(tmp, 8, 16)
	clean_combine_graph$heights[7] <- unit(0.5, "in") # Reduce motif panel height

	pdf(plot_name, width=4+0.5*len_motif, height=9)
	grid.draw(clean_combine_graph)
	dev.off()
}

correct.prediction.bipartite.motif <- function(motif, predicted_mod_pos){
	motif_position <- seq(1,nchar(motif))
	motif_N_position <- grepl("N",strsplit(motif,"")[[1]])
	distance_to_prediction <- abs(motif_position[!motif_N_position] - (predicted_mod_pos + 1) ) # predicted_mod_pos is 0-based
	alternative_positions <- which(distance_to_prediction==min(distance_to_prediction))
	if(length(alternative_positions)==1){ # Only one best alternative
		return(alternative_positions - 1) # Convert to 0-based
	}else if(length(alternative_positions)>1){ # Multiple alternatives
		if(length(motif_N_position[motif_N_position==TRUE])==1){ # If one N
			return(predicted_mod_pos) # Do not change prediction
		}else{ # Multiple N with prediction centered
			return(NA) # Do not propose prediction
		}
	}
}

find.signature.center <- function(final_stat_data, motif_summary, genome, nbCPU, seq_params, iupac_nc, signal_margin, signal_len, base_name, path_output=NA){
	expected_signal_left <- -6
	expected_signal_right <- -1
	min_cov <- 0 # TODO Can tweak

	results <- foreach(idx_motif=seq(1,nrow(motif_summary)), .combine=rbind) %do% {
		motif <- motif_summary$motif[idx_motif]
		print_message(paste0("  Process ",motif))
		mod_pos <- motif_summary$mod_pos[idx_motif]
		discovered_motifs <- motif_summary$motif[-c(idx_motif)] # Every discovered motifs except current

		# Signal with motif start reference frame
		print_message(paste0("    Tag ",motif," occurrences"))
		tagged_final_stat_data <- tag.motifs(motif, discovered_motifs, final_stat_data, genome, iupac_nc, expected_signal_left, expected_signal_right, signal_margin, min_cov)
		# Align signal with sequence
		# Nanopolish used 6-mer signal at -1 from 6-mer, so 0 with motif start reference frame
		# Signal peaks seems at 3/6 of 6-mer (-4 to 0)
		tagged_final_stat_data$distance <- tagged_final_stat_data$distance + 4

		print_message(paste0("    Score ",motif," modified position"))
		motif_score <- tagged_final_stat_data %>% # TODO could maybe improve summary function
			group_by(distance) %>%
			summarize(score_position=abs(mean(mean_diff, na.rm=TRUE)), .groups="drop_last") %>%
			mutate(score_window=rollapplyr(score_position, 5, mean, partial=TRUE, fill=NA, align="center"))
		signature_center <- motif_score %>%
			filter(max(score_window)==score_window)
		# Motif reference frame, 0-based (e.g. Gm6ATC -> 0.1.2.3 and m6A=1)
		modified_base_prediction <- cbind(motif_summary[idx_motif,], signature_center) %>%
			mutate(predicted_mod_pos=distance) %>% # Distance is relative to motif start; 0-based
			mutate(len_motif=nchar(motif)) %>%
			mutate(possible_mod_pos=ifelse(predicted_mod_pos+1>len_motif, len_motif-1, ifelse(predicted_mod_pos<0, 0, ifelse(strsplit(motif,"")[[1]][predicted_mod_pos+1]=="N", correct.prediction.bipartite.motif(motif, predicted_mod_pos), predicted_mod_pos)))) # Predicted further than motif | before motif | degenerated (N)

		plot_name <- paste0("Signature_center_detection_",base_name,"_",motif,"_v1.pdf")
		if(!is.na(path_output)){
			plot_name <- paste0(path_output,plot_name)
		}
		draw.signature.center.detection(tagged_final_stat_data, motif_score, modified_base_prediction, motif, expected_signal_left, expected_signal_right, signal_margin, signal_len, mod_pos, plot_name)

		return(modified_base_prediction)
	}

	return(results)
}

#  ___  ___      _   _  __       _____ _                          _            _          _   _             
#  |  \/  |     | | (_)/ _|     /  __ \ |                        | |          (_)        | | (_)            
#  | .  . | ___ | |_ _| |_ ___  | /  \/ |__   __ _ _ __ __ _  ___| |_ ___ _ __ _ ______ _| |_ _  ___  _ __  
#  | |\/| |/ _ \| __| |  _/ __| | |   | '_ \ / _` | '__/ _` |/ __| __/ _ \ '__| |_  / _` | __| |/ _ \| '_ \ 
#  | |  | | (_) | |_| | | \__ \ | \__/\ | | | (_| | | | (_| | (__| ||  __/ |  | |/ / (_| | |_| | (_) | | | |
#  \_|  |_/\___/ \__|_|_| |___/  \____/_| |_|\__,_|_|  \__,_|\___|\__\___|_|  |_/___\__,_|\__|_|\___/|_| |_|
#                                                                                                           
#                                                                                                           

characterize.signature <- function(classification_data, base_name){
	# Retrieve color palette
	tmp <- unique(subset(attr(classification_data, "annotation_mod"), select=c("label","col_label")))
	color_mod <- as.character(tmp$col_label)
	names(color_mod) <- as.character(tmp$label)

	# Prepare data
	classification_data$id <- rownames(classification_data)
	rownames(classification_data) <- NULL
	classification_data$motif <- attr(classification_data, "annotation_motif")$label
	classification_data$mod <- attr(classification_data, "annotation_mod")$label
	classification_data$dir <- attr(classification_data, "annotation_dir")$label
	classification_data$strain <- attr(classification_data, "annotation_strain")$label
	classification_data <- classification_data %>%
		gather(distance, mean_diff, -c(id, motif, mod, dir, strain))

	classification_data$distance <- as.numeric(classification_data$distance) + 4 # Center third 6-mer position on modified base
	ordered_distances <- as.character(sort(as.numeric(unique(classification_data$distance))))
	classification_data$distance <- ordered(classification_data$distance, levels=ordered_distances)

	# Plotting
	gp_overall <- ggplot(classification_data, aes(x=as.factor(distance), y=mean_diff)) +
		geom_boxplot(outlier.shape=NA) +
		geom_hline(yintercept=0, col="red") +
		labs(title=paste0("Overall motif signature")) +
		labs(y="Mean current differences (pA)", x="Distance from methylated base") +
		coord_cartesian(ylim=c(-10,10), expand=FALSE) +
		theme_bw() +
		theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold"))
	pdf(paste0("Signature_Statistics_Overall_",base_name,"_v1.pdf"), width=5, height=3)
	print(gp_overall)
	dev.off()

	stat_names <- c(
		abs_mean="Average\nabsolute",
		sd="Standard\ndeviation"
	)

	classification_data_summary_all <- classification_data %>%
		group_by(distance) %>%
		summarize(abs_mean=abs(mean(mean_diff,na.rm=TRUE)), sd=sd(mean_diff,na.rm=TRUE), .groups="drop_last") %>%
		gather(stat, value, -c(distance))

	gp_overall_summary <- ggplot(classification_data_summary_all) +
		geom_bar(aes(x=distance, y=value), stat="identity", position=position_dodge()) +
		facet_grid(stat~., scales="free_y", labeller=labeller(stat=stat_names)) +
		labs(title=paste0("Overall motif signature characteristics")) +
		labs(x="Distance from methylated base", y="Mean current differences (pA)") +
		theme_bw() +
		theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold")) +
		theme(strip.background=element_rect(fill="white"), strip.text=element_text(size=10))
	pdf(paste0("Signature_Statistics_Overall_Summary_",base_name,"_v1.pdf"), width=6, height=3)
	print(gp_overall_summary)
	dev.off()

	gp_overall_ByModType <- ggplot(classification_data, aes(x=as.factor(distance), y=mean_diff)) +
		geom_boxplot(outlier.shape=NA) +
		stat_summary(fun.y=mean, geom="point", shape=16, size=1) +
		facet_grid(mod~.) +
		labs(title=paste0("Motif signatures by methylation type")) +
		labs(x="Distance from methylated base", y="Mean current differences (pA)") +
		geom_hline(yintercept=0, col="red") +
		coord_cartesian(ylim=c(-10,10), expand=FALSE) +
		theme_bw() +
		theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold")) +
		theme(strip.background=element_rect(fill="white"), strip.text=element_text(size=10)) +
		theme(panel.spacing=unit(0.5, "lines"))
	pdf(paste0("Signature_Statistics_ByModType_",base_name,"_v1.pdf"), width=6, height=5)
	print(gp_overall_ByModType)
	dev.off()

	classification_data_summary_mod <- classification_data %>%
		group_by(mod, distance) %>%
		summarize(abs_mean=abs(mean(mean_diff,na.rm=TRUE)), sd=sd(mean_diff,na.rm=TRUE), .groups="drop_last") %>%
		gather(stat, value, -c(mod, distance))

	gp_ByModType_summary <- ggplot(classification_data_summary_mod) +
		geom_bar(aes(x=distance, y=value, fill=mod), stat="identity", position=position_dodge()) +
		facet_grid(stat~., scales="free_y", labeller=labeller(stat=stat_names)) +
		scale_fill_manual(values=color_mod) +
		labs(title=paste0("Motif signature characteristics by methylation type")) +
		labs(x="Distance from methylated base", y="Mean current differences (pA)") +
		labs(fill="Methylation\nTypes") +
		theme_bw() +
		theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold")) +
		theme(strip.background=element_rect(fill="white"), strip.text=element_text(size=10))
	pdf(paste0("Signature_Statistics_ByModType_Summary_",base_name,"_v1.pdf"), width=7, height=5)
	print(gp_ByModType_summary)
	dev.off()

	classification_data_summary_motif <- classification_data %>%
		group_by(motif, distance) %>%
		summarize(abs_mean=abs(mean(mean_diff,na.rm=TRUE)), sd=sd(mean_diff,na.rm=TRUE), .groups="drop_last") %>%
		gather(stat, value, -c(motif, distance))

	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	motif_col <- myPalette(length(unique(classification_data_summary_motif$motif)))
	names(motif_col) <- unique(classification_data_summary_motif$motif)
	gp_ByMotif_summary <- ggplot(classification_data_summary_motif) +
		geom_bar(aes(x=distance, y=value, fill=motif), stat="identity", position=position_dodge()) +
		facet_grid(stat~., scales="free_y", labeller=labeller(stat=stat_names)) +
		scale_fill_manual(values=motif_col) +
		labs(title=paste0("Signature statistics by motif")) +
		labs(x="Distance from methylated base", y="Mean current differences (pA)") +
		guides(fill=FALSE) +
		theme_bw() +
		theme(axis.text=element_text(size=8), axis.title=element_text(size=10,face="bold")) +
		theme(strip.background=element_rect(fill="white"), strip.text=element_text(size=10))
	pdf(paste0("Signature_Statistics_ByMotif_Summary_",base_name,"_v1.pdf"), width=9, height=5)
	print(gp_ByMotif_summary)
	dev.off()
}

#  ___  ___      _   _  __       _____ _               _  __ _           _   _             
#  |  \/  |     | | (_)/ _|     /  __ \ |             (_)/ _(_)         | | (_)            
#  | .  . | ___ | |_ _| |_ ___  | /  \/ | __ _ ___ ___ _| |_ _  ___ __ _| |_ _  ___  _ __  
#  | |\/| |/ _ \| __| |  _/ __| | |   | |/ _` / __/ __| |  _| |/ __/ _` | __| |/ _ \| '_ \ 
#  | |  | | (_) | |_| | | \__ \ | \__/\ | (_| \__ \__ \ | | | | (_| (_| | |_| | (_) | | | |
#  \_|  |_/\___/ \__|_|_| |___/  \____/_|\__,_|___/___/_|_| |_|\___\__,_|\__|_|\___/|_| |_|
#                                                                                          
#                                                                                          

pca.plot <- function(pca_data, annotation, PCx="PC1", PCy="PC2") {
	data <- data.frame(id=row.names(pca_data$x), pca_data$x)
	data <- merge(data, annotation, by.x="id", by.y="id")

	labels <- unique(subset(data, select=c("label","col_label")))
	labels[] <- lapply(labels, as.character) # convert to character

	gp <- ggplot(data, aes_string(x=PCx, y=PCy, colour="label")) +
		geom_point(alpha=1, pch=46) +
		geom_hline(aes(yintercept=0), size=0.3) +
		geom_vline(aes(xintercept=0), size=0.3) +
		coord_equal() +
		scale_colour_manual(name="Colors:", breaks=labels$label, values=labels$col_label, labels=labels$label) +
		labs(title=paste0("Biplot projection on ",PCx," and ",PCy)) +
		theme(legend.position="right")
	print(gp)
}

pca.3D.plot <- function(pca_data, annotation, noAlpha=FALSE, PCx="PC1", PCy="PC2", PCz="PC3"){
	data <- data.frame(id=row.names(pca_data$x), pca_data$x)
	data <- merge(data, annotation, by.x="id", by.y="id")

	labels <- unique(subset(data, select=c("label","col_label")))
	labels[] <- lapply(labels, as.character) # convert to character

	if(!noAlpha){
		alpha_val <- 0.3
	}else{
		alpha_val <- 1
	}

	plot3d(data[,PCx],data[,PCy],data[,PCz], xlab=PCx, ylab=PCy, zlab=PCz, col=data$col_label, alpha=alpha_val, size=0.6, type="s")
	legend3d("topright", legend=labels$label, pch=15, col=labels$col_label, cex=1.5)
}

gglegend <- function(gp){ 
	tmp <- ggplot_gtable(ggplot_build(gp)) 
	leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
	
	return(tmp$grobs[[leg]])
}

tsne.plot <- function(tsne_data, signature_data, annotation, excluded_motifs, with_legend) {
	data <- as.data.frame(tsne_data$Y)
	colnames(data) <- c("tsne1","tsne2")
	data$id <- rownames(signature_data)
	data <- merge(data, annotation, by.x="id", by.y="id")

	labels <- subset(data[!duplicated(data$label),], select=c("label","col_label"))
	labels[] <- lapply(labels, as.character) # convert to character
	color_palette <- labels$col_label
	names(color_palette) <- labels$label

	if(!is.na(excluded_motifs)){
		data <- subset(data, !grepl(excluded_motifs, id))
	}

	gp <- ggplot(data) +
		geom_point(aes(x=tsne1, y=tsne2, colour=as.character(label)), pch=46, alpha=1) +
		# geom_point(aes(x=tsne1, y=tsne2, colour=as.character(label)), pch=16, alpha=1) +
		scale_colour_manual(name="Colors:", values=color_palette) +
		labs(title="Projection using t-SNE") +
		labs(x="Projection vector X", y="Projection vector Y") +
		theme_bw() +
		theme(panel.border=element_blank(), panel.grid.major=element_blank()) +
		theme(panel.grid.minor=element_blank(), axis.line=element_line(colour="black")) +
		theme(legend.position="bottom") +
		# theme(legend.position="none") +
		# coord_cartesian(xlim=c(-13,-7), ylim=c(-5,5)) + 
		guides(colour=guide_legend(override.aes=list(shape=15, size=4)))

	if(with_legend==0){
		gp <- gp +
			theme(legend.position="none")
		print(gp)
	}else if(with_legend==1){
		print(gp)
	}else if(with_legend==2){
		gp_with_legend <- gp
		gp <- gp +
			theme(legend.position="none")
		print(gp)
		grid.newpage()
		gp_legend <- gglegend(gp_with_legend)
		grid.draw(gp_legend)
	}
}

tsne.3D.plot <- function(tsne_3d_data, signature_data, annotation, noAlpha=FALSE, PCx="tsne1", PCy="tsne2", PCz="tsne3"){
	data <- as.data.frame(tsne_3d_data$Y)
	colnames(data) <- paste0("tsne",seq(1,ncol(tsne_3d_data$Y)))
	data$id <- rownames(signature_data)
	data <- merge(data, annotation, by.x="id", by.y="id")

	if(!noAlpha){
		alpha_val <- 0.3
	}else{
		alpha_val <- 1
	}

	labels <- unique(subset(data, select=c("label","col_label")))
	labels[] <- lapply(labels, as.character) # convert to character

	plot3d(data[,PCx],data[,PCy],data[,PCz], xlab=PCx, ylab=PCy, zlab=PCz, col=data$col_label, alpha=alpha_val, size=0.6, type="s") #, size=7
	legend3d("topright", legend=labels$label, pch=15, col=labels$col_label, cex=1.5)
}

remove.excluded.motifs <- function(classification_data, excluded_motifs){
	motif_annotation <- attr(classification_data, "annotation_motif")
	row_to_remove <- grepl(excluded_motifs, motif_annotation$id)

	classification_data <- classification_data[!row_to_remove,]
	attr(classification_data, "annotation_motif") <- attr(classification_data, "annotation_motif")[!row_to_remove,]
	attr(classification_data, "annotation_mod") <- attr(classification_data, "annotation_mod")[!row_to_remove,]
	attr(classification_data, "annotation_dir") <- attr(classification_data, "annotation_dir")[!row_to_remove,]
	if("annotation_strain" %in% names(attributes(classification_data))){
		attr(classification_data, "annotation_strain") <- attr(classification_data, "annotation_strain")[!row_to_remove,]
	}

	return(classification_data)
}

prepare.classification.data <- function(methylation_signal, strain_id, motif_summary, genome, min_cov, filter_iso, iupac_nc, nbCPU, data_type="train"){
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 8 # Huge margin needed for offseted windows
				
	motifs_signature <- extract.motifs.signature(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU)
	if(data_type=="train"){
		motifs_signature <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,mean_diff,distance)) %>%
			spread(distance, mean_diff)
	}else if(data_type=="classify"){
		motifs_signature <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,mean_diff,distance)) %>%
			spread(distance, mean_diff)
	}
	rownames(motifs_signature) <- motifs_signature$id
	motifs_signature$id <- NULL
	motifs_signature <- motifs_signature[complete.cases(motifs_signature),] # TODO check if can be handle better

	return(motifs_signature)
}

prepare.classification.data.basic_group <- function(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, filter_iso, iupac_nc, nbCPU, data_type="train"){
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 8 # Huge margin needed for offseted windows
	ref_gene<- readDNAStringSet(genome)
	#ref_gene_rev<- readDNAStringSet(genome_rev)
	ref_gene_rev<-complement(ref_gene)
	ref_gene<-toString(ref_gene[[1]])
	ref_gene_rev<-toString(ref_gene_rev[[1]])
	motifs_signature <- extract.motifs.signature(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU)
	motifs_signature<-motifs_signature$methylation_at_motifs
	print(dim(motifs_signature))
	if(data_type=="train"){ #just id is different
		motifs_signature.nat <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,nat,distance)) %>%
		spread(distance, nat)

		motifs_signature.wga <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,wga,distance)) %>%
		spread(distance, wga)

		motifs_signature.mean_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_diff,distance)) %>%
		spread(distance, mean_diff)

		motifs_signature.sd_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,sd_diff,distance)) %>%
		spread(distance, sd_diff)

		motifs_signature.mean_nat <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_nat,distance)) %>%
		spread(distance,mean_nat)

		motifs_signature.mean_wga <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_wga,distance)) %>%
		spread(distance,mean_wga)

		motifs_signature.sd_nat <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,sd_nat,distance)) %>%
		spread(distance, sd_nat)

		motifs_signature.sd_wga <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,sd_wga,distance)) %>%
		spread(distance,sd_wga)

		motifs_signature.basic_group <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,pos_signal,dir,distance)) 

		for(i in 1:length(motifs_signature.basic_group$pos_signal)){
			if(motifs_signature.basic_group$dir[i]=="fwd"){
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene,motifs_signature.basic_group$pos_signal[i],motifs_signature.basic_group$pos_signal[i])
			}
			else{
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene_rev,motifs_signature.basic_group$pos_signal[i]+7,motifs_signature.basic_group$pos_signal[i]+7)
			} #lapply(motifs_signature.basic_group$basic_group,function(x) )
		}
		motifs_signature.basic_group <- motifs_signature.basic_group %>%
		dplyr::select(c(id,basic_group,distance)) %>%
		spread(distance,basic_group)

	}else if(data_type=="classify"){
		motifs_signature.basic_group <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,pos_signal,dir,distance))

		for(i in 1:length(motifs_signature.basic_group$pos_signal)){
			if(motifs_signature.basic_group$dir[i]=="fwd"){
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene,motifs_signature.basic_group$pos_signal[i],motifs_signature.basic_group$pos_signal[i])
			}
			else{
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene_rev,motifs_signature.basic_group$pos_signal[i]+7,motifs_signature.basic_group$pos_signal[i]+7)
			}
		}	
		motifs_signature.basic_group <- motifs_signature.basic_group %>%
		dplyr::select(c(id,basic_group,distance)) %>%
		spread(distance,basic_group)

		motifs_signature.mean_nat <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_nat,distance)) %>%
		spread(distance, mean_nat)

		motifs_signature.mean_wga <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_wga,distance)) %>%
		spread(distance, mean_wga)

		motifs_signature.sd_nat <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,sd_nat,distance)) %>%
		spread(distance,sd_nat)

		motifs_signature.sd_wga <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,sd_wga,distance)) %>%
		spread(distance,sd_wga)

		motifs_signature.mean_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_diff,distance)) %>%
		spread(distance, mean_diff)

		motifs_signature.sd_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,sd_diff,distance)) %>%
		spread(distance, sd_diff)

		motifs_signature.nat<- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,nat,distance)) %>%
		spread(distance, nat)

		motifs_signature.wga <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,wga,distance)) %>%
		spread(distance, wga)
	}
	rownames(motifs_signature.mean_diff) <- motifs_signature.mean_diff$id
	motifs_signature.mean_diff$id <- NULL
	motifs_signature.complete.cases<-complete.cases(motifs_signature.mean_diff)
	motifs_signature.mean_diff <- motifs_signature.mean_diff[motifs_signature.complete.cases,] # TODO check if can be handle better

	rownames(motifs_signature.nat) <- motifs_signature.nat$id
	motifs_signature.nat$id <- NULL
	motifs_signature.nat <- motifs_signature.nat[motifs_signature.complete.cases,] # TODO check if can be handle better

	rownames(motifs_signature.wga) <- motifs_signature.wga$id
	motifs_signature.wga$id <- NULL
	motifs_signature.wga <- motifs_signature.wga[motifs_signature.complete.cases,] # TODO check if can be handle better

	rownames(motifs_signature.basic_group) <- motifs_signature.basic_group$id
	motifs_signature.basic_group$id <- NULL
	motifs_signature.basic_group <- motifs_signature.basic_group[motifs_signature.complete.cases,] # TODO check if can be handle better

	rownames(motifs_signature.sd_diff) <- motifs_signature.sd_diff$id
	motifs_signature.sd_diff$id <- NULL
	motifs_signature.sd_diff <- motifs_signature.sd_diff[motifs_signature.complete.cases,]

	rownames(motifs_signature.mean_nat) <- motifs_signature.mean_nat$id
	motifs_signature.mean_nat$id <- NULL
	motifs_signature.mean_nat <- motifs_signature.mean_nat[motifs_signature.complete.cases,]
	rownames(motifs_signature.mean_wga) <- motifs_signature.mean_wga$id
	motifs_signature.mean_wga$id <- NULL
	motifs_signature.mean_wga <- motifs_signature.mean_wga[motifs_signature.complete.cases,]
	rownames(motifs_signature.sd_nat) <- motifs_signature.sd_nat$id
	motifs_signature.sd_nat$id <- NULL
	motifs_signature.sd_nat <- motifs_signature.sd_nat[motifs_signature.complete.cases,]
	rownames(motifs_signature.sd_wga) <- motifs_signature.sd_wga$id
	motifs_signature.sd_wga$id <- NULL
	motifs_signature.sd_wga <- motifs_signature.sd_wga[motifs_signature.complete.cases,]
	return(list(motifs_signature.mean_diff=motifs_signature.mean_diff,motifs_signature.sd_diff=motifs_signature.sd_diff,motifs_signature.basic_group=motifs_signature.basic_group,motifs_signature.nat=motifs_signature.nat,motifs_signature.wga=motifs_signature.wga,
	motifs_signature.mean_nat=motifs_signature.mean_nat,motifs_signature.mean_wga=motifs_signature.mean_wga,motifs_signature.sd_nat=motifs_signature.sd_nat,motifs_signature.sd_wga=motifs_signature.sd_wga))
}
basic_function <- function(x,ref_gene) {
  return(substr(ref_gene,x,x))
}
prepare.classification.data.basic_group.time <- function(methylation_signal,basecall, strain_id, motif_summary, genome,genome_rev, min_cov, filter_iso, iupac_nc, nbCPU, data_type="train"){
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 8 # Huge margin needed for offseted windows
	ref_gene<- readDNAStringSet(genome)
	#ref_gene_rev<- readDNAStringSet(genome_rev)
	ref_gene_rev<-complement(ref_gene)
	ref_gene<-toString(ref_gene[[1]])
	ref_gene_rev<-toString(ref_gene_rev[[1]])
	motifs_signature <- extract.motifs.signature(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU)
	motifs_signature<-motifs_signature$methylation_at_motifs
	print("bbbbbbb")
	if(data_type=="train"){ #just id is different
		motifs_signature.t_test_pval <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,t_test_pval,distance)) %>%
		spread(distance,t_test_pval)

		motifs_signature.u_test_pval <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,u_test_pval,distance)) %>%
		spread(distance, u_test_pval)

		motifs_signature.mean_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_diff,distance)) %>%
		spread(distance, mean_diff)

		motifs_signature.std_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,std_diff,distance)) %>%
		spread(distance,std_diff)

		motifs_signature.time_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,time_diff,distance)) %>%
		spread(distance,time_diff)

		motifs_signature.basic_group <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			#mutate(basic_group=ifelse(dir=="fwd",substr(ref_gene,pos_signal,pos_signal),substr(ref_gene_rev,pos_signal,pos_signal))) %>%
			dplyr::select(c(id,pos_signal,dir,distance)) 
			#dplyr::select(c(id,basic_group,distance)) %>%
			#spread(distance,basic_group)

		motifs_signature.basic_group$basic_group<-sapply(motifs_signature.basic_group$pos_signal,basic_function,ref_gene=ref_gene)
		motifs_signature.basic_group$basic_group[which(motifs_signature.basic_group$dir=="rev")]<-sapply(motifs_signature.basic_group$pos_signal[which(motifs_signature.basic_group$dir=="rev")]+7,basic_function,ref_gene=ref_gene_rev)
		#motifs_signature.basic_group$basic_group<-lapply(motifs_signature.basic_group$pos_signal,function(x) ifelse(motifs_signature.basic_group$dir[x]=="fwd",substr(ref_gene,x,x),substr(ref_gene_rev,x,x)) )
		motifs_signature.basic_group <- motifs_signature.basic_group %>%
		dplyr::select(c(id,basic_group,distance)) %>%
		spread(distance,basic_group)

		basecall<-subset(basecall, select = -X)
		motifs_signature.merged_df <- merge(motifs_signature, basecall, by.x ="pos_signal" ,by.y = "position")

		motifs_signature.avg_basequality <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,avg_basequality,distance)) %>%
		spread(distance,avg_basequality)

		motifs_signature.base_num_fraction <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,base_num_fraction,distance)) %>%
		spread(distance,base_num_fraction)

		motifs_signature.g_num <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,g_num,distance)) %>%
		spread(distance,g_num)

		motifs_signature.g_value <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,g_value,distance)) %>%
		spread(distance,g_value)

	}else if(data_type=="classify"){
		basecall<-subset(basecall, select = -X)
		motifs_signature.merged_df <- merge(motifs_signature, basecall, by.x ="pos_signal" ,by.y = "position")

		motifs_signature.avg_basequality <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,avg_basequality,distance)) %>%
		spread(distance,avg_basequality)

		motifs_signature.base_num_fraction <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,base_num_fraction,distance)) %>%
		spread(distance,base_num_fraction)

		motifs_signature.g_num <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,g_num,distance)) %>%
		spread(distance,g_num)

		motifs_signature.g_value <- motifs_signature.merged_df %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,g_value,distance)) %>%
		spread(distance,g_value)
		
		motifs_signature.mean_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_diff,distance)) %>%
		spread(distance, mean_diff)

		motifs_signature.std_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,std_diff,distance)) %>%
		spread(distance, std_diff)

		motifs_signature.time_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,time_diff,distance)) %>%
		spread(distance, time_diff)

		motifs_signature.t_test_pval<- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,t_test_pval,distance)) %>%
		spread(distance,t_test_pval)

		motifs_signature.u_test_pval <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,u_test_pval,distance)) %>%
		spread(distance,u_test_pval)

		motifs_signature.basic_group <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,pos_signal,dir,distance))

		motifs_signature.basic_group$basic_group<-sapply(motifs_signature.basic_group$pos_signal,basic_function,ref_gene=ref_gene)
		motifs_signature.basic_group$basic_group[which(motifs_signature.basic_group$dir=="rev")]<-sapply(motifs_signature.basic_group$pos_signal[which(motifs_signature.basic_group$dir=="rev")]+7,basic_function,ref_gene=ref_gene_rev)
		#motifs_signature.basic_group$basic_group<-lapply(motifs_signature.basic_group$pos_signal,function(x) ifelse(motifs_signature.basic_group$dir[x]=="fwd",substr(ref_gene,x,x),substr(ref_gene_rev,x,x)) )
		motifs_signature.basic_group <- motifs_signature.basic_group %>%
		dplyr::select(c(id,basic_group,distance)) %>%
		spread(distance,basic_group)
	}
	rownames(motifs_signature.mean_diff) <- motifs_signature.mean_diff$id
	motifs_signature.mean_diff$id <- NULL
	motifs_signature.complete.cases<-complete.cases(motifs_signature.mean_diff)
	motifs_signature.mean_diff <- motifs_signature.mean_diff[motifs_signature.complete.cases,] # TODO check if can be handle better
	print(dim(motifs_signature.mean_diff))
	rownames(motifs_signature.t_test_pval) <- motifs_signature.t_test_pval$id
	motifs_signature.t_test_pval$id <- NULL
	motifs_signature.t_test_pval <- motifs_signature.t_test_pval[motifs_signature.complete.cases,] # TODO check if can be handle better

	rownames(motifs_signature.u_test_pval) <- motifs_signature.u_test_pval$id
	motifs_signature.u_test_pval$id <- NULL
	motifs_signature.u_test_pval <- motifs_signature.u_test_pval[motifs_signature.complete.cases,] # TODO check if can be handle better

	rownames(motifs_signature.basic_group) <- motifs_signature.basic_group$id
	motifs_signature.basic_group$id <- NULL
	motifs_signature.basic_group <- motifs_signature.basic_group[motifs_signature.complete.cases,] # TODO check if can be handle better

	rownames(motifs_signature.std_diff) <- motifs_signature.std_diff$id
	motifs_signature.std_diff$id <- NULL
	motifs_signature.std_diff <- motifs_signature.std_diff[motifs_signature.complete.cases,]

	rownames(motifs_signature.time_diff) <- motifs_signature.time_diff$id
	motifs_signature.time_diff$id <- NULL
	motifs_signature.time_diff <- motifs_signature.time_diff[motifs_signature.complete.cases,]

	rownames(motifs_signature.avg_basequality) <- motifs_signature.avg_basequality$id
	motifs_signature.avg_basequality$id <- NULL
	motifs_signature.avg_basequality <- motifs_signature.avg_basequality[motifs_signature.complete.cases,]

	rownames(motifs_signature.base_num_fraction) <- motifs_signature.base_num_fraction$id
	motifs_signature.base_num_fraction$id <- NULL
	motifs_signature.base_num_fraction <- motifs_signature.base_num_fraction[motifs_signature.complete.cases,]

	rownames(motifs_signature.g_num) <- motifs_signature.g_num$id
	motifs_signature.g_num$id <- NULL
	motifs_signature.g_num <- motifs_signature.g_num[motifs_signature.complete.cases,]

	rownames(motifs_signature.g_value) <- motifs_signature.g_value$id
	motifs_signature.g_value$id <- NULL
	motifs_signature.g_value <- motifs_signature.g_value[motifs_signature.complete.cases,]

	return(list(motifs_signature.mean_diff=motifs_signature.mean_diff,motifs_signature.std_diff=motifs_signature.std_diff,motifs_signature.time_diff=motifs_signature.time_diff,motifs_signature.basic_group=motifs_signature.basic_group,motifs_signature.t_test_pval=motifs_signature.t_test_pval,motifs_signature.u_test_pval=motifs_signature.u_test_pval,
	motifs_signature.avg_basequality=motifs_signature.avg_basequality,motifs_signature.base_num_fraction=motifs_signature.base_num_fraction,motifs_signature.g_value=motifs_signature.g_value,
	motifs_signature.g_num=motifs_signature.g_num))
}

prepare.classification.data.basic_group.long <- function(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, filter_iso, iupac_nc, nbCPU, data_type="train"){
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -512
	expected_signal_right <- 512
	signal_margin <- 16 # Huge margin needed for offseted windows
	ref_gene<- readDNAStringSet(genome)
	#ref_gene_rev<- readDNAStringSet(genome_rev)
	ref_gene_rev<-complement(ref_gene)
	ref_gene<-toString(ref_gene[[1]])
	ref_gene_rev<-toString(ref_gene_rev[[1]])
	motifs_signature <- extract.motifs.signature(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU)
	motifs_signature<-motifs_signature$methylation_at_motifs

	if(data_type=="train"){ #just id is different
		motifs_signature.mean_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_diff,distance)) %>%
		spread(distance, mean_diff)

	}else if(data_type=="classify"){

		motifs_signature.mean_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_diff,distance)) %>%
		spread(distance, mean_diff)
	}
	rownames(motifs_signature.mean_diff) <- motifs_signature.mean_diff$id
	motifs_signature.mean_diff$id <- NULL
	motifs_signature.complete.cases<-complete.cases(motifs_signature.mean_diff)
	motifs_signature.mean_diff <- motifs_signature.mean_diff[motifs_signature.complete.cases,] # TODO check if can be handle better

	return(list(motifs_signature.mean_diff=motifs_signature.mean_diff))
}


prepare.classification.data.basic_group.array <- function(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, filter_iso, iupac_nc, nbCPU, data_type="train"){
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 8 # Huge margin needed for offseted windows
	ref_gene<- readDNAStringSet(genome)
	#ref_gene_rev<- readDNAStringSet(genome_rev)
	ref_gene_rev<-complement(ref_gene)
	ref_gene<-toString(ref_gene[[1]])
	ref_gene_rev<-toString(ref_gene_rev[[1]])
	motifs_signature <- extract.motifs.signature(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU)
	methylation_at_motifs_list<-motifs_signature$methylation_at_motifs_list
	motifs_signature<-motifs_signature$methylation_at_motifs
	motifs_signature<-motifs_signature %>%
			mutate(label=motif_summary$mod_type[match(motif,motif_summary$motif)])
	#print(sum(is.na(methylation_at_motifs_list$nat)))
	tem_dataframe_nat=as.data.frame(methylation_at_motifs_list$nat)
	dim_tem_dataframe_nat=dim(tem_dataframe_nat)
	nat=array(as.matrix(tem_dataframe_nat),c(dim_tem_dataframe_nat[1],22,dim_tem_dataframe_nat[2]%/%22))
	tem_dataframe_wga=as.data.frame(methylation_at_motifs_list$wga)
	dim_tem_dataframe_wga=dim(tem_dataframe_wga)
	wga=array(as.matrix(tem_dataframe_wga),c(dim_tem_dataframe_wga[1],22,dim_tem_dataframe_wga[2]%/%22))	
	if(data_type=="train"){ #just id is different
		dim_tem_dataframe_mean_diff=dim(motifs_signature)
		motifs_signature.mean_diff=array(motifs_signature$mean_diff,c(22,dim_tem_dataframe_mean_diff[1]%/%22))
		#sink("test/motifs_signature.mean_diff_digit.txt")
		#print(dim(motifs_signature.mean_diff))
		#print(motifs_signature.mean_diff)
		#sink()
		#print("lll")
		motifs_signature.label=array(as.matrix(motifs_signature$label),c(22,dim_tem_dataframe_mean_diff[1]%/%22))
		dim_tem_dataframe_sd_diff=dim(motifs_signature)
		motifs_signature.sd_diff=array(motifs_signature$sd_diff,c(22,dim_tem_dataframe_sd_diff[1]%/%22))
		motifs_signature.basic_group <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			#mutate(basic_group=substr(ref_gene,pos_signal,pos_signal)) %>%
			dplyr::select(c(id,pos_signal,dir,distance)) #%>%
			#spread(distance,basic_group)
		#motifs_signature$basic_group <- substr(ref_gene,as.integer(motifs_signature$pos_signal),as.integer(motifs_signature$pos_signal))
		for(i in 1:length(motifs_signature.basic_group$pos_signal)){
			if(motifs_signature.basic_group$dir[i]=="fwd"){
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene,motifs_signature.basic_group$pos_signal[i],motifs_signature.basic_group$pos_signal[i])
			}
			else{
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene_rev,motifs_signature.basic_group$pos_signal[i]+7,motifs_signature.basic_group$pos_signal[i]+7)
			} #lapply(motifs_signature.basic_group$basic_group,function(x) )
		}
		dim_tem_dataframe_basic_group=dim(motifs_signature.basic_group)
		motifs_signature.basic_group=array(motifs_signature.basic_group$basic_group,c(22,dim_tem_dataframe_basic_group[1]%/%22))

	}else if(data_type=="classify"){
		motifs_signature.basic_group <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,pos_signal,dir,distance))

		for(i in 1:length(motifs_signature.basic_group$pos_signal)){
			if(motifs_signature.basic_group$dir[i]=="fwd"){
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene,motifs_signature.basic_group$pos_signal[i],motifs_signature.basic_group$pos_signal[i])
			}
			else{
				motifs_signature.basic_group$basic_group[i]<-substr(ref_gene_rev,motifs_signature.basic_group$pos_signal[i]+7,motifs_signature.basic_group$pos_signal[i]+7)
			}
		}	
		motifs_signature.basic_group <- motifs_signature.basic_group %>%
		dplyr::select(c(id,basic_group,distance)) %>%
		spread(distance,basic_group)

		motifs_signature.mean_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,mean_diff,distance)) %>%
		spread(distance, mean_diff)

		motifs_signature.sd_diff <- motifs_signature %>%
		mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
		dplyr::select(c(id,sd_diff,distance)) %>%
		spread(distance, sd_diff)
	}
	motifs_signature.mean_diff_data_frame=as.data.frame(t(motifs_signature.mean_diff))
	motifs_signature.complete.cases<-complete.cases(motifs_signature.mean_diff_data_frame)

	motifs_signature.mean_diff <- motifs_signature.mean_diff[,motifs_signature.complete.cases] 
	motifs_signature.basic_group <- motifs_signature.basic_group[,motifs_signature.complete.cases] 
	motifs_signature.sd_diff <- motifs_signature.sd_diff[,motifs_signature.complete.cases] 
	motifs_signature.label <- motifs_signature.label[1,motifs_signature.complete.cases]
	nat<-nat[,,motifs_signature.complete.cases]
	wga<-wga[,,motifs_signature.complete.cases]

	#motifs_signature.label_data_frame=as.data.frame(motifs_signature.label)
	#colnames(motifs_signature.label_data_frame) <- c("label")
	#sample_num<- data.frame(label=NA,train_idx=NA)
	#max_sample_size <- min(table(motifs_signature.label_data_frame[,1])) # Find smallest category

	#unique_label=unique(motifs_signature.label_data_frame[,1])
	# filter 4mc 5mc 6ma to same number
	#result<-foreach(label_type=unique_label, .combine=rbind) %do% {
	#	subtype_training_data_label <- subset(motifs_signature.label_data_frame, label==label_type)
	#	rowname_subtype_training_data_label<-rownames(subtype_training_data_label)
	#	train_idx <- sample(rowname_subtype_training_data_label, size=max_sample_size, replace=FALSE) # Subsample every type to smallest category
    #    new_row <- data.frame(label=label_type,train_idx=train_idx)
    #    sample_num <- rbind(sample_num, new_row)
	#}
	#sample_num<-sample_num[-1,]

	#sample_num$train_idx <- as.integer(sample_num$train_idx) 
	#motifs_signature.label <- motifs_signature.label[sample_num$train_idx]
	#motifs_signature.mean_diff <- motifs_signature.mean_diff[,sample_num$train_idx]
	#motifs_signature.basic_group <- motifs_signature.basic_group[,sample_num$train_idx]
	#motifs_signature.sd_diff<-motifs_signature.sd_diff[,sample_num$train_idx]
	#nat<-nat[,,sample_num$train_idx]
	#wga<-wga[,,sample_num$train_idx]
	motifs_signature.label<-array(motifs_signature.label,c(1,length(motifs_signature.label)))
	return(list(motifs_signature.mean_diff=motifs_signature.mean_diff,motifs_signature.sd_diff=motifs_signature.sd_diff,motifs_signature.basic_group=motifs_signature.basic_group,motifs_signature.nat=nat,motifs_signature.wga=wga,motifs_signature.label=motifs_signature.label))
}

prepare.annotation.data <- function(motifs_signature, ann_type, motif_summary=NA){
	if(ann_type=="motif"){
		motif_summary$old_motif <- motif_summary$motif # TODO modify motif_summary
		motif_summary$motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,2]) # Motif
		annotation <- merge(annotation, subset(motif_summary, select=c(motif, col_motif)), by.x=c("label"), by.y=c("motif"))
		names(annotation)[names(annotation)=='col_motif'] <- 'col_label'
		annotation$col_label <- as.factor(annotation$col_label)
	}else if(ann_type=="mod_type"){
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,3]) # Modification types
		annotation <- merge(annotation, data.frame(mod_type=c("4mC","5mC","6mA"), col_label=c("#36D2A0","#444CF0","#FF5733")), by.x=c("label"), by.y=c("mod_type"))
		annotation$col_label <- as.factor(annotation$col_label)
	}else if(ann_type=="dir"){
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,6]) # Strand types
		annotation <- merge(annotation, data.frame(dir=c("fwd","rev"), col_label=c("#444CF0","#FF5733")), by.x=c("label"), by.y=c("dir"))
		annotation$col_label <- as.factor(annotation$col_label)
	}else if(ann_type=="classify"){
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,2]) # Motif
		annotation <- merge(annotation, subset(motif_summary, select=c(motif, col_motif)), by.x=c("label"), by.y=c("motif"))
		names(annotation)[names(annotation)=='col_motif'] <- 'col_label'
		annotation$col_label <- as.factor(annotation$col_label)
	}else if(ann_type=="train_label"){
		base_df <- factor(c("4mC_-10","4mC_-11","4mC_-12","4mC_-6","4mC_-7","4mC_-8","4mC_-9"
		,"5mC_-10","5mC_-11","5mC_-12","5mC_-6","5mC_-7","5mC_-8","5mC_-9"
		,"6mA_-10","6mA_-11","6mA_-12","6mA_-6","6mA_-7","6mA_-8","6mA_-9"))
		annotation <- data.frame(id=rownames(motifs_signature), motif=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,2]) # Motif
		annotation <- merge(annotation, subset(motif_summary, select=c(explicit_motif,train_label,col_motif)), by.x=c("motif"), by.y=c("explicit_motif"))
		names(annotation)[names(annotation)=='col_motif'] <- 'col_label'
		annotation$col_label <- as.factor(annotation$col_label)
		annotation$train_label <- factor(annotation$train_label,levels=levels(base_df))
	}
	annotation <- annotation[match(rownames(motifs_signature), annotation$id),]
	return(annotation)
}

prepare.meta.classification.data <- function(methylation_signal, strain_id, motif_summary, genome, min_cov, keepIsolated, iupac_nc, nbCPU){
	print("Prepare motifs signatures.")
	motifs_signature <- prepare.classification.data(methylation_signal, strain_id, motif_summary, genome, min_cov, keepIsolated, iupac_nc, nbCPU)
	not_dup <- !(duplicated(motifs_signature) | duplicated(motifs_signature, fromLast=TRUE))
	unique_motifs_signature <- motifs_signature[not_dup,] # Need to remove duplicate when same base in two motifs for tsne
	print("Prepare motifs signatures annotation.")
	attr(unique_motifs_signature, "annotation_motif") <- prepare.annotation.data(motifs_signature, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature, "annotation_mod") <- prepare.annotation.data(motifs_signature, "mod_type")
	attr(unique_motifs_signature, "annotation_dir") <- prepare.annotation.data(motifs_signature, "dir")
	return(unique_motifs_signature)
}

prepare.meta.classification.data.basic_group <- function(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU){
	print("Prepare motifs signatures.")
	motifs_signature <- prepare.classification.data.basic_group(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU)
	#not_dup <- !(duplicated(motifs_signature$motifs_signature) | duplicated(motifs_signature$motifs_signature, fromLast=TRUE))
	#unique_motifs_signature <- motifs_signature$motifs_signature[not_dup,] # Need to remove duplicate when same base in two motifs for tsne
	#unique_motifs_signature.basic_group <- motifs_signature$motifs_signature.basic_group[not_dup,]
	unique_motifs_signature.mean_diff<-motifs_signature$motifs_signature.mean_diff
	unique_motifs_signature.sd_diff<-motifs_signature$motifs_signature.sd_diff
	unique_motifs_signature.basic_group<-motifs_signature$motifs_signature.basic_group
	unique_motifs_signature.nat<-motifs_signature$motifs_signature.nat
	unique_motifs_signature.wga<-motifs_signature$motifs_signature.wga
	unique_motifs_signature.mean_nat<-motifs_signature$motifs_signature.mean_nat
	unique_motifs_signature.mean_wga<-motifs_signature$motifs_signature.mean_wga
	unique_motifs_signature.sd_nat<-motifs_signature$motifs_signature.sd_nat
	unique_motifs_signature.sd_wga<-motifs_signature$motifs_signature.sd_wga
	print("Prepare motifs signatures annotation.")
	attr(unique_motifs_signature.mean_diff, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.mean_diff, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.mean_diff, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.sd_diff, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.sd_diff, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.sd_diff, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.basic_group, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.basic_group, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.basic_group, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")
	
	attr(unique_motifs_signature.nat, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.nat, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.nat, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.wga, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.wga, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.wga, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.mean_nat, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.mean_nat, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.mean_nat, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.mean_wga, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.mean_wga, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.mean_wga, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.sd_nat, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.sd_nat, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.sd_nat, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.sd_wga, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.sd_wga, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.sd_wga, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")
	rm(motifs_signature)
	return(list(unique_motifs_signature.mean_diff=unique_motifs_signature.mean_diff,unique_motifs_signature.sd_diff=unique_motifs_signature.sd_diff,unique_motifs_signature.basic_group=unique_motifs_signature.basic_group
	,unique_motifs_signature.nat=unique_motifs_signature.nat,unique_motifs_signature.wga=unique_motifs_signature.wga,unique_motifs_signature.mean_nat=unique_motifs_signature.mean_nat,unique_motifs_signature.mean_wga=unique_motifs_signature.mean_wga
	,unique_motifs_signature.sd_nat=unique_motifs_signature.sd_nat,unique_motifs_signature.sd_wga=unique_motifs_signature.sd_wga))
}

prepare.meta.classification.data.basic_group.long <- function(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU){
	print("Prepare motifs signatures.")
	motifs_signature <- prepare.classification.data.basic_group.long(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU)
	#not_dup <- !(duplicated(motifs_signature$motifs_signature) | duplicated(motifs_signature$motifs_signature, fromLast=TRUE))
	#unique_motifs_signature <- motifs_signature$motifs_signature[not_dup,] # Need to remove duplicate when same base in two motifs for tsne
	#unique_motifs_signature.basic_group <- motifs_signature$motifs_signature.basic_group[not_dup,]
	unique_motifs_signature.mean_diff<-motifs_signature$motifs_signature.mean_diff
	print("Prepare motifs signatures annotation.")
	attr(unique_motifs_signature.mean_diff, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.mean_diff, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.mean_diff, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	rm(motifs_signature)
	return(list(unique_motifs_signature.mean_diff=unique_motifs_signature.mean_diff))
}

prepare.meta.classification.data.basic_group.time <- function(methylation_signal,basecall, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU){
	print("Prepare motifs signatures.")
	motifs_signature <- prepare.classification.data.basic_group.time(methylation_signal,basecall, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU)
	#not_dup <- !(duplicated(motifs_signature$motifs_signature) | duplicated(motifs_signature$motifs_signature, fromLast=TRUE))
	#unique_motifs_signature <- motifs_signature$motifs_signature[not_dup,] # Need to remove duplicate when same base in two motifs for tsne
	#unique_motifs_signature.basic_group <- motifs_signature$motifs_signature.basic_group[not_dup,]
	unique_motifs_signature.mean_diff<-motifs_signature$motifs_signature.mean_diff
	unique_motifs_signature.t_test_pval<-motifs_signature$motifs_signature.t_test_pval
	unique_motifs_signature.u_test_pval<-motifs_signature$motifs_signature.u_test_pval
	unique_motifs_signature.std_diff<-motifs_signature$motifs_signature.std_diff
	unique_motifs_signature.time_diff<-motifs_signature$motifs_signature.time_diff
	unique_motifs_signature.avg_basequality<-motifs_signature$motifs_signature.avg_basequality
	unique_motifs_signature.base_num_fraction<-motifs_signature$motifs_signature.base_num_fraction
	unique_motifs_signature.g_num<-motifs_signature$motifs_signature.g_num
	unique_motifs_signature.g_value<-motifs_signature$motifs_signature.g_value
	unique_motifs_signature.basic_group<-motifs_signature$motifs_signature.basic_group

	print("Prepare motifs signatures annotation.")
	attr(unique_motifs_signature.mean_diff, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.mean_diff, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.mean_diff, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.std_diff, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.std_diff, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.std_diff, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.time_diff, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.time_diff, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.time_diff, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.basic_group, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.basic_group, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.basic_group, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")
	
	attr(unique_motifs_signature.t_test_pval, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.t_test_pval, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.t_test_pval, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.u_test_pval, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.u_test_pval, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.u_test_pval, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.avg_basequality, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.avg_basequality, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.avg_basequality, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.base_num_fraction, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.base_num_fraction, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.base_num_fraction, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.g_num, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.g_num, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.g_num, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")

	attr(unique_motifs_signature.g_value, "annotation_motif") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "motif", motif_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.g_value, "annotation_mod") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "mod_type")
	attr(unique_motifs_signature.g_value, "annotation_dir") <- prepare.annotation.data(motifs_signature$motifs_signature.mean_diff, "dir")
	
	rm(motifs_signature)
	return(list(unique_motifs_signature.mean_diff=unique_motifs_signature.mean_diff,unique_motifs_signature.t_test_pval=unique_motifs_signature.t_test_pval,unique_motifs_signature.u_test_pval=unique_motifs_signature.u_test_pval
	,unique_motifs_signature.std_diff=unique_motifs_signature.std_diff,unique_motifs_signature.time_diff=unique_motifs_signature.time_diff,unique_motifs_signature.avg_basequality=unique_motifs_signature.avg_basequality
	,unique_motifs_signature.basic_group=unique_motifs_signature.basic_group,unique_motifs_signature.base_num_fraction=unique_motifs_signature.base_num_fraction,unique_motifs_signature.g_num=unique_motifs_signature.g_num,unique_motifs_signature.g_value=unique_motifs_signature.g_value))
}

prepare.meta.classification.data.basic_group.array <- function(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU){
	print("Prepare motifs signatures.")
	motifs_signature <- prepare.classification.data.basic_group.array(methylation_signal, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU)
	print("Prepare motifs signatures annotation.")
	return(list(unique_motifs_signature.mean_diff=motifs_signature$motifs_signature.mean_diff,unique_motifs_signature.sd_diff=motifs_signature$motifs_signature.sd_diff,unique_motifs_signature.basic_group=motifs_signature$motifs_signature.basic_group,
	unique_motifs_signature.nat=motifs_signature$motifs_signature.nat,unique_motifs_signature.wga=motifs_signature$motifs_signature.wga,unique_motifs_signature.label=motifs_signature$motifs_signature.label))
}

merge.classification.features <- function(list_classification_data, merged_classification_data){
	stifle <- foreach(feature=c("annotation_motif","annotation_mod","annotation_dir","annotation_strain")) %do% {
		attr(merged_classification_data, feature) <- foreach(idx_classification_data=seq(1,length(list_classification_data)), .combine=rbind) %do% {

			return(attr(list_classification_data[[idx_classification_data]], feature))
		}
	}

	return(merged_classification_data)
}

merge.classification.data <- function(list_classification_data, list_strain_name){
	list_classification_data <- foreach(idx_classification_data=seq(1,length(list_classification_data))) %do% {
		classification_data <- list_classification_data[[idx_classification_data]]
		# Create strain annotation
		attr(classification_data, "annotation_strain") <- data.frame(
			label=as.factor(names(list_strain_name[idx_classification_data])),
			id=rownames(classification_data),
			col_label=as.factor(list_strain_name[[idx_classification_data]])
		)

		return(classification_data)
	}
	merged_classification_data <- rbindlist(list_classification_data)
	merged_classification_data <- as.data.frame(merged_classification_data)
	rownames(merged_classification_data) <- foreach(idx_classification_data=seq(1,length(list_classification_data)), .combine=c) %do% {

		return(rownames(list_classification_data[[idx_classification_data]]))
	}
	merged_classification_data <- merge.classification.features(list_classification_data, merged_classification_data)

	return(merged_classification_data)
}

draw.signature.clusters <- function(tsne_data, classification_data, base_name, path_output, excluded_motifs=NA, with_legend=0){
	pdf(paste0(path_output,"tsne_motif_",base_name,"_v12.pdf"), width=10, height=10)
	tsne.plot(tsne_data, classification_data, attr(classification_data, "annotation_motif"), excluded_motifs, with_legend)
	dev.off()
	pdf(paste0(path_output,"tsne_mod_type_",base_name,"_v12.pdf"), width=10, height=10)
	tsne.plot(tsne_data, classification_data, attr(classification_data, "annotation_mod"), excluded_motifs, with_legend)
	dev.off()
	pdf(paste0(path_output,"tsne_dir_type_",base_name,"_v12.pdf"), width=10, height=10)
	tsne.plot(tsne_data, classification_data, attr(classification_data, "annotation_dir"), excluded_motifs, with_legend)
	dev.off()
	if("annotation_strain" %in% names(attributes(classification_data))){
		pdf(paste0(path_output,"tsne_strain_",base_name,"_v12.pdf"), width=10, height=10)
		tsne.plot(tsne_data, classification_data, attr(classification_data, "annotation_strain"), excluded_motifs, with_legend)
		dev.off()
	}
}

draw.signature.heatmap <- function(motifs_signature, motif_summary, selected_pos, hc_pos, base_name, nbCPU){
	heatmap_data <- subset(motifs_signature, select=as.character(selected_pos))
	infos <- do.call(rbind,strsplit(rownames(heatmap_data),"_"))
	heatmap_data$motif <- as.factor(infos[,2])
	heatmap_data$mod_type <- as.factor(infos[,3])
	heatmap_data$id <- as.factor(paste0(infos[,4],"_",infos[,5],"_",infos[,6])) # contig modified compare reference, 7 not used
	heatmap_data <- melt(heatmap_data, id=c("motif","mod_type","id"), variable_name="distance")
	colnames(heatmap_data)[which(names(heatmap_data)=="value")] <- "mean_diff"

	if(!dir.exists(paste0("tmp_",base_name))){
		dir.create(paste0("tmp_",base_name))
	}

	# Remove motifs hidden for clustering
	motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(heatmap_data$motif),]
	# gp_heights <- heatmap_data %>%
	# 	group_by(motif) %>%
	# 	summarize(height=6 + (n()/5000), .groups="drop_last")
	# motif_summary <- merge(motif_summary, gp_heights, by.x="explicit_motif", by.y="motif")

	registerDoMC(nbCPU)
	list_plots <- foreach(idx_selected_motif=seq(1,nrow(motif_summary)), .inorder=TRUE) %dopar% {
		motif <- motif_summary$motif[idx_selected_motif]
		mod_pos <- motif_summary$mod_pos[idx_selected_motif]
		mod_type <- motif_summary$mod_type[idx_selected_motif]
		selected_motif <- paste0(substr(motif,1,mod_pos-1), as.character(mod_type), substr(motif,mod_pos+1,nchar(motif)))

		# # Adapt height to motif frequency
		subset_heatmap_data <- subset(heatmap_data, motif==selected_motif)
		gp_height <- 6 + nrow(subset_heatmap_data)/5000

		# Hierarchical Clustering
		dendro_data <- reshape::cast(subset(subset_heatmap_data, distance %in% hc_pos), id~distance, value="mean_diff")
		# dd <- dist(scale(dendro_data), method="euclidean") # TODO I think scaling not necessary because difference are standardized
		dd <- dist(dendro_data, method="euclidean")
		hc <- hclust(dd, method="ward.D2")

		# Prepare dendrogram
		ddata <- dendro_data(hc, type="rectangle") # Long
		gd_plot <- ggplot(segment(ddata)) + # ddata$segments
			geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
			coord_flip() +
			scale_x_continuous(expand=c(0,0)) +
			scale_y_reverse(expand=c(0,0)) +
			theme_bw() + theme(panel.grid=element_blank(), panel.border=element_blank()) +
			theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

		# Reorder by "classes"
		subset_heatmap_data$id <- ordered(subset_heatmap_data$id, levels=hc$label[hc$order])

		# Adapt color scale to motif, keep 0 centered
		percentile <- ecdf(subset_heatmap_data$mean_diff)
		q0 <- percentile(0)
		qn <- quantile(subset_heatmap_data$mean_diff, c(0.01, q0, 0.99), na.rm = TRUE)
		qn_range <- rescale(c(qn, range(subset_heatmap_data$mean_diff)))
		values_fill_scale <- c(0, seq(qn_range[1], qn_range[2], length.out=49), seq(qn_range[2], qn_range[3], length.out=49), 1)
		myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

		# Prepare signature heatmap
		gr_plot <- ggplot(subset_heatmap_data) +
			geom_raster(aes(x=distance, y=id, fill=mean_diff)) +
			scale_fill_gradientn(colours=myPalette(100), values=values_fill_scale) +
			theme_bw() +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(panel.grid=element_blank(), panel.border=element_blank()) +
			labs(x="Relative distance to modified base", fill="Mean pA\ndifferences")

		# return(list(gd=gd_plot, gr=gr_plot))

		# Combine both
		pdf(paste0("tmp_",base_name,"/Signature_heatmap_",motif,"_",base_name,"_v1.pdf"), height=gp_height)
		ggarrange(gd_plot, gr_plot, ncol=2, newpage=FALSE, top=paste0("Signature heatmap for ",selected_motif," motif"))
		dev.off()

		return(NA)
	}
	registerDoSEQ()

	pdf_list <- paste0("tmp_",base_name,"/Signature_heatmap_",motif_summary$motif,"_",base_name,"_v1.pdf", collapse=" ")
	merging_cmd <- paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dColorConversionStrategy=/LeaveColorUnchanged -dEncodeColorImages=false -dEncodeGrayImages=false -dEncodeMonoImages=false -sOutputFile=","Signature_heatmaps_",base_name,"_v1.pdf"," ",pdf_list, collapse=" ")
	system(merging_cmd)

	# Remove temporary directory
	unlink(paste0("tmp_",base_name), recursive=TRUE)
}

draw.signature.heatmap.detailed <- function(motifs_signature, motif_summary, genome, selected_pos, hc_pos, base_name, nbCPU){
	heatmap_data <- subset(motifs_signature, select=as.character(selected_pos))
	infos <- do.call(rbind,strsplit(rownames(heatmap_data),"_"))
	heatmap_data$motif <- as.factor(infos[,2])
	heatmap_data$mod_type <- as.factor(infos[,3])
	heatmap_data$id <- as.factor(paste0(infos[,4],"_",infos[,5],"_",infos[,6])) # contig modified compare reference, 7 not used
	heatmap_data <- melt(heatmap_data, id=c("motif","mod_type","id"), variable_name="distance")
	colnames(heatmap_data)[which(names(heatmap_data)=="value")] <- "mean_diff"

	if(!dir.exists(paste0("tmp_",base_name))){
		dir.create(paste0("tmp_",base_name))
	}

	# Remove motifs hidden for clustering
	motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(heatmap_data$motif),]

	registerDoMC(nbCPU)
	list_plots <- foreach(idx_selected_motif=seq(1,nrow(motif_summary)), .inorder=TRUE) %dopar% {
		motif <- motif_summary$motif[idx_selected_motif]
		mod_pos <- motif_summary$mod_pos[idx_selected_motif]
		mod_type <- motif_summary$mod_type[idx_selected_motif]
		selected_motif <- paste0(substr(motif,1,mod_pos-1), as.character(mod_type), substr(motif,mod_pos+1,nchar(motif)))

		# # Adapt height to motif frequency
		subset_heatmap_data <- subset(heatmap_data, motif==selected_motif)
		gp_height <- 6 + nrow(subset_heatmap_data)/5000

		# Hierarchical Clustering
		dendro_data <- reshape::cast(subset(subset_heatmap_data, distance %in% hc_pos), id~distance, value="mean_diff")
		# dd <- dist(scale(dendro_data), method="euclidean") # TODO I think scaling not necessary because difference are standardized
		dd <- dist(dendro_data, method="euclidean")
		hc <- hclust(dd, method="ward.D2")

		# Labelling Sub-Clusters
		clustering_info <- pamk(dendro_data) # Gave different clustering groups than hclust

		# Analyzing Sub-Clusters Sequences
		clustering_id <- as.data.frame(as.vector(cutree(hc, k=clustering_info$nc)))
		colnames(clustering_id) <- "cluster_id"
		clustering_id$id <- names(cutree(hc, k=clustering_info$nc))
		rownames(clustering_id) <- NULL
		# clustering_id <- as.data.frame(clustering_info$pamobject$clustering)
		# colnames(clustering_id) <- "cluster_id"
		# clustering_id$id <- rownames(clustering_id)
		# rownames(clustering_id) <- NULL

		list_left_bound_size <- 10
		list_right_bound_size <- 10
		clustering_id <- clustering_id %>%
			separate(id, c("contig", "position", "strand"), "_", remove=FALSE, convert=TRUE) %>%
			mutate(left_bound=ifelse(strand=="fwd", position-list_left_bound_size, position-list_right_bound_size)) %>%
			mutate(right_bound=ifelse(strand=="fwd", position+list_right_bound_size, position+list_left_bound_size))
		clustering_id <- clustering_id[match(hc$label[hc$order], clustering_id$id),] # Reorder according to hclust

		g_seq <- readDNAStringSet(genome) 
		list_contig_name <- gsub("_","",str_split(names(g_seq)," ", simplify=TRUE)[,1])
		gr <- GRanges(
			seqnames=Rle(names(g_seq)[match(clustering_id$contig, list_contig_name)]),
			ranges=IRanges(start=clustering_id$left_bound, end=clustering_id$right_bound),
			strand=ifelse(clustering_id$strand=="fwd","+","-") # Is strand specific
		)
		sequences <- getSeq(g_seq, gr)
		names(sequences) <- clustering_id$id

		logo_data <- foreach(chosen_cluster_id=rev(unique(clustering_id$cluster_id))) %do% { # Put in reverse order
			initiation_seqs <- sequences[names(sequences) %in% subset(clustering_id, cluster_id==chosen_cluster_id)$id]
			
			return(as.character(initiation_seqs))
		}
		names(logo_data) <- paste0("Cluster ",seq(1,clustering_info$nc)) # Use natural order not rev(unique(clustering_id$cluster_id))
		gp_logo <- ggplot() +
			geom_logo(logo_data) +
			theme_logo() + 
			facet_wrap(~seq_group, ncol=1, scales='free_x') #, nrow=8
		gp_empty <- ggplot() +
			geom_blank() +
			theme_bw() +
			theme(panel.border=element_blank())

		# Prepare dendrogram
		ddata <- as.dendrogram(hc) %>%
			set("branches_lwd", 0.3) %>%
			set("branches_k_color", k=clustering_info$nc) %>%
			as.ggdend()
		gd_plot <- ggplot(ddata$segments) + # ddata$segments or segment(ddata)
			geom_segment(aes(x=x, y=y, xend=xend, yend=yend, col=col)) + #, lwd=lwd
			coord_flip() + # Vertical dendrogram
			scale_x_continuous(expand=c(0,0)) +
			scale_y_reverse(expand=c(0,0)) + # Reverse from left to right
			theme_bw() + theme(panel.grid=element_blank(), panel.border=element_blank()) +
			theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
			guides(col=FALSE)

		# Reorder by "classes"
		subset_heatmap_data$id <- ordered(subset_heatmap_data$id, levels=hc$label[hc$order])

		# Adapt color scale to motif, keep 0 centered
		percentile <- ecdf(subset_heatmap_data$mean_diff)
		q0 <- percentile(0)
		qn <- quantile(subset_heatmap_data$mean_diff, c(0.01, q0, 0.99), na.rm = TRUE)
		qn_range <- rescale(c(qn, range(subset_heatmap_data$mean_diff)))
		values_fill_scale <- c(0, seq(qn_range[1], qn_range[2], length.out=49), seq(qn_range[2], qn_range[3], length.out=49), 1)
		myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

		# Prepare signature heatmap
		gr_plot <- ggplot(subset_heatmap_data) +
			geom_raster(aes(x=distance, y=id, fill=mean_diff)) +
			# geom_point(data=subset(subset_heatmap_data, id %in% res1$id), aes(x=distance, y=id)) +
			scale_fill_gradientn(colours=myPalette(100), values=values_fill_scale) +
			theme_bw() +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(panel.grid=element_blank(), panel.border=element_blank()) +
			labs(x="Relative distance to modified base", fill="Mean pA\ndifferences")

		# Combine both
		pdf(file=NULL) # Avoid opening empty window
		dendro_cluster <- ggarrange(gd_plot, gr_plot, ncol=2, newpage=FALSE, top=paste0("Signature heatmap for ",selected_motif," motif"), widths=c(3,3))
		logo_cluster <- ggarrange(gp_logo, gp_empty, ncol=1, newpage=FALSE, top=paste0("Sub-clusters ",selected_motif," motif"), heights=c(clustering_info$nc, round(gp_height)-clustering_info$nc))
		dev.off()
		pdf(paste0("tmp_",base_name,"/Signature_heatmap_detailed_",motif,"_",base_name,"_v1.pdf"), height=gp_height, width=9)
		grid.arrange(dendro_cluster, logo_cluster, ncol=2, widths=c(6,6))
		dev.off()

		return(NA)
	}
	registerDoSEQ()

	pdf_list <- paste0("tmp_",base_name,"/Signature_heatmap_detailed_",motif_summary$motif,"_",base_name,"_v1.pdf", collapse=" ")
	merging_cmd <- paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dColorConversionStrategy=/LeaveColorUnchanged -dEncodeColorImages=false -dEncodeGrayImages=false -dEncodeMonoImages=false -sOutputFile=","Signature_heatmaps_detailed_",base_name,"_v1.pdf"," ",pdf_list, collapse=" ")
	system(merging_cmd)

	# Remove temporary directory
	unlink(paste0("tmp_",base_name), recursive=TRUE)
}

gg_color_hue <- function(n){
	hues <- seq(15, 375, length=n + 1)

	return(hcl(h=hues, l=65, c=100)[1:n])
}

draw.tsne.signature.heatmap.detailed <- function(motifs_signature, motif_summary, tsne_data, genome, selected_pos, hc_pos, base_name, nbCPU, clustering_method="kmeans"){
	tsne_data <- as.data.frame(tsne_data$Y)
	tsne_data <- cbind(tsne_data, attr(motifs_signature, "annotation_motif"))
	tsne_data <- tsne_data %>%
		separate(id, into=c("genome","motif","mod_type","contig_name","pos","strand","type"), sep="_") %>%
		unite(id, contig_name, pos, strand, sep="_")

	heatmap_data <- subset(motifs_signature, select=as.character(selected_pos))
	infos <- do.call(rbind,strsplit(rownames(heatmap_data),"_"))
	heatmap_data$motif <- as.factor(infos[,2])
	heatmap_data$mod_type <- as.factor(infos[,3])
	heatmap_data$id <- as.factor(paste0(infos[,4],"_",infos[,5],"_",infos[,6])) # contig modified compare reference, 7 not used
	heatmap_data <- melt(heatmap_data, id=c("motif","mod_type","id"), variable_name="distance")
	colnames(heatmap_data)[which(names(heatmap_data)=="value")] <- "mean_diff"

	if(!dir.exists(paste0("tmp_",base_name))){
		dir.create(paste0("tmp_",base_name))
	}

	# Remove motifs hidden for clustering
	motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(heatmap_data$motif),]

	list_plots <- foreach(idx_selected_motif=seq(1,nrow(motif_summary)), .inorder=TRUE) %do% {
		motif <- motif_summary$motif[idx_selected_motif]
		mod_pos <- motif_summary$mod_pos[idx_selected_motif]
		mod_type <- motif_summary$mod_type[idx_selected_motif]
		selected_motif <- paste0(substr(motif,1,mod_pos-1), as.character(mod_type), substr(motif,mod_pos+1,nchar(motif)))

		# # Adapt height to motif frequency
		subset_heatmap_data <- subset(heatmap_data, motif==selected_motif)
		gp_height <- 6 + nrow(subset_heatmap_data)/5000

		# Hierarchical Clustering
		dendro_data <- reshape::cast(subset(subset_heatmap_data, distance %in% hc_pos), id~distance, value="mean_diff")
		# dd <- dist(scale(dendro_data), method="euclidean") # TODO I think scaling not necessary because difference are standardized
		dd <- dist(dendro_data, method="euclidean")
		hc <- hclust(dd, method="ward.D2")

		# Labelling Sub-Clusters
		clustering_info <- pamk(dendro_data) # Gave different clustering groups than hclust

		# Analyzing Sub-Clusters Sequences
		clustering_id <- as.data.frame(as.vector(cutree(hc, k=clustering_info$nc)))
		colnames(clustering_id) <- "cluster_id"
		clustering_id$id <- names(cutree(hc, k=clustering_info$nc))
		rownames(clustering_id) <- NULL

		tsne_data_tmp <- subset(tsne_data, label==selected_motif)
		tsne_data_tmp <- merge(tsne_data_tmp, clustering_id, by="id")
		gp_tsne <- ggplot(tsne_data_tmp) +
			geom_point(aes(x=V1, y=V2, col=as.factor(cluster_id))) +
			labs(col="Clusters") +
			theme(legend.justification=c(1,0), legend.position=c(1,0))

		list_left_bound_size <- 10
		list_right_bound_size <- 10
		clustering_id <- clustering_id %>%
			separate(id, c("contig", "position", "strand"), "_", remove=FALSE, convert=TRUE) %>%
			mutate(left_bound=ifelse(strand=="fwd", position-list_left_bound_size, position-list_right_bound_size)) %>%
			mutate(right_bound=ifelse(strand=="fwd", position+list_right_bound_size, position+list_left_bound_size))
		clustering_id <- clustering_id[match(hc$label[hc$order], clustering_id$id),] # Reorder according to hclust

		g_seq <- readDNAStringSet(genome) 
		list_contig_name <- gsub("_","",str_split(names(g_seq)," ", simplify=TRUE)[,1])
		gr <- GRanges(
			seqnames=Rle(names(g_seq)[match(clustering_id$contig, list_contig_name)]),
			ranges=IRanges(start=clustering_id$left_bound, end=clustering_id$right_bound),
			strand=ifelse(clustering_id$strand=="fwd","+","-") # Is strand specific
		)
		sequences <- getSeq(g_seq, gr)
		names(sequences) <- clustering_id$id

		logo_data <- foreach(chosen_cluster_id=rev(unique(clustering_id$cluster_id))) %do% { # Put in reverse order | hc order?
			initiation_seqs <- sequences[names(sequences) %in% subset(clustering_id, cluster_id==chosen_cluster_id)$id]
			
			return(as.character(initiation_seqs))
		}
		names(logo_data) <- paste0("Cluster ",seq(1,clustering_info$nc)) # Use natural order not rev(unique(clustering_id$cluster_id))
		gp_logo <- ggplot() +
			geom_logo(logo_data) +
			theme_logo() + 
			facet_wrap(~seq_group, ncol=1, scales='free_x') #, nrow=8
		gp_empty <- ggplot() +
			geom_blank() +
			theme_bw() +
			theme(panel.border=element_blank())

		# Prepare dendrogram
		ddata <- as.dendrogram(hc) %>%
			set("branches_lwd", 0.3) %>%
			set("branches_k_color", k=clustering_info$nc) %>%
			as.ggdend()
		gd_plot <- ggplot(ddata$segments) + # ddata$segments or segment(ddata)
			geom_segment(aes(x=x, y=y, xend=xend, yend=yend, col=col)) + #, lwd=lwd
			coord_flip() + # Vertical dendrogram
			scale_x_continuous(expand=c(0,0)) +
			scale_y_reverse(expand=c(0,0)) + # Reverse from left to right
			theme_bw() + theme(panel.grid=element_blank(), panel.border=element_blank()) +
			theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
			guides(col=FALSE)

		# Reorder by "classes"
		subset_heatmap_data$id <- ordered(subset_heatmap_data$id, levels=hc$label[hc$order])

		# Adapt color scale to motif, keep 0 centered
		percentile <- ecdf(subset_heatmap_data$mean_diff)
		q0 <- percentile(0)
		qn <- quantile(subset_heatmap_data$mean_diff, c(0.01, q0, 0.99), na.rm = TRUE)
		qn_range <- rescale(c(qn, range(subset_heatmap_data$mean_diff)))
		values_fill_scale <- c(0, seq(qn_range[1], qn_range[2], length.out=49), seq(qn_range[2], qn_range[3], length.out=49), 1)
		myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

		# Prepare signature heatmap
		gr_plot <- ggplot(subset_heatmap_data) +
			geom_raster(aes(x=distance, y=id, fill=mean_diff)) +
			# geom_point(data=subset(subset_heatmap_data, id %in% res1$id), aes(x=distance, y=id)) +
			scale_fill_gradientn(colours=myPalette(100), values=values_fill_scale) +
			theme_bw() +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(panel.grid=element_blank(), panel.border=element_blank()) +
			labs(x="Relative distance to modified base", fill="Mean pA\ndifferences")

		# Combine both
		pdf(file=NULL) # Avoid opening empty window
		dendro_cluster <- ggarrange(gd_plot, gr_plot, ncol=2, newpage=FALSE, top=paste0("Signature heatmap for ",selected_motif," motif"), widths=c(3,3))
		logo_cluster <- ggarrange(gp_logo, gp_empty, ncol=1, newpage=FALSE, top=paste0("Sub-clusters ",selected_motif," motif"), heights=c(clustering_info$nc, round(gp_height)-clustering_info$nc))
		tsne_cluster <- ggarrange(gp_tsne, gp_empty, ncol=1, newpage=FALSE, top=paste0("Sub-clusters ",selected_motif," motif"), heights=c(min(5, round(gp_height)), round(gp_height)-min(5, round(gp_height))))
		dev.off()
		pdf(paste0("tmp_",base_name,"/Signature_heatmap_detailed_",motif,"_",base_name,"_v1.pdf"), height=gp_height, width=12)
		grid.arrange(dendro_cluster, logo_cluster, tsne_cluster, ncol=3, widths=c(6,6,9))
		dev.off()

		tsne_data_tmp <- subset(tsne_data, label==selected_motif)
		if(clustering_method=="kmeans"){
			# Elbow + manual
			# Determine number of clusters
			wss <- (nrow(tsne_data_tmp)-1)*sum(apply(subset(tsne_data_tmp, select=c("V1","V2")),2,var))
			for (i in 2:15) wss[i] <- sum(kmeans(subset(tsne_data_tmp, select=c("V1","V2")), centers=i)$withinss)
			dev.new()
			plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
			choice_k <- as.integer(readline(prompt="How many clusters? "))
			dev.off()
			fit <- kmeans(subset(tsne_data_tmp, select=c("V1","V2")), centers=choice_k, nstart=50, iter.max=500)
			tsne_data_tmp$cluster_id2 <- as.vector(fit$cluster)
		}else if(clustering_method=="Mclust"){
			# M clust
			# Worst than elbow + manual?
			library(mclust)
			tsne_data_Mclust <- Mclust(as.matrix(subset(tsne_data_tmp, select=c("V1","V2"))), G=1:15, modelNames=mclust.options("emModelNames"))
			dev.new()
			plot(tsne_data_Mclust, what="BIC")
			choice_k <- as.integer(readline(prompt="How many clusters? "))
			dev.off()
			fit <- Mclust(as.matrix(subset(tsne_data_tmp, select=c("V1","V2"))), G=choice_k, modelNames=mclust.options("emModelNames"))
			tsne_data_tmp$cluster_id2 <- as.vector(fit$classification)
		}

		gp_tsne <- ggplot(tsne_data_tmp) +
			geom_point(aes(x=V1, y=V2, col=as.factor(cluster_id2))) +
			labs(col="Clusters") +
			theme(legend.justification=c(1,0), legend.position=c(1,0))

		# Regenerate clean version clustering_id
		clustering_id <- as.data.frame(as.vector(cutree(hc, k=clustering_info$nc)))
		colnames(clustering_id) <- "cluster_id"
		clustering_id$id <- names(cutree(hc, k=clustering_info$nc))
		rownames(clustering_id) <- NULL
		clustering_id$cluster_id <- tsne_data_tmp$cluster_id2[match(clustering_id$id, tsne_data_tmp$id)]
		clustering_info$nc <- choice_k

		subset_heatmap_data$cluster_id <- paste0("Cluster ",clustering_id$cluster_id[match(subset_heatmap_data$id, clustering_id$id)])

		# Hierarchical Clustering
		list_cluster_color <- gg_color_hue(choice_k)
		list_hc <- foreach(current_cluster=seq(1,choice_k)) %do% {
			cluster_subset_heatmap_data <- subset(subset_heatmap_data, id %in% clustering_id$id[clustering_id$cluster_id==current_cluster])

			cluster_dendro_data <- reshape::cast(subset(cluster_subset_heatmap_data, distance %in% hc_pos), id~distance, value="mean_diff")
			# cluster_dd <- dist(scale(cluster_dendro_data), method="euclidean") # TODO I think scaling not necessary because difference are standardized
			cluster_dd <- dist(cluster_dendro_data, method="euclidean")
			cluster_hc <- hclust(cluster_dd, method="ward.D2")

			return(cluster_hc)
		}

		ordered_id <- unlist(sapply(list_hc,function(x) return(x$label[x$order])))

		# Reorder by "classes"
		subset_heatmap_data$id <- ordered(subset_heatmap_data$id, levels=ordered_id)

		# Adapt color scale to motif, keep 0 centered
		percentile <- ecdf(subset_heatmap_data$mean_diff)
		q0 <- percentile(0)
		qn <- quantile(subset_heatmap_data$mean_diff, c(0.01, q0, 0.99), na.rm = TRUE)
		qn_range <- rescale(c(qn, range(subset_heatmap_data$mean_diff)))
		values_fill_scale <- c(0, seq(qn_range[1], qn_range[2], length.out=49), seq(qn_range[2], qn_range[3], length.out=49), 1)
		myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

		# Prepare signature heatmap
		gr_plot <- ggplot(subset_heatmap_data) +
			geom_raster(aes(x=distance, y=id, fill=mean_diff)) +
			scale_fill_gradientn(colours=myPalette(100), values=values_fill_scale) +
			facet_grid(cluster_id~., scales="free_y", space="free_y") +
			theme_bw() +
			theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			theme(panel.grid=element_blank(), panel.border=element_blank()) +
			labs(x="Relative distance to modified base", fill="Mean pA\ndifferences")

		list_left_bound_size <- 10
		list_right_bound_size <- 10
		clustering_id <- clustering_id %>%
			separate(id, c("contig", "position", "strand"), "_", remove=FALSE, convert=TRUE) %>%
			mutate(left_bound=ifelse(strand=="fwd", position-list_left_bound_size, position-list_right_bound_size)) %>%
			mutate(right_bound=ifelse(strand=="fwd", position+list_right_bound_size, position+list_left_bound_size))

		g_seq <- readDNAStringSet(genome) 
		list_contig_name <- gsub("_","",str_split(names(g_seq)," ", simplify=TRUE)[,1])
		gr <- GRanges(
			seqnames=Rle(names(g_seq)[match(clustering_id$contig, list_contig_name)]),
			ranges=IRanges(start=clustering_id$left_bound, end=clustering_id$right_bound),
			strand=ifelse(clustering_id$strand=="fwd","+","-") # Is strand specific
		)
		sequences <- getSeq(g_seq, gr)
		names(sequences) <- clustering_id$id

		logo_data <- foreach(chosen_cluster_id=seq(1,clustering_info$nc)) %do% {
			initiation_seqs <- sequences[names(sequences) %in% subset(clustering_id, cluster_id==chosen_cluster_id)$id]
			
			return(as.character(initiation_seqs))
		}
		names(logo_data) <- paste0("Cluster ",seq(1,clustering_info$nc)) # Same order as in foreach
		gp_logo <- ggplot() +
			geom_logo(logo_data) +
			theme_logo() + 
			facet_wrap(~seq_group, ncol=1, scales='free_x') #, nrow=8
		gp_empty <- ggplot() +
			geom_blank() +
			theme_bw() +
			theme(panel.border=element_blank())

		# Combine both
		pdf(file=NULL) # Avoid opening empty window
		dendro_cluster <- ggarrange(gr_plot, ncol=1, newpage=FALSE, top=paste0("Signature heatmap for ",selected_motif," motif"), widths=c(6))
		logo_cluster <- ggarrange(gp_logo, gp_empty, ncol=1, newpage=FALSE, top=paste0("Sub-clusters ",selected_motif," motif"), heights=c(clustering_info$nc, round(gp_height)-clustering_info$nc))
		tsne_cluster <- ggarrange(gp_tsne, gp_empty, ncol=1, newpage=FALSE, top=paste0("Sub-clusters ",selected_motif," motif"), heights=c(min(5, round(gp_height)), round(gp_height)-min(5, round(gp_height))))
		dev.off()
		pdf(paste0("tmp_",base_name,"/Signature_heatmap_detailed2_",motif,"_",base_name,"_v1.pdf"), height=gp_height, width=12)
		grid.arrange(dendro_cluster, logo_cluster, tsne_cluster, ncol=3, widths=c(6,6,9))
		dev.off()

		return(NA)
	}

	pdf_list <- paste0("tmp_",base_name,"/Signature_heatmap_detailed_",motif_summary$motif,"_",base_name,"_v1.pdf", collapse=" ")
	merging_cmd <- paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dColorConversionStrategy=/LeaveColorUnchanged -dEncodeColorImages=false -dEncodeGrayImages=false -dEncodeMonoImages=false -sOutputFile=","Signature_heatmaps_detailed_",base_name,"_v1.pdf"," ",pdf_list, collapse=" ")
	system(merging_cmd)

	pdf_list <- paste0("tmp_",base_name,"/Signature_heatmap_detailed2_",motif_summary$motif,"_",base_name,"_v1.pdf", collapse=" ")
	merging_cmd <- paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dColorConversionStrategy=/LeaveColorUnchanged -dEncodeColorImages=false -dEncodeGrayImages=false -dEncodeMonoImages=false -sOutputFile=","Signature_heatmaps_detailed2_",base_name,"_v1.pdf"," ",pdf_list, collapse=" ")
	system(merging_cmd)

	# Remove temporary directory
	unlink(paste0("tmp_",base_name), recursive=TRUE)
}

######
## Classifier of DNA modification and motifs
######

# Linear Classification in R from http://machinelearningmastery.com/linear-classification-in-r/
logistic.regression <- function(training_data, testing_data){
	# Logistic Regression
	library(VGAM)
	fit <- vglm(label~., family=multinomial, data=training_data) # Warnings
	probabilities <- predict(fit, testing_data, type="response")
	predictions <- apply(probabilities, 1, which.max)
	predictions <- levels(training_data$label)[predictions]

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- factor(as.factor(predictions), levels=levels(training_data$label))
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

linear.discriminant.analysis <- function(training_data, testing_data){
	# Linear Discriminant Analysis
	library(MASS)
	fit <- lda(label~., data=training_data)
	predictions <- predict(fit, testing_data)$class

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

partial.least.square.discriminant.analysis <- function(training_data, testing_data){
	# Partial Least Squares Discriminant Analysis
	library(caret)
	fit <- plsda(subset(training_data, select=-c(label)), training_data$label, probMethod="Bayes")
	predictions <- predict(fit, subset(testing_data, select=-c(label)))

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

# Non-Linear Classification in R from http://machinelearningmastery.com/non-linear-classification-in-r/
mixture.discriminant.analysis <- function(training_data, testing_data){
	# Mixture Discriminant Analysis
	library(mda)
	# fit <- mda(label~., data=training_data) # Default
	fit <- mda(label~., data=training_data, subclasses=5, iter=50)
	predictions <- predict(fit, testing_data)

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

quadratic.discriminant.analysis <- function(training_data, testing_data){
	# Quadratic Discriminant Analysis
	library(MASS)
	fit <- qda(label~., data=training_data)
	predictions <- predict(fit, testing_data)$class

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

regularized.discriminant.analysis <- function(training_data, testing_data){
	# Regularized Discriminant Analysis
	library(klaR)
	# fit <- rda(subset(training_data, select=-c(label)), training_data$label, gamma=0.05, lambda=0.01) # Default used
	fit <- rda(subset(training_data, select=-c(label)), training_data$label, gamma=0.08, lambda=0.08)
	predictions <- predict(fit, subset(testing_data, select=-c(label)))$class

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- factor(predictions, levels=levels(training_data$label))
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

neural.network <- function(training_data, testing_data,motif_name){
	# Neural Network
	library(nnet)
	# fit <- nnet(label~., data=training_data, size=4, decay=0.0001, maxit=500, trace=FALSE) # Default used
	fit <- nnet(label~., data=training_data, size=250, decay=0.00001, maxit=500, MaxNWts=40000, trace=FALSE)
	predictions <- predict(fit, subset(testing_data, select=-c(label)), type="class")
	testing_data$prediction <- factor(predictions, levels=levels(training_data$label))
	saveRDS(fit, file = paste("my_train_models/loocv/",motif_name, "_model_balance.RDS", sep = ""))
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")
	return(list(results=res, model=fit))
}

neural.network.notest <- function(training_data, testing_data,motif_name){
	# Neural Network
	library(nnet)
	# fit <- nnet(label~., data=training_data, size=4, decay=0.0001, maxit=500, trace=FALSE) # Default used
	fit <- nnet(label~., data=training_data, size=250, decay=0.00001, maxit=500, MaxNWts=40000, trace=FALSE)
	#saveRDS(fit, file = "my_train_models/loocv/All46_model_balance.RDS")
	return(list(results=NA, model=fit))
}

neural.network.predict <- function(training_data, testing_data,motif_name,possible_mod_pos){
	# Neural Network
	library(nnet)
	path_output <- "my_predict/loocv/"
	#fit <- readRDS(file=paste("my_train_models/loocv/", "All46_model.RDS", sep = ""))
	fit <- readRDS(file=paste("my_train_models/loocv/",motif_name, "_model_balance.RDS", sep = ""))
	# fit <- nnet(label~., data=training_data, size=4, decay=0.0001, maxit=500, trace=FALSE) # Default used
	#fit <- nnet(label~., data=training_data, size=250, decay=0.00001, maxit=500, MaxNWts=40000, trace=FALSE)
	predictions <- predict(fit, subset(testing_data, select=-c(label)), type="class")
	# Tidy prediction & summarize for motif
	prediction_summary <- data.frame(table(predictions))
	colnames(prediction_summary) <- c("label","count")

	# Generate complete data.frame
	base_df <- expand.grid(mod=c("4mC","5mC","6mA"), pos=seq(-12,-6)) %>%
		mutate(label=paste0(mod,"_",pos)) %>%
		arrange(label) %>%
		dplyr::select(label)
	
	prediction_summary <- merge(prediction_summary, base_df, by=c("label"), all=TRUE) %>% # Add combination potentially missing
		mutate(count=ifelse(is.na(count),0,count)) %>% # Put 0 when missing
		mutate(motif=motif_name) %>%
		mutate(mod_type=str_split(label,"_",simplify=TRUE)[,1]) %>%
		mutate(mod_pos=possible_mod_pos - (as.integer(str_split(label,"_",simplify=TRUE)[,2]) - (-9))) %>% # Leave it 0-based
		mutate(score=(count/sum(count))*100)
	#draw.classification.results(prediction_summary, paste0(motif_name,"_predict"), path_output)
	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- factor(predictions, levels=levels(training_data$label))
	#saveRDS(fit, file = paste("my_train_models/loocv/",motif_name, "_model.RDS", sep = ""))
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")
	return(list(results=res, model=fit,prediction_summary=prediction_summary))
}

flexible.discriminant.analysis <- function(training_data, testing_data){
	# Flexible Discriminant Analysis
	library(mda)
	library(earth)
	# fit <- fda(label~., data=training_data) # Default used
	fit <- fda(label~., data=training_data, method=earth, degree=1, nprune=21) # No pruning
	predictions <- predict(fit, testing_data)

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

support.vector.machine <- function(training_data, testing_data){
	# Support Vector Machine
	library(kernlab)
	fit <- ksvm(label~., data=training_data)
	predictions <- predict(fit, testing_data, type="response")

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

# Maybe all feature on same scale
knearest.neighbors <- function(training_data, testing_data){
	# k-Nearest Neighbors
	library(caret)
	# fit <- knn3(subset(training_data, select=-c(label)), training_data$label, k=5) # Default used
	fit <- knn3(subset(training_data, select=-c(label)), training_data$label, k=10)
	predictions <- predict(fit, subset(testing_data, select=-c(label)), type="class")

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

naive.bayes <- function(training_data, testing_data){
	# Naive Bayes
	library(klaR)
	# fit <- NaiveBayes(label~., data=training_data) # Likely  fL=0, usekernel=FALSE, adjust=1
	fit <- NaiveBayes(label~., data=training_data, fL=0, usekernel=TRUE, adjust=1.00) # Tuned
	predictions <- predict(fit, subset(testing_data, select=-label))

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- factor(as.vector(predictions$class), levels=levels(testing_data$label))
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

random.forest <- function(training_data, testing_data){
	# Random Forest
	library(randomForest)
	
	#Error if original colnames
	colnames(training_data) <- c("label",paste0("pos",seq(1,ncol(training_data)-1)))
	colnames(testing_data) <- c("label",paste0("pos",seq(1,ncol(testing_data)-1)))

	# fit <- randomForest(label~., data=training_data) # Default used should be mtry=3 (floor(sqrt(12))), ntree=500
	fit <- randomForest(label~., data=training_data, mtry=4, ntree=500)
	predictions <- predict(fit, testing_data)

	# res <- table(predictions, testing_data$label)
	testing_data$prediction <- predictions
	res <- confusionMatrix(testing_data$prediction, testing_data$label, mode="everything")

	return(list(results=res, model=fit))
}

prepare.classifier.input <- function(classification_data, annotation){
	processed_classification_data <- classification_data
	processed_classification_data$id <- rownames(processed_classification_data)
	processed_classification_data <- merge(processed_classification_data, annotation, by.x="id", by.y="id")
	processed_classification_data$motif <- do.call(rbind,strsplit(processed_classification_data$id,"_"))[,2]
	processed_classification_data <- subset(processed_classification_data, select=-c(id, col_label))

	return(processed_classification_data)
}

generate.offsetted.dataset <- function(processed_classification_data, length_vector, selected_starts, background_type="multi"){
	offseted_classification_data <- foreach(vector_start=selected_starts, .combine=rbind) %do% { # -6 to start at expected signal start
		tmp_col <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		tmp <- subset(processed_classification_data, select=c("label","motif",tmp_col))
		colnames(tmp)[which(names(tmp) %in% tmp_col)] <- seq(1,length_vector)
		tmp$label <- as.factor(paste0(tmp$label,"_",vector_start))

		return(tmp)
	}

	if(background_type=="single"){
		offseted_classification_data <- subset(offseted_classification_data, !label %in% paste0("none_",c(-12,-11,-10,-8,-7,-6)))
		offseted_classification_data$label <- droplevels(offseted_classification_data$label)
	}else if(background_type=="multi"){
		# Keep everything
	}else{
		# No background used
	}

	return(offseted_classification_data)
}

generate.offsetted.dataset.22 <- function(processed_classification_data, length_vector, selected_starts, background_type="multi"){
	tmp_col <- as.character(seq(-14, 7))
	tmp <- subset(processed_classification_data, select=c("label","motif",tmp_col))
	colnames(tmp)[which(names(tmp) %in% tmp_col)] <- seq(1,22)
	tmp$label <- as.factor(paste0(tmp$label,"_","-14"))
	return(tmp)
}

generate.offsetted.dataset.array.label <- function(processed_classification_data, length_vector, selected_starts, background_type="multi"){
	offseted_classification_data <- foreach(vector_start=seq(-12+14+1,-6+14+1), .combine=cbind) %do% { # -6 to start at expected signal start
	#tmp_col <- seq(vector_start, vector_start + (12 - 1))
  	tmp<-apply(processed_classification_data,MARGIN=c(1,2), function(x) paste0(x,"_",vector_start-15))
	return(tmp)
	}
	return(offseted_classification_data)
}

generate.offsetted.dataset.array <- function(processed_classification_data, length_vector, selected_starts, background_type="multi"){
	offseted_classification_data <- foreach(vector_start=seq(-12+14+1,-6+14+1), .combine=cbind) %do% { # -6 to start at expected signal start
	tmp_col <- seq(vector_start, vector_start + (12 - 1))
  	#tmp<-apply(processed_classification_data[,tmp_col],MARGIN=c(1,2), function(x) paste0(x,"_",vector_start))
	return(processed_classification_data[tmp_col,])
	}
	return(offseted_classification_data)
}

generate.offsetted.dataset.array.electric <- function(processed_classification_data, length_vector, selected_starts, background_type="multi"){
	len_data=dim(processed_classification_data)
	offseted_classification_data<-array(c(processed_classification_data[,seq(3,14),],processed_classification_data[,seq(4,15),],processed_classification_data[,seq(5,16),],
	processed_classification_data[,seq(6,17),],processed_classification_data[,seq(7,18),],
	processed_classification_data[,seq(8,19),],processed_classification_data[,seq(9,20),]),c(128,12,len_data[3]*7))

	#offseted_classification_data <- foreach(vector_start=seq(-12+14+1,-6+14+1), .combine=cbind) %do% { # -6 to start at expected signal start
	#tmp_col <- seq(vector_start, vector_start + (12 - 1))
  	#tmp<-apply(processed_classification_data[,tmp_col],MARGIN=c(1,2), function(x) paste0(x,"_",vector_start))
	#return(processed_classification_data[,tmp_col,])
	#}
	return(offseted_classification_data)
}

remove.impossible.exclusion <- function(excluded_motifs, motif_summary){
	removed_all_modType_occ <- foreach(idx=seq(1,nrow(excluded_motifs)), .combine=rbind) %do% {
		res <- any(table(motif_summary$mod_type[match(excluded_motifs[idx,], table=motif_summary$explicit_motif)]) >= table(motif_summary$mod_type))
	
		return(res)
	}
	excluded_motifs <- as.matrix(excluded_motifs[!removed_all_modType_occ,]) # Keep structure with only one motif

	return(excluded_motifs)
}

balancing.classifier.classes <- function(training_data){
	# Avoid unbalanced classes
	max_sample_size <- min(table(training_data$label)) # Find smallest category
	balanced_training_data <- foreach(label_type=unique(training_data$label), .combine=rbind) %do% {
		subtype_training_data <- subset(training_data, label==label_type)
		train_idx <- sample(seq_len(nrow(subtype_training_data)), size=max_sample_size, replace=FALSE) # Subsample every type to smallest category
		
		return(subtype_training_data[train_idx,])
	}

	return(balanced_training_data)
}

balancing.classifier.classes.basic_group <- function(training_data.basic_group,training_data,training_data.sd_diff,training_data.nat,training_data.wga,training_data.mean_nat,training_data.mean_wga,training_data.sd_nat,training_data.sd_wga){
	# Avoid unbalanced classes
	sample_num<- data.frame(label=NA,train_idx=NA)
	max_sample_size <- min(table(training_data$label)) # Find smallest category
	unique_label=unique(training_data$label)
	balanced_training_data <- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data <- subset(training_data, label==label_type)
		train_idx <- sample(seq_len(nrow(subtype_training_data)), size=max_sample_size, replace=FALSE) # Subsample every type to smallest category
        new_row <- data.frame(label=label_type,train_idx=train_idx)
        sample_num <- rbind(sample_num, new_row)
		return(subtype_training_data[train_idx,])
	}

	sample_num<-sample_num[-1,]

	balanced_training_data.basic_group <- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.basic_group <- subset(training_data.basic_group, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.basic_group[train_idx,])
	}

	balanced_training_data.sd_diff<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.sd_diff <- subset(training_data.sd_diff, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.sd_diff[train_idx,])
	}

	balanced_training_data.nat<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.nat <- subset(training_data.nat, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.nat[train_idx,])
	}

	balanced_training_data.wga<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.wga <- subset(training_data.wga, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.wga[train_idx,])
	}

	balanced_training_data.mean_nat<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.mean_nat <- subset(training_data.mean_nat, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.mean_nat[train_idx,])
	}

	balanced_training_data.mean_wga<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.mean_wga <- subset(training_data.mean_wga, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.mean_wga[train_idx,])
	}

		balanced_training_data.sd_nat<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.sd_nat <- subset(training_data.sd_nat, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.sd_nat[train_idx,])
	}

	balanced_training_data.sd_wga<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.sd_wga <- subset(training_data.sd_wga, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.sd_wga[train_idx,])
	}
	return(list(balanced_training_data=balanced_training_data,balanced_training_data.sd_diff=balanced_training_data.sd_diff,balanced_training_data.basic_group=balanced_training_data.basic_group,balanced_training_data.nat=balanced_training_data.nat,balanced_training_data.wga=balanced_training_data.wga,
	balanced_training_data.mean_nat=balanced_training_data.mean_nat,balanced_training_data.mean_wga=balanced_training_data.mean_wga,balanced_training_data.sd_nat=balanced_training_data.sd_nat,balanced_training_data.sd_wga=balanced_training_data.sd_wga))
}

balancing.classifier.classes.basic_group.time <- function(training_data.basic_group,training_data,training_data.std_diff,training_data.time_diff,training_data.t_test_pval,training_data.u_test_pval,training_data.avg_basequality,training_data.base_num_fraction,training_data.g_num,training_data.g_value){
	# Avoid unbalanced classes
	sample_num<- data.frame(label=NA,train_idx=NA)
	max_sample_size <- max(table(training_data$label)) # Find smallest category min
	unique_label=unique(training_data$label)
	balanced_training_data <- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data <- subset(training_data, label==label_type)
		train_idx <- sample(seq_len(nrow(subtype_training_data)), size=max_sample_size, replace=TRUE ) #FALSE Subsample every type to smallest category
        new_row <- data.frame(label=label_type,train_idx=train_idx)
        sample_num <- rbind(sample_num, new_row)
		return(subtype_training_data[train_idx,])
	}

	sample_num<-sample_num[-1,]

	balanced_training_data.basic_group <- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.basic_group <- subset(training_data.basic_group, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.basic_group[train_idx,])
	}

	balanced_training_data.std_diff<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.std_diff <- subset(training_data.std_diff, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.std_diff[train_idx,])
	}

	balanced_training_data.time_diff<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.time_diff <- subset(training_data.time_diff, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.time_diff[train_idx,])
	}

	balanced_training_data.t_test_pval<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.t_test_pval <- subset(training_data.t_test_pval, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.t_test_pval[train_idx,])
	}

	balanced_training_data.u_test_pval<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.u_test_pval <- subset(training_data.u_test_pval, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.u_test_pval[train_idx,])
	}

	balanced_training_data.avg_basequality<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.avg_basequality <- subset(training_data.avg_basequality, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.avg_basequality[train_idx,])
	}

	balanced_training_data.base_num_fraction<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.base_num_fraction <- subset(training_data.base_num_fraction, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.base_num_fraction[train_idx,])
	}

		balanced_training_data.g_num<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.g_num <- subset(training_data.g_num, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.g_num[train_idx,])
	}

	balanced_training_data.g_value<- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data.g_value <- subset(training_data.g_value, label==label_type)
		train_idx <- subset(sample_num, label==label_type)$train_idx
		return(subtype_training_data.g_value[train_idx,])
	}
	return(list(balanced_training_data=balanced_training_data,balanced_training_data.std_diff=balanced_training_data.std_diff,balanced_training_data.time_diff=balanced_training_data.time_diff,balanced_training_data.basic_group=balanced_training_data.basic_group,balanced_training_data.t_test_pval=balanced_training_data.t_test_pval,balanced_training_data.u_test_pval=balanced_training_data.u_test_pval,
	balanced_training_data.avg_basequality=balanced_training_data.avg_basequality,balanced_training_data.base_num_fraction=balanced_training_data.base_num_fraction,balanced_training_data.g_num=balanced_training_data.g_num,balanced_training_data.g_value=balanced_training_data.g_value))
}

balancing.classifier.classes.basic_group.long <- function(training_data){
	# Avoid unbalanced classes
	sample_num<- data.frame(label=NA,train_idx=NA)
	max_sample_size <- min(table(training_data$label)) # Find smallest category
	unique_label=unique(training_data$label)
	balanced_training_data <- foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data <- subset(training_data, label==label_type)
		train_idx <- sample(seq_len(nrow(subtype_training_data)), size=max_sample_size, replace=FALSE) # Subsample every type to smallest category
        new_row <- data.frame(label=label_type,train_idx=train_idx)
        sample_num <- rbind(sample_num, new_row)
		return(subtype_training_data[train_idx,])
	}

	sample_num<-sample_num[-1,]

	return(list(balanced_training_data=balanced_training_data))
}

prepare.classifier.dataset <- function(processed_classification_data, excluded_motifs, idx, balancing, motif_center_summary){
	# Select motif for training
	training_data <- subset(processed_classification_data, ! motif %in% excluded_motifs[idx,])
	training_data$label <- droplevels(training_data$label)
	if(balancing==1){ # Avoid unbalanced classes
		training_data <- balancing.classifier.classes(training_data)
	}
	training_data <- subset(training_data, select=-c(motif))
	testing_data <- subset(processed_classification_data, motif %in% excluded_motifs[idx,]) # Keep every occurrences
	testing_data <- subset(testing_data, select=-c(motif))

	if(!any(is.na(motif_center_summary))){ # If provided, only test predicted position
		predicted_start <- motif_center_summary$distance_prediction[motif_center_summary$explicit_motif %in% excluded_motifs[idx,]] - 9
		testing_data <- testing_data %>%
			filter(as.integer(str_split(label,"_",simplify=TRUE)[,2])==predicted_start)
		# testing_data <- droplevels(testing_data)
	}
	return(list(training_data=training_data,testing_data=testing_data))
}

prepare.classifier.dataset.notest <- function(processed_classification_data, excluded_motifs, idx, balancing, motif_center_summary){
	# Select motif for training
	#training_data <- subset(processed_classification_data, ! motif %in% excluded_motifs[idx,])
	training_data <- processed_classification_data
	#training_data$label <- droplevels(training_data$label)
	if(balancing==1){ # Avoid unbalanced classes
		training_data <- balancing.classifier.classes(training_data)
	}
	training_data <- subset(training_data, select=-c(motif))
	#testing_data <- subset(processed_classification_data, motif %in% excluded_motifs[idx,]) # Keep every occurrences
	#testing_data <- subset(testing_data, select=-c(motif))

	if(!any(is.na(motif_center_summary))){ # If provided, only test predicted position
		predicted_start <- motif_center_summary$distance_prediction[motif_center_summary$explicit_motif %in% excluded_motifs[idx,]] - 9
		testing_data <- testing_data %>%
			filter(as.integer(str_split(label,"_",simplify=TRUE)[,2])==predicted_start)
		# testing_data <- droplevels(testing_data)
	}
	return(list(training_data=training_data,testing_data=NA))
}

prepare.classifier.dataset.notest.basic_group <- function(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.sd_diff,processed_classification_data.nat,processed_classification_data.wga,processed_classification_data.mean_nat,processed_classification_data.mean_wga,processed_classification_data.sd_nat,processed_classification_data.sd_wga, excluded_motifs, idx, balancing, motif_center_summary){
	# Select motif for training
	#training_data <- subset(processed_classification_data, ! motif %in% excluded_motifs[idx,])
	#training_data <- processed_classification_data
	#training_data.basic_group <- processed_classification_data.basic_group
	#training_data$label <- droplevels(training_data$label)
	if(balancing==1){ # Avoid unbalanced classes
		training_data <- balancing.classifier.classes.basic_group(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.sd_diff,processed_classification_data.nat,processed_classification_data.wga,processed_classification_data.mean_nat,processed_classification_data.mean_wga,processed_classification_data.sd_nat,processed_classification_data.sd_wga)
	}

	training_data_balanced<-training_data$balanced_training_data
	training_data_balanced.basic_group<-training_data$balanced_training_data.basic_group
	training_data_balanced.sd_diff<-training_data$balanced_training_data.sd_diff
	training_data_balanced.nat<-training_data$balanced_training_data.nat
	training_data_balanced.wga<-training_data$balanced_training_data.wga
	training_data_balanced.mean_nat<-training_data$balanced_training_data.mean_nat
	training_data_balanced.mean_wga<-training_data$balanced_training_data.mean_wga
	training_data_balanced.sd_nat<-training_data$balanced_training_data.sd_nat
	training_data_balanced.sd_wga<-training_data$balanced_training_data.sd_wga

	#training_data_balanced <- subset(training_data$balanced_training_data, select=-c(motif))
	#training_data_balanced.basic_group <- subset(training_data$balanced_training_data.basic_group, select=-c(motif))
	#training_data_balanced.sd_diff <- subset(training_data$balanced_training_data.sd_diff, select=-c(motif))
	#training_data_balanced.nat <- subset(training_data$balanced_training_data.nat, select=-c(motif))
	#training_data_balanced.wga <- subset(training_data$balanced_training_data.wga, select=-c(motif))
	#training_data_balanced.mean_nat <- subset(training_data$balanced_training_data.mean_nat, select=-c(motif))
	#training_data_balanced.mean_wga <- subset(training_data$balanced_training_data.mean_wga, select=-c(motif))
	#training_data_balanced.sd_nat <- subset(training_data$balanced_training_data.sd_nat, select=-c(motif))
	#training_data_balanced.sd_wga <- subset(training_data$balanced_training_data.sd_wga, select=-c(motif))

	#testing_data <- subset(processed_classification_data, motif %in% excluded_motifs[idx,]) # Keep every occurrences
	#testing_data <- subset(testing_data, select=-c(motif))

	if(!any(is.na(motif_center_summary))){ # If provided, only test predicted position
		predicted_start <- motif_center_summary$distance_prediction[motif_center_summary$explicit_motif %in% excluded_motifs[idx,]] - 9
		testing_data <- testing_data %>%
			filter(as.integer(str_split(label,"_",simplify=TRUE)[,2])==predicted_start)
		# testing_data <- droplevels(testing_data)
	}
	return(list(training_data=training_data_balanced,training_data.sd_diff=training_data_balanced.sd_diff,training_data.basic_group=training_data_balanced.basic_group,training_data.nat=training_data_balanced.nat,training_data.wga=training_data_balanced.wga,
	training_data.mean_nat=training_data_balanced.mean_nat,training_data.mean_wga=training_data_balanced.mean_wga,training_data.sd_nat=training_data_balanced.sd_nat,training_data.sd_wga=training_data_balanced.sd_wga))
}

prepare.classifier.dataset.notest.basic_group.time <- function(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.std_diff,processed_classification_data.time_diff,processed_classification_data.t_test_pval,processed_classification_data.u_test_pval,processed_classification_data.avg_basequality,processed_classification_data.base_num_fraction,processed_classification_data.g_num,processed_classification_data.g_value, excluded_motifs, idx, balancing, motif_center_summary){
	# Select motif for training
	#training_data <- subset(processed_classification_data, ! motif %in% excluded_motifs[idx,])
	#training_data <- processed_classification_data
	#training_data.basic_group <- processed_classification_data.basic_group
	#training_data$label <- droplevels(training_data$label)
	if(balancing==2){ # Avoid unbalanced classes
		training_data <- balancing.classifier.classes.basic_group.time(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.std_diff,processed_classification_data.time_diff,processed_classification_data.t_test_pval,processed_classification_data.u_test_pval,processed_classification_data.avg_basequality,processed_classification_data.base_num_fraction,processed_classification_data.g_num,processed_classification_data.g_value)
	}

	training_data_balanced<-processed_classification_data
	training_data_balanced.basic_group<-processed_classification_data.basic_group
	training_data_balanced.std_diff<-processed_classification_data.std_diff
	training_data_balanced.time_diff<-processed_classification_data.time_diff
	training_data_balanced.t_test_pval<-processed_classification_data.t_test_pval
	training_data_balanced.u_test_pval<-processed_classification_data.u_test_pval
	training_data_balanced.avg_basequality<-processed_classification_data.avg_basequality
	training_data_balanced.base_num_fraction<-processed_classification_data.base_num_fraction
	training_data_balanced.g_num<-processed_classification_data.g_num
	training_data_balanced.g_value<-processed_classification_data.g_value

	if(!any(is.na(motif_center_summary))){ # If provided, only test predicted position
		predicted_start <- motif_center_summary$distance_prediction[motif_center_summary$explicit_motif %in% excluded_motifs[idx,]] - 9
		testing_data <- testing_data %>%
			filter(as.integer(str_split(label,"_",simplify=TRUE)[,2])==predicted_start)
		# testing_data <- droplevels(testing_data)
	}
	return(list(training_data=training_data_balanced,training_data.std_diff=training_data_balanced.std_diff,training_data.time_diff=training_data_balanced.time_diff,training_data.basic_group=training_data_balanced.basic_group,training_data.t_test_pval=training_data_balanced.t_test_pval,training_data.u_test_pval=training_data_balanced.u_test_pval,
	training_data.avg_basequality=training_data_balanced.avg_basequality,training_data.base_num_fraction=training_data_balanced.base_num_fraction,training_data.g_num=training_data_balanced.g_num,training_data.g_value=training_data_balanced.g_value))
}

prepare.classifier.dataset.notest.basic_group.long <- function(processed_classification_data, excluded_motifs, idx, balancing, motif_center_summary){
	# Select motif for training
	#training_data <- subset(processed_classification_data, ! motif %in% excluded_motifs[idx,])
	#training_data <- processed_classification_data
	#training_data.basic_group <- processed_classification_data.basic_group
	#training_data$label <- droplevels(training_data$label)
	if(balancing==1){ # Avoid unbalanced classes
		training_data <- balancing.classifier.classes.basic_group.long(processed_classification_data)
	}

	training_data_balanced<-training_data$balanced_training_data

	if(!any(is.na(motif_center_summary))){ # If provided, only test predicted position
		predicted_start <- motif_center_summary$distance_prediction[motif_center_summary$explicit_motif %in% excluded_motifs[idx,]] - 9
		testing_data <- testing_data %>%
			filter(as.integer(str_split(label,"_",simplify=TRUE)[,2])==predicted_start)
		# testing_data <- droplevels(testing_data)
	}
	return(list(training_data=training_data_balanced))
}

prepare.classifier.dataset.notest.basic_group.array <- function(processed_classification_data,processed_classification_data.sd_diff,processed_classification_data.nat,processed_classification_data.wga, excluded_motifs, idx, balancing, motif_center_summary){
	# Select motif for training
	#training_data <- subset(processed_classification_data, ! motif %in% excluded_motifs[idx,])
	#training_data <- processed_classification_data
	#training_data.basic_group <- processed_classification_data.basic_group
	#training_data$label <- droplevels(training_data$label)
	if(balancing==1){ # Avoid unbalanced classes
		training_data <- balancing.classifier.classes.basic_group(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.sd_diff,processed_classification_data.nat,processed_classification_data.wga)
	}
	training_data_balanced <- subset(training_data$balanced_training_data, select=-c(motif))


	return(list(training_data=training_data_balanced,training_data.sd_diff=training_data_balanced.sd_diff,training_data.basic_group=training_data_balanced.basic_group))
}

evaluate.classifier <- function(training_data, testing_data, list_classifier,motif_name){
	classifier_results <- foreach(classifier_name=list_classifier) %do% {
		classifier <- get(classifier_name)
		tmp_classifier_data <- classifier(training_data, testing_data,motif_name)

		tmp_classifier_data$model <- NULL # TODO crash if return model because too large

		return(tmp_classifier_data)
	}

	return(classifier_results)
}
evaluate.classifier.predict <- function(training_data, testing_data, list_classifier,motif_name,motif_possible_mod_pos){
	classifier_name=list_classifier[1]
	classifier <- get(classifier_name)
	tmp_classifier_data <- classifier(training_data,testing_data,motif_name,motif_possible_mod_pos)
	tmp_classifier_data$model <- NULL # TODO crash if return model because too large
	return(tmp_classifier_data$prediction_summary)
}
tidy.classifier.results <- function(classifier_results, list_classifier, excluded_motifs, idx){
	classifier_summary <- NA
	classifier_summary <- foreach(classifier_result=classifier_results, .combine=rbind) %do%{
		tmp_summary <- data.frame(
			accuracy=classifier_result$results$overall[["Accuracy"]]
		)

		return(tmp_summary)
	}
	classifier_summary$classifier <- list_classifier
	classifier_summary$excluded_motifs <- as.factor(paste0(excluded_motifs[idx,], collapse="/"))

	classifier_details <- NA
	classifier_details <- foreach(idx_classifier=seq(1,length(list_classifier)), .combine=rbind) %do%{
		contingency_table <- as.data.frame(classifier_results[[idx_classifier]]$results$table)
		contingency_table$classifier <- list_classifier[idx_classifier]

		return(contingency_table)
	}
	classifier_details$excluded_motifs <- as.factor(paste0(excluded_motifs[idx,], collapse="/"))
	classifier_details$Prediction <- ordered(classifier_details$Prediction, levels=levels(classifier_details$Reference))
	classifier_details$Reference <- ordered(classifier_details$Reference, levels=levels(classifier_details$Reference))

	# classifier_results$model <- foreach(prediction=classifier_results) %do%{
	# 	return(prediction$model)
	# }
	# names(classifier_results$model) <- list_classifier

	return(list(summary=classifier_summary, details=classifier_details))
}

tidy.meta_classifier.data <- function(meta_classifier_data){
	print(paste0(Sys.time()," Tidy classifiers results."))
	meta_classifier_summary <- foreach(classifier_data=meta_classifier_data, .combine=rbind) %do%{
		classifier_data$summary$classifier <- as.factor(classifier_data$summary$classifier)

		return(classifier_data$summary)
	}
	meta_classifier_details <- foreach(classifier_data=meta_classifier_data, .combine=rbind) %do%{
		classifier_data$details$classifier <- as.factor(classifier_data$details$classifier)

		return(classifier_data$details)
	}

	return(list(summary=meta_classifier_summary, details=meta_classifier_details))
}

evaluate.classifiers.performance <- function(classification_data, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)

	processed_classification_data <- prepare.classifier.input(classification_data, attr(classification_data, annotation_type))
	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
	}	
	# Generate exclusion list
	motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences
	registerDoMC(16)
	meta_classifier_data <- foreach(idx=seq(1,nrow(excluded_motifs)), .options.multicore=list(preschedule=TRUE)) %dopar% { #nrow(excluded_motifs)
		# Split datsets in training/testing
		classifier_datasets <- prepare.classifier.dataset(processed_classification_data, excluded_motifs, idx, balancing, motif_center_summary)
		# Train and test classifier
		classifier_results <- evaluate.classifier(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier)
		# Tidy classifier data
		tmp_meta_classifier_data <- tidy.classifier.results(classifier_results, list_classifier, excluded_motifs, idx)

		return(tmp_meta_classifier_data)
	}
	registerDoSEQ()

	# Tidy meta classifier data
	performance_classifier <- tidy.meta_classifier.data(meta_classifier_data)
	performance_classifier$summary$classifier <- factor(performance_classifier$summary$classifier, levels=list_classifier)
	performance_classifier$details$classifier <- factor(performance_classifier$details$classifier, levels=list_classifier)
	
	print(paste0(Sys.time()," Done."))

	return(performance_classifier)
}

evaluate.classifiers.performance.LOOCV <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)

	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data, attr(BA_classification_data, annotation_type))
	BF_processed_classification_data <- prepare.classifier.input(BF_classification_data, attr(BF_classification_data, annotation_type))
	CP_processed_classification_data <- prepare.classifier.input(CP_classification_data, attr(CP_classification_data, annotation_type))
	EC_processed_classification_data <- prepare.classifier.input(EC_classification_data, attr(EC_classification_data, annotation_type))
	HP_processed_classification_data <- prepare.classifier.input(HP_classification_data, attr(HP_classification_data, annotation_type))
	MH_processed_classification_data <- prepare.classifier.input(MH_classification_data, attr(MH_classification_data, annotation_type))
	NG_processed_classification_data <- prepare.classifier.input(NG_classification_data, attr(NG_classification_data, annotation_type))

	processed_classification_data<-rbind(BA_processed_classification_data,BF_processed_classification_data,CP_processed_classification_data,EC_processed_classification_data,
	HP_processed_classification_data,MH_processed_classification_data,NG_processed_classification_data)
	motifs_possible_mod_pos<-c(2,2,3,0,1,4,4,4,5,1,2,1,2,0,5,2,3,2,4,2,2,1,1,2,3,2,4,2,4,2,3,5,3,0,2,4,6,3,2,1,2,2,4,4,2,2)
	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
		#processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
	}	
	# Generate exclusion list
	#下面的源代码没注释
	#motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	#classification_data_label<-rbind(attr(BA_classification_data, "annotation_motif")$label)
	#motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences

	#train ALL
	#classifier_datasets <- prepare.classifier.dataset.notest(processed_classification_data, excluded_motifs, 1, balancing, motif_center_summary)
	#meta_classifier_data <- evaluate.classifier(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,NA)
	#return(excluded_motifs)
	registerDoMC(nbCPU)
	meta_classifier_data <- foreach(idx=seq(1,nrow(excluded_motifs)), .options.multicore=list(preschedule=TRUE)) %dopar% { #nrow(excluded_motifs)
		# Split datsets in training/testing
		classifier_datasets <- prepare.classifier.dataset(processed_classification_data, excluded_motifs, idx, balancing, motif_center_summary)
		# Train and test classifier
		print(excluded_motifs[idx,])
		print(motifs_possible_mod_pos[idx])
		classifier_results <- evaluate.classifier(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,excluded_motifs[idx,])
		# Tidy classifier data
		tmp_meta_classifier_data <- tidy.classifier.results(classifier_results, list_classifier, excluded_motifs, idx)

		return(tmp_meta_classifier_data)
	}
	registerDoSEQ()
	# Tidy meta classifier data
	performance_classifier <- tidy.meta_classifier.data(meta_classifier_data)
	performance_classifier$summary$classifier <- factor(performance_classifier$summary$classifier, levels=list_classifier)
	performance_classifier$details$classifier <- factor(performance_classifier$details$classifier, levels=list_classifier)
	
	print(paste0(Sys.time()," Done."))
	return(performance_classifier)
}
evaluate.classifiers.performance.basic_group.ALL <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)
	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.mean_diff, attr(BA_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	BF_processed_classification_data <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.mean_diff, attr(BF_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	CP_processed_classification_data <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.mean_diff, attr(CP_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	EC_processed_classification_data <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.mean_diff, attr(EC_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	HP_processed_classification_data <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.mean_diff, attr(HP_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	MH_processed_classification_data <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.mean_diff, attr(MH_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	NG_processed_classification_data <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.mean_diff, attr(NG_classification_data$unique_motifs_signature.mean_diff, annotation_type))

	processed_classification_data<-rbind(BA_processed_classification_data,BF_processed_classification_data,CP_processed_classification_data,EC_processed_classification_data,
	HP_processed_classification_data,MH_processed_classification_data,NG_processed_classification_data)

	BA_processed_classification_data.sd_diff <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.sd_diff, attr(BA_classification_data$unique_motifs_signature.sd_diff, annotation_type))
	BF_processed_classification_data.sd_diff <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.sd_diff, attr(BF_classification_data$unique_motifs_signature.sd_diff, annotation_type))
	CP_processed_classification_data.sd_diff <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.sd_diff, attr(CP_classification_data$unique_motifs_signature.sd_diff, annotation_type))
	EC_processed_classification_data.sd_diff <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.sd_diff, attr(EC_classification_data$unique_motifs_signature.sd_diff, annotation_type))
	HP_processed_classification_data.sd_diff <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.sd_diff, attr(HP_classification_data$unique_motifs_signature.sd_diff, annotation_type))
	MH_processed_classification_data.sd_diff <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.sd_diff, attr(MH_classification_data$unique_motifs_signature.sd_diff, annotation_type))
	NG_processed_classification_data.sd_diff <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.sd_diff, attr(NG_classification_data$unique_motifs_signature.sd_diff, annotation_type))

	processed_classification_data.sd_diff<-rbind(BA_processed_classification_data.sd_diff,BF_processed_classification_data.sd_diff,CP_processed_classification_data.sd_diff,EC_processed_classification_data.sd_diff,
	HP_processed_classification_data.sd_diff,MH_processed_classification_data.sd_diff,NG_processed_classification_data.sd_diff)

	BA_processed_classification_data.basic_group <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.basic_group, attr(BA_classification_data$unique_motifs_signature.basic_group, annotation_type))
	BF_processed_classification_data.basic_group <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.basic_group, attr(BF_classification_data$unique_motifs_signature.basic_group, annotation_type))
	CP_processed_classification_data.basic_group <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.basic_group, attr(CP_classification_data$unique_motifs_signature.basic_group, annotation_type))
	EC_processed_classification_data.basic_group <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.basic_group, attr(EC_classification_data$unique_motifs_signature.basic_group, annotation_type))
	HP_processed_classification_data.basic_group <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.basic_group, attr(HP_classification_data$unique_motifs_signature.basic_group, annotation_type))
	MH_processed_classification_data.basic_group <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.basic_group, attr(MH_classification_data$unique_motifs_signature.basic_group, annotation_type))
	NG_processed_classification_data.basic_group <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.basic_group, attr(NG_classification_data$unique_motifs_signature.basic_group, annotation_type))


	processed_classification_data.basic_group<-rbind(BA_processed_classification_data.basic_group,BF_processed_classification_data.basic_group,CP_processed_classification_data.basic_group,EC_processed_classification_data.basic_group,
	HP_processed_classification_data.basic_group,MH_processed_classification_data.basic_group,NG_processed_classification_data.basic_group)

	BA_processed_classification_data.nat <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.nat, attr(BA_classification_data$unique_motifs_signature.nat, annotation_type))
	BF_processed_classification_data.nat <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.nat, attr(BF_classification_data$unique_motifs_signature.nat, annotation_type))
	CP_processed_classification_data.nat <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.nat, attr(CP_classification_data$unique_motifs_signature.nat, annotation_type))
	EC_processed_classification_data.nat <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.nat, attr(EC_classification_data$unique_motifs_signature.nat, annotation_type))
	HP_processed_classification_data.nat <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.nat, attr(HP_classification_data$unique_motifs_signature.nat, annotation_type))
	MH_processed_classification_data.nat <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.nat, attr(MH_classification_data$unique_motifs_signature.nat, annotation_type))
	NG_processed_classification_data.nat <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.nat, attr(NG_classification_data$unique_motifs_signature.nat, annotation_type))


	processed_classification_data.nat<-rbind(BA_processed_classification_data.nat,BF_processed_classification_data.nat,CP_processed_classification_data.nat,EC_processed_classification_data.nat,
	HP_processed_classification_data.nat,MH_processed_classification_data.nat,NG_processed_classification_data.nat)

	BA_processed_classification_data.wga <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.wga, attr(BA_classification_data$unique_motifs_signature.wga, annotation_type))
	BF_processed_classification_data.wga <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.wga, attr(BF_classification_data$unique_motifs_signature.wga, annotation_type))
	CP_processed_classification_data.wga <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.wga, attr(CP_classification_data$unique_motifs_signature.wga, annotation_type))
	EC_processed_classification_data.wga <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.wga, attr(EC_classification_data$unique_motifs_signature.wga, annotation_type))
	HP_processed_classification_data.wga <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.wga, attr(HP_classification_data$unique_motifs_signature.wga, annotation_type))
	MH_processed_classification_data.wga <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.wga, attr(MH_classification_data$unique_motifs_signature.wga, annotation_type))
	NG_processed_classification_data.wga <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.wga, attr(NG_classification_data$unique_motifs_signature.wga, annotation_type))


	processed_classification_data.wga<-rbind(BA_processed_classification_data.wga,BF_processed_classification_data.wga,CP_processed_classification_data.wga,EC_processed_classification_data.wga,
	HP_processed_classification_data.wga,MH_processed_classification_data.wga,NG_processed_classification_data.wga)

	BA_processed_classification_data.mean_nat <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.mean_nat, attr(BA_classification_data$unique_motifs_signature.mean_nat, annotation_type))
	BF_processed_classification_data.mean_nat <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.mean_nat, attr(BF_classification_data$unique_motifs_signature.mean_nat, annotation_type))
	CP_processed_classification_data.mean_nat <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.mean_nat, attr(CP_classification_data$unique_motifs_signature.mean_nat, annotation_type))
	EC_processed_classification_data.mean_nat <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.mean_nat, attr(EC_classification_data$unique_motifs_signature.mean_nat, annotation_type))
	HP_processed_classification_data.mean_nat <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.mean_nat, attr(HP_classification_data$unique_motifs_signature.mean_nat, annotation_type))
	MH_processed_classification_data.mean_nat <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.mean_nat, attr(MH_classification_data$unique_motifs_signature.mean_nat, annotation_type))
	NG_processed_classification_data.mean_nat <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.mean_nat, attr(NG_classification_data$unique_motifs_signature.mean_nat, annotation_type))


	processed_classification_data.mean_nat<-rbind(BA_processed_classification_data.mean_nat,BF_processed_classification_data.mean_nat,CP_processed_classification_data.mean_nat,EC_processed_classification_data.mean_nat,
	HP_processed_classification_data.mean_nat,MH_processed_classification_data.mean_nat,NG_processed_classification_data.mean_nat)

	BA_processed_classification_data.mean_wga <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.mean_wga, attr(BA_classification_data$unique_motifs_signature.mean_wga, annotation_type))
	BF_processed_classification_data.mean_wga <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.mean_wga, attr(BF_classification_data$unique_motifs_signature.mean_wga, annotation_type))
	CP_processed_classification_data.mean_wga <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.mean_wga, attr(CP_classification_data$unique_motifs_signature.mean_wga, annotation_type))
	EC_processed_classification_data.mean_wga <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.mean_wga, attr(EC_classification_data$unique_motifs_signature.mean_wga, annotation_type))
	HP_processed_classification_data.mean_wga <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.mean_wga, attr(HP_classification_data$unique_motifs_signature.mean_wga, annotation_type))
	MH_processed_classification_data.mean_wga <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.mean_wga, attr(MH_classification_data$unique_motifs_signature.mean_wga, annotation_type))
	NG_processed_classification_data.mean_wga <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.mean_wga, attr(NG_classification_data$unique_motifs_signature.mean_wga, annotation_type))


	processed_classification_data.mean_wga<-rbind(BA_processed_classification_data.mean_wga,BF_processed_classification_data.mean_wga,CP_processed_classification_data.mean_wga,EC_processed_classification_data.mean_wga,
	HP_processed_classification_data.mean_wga,MH_processed_classification_data.mean_wga,NG_processed_classification_data.mean_wga)

	BA_processed_classification_data.sd_nat <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.sd_nat, attr(BA_classification_data$unique_motifs_signature.sd_nat, annotation_type))
	BF_processed_classification_data.sd_nat <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.sd_nat, attr(BF_classification_data$unique_motifs_signature.sd_nat, annotation_type))
	CP_processed_classification_data.sd_nat <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.sd_nat, attr(CP_classification_data$unique_motifs_signature.sd_nat, annotation_type))
	EC_processed_classification_data.sd_nat <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.sd_nat, attr(EC_classification_data$unique_motifs_signature.sd_nat, annotation_type))
	HP_processed_classification_data.sd_nat <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.sd_nat, attr(HP_classification_data$unique_motifs_signature.sd_nat, annotation_type))
	MH_processed_classification_data.sd_nat <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.sd_nat, attr(MH_classification_data$unique_motifs_signature.sd_nat, annotation_type))
	NG_processed_classification_data.sd_nat <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.sd_nat, attr(NG_classification_data$unique_motifs_signature.sd_nat, annotation_type))


	processed_classification_data.sd_nat<-rbind(BA_processed_classification_data.sd_nat,BF_processed_classification_data.sd_nat,CP_processed_classification_data.sd_nat,EC_processed_classification_data.sd_nat,
	HP_processed_classification_data.sd_nat,MH_processed_classification_data.sd_nat,NG_processed_classification_data.sd_nat)

	BA_processed_classification_data.sd_wga <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.sd_wga, attr(BA_classification_data$unique_motifs_signature.sd_wga, annotation_type))
	BF_processed_classification_data.sd_wga <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.sd_wga, attr(BF_classification_data$unique_motifs_signature.sd_wga, annotation_type))
	CP_processed_classification_data.sd_wga <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.sd_wga, attr(CP_classification_data$unique_motifs_signature.sd_wga, annotation_type))
	EC_processed_classification_data.sd_wga <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.sd_wga, attr(EC_classification_data$unique_motifs_signature.sd_wga, annotation_type))
	HP_processed_classification_data.sd_wga <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.sd_wga, attr(HP_classification_data$unique_motifs_signature.sd_wga, annotation_type))
	MH_processed_classification_data.sd_wga <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.sd_wga, attr(MH_classification_data$unique_motifs_signature.sd_wga, annotation_type))
	NG_processed_classification_data.sd_wga <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.sd_wga, attr(NG_classification_data$unique_motifs_signature.sd_wga, annotation_type))


	processed_classification_data.sd_wga<-rbind(BA_processed_classification_data.sd_wga,BF_processed_classification_data.sd_wga,CP_processed_classification_data.sd_wga,EC_processed_classification_data.sd_wga,
	HP_processed_classification_data.sd_wga,MH_processed_classification_data.sd_wga,NG_processed_classification_data.sd_wga)
	motifs_possible_mod_pos<-c(2,2,3,0,1,4,4,4,5,1,2,1,2,0,5,2,3,2,4,2,2,1,1,2,3,2,4,2,4,2,3,5,3,0,2,4,6,3,2,1,2,2,4,4,2,2)
	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data.basic_group <- generate.offsetted.dataset(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset(processed_classification_data.wga, length_vector, selected_starts, background_type)
		processed_classification_data.mean_nat <- generate.offsetted.dataset(processed_classification_data.mean_nat, length_vector, selected_starts, background_type)
		processed_classification_data.mean_wga <- generate.offsetted.dataset(processed_classification_data.mean_wga, length_vector, selected_starts, background_type)
		processed_classification_data.sd_nat <- generate.offsetted.dataset(processed_classification_data.sd_nat, length_vector, selected_starts, background_type)
		processed_classification_data.sd_wga <- generate.offsetted.dataset(processed_classification_data.sd_wga, length_vector, selected_starts, background_type)
		#processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data.basic_group <- generate.offsetted.dataset(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset(processed_classification_data.wga, length_vector, selected_starts, background_type)
		processed_classification_data.mean_nat <- generate.offsetted.dataset(processed_classification_data.mean_nat, length_vector, selected_starts, background_type)
		processed_classification_data.mean_wga <- generate.offsetted.dataset(processed_classification_data.mean_wga, length_vector, selected_starts, background_type)
		processed_classification_data.sd_nat <- generate.offsetted.dataset(processed_classification_data.sd_nat, length_vector, selected_starts, background_type)
		processed_classification_data.sd_wga <- generate.offsetted.dataset(processed_classification_data.sd_wga, length_vector, selected_starts, background_type)
	}	
	# Generate exclusion list
	#下面的源代码没注释
	#motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	#classification_data_label<-rbind(attr(BA_classification_data, "annotation_motif")$label)
	#motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences

	#train ALL
	classifier_datasets <- prepare.classifier.dataset.notest.basic_group(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.sd_diff,processed_classification_data.nat,processed_classification_data.wga,processed_classification_data.mean_nat,processed_classification_data.mean_wga,processed_classification_data.sd_nat,processed_classification_data.sd_wga, excluded_motifs, 1, balancing, motif_center_summary)
	saveRDS(classifier_datasets$training_data, file = "my_train_models/all_46_basic_group/All46_train_data_balance_1_290.rds")
	saveRDS(classifier_datasets$training_data.sd_diff, file = "my_train_models/all_46_basic_group/All46_train_data_balance_sd_diff_1_290.rds")
	saveRDS(classifier_datasets$training_data.basic_group, file = "my_train_models/all_46_basic_group/All46_train_data_balance_basic_group_1_290.rds")
	saveRDS(classifier_datasets$training_data.nat, file = "my_train_models/all_46_basic_group/All46_train_data_balance_nat_1_290.rds")
	saveRDS(classifier_datasets$training_data.wga, file = "my_train_models/all_46_basic_group/All46_train_data_balance_wga_1_290.rds")
	saveRDS(classifier_datasets$training_data.mean_nat, file = "my_train_models/all_46_basic_group/All46_train_data_balance_mean_nat_1_290.rds")
	saveRDS(classifier_datasets$training_data.mean_wga, file = "my_train_models/all_46_basic_group/All46_train_data_balance_mean_wga_1_290.rds")
	saveRDS(classifier_datasets$training_data.sd_nat, file = "my_train_models/all_46_basic_group/All46_train_data_balance_sd_nat_1_290.rds")
	saveRDS(classifier_datasets$training_data.sd_wga, file = "my_train_models/all_46_basic_group/All46_train_data_balance_sd_wga_1_290.rds")
	#meta_classifier_data <- evaluate.classifier(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,NA)
	return(excluded_motifs)
}

evaluate.classifiers.performance.basic_group.ALL.long <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)
	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.mean_diff, attr(BA_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	BF_processed_classification_data <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.mean_diff, attr(BF_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	CP_processed_classification_data <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.mean_diff, attr(CP_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	EC_processed_classification_data <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.mean_diff, attr(EC_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	HP_processed_classification_data <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.mean_diff, attr(HP_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	MH_processed_classification_data <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.mean_diff, attr(MH_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	NG_processed_classification_data <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.mean_diff, attr(NG_classification_data$unique_motifs_signature.mean_diff, annotation_type))

	processed_classification_data<-rbind(BA_processed_classification_data,BF_processed_classification_data,CP_processed_classification_data,EC_processed_classification_data,
	HP_processed_classification_data,MH_processed_classification_data,NG_processed_classification_data)

	motifs_possible_mod_pos<-c(2,2,3,0,1,4,4,4,5,1,2,1,2,0,5,2,3,2,4,2,2,1,1,2,3,2,4,2,4,2,3,5,3,0,2,4,6,3,2,1,2,2,4,4,2,2)
	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
	}else{
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
	}	
	# Generate exclusion list
	#下面的源代码没注释
	#motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	#classification_data_label<-rbind(attr(BA_classification_data, "annotation_motif")$label)
	#motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences

	#train ALL
	classifier_datasets <- prepare.classifier.dataset.notest.basic_group.long(processed_classification_data, excluded_motifs, 1, balancing, motif_center_summary)
	saveRDS(classifier_datasets$training_data, file = "my_train_models/all_46_basic_group/All46_train_data_balance_long_1_291.rds")

	return(excluded_motifs)
}

evaluate.classifiers.performance.basic_group.ALL.time <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)
	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.mean_diff, attr(BA_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	BF_processed_classification_data <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.mean_diff, attr(BF_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	CP_processed_classification_data <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.mean_diff, attr(CP_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	EC_processed_classification_data <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.mean_diff, attr(EC_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	HP_processed_classification_data <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.mean_diff, attr(HP_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	MH_processed_classification_data <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.mean_diff, attr(MH_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	NG_processed_classification_data <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.mean_diff, attr(NG_classification_data$unique_motifs_signature.mean_diff, annotation_type))

	processed_classification_data<-rbind(BA_processed_classification_data,BF_processed_classification_data,CP_processed_classification_data,EC_processed_classification_data,
	HP_processed_classification_data,MH_processed_classification_data,NG_processed_classification_data)

	BA_processed_classification_data.std_diff <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.std_diff, attr(BA_classification_data$unique_motifs_signature.std_diff, annotation_type))
	BF_processed_classification_data.std_diff <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.std_diff, attr(BF_classification_data$unique_motifs_signature.std_diff, annotation_type))
	CP_processed_classification_data.std_diff <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.std_diff, attr(CP_classification_data$unique_motifs_signature.std_diff, annotation_type))
	EC_processed_classification_data.std_diff <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.std_diff, attr(EC_classification_data$unique_motifs_signature.std_diff, annotation_type))
	HP_processed_classification_data.std_diff <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.std_diff, attr(HP_classification_data$unique_motifs_signature.std_diff, annotation_type))
	MH_processed_classification_data.std_diff <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.std_diff, attr(MH_classification_data$unique_motifs_signature.std_diff, annotation_type))
	NG_processed_classification_data.std_diff <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.std_diff, attr(NG_classification_data$unique_motifs_signature.std_diff, annotation_type))

	processed_classification_data.std_diff<-rbind(BA_processed_classification_data.std_diff,BF_processed_classification_data.std_diff,CP_processed_classification_data.std_diff,EC_processed_classification_data.std_diff,
	HP_processed_classification_data.std_diff,MH_processed_classification_data.std_diff,NG_processed_classification_data.std_diff)

	BA_processed_classification_data.time_diff <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.time_diff, attr(BA_classification_data$unique_motifs_signature.time_diff, annotation_type))
	BF_processed_classification_data.time_diff <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.time_diff, attr(BF_classification_data$unique_motifs_signature.time_diff, annotation_type))
	CP_processed_classification_data.time_diff <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.time_diff, attr(CP_classification_data$unique_motifs_signature.time_diff, annotation_type))
	EC_processed_classification_data.time_diff <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.time_diff, attr(EC_classification_data$unique_motifs_signature.time_diff, annotation_type))
	HP_processed_classification_data.time_diff <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.time_diff, attr(HP_classification_data$unique_motifs_signature.time_diff, annotation_type))
	MH_processed_classification_data.time_diff <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.time_diff, attr(MH_classification_data$unique_motifs_signature.time_diff, annotation_type))
	NG_processed_classification_data.time_diff <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.time_diff, attr(NG_classification_data$unique_motifs_signature.time_diff, annotation_type))

	processed_classification_data.time_diff<-rbind(BA_processed_classification_data.time_diff,BF_processed_classification_data.time_diff,CP_processed_classification_data.time_diff,EC_processed_classification_data.time_diff,
	HP_processed_classification_data.time_diff,MH_processed_classification_data.time_diff,NG_processed_classification_data.time_diff)

	BA_processed_classification_data.t_test_pval <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.t_test_pval, attr(BA_classification_data$unique_motifs_signature.t_test_pval, annotation_type))
	BF_processed_classification_data.t_test_pval <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.t_test_pval, attr(BF_classification_data$unique_motifs_signature.t_test_pval, annotation_type))
	CP_processed_classification_data.t_test_pval <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.t_test_pval, attr(CP_classification_data$unique_motifs_signature.t_test_pval, annotation_type))
	EC_processed_classification_data.t_test_pval <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.t_test_pval, attr(EC_classification_data$unique_motifs_signature.t_test_pval, annotation_type))
	HP_processed_classification_data.t_test_pval <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.t_test_pval, attr(HP_classification_data$unique_motifs_signature.t_test_pval, annotation_type))
	MH_processed_classification_data.t_test_pval <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.t_test_pval, attr(MH_classification_data$unique_motifs_signature.t_test_pval, annotation_type))
	NG_processed_classification_data.t_test_pval <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.t_test_pval, attr(NG_classification_data$unique_motifs_signature.t_test_pval, annotation_type))

	processed_classification_data.t_test_pval<-rbind(BA_processed_classification_data.t_test_pval,BF_processed_classification_data.t_test_pval,CP_processed_classification_data.t_test_pval,EC_processed_classification_data.t_test_pval,
	HP_processed_classification_data.t_test_pval,MH_processed_classification_data.t_test_pval,NG_processed_classification_data.t_test_pval)
	
	BA_processed_classification_data.u_test_pval <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.u_test_pval, attr(BA_classification_data$unique_motifs_signature.u_test_pval, annotation_type))
	BF_processed_classification_data.u_test_pval <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.u_test_pval, attr(BF_classification_data$unique_motifs_signature.u_test_pval, annotation_type))
	CP_processed_classification_data.u_test_pval <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.u_test_pval, attr(CP_classification_data$unique_motifs_signature.u_test_pval, annotation_type))
	EC_processed_classification_data.u_test_pval <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.u_test_pval, attr(EC_classification_data$unique_motifs_signature.u_test_pval, annotation_type))
	HP_processed_classification_data.u_test_pval <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.u_test_pval, attr(HP_classification_data$unique_motifs_signature.u_test_pval, annotation_type))
	MH_processed_classification_data.u_test_pval <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.u_test_pval, attr(MH_classification_data$unique_motifs_signature.u_test_pval, annotation_type))
	NG_processed_classification_data.u_test_pval <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.u_test_pval, attr(NG_classification_data$unique_motifs_signature.u_test_pval, annotation_type))

	processed_classification_data.u_test_pval<-rbind(BA_processed_classification_data.u_test_pval,BF_processed_classification_data.u_test_pval,CP_processed_classification_data.u_test_pval,EC_processed_classification_data.u_test_pval,
	HP_processed_classification_data.u_test_pval,MH_processed_classification_data.u_test_pval,NG_processed_classification_data.u_test_pval)
	
	BA_processed_classification_data.basic_group <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.basic_group, attr(BA_classification_data$unique_motifs_signature.basic_group, annotation_type))
	BF_processed_classification_data.basic_group <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.basic_group, attr(BF_classification_data$unique_motifs_signature.basic_group, annotation_type))
	CP_processed_classification_data.basic_group <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.basic_group, attr(CP_classification_data$unique_motifs_signature.basic_group, annotation_type))
	EC_processed_classification_data.basic_group <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.basic_group, attr(EC_classification_data$unique_motifs_signature.basic_group, annotation_type))
	HP_processed_classification_data.basic_group <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.basic_group, attr(HP_classification_data$unique_motifs_signature.basic_group, annotation_type))
	MH_processed_classification_data.basic_group <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.basic_group, attr(MH_classification_data$unique_motifs_signature.basic_group, annotation_type))
	NG_processed_classification_data.basic_group <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.basic_group, attr(NG_classification_data$unique_motifs_signature.basic_group, annotation_type))

	processed_classification_data.basic_group<-rbind(BA_processed_classification_data.basic_group,BF_processed_classification_data.basic_group,CP_processed_classification_data.basic_group,EC_processed_classification_data.basic_group,
	HP_processed_classification_data.basic_group,MH_processed_classification_data.basic_group,NG_processed_classification_data.basic_group)

	BA_processed_classification_data.avg_basequality <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.avg_basequality, attr(BA_classification_data$unique_motifs_signature.avg_basequality, annotation_type))
	BF_processed_classification_data.avg_basequality <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.avg_basequality, attr(BF_classification_data$unique_motifs_signature.avg_basequality, annotation_type))
	CP_processed_classification_data.avg_basequality <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.avg_basequality, attr(CP_classification_data$unique_motifs_signature.avg_basequality, annotation_type))
	EC_processed_classification_data.avg_basequality <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.avg_basequality, attr(EC_classification_data$unique_motifs_signature.avg_basequality, annotation_type))
	HP_processed_classification_data.avg_basequality <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.avg_basequality, attr(HP_classification_data$unique_motifs_signature.avg_basequality, annotation_type))
	MH_processed_classification_data.avg_basequality <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.avg_basequality, attr(MH_classification_data$unique_motifs_signature.avg_basequality, annotation_type))
	NG_processed_classification_data.avg_basequality <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.avg_basequality, attr(NG_classification_data$unique_motifs_signature.avg_basequality, annotation_type))

	processed_classification_data.avg_basequality<-rbind(BA_processed_classification_data.avg_basequality,BF_processed_classification_data.avg_basequality,CP_processed_classification_data.avg_basequality,EC_processed_classification_data.avg_basequality,
	HP_processed_classification_data.avg_basequality,MH_processed_classification_data.avg_basequality,NG_processed_classification_data.avg_basequality)

	BA_processed_classification_data.base_num_fraction <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.base_num_fraction, attr(BA_classification_data$unique_motifs_signature.base_num_fraction, annotation_type))
	BF_processed_classification_data.base_num_fraction <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.base_num_fraction, attr(BF_classification_data$unique_motifs_signature.base_num_fraction, annotation_type))
	CP_processed_classification_data.base_num_fraction <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.base_num_fraction, attr(CP_classification_data$unique_motifs_signature.base_num_fraction, annotation_type))
	EC_processed_classification_data.base_num_fraction <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.base_num_fraction, attr(EC_classification_data$unique_motifs_signature.base_num_fraction, annotation_type))
	HP_processed_classification_data.base_num_fraction <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.base_num_fraction, attr(HP_classification_data$unique_motifs_signature.base_num_fraction, annotation_type))
	MH_processed_classification_data.base_num_fraction <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.base_num_fraction, attr(MH_classification_data$unique_motifs_signature.base_num_fraction, annotation_type))
	NG_processed_classification_data.base_num_fraction <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.base_num_fraction, attr(NG_classification_data$unique_motifs_signature.base_num_fraction, annotation_type))

	processed_classification_data.base_num_fraction<-rbind(BA_processed_classification_data.base_num_fraction,BF_processed_classification_data.base_num_fraction,CP_processed_classification_data.base_num_fraction,EC_processed_classification_data.base_num_fraction,
	HP_processed_classification_data.base_num_fraction,MH_processed_classification_data.base_num_fraction,NG_processed_classification_data.base_num_fraction)

	BA_processed_classification_data.g_num <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.g_num, attr(BA_classification_data$unique_motifs_signature.g_num, annotation_type))
	BF_processed_classification_data.g_num <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.g_num, attr(BF_classification_data$unique_motifs_signature.g_num, annotation_type))
	CP_processed_classification_data.g_num <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.g_num, attr(CP_classification_data$unique_motifs_signature.g_num, annotation_type))
	EC_processed_classification_data.g_num <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.g_num, attr(EC_classification_data$unique_motifs_signature.g_num, annotation_type))
	HP_processed_classification_data.g_num <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.g_num, attr(HP_classification_data$unique_motifs_signature.g_num, annotation_type))
	MH_processed_classification_data.g_num <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.g_num, attr(MH_classification_data$unique_motifs_signature.g_num, annotation_type))
	NG_processed_classification_data.g_num <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.g_num, attr(NG_classification_data$unique_motifs_signature.g_num, annotation_type))

	processed_classification_data.g_num<-rbind(BA_processed_classification_data.g_num,BF_processed_classification_data.g_num,CP_processed_classification_data.g_num,EC_processed_classification_data.g_num,
	HP_processed_classification_data.g_num,MH_processed_classification_data.g_num,NG_processed_classification_data.g_num)

	BA_processed_classification_data.g_value <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.g_value, attr(BA_classification_data$unique_motifs_signature.g_value, annotation_type))
	BF_processed_classification_data.g_value <- prepare.classifier.input(BF_classification_data$unique_motifs_signature.g_value, attr(BF_classification_data$unique_motifs_signature.g_value, annotation_type))
	CP_processed_classification_data.g_value <- prepare.classifier.input(CP_classification_data$unique_motifs_signature.g_value, attr(CP_classification_data$unique_motifs_signature.g_value, annotation_type))
	EC_processed_classification_data.g_value <- prepare.classifier.input(EC_classification_data$unique_motifs_signature.g_value, attr(EC_classification_data$unique_motifs_signature.g_value, annotation_type))
	HP_processed_classification_data.g_value <- prepare.classifier.input(HP_classification_data$unique_motifs_signature.g_value, attr(HP_classification_data$unique_motifs_signature.g_value, annotation_type))
	MH_processed_classification_data.g_value <- prepare.classifier.input(MH_classification_data$unique_motifs_signature.g_value, attr(MH_classification_data$unique_motifs_signature.g_value, annotation_type))
	NG_processed_classification_data.g_value <- prepare.classifier.input(NG_classification_data$unique_motifs_signature.g_value, attr(NG_classification_data$unique_motifs_signature.g_value, annotation_type))

	processed_classification_data.g_value<-rbind(BA_processed_classification_data.g_value,BF_processed_classification_data.g_value,CP_processed_classification_data.g_value,EC_processed_classification_data.g_value,
	HP_processed_classification_data.g_value,MH_processed_classification_data.g_value,NG_processed_classification_data.g_value)

	motifs_possible_mod_pos<-c(2,2,3,0,1,4,4,4,5,1,2,1,2,0,5,2,3,2,4,2,2,1,1,2,3,2,4,2,4,2,3,5,3,0,2,4,6,3,2,1,2,2,4,4,2,2)

	processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
	processed_classification_data.basic_group <- generate.offsetted.dataset(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
	processed_classification_data.time_diff <- generate.offsetted.dataset(processed_classification_data.time_diff, length_vector, selected_starts, background_type)
	processed_classification_data.std_diff <- generate.offsetted.dataset(processed_classification_data.std_diff, length_vector, selected_starts, background_type)
	processed_classification_data.t_test_pval <- generate.offsetted.dataset(processed_classification_data.t_test_pval, length_vector, selected_starts, background_type)
	processed_classification_data.u_test_pval <- generate.offsetted.dataset(processed_classification_data.u_test_pval, length_vector, selected_starts, background_type)
	processed_classification_data.avg_basequality <- generate.offsetted.dataset(processed_classification_data.avg_basequality, length_vector, selected_starts, background_type)
	processed_classification_data.base_num_fraction <- generate.offsetted.dataset(processed_classification_data.base_num_fraction, length_vector, selected_starts, background_type)
	processed_classification_data.g_num <- generate.offsetted.dataset(processed_classification_data.g_num, length_vector, selected_starts, background_type)
	processed_classification_data.g_value <- generate.offsetted.dataset(processed_classification_data.g_value, length_vector, selected_starts, background_type)

	# Generate exclusion list
	#下面的源代码没注释
	#motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	#classification_data_label<-rbind(attr(BA_classification_data, "annotation_motif")$label)
	#motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences

	#train ALL
	classifier_datasets <- prepare.classifier.dataset.notest.basic_group.time(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.std_diff,processed_classification_data.time_diff,processed_classification_data.t_test_pval,processed_classification_data.u_test_pval,processed_classification_data.avg_basequality,processed_classification_data.base_num_fraction,processed_classification_data.g_num,processed_classification_data.g_value, excluded_motifs, 1, balancing, motif_center_summary)
	saveRDS(classifier_datasets$training_data, file = "my_train_models/all_46_basic_group/All46_train_data_balance_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.std_diff, file = "my_train_models/all_46_basic_group/All46_train_data_balance_std_diff_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.time_diff, file = "my_train_models/all_46_basic_group/All46_train_data_balance_time_diff_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.basic_group, file = "my_train_models/all_46_basic_group/All46_train_data_balance_basic_group_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.t_test_pval, file = "my_train_models/all_46_basic_group/All46_train_data_balance_t_test_pval_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.u_test_pval, file = "my_train_models/all_46_basic_group/All46_train_data_balance_u_test_pval_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.avg_basequality, file = "my_train_models/all_46_basic_group/All46_train_data_balance_avg_basequality_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.base_num_fraction, file = "my_train_models/all_46_basic_group/All46_train_data_balance_base_num_fraction_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.g_num, file = "my_train_models/all_46_basic_group/All46_train_data_balance_g_num_time_12_code_filter_iso_offset_wo_balance.rds")
	saveRDS(classifier_datasets$training_data.g_value, file = "my_train_models/all_46_basic_group/All46_train_data_balance_g_value_time_12_code_filter_iso_offset_wo_balance.rds")
	return(excluded_motifs)
}
evaluate.classifiers.performance.basic_group.single <- function(BA_classification_data, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)
	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.mean_diff, attr(BA_classification_data$unique_motifs_signature.mean_diff, annotation_type))
	processed_classification_data<-BA_processed_classification_data #add motif attribute

	BA_processed_classification_data.sd_diff <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.sd_diff, attr(BA_classification_data$unique_motifs_signature.sd_diff, annotation_type))
	processed_classification_data.sd_diff<-BA_processed_classification_data.sd_diff

	BA_processed_classification_data.nat <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.nat, attr(BA_classification_data$unique_motifs_signature.nat, annotation_type))
	processed_classification_data.nat<-BA_processed_classification_data.nat

	BA_processed_classification_data.wga <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.wga, attr(BA_classification_data$unique_motifs_signature.wga, annotation_type))
	processed_classification_data.wga<-BA_processed_classification_data.wga

	BA_processed_classification_data.basic_group <- prepare.classifier.input(BA_classification_data$unique_motifs_signature.basic_group, attr(BA_classification_data$unique_motifs_signature.basic_group, annotation_type))
	processed_classification_data.basic_group<-BA_processed_classification_data.basic_group

	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data.basic_group <- generate.offsetted.dataset(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset(processed_classification_data.wga, length_vector, selected_starts, background_type)
		#processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data.basic_group <- generate.offsetted.dataset(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset(processed_classification_data.wga, length_vector, selected_starts, background_type)
	}	
	# Generate exclusion list
	#下面的源代码没注释
	#motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	#classification_data_label<-rbind(attr(BA_classification_data, "annotation_motif")$label)
	#motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences

	#train ALL
	classifier_datasets <- prepare.classifier.dataset.notest.basic_group(processed_classification_data.basic_group,processed_classification_data,processed_classification_data.sd_diff,processed_classification_data.nat,processed_classification_data.wga, excluded_motifs, 1, balancing, motif_center_summary)
	saveRDS(classifier_datasets$training_data, file = "my_train_models/all_46_basic_group/BA_train_data_balance_1_17.rds")
	saveRDS(classifier_datasets$training_data.sd_diff, file = "my_train_models/all_46_basic_group/BA_train_data_balance_sd_diff_1_17.rds")
	saveRDS(classifier_datasets$training_data.basic_group, file = "my_train_models/all_46_basic_group/BA_train_data_balance_basic_group_1_17.rds")
	saveRDS(classifier_datasets$training_data.nat, file = "my_train_models/all_46_basic_group/BA_train_data_balance_nat_1_17.rds")
	saveRDS(classifier_datasets$training_data.wga, file = "my_train_models/all_46_basic_group/BA_train_data_balance_wga_1_17.rds")
	#meta_classifier_data <- evaluate.classifier(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,NA)
	return(excluded_motifs)
}

evaluate.classifiers.performance.basic_group.single.array <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)
	#BA_processed_classification_data <- BA_classification_data$unique_motifs_signature.mean_diff
	processed_classification_data<-cbind(BA_classification_data$unique_motifs_signature.mean_diff,BF_classification_data$unique_motifs_signature.mean_diff,CP_classification_data$unique_motifs_signature.mean_diff,
	EC_classification_data$unique_motifs_signature.mean_diff,HP_classification_data$unique_motifs_signature.mean_diff,
	MH_classification_data$unique_motifs_signature.mean_diff,NG_classification_data$unique_motifs_signature.mean_diff)

	#BA_processed_classification_data.sd_diff <- BA_classification_data$unique_motifs_signature.sd_diff
	processed_classification_data.sd_diff<-cbind(BA_classification_data$unique_motifs_signature.sd_diff,BF_classification_data$unique_motifs_signature.sd_diff,CP_classification_data$unique_motifs_signature.sd_diff,
	EC_classification_data$unique_motifs_signature.sd_diff,HP_classification_data$unique_motifs_signature.sd_diff,
	MH_classification_data$unique_motifs_signature.sd_diff,NG_classification_data$unique_motifs_signature.sd_diff)

	#BA_processed_classification_data.nat <- BA_classification_data$unique_motifs_signature.nat
	dim_BA=dim(BA_classification_data$unique_motifs_signature.nat)
	dim_BF=dim(BF_classification_data$unique_motifs_signature.nat)
	dim_CP=dim(CP_classification_data$unique_motifs_signature.nat)
	dim_EC=dim(EC_classification_data$unique_motifs_signature.nat)
	dim_HP=dim(HP_classification_data$unique_motifs_signature.nat)
	dim_MH=dim(MH_classification_data$unique_motifs_signature.nat)
	dim_NG=dim(NG_classification_data$unique_motifs_signature.nat)
	dim_All=dim_BA[3]+dim_BF[3]+dim_CP[3]+dim_EC[3]+dim_HP[3]+dim_MH[3]+dim_NG[3]
	processed_classification_data.nat<-array(c(BA_classification_data$unique_motifs_signature.nat,BF_classification_data$unique_motifs_signature.nat,CP_classification_data$unique_motifs_signature.nat,
	EC_classification_data$unique_motifs_signature.nat,HP_classification_data$unique_motifs_signature.nat,
	MH_classification_data$unique_motifs_signature.nat,NG_classification_data$unique_motifs_signature.nat),c(128,22,dim_All))

	#BA_processed_classification_data.wga <- BA_classification_data$unique_motifs_signature.wga
	processed_classification_data.wga<-array(c(BA_classification_data$unique_motifs_signature.wga,BF_classification_data$unique_motifs_signature.wga,CP_classification_data$unique_motifs_signature.wga,
	EC_classification_data$unique_motifs_signature.wga,HP_classification_data$unique_motifs_signature.wga,
	MH_classification_data$unique_motifs_signature.wga,NG_classification_data$unique_motifs_signature.wga),c(128,22,dim_All))

	#BA_processed_classification_data.basic_group <- BA_classification_data$unique_motifs_signature.basic_group
	processed_classification_data.basic_group<-cbind(BA_classification_data$unique_motifs_signature.basic_group,BF_classification_data$unique_motifs_signature.basic_group,CP_classification_data$unique_motifs_signature.basic_group,
	EC_classification_data$unique_motifs_signature.basic_group,HP_classification_data$unique_motifs_signature.basic_group,
	MH_classification_data$unique_motifs_signature.basic_group,NG_classification_data$unique_motifs_signature.basic_group)

	#BA_processed_classification_data.label <- BA_classification_data$unique_motifs_signature.label
	processed_classification_data.label<-cbind(BA_classification_data$unique_motifs_signature.label,BF_classification_data$unique_motifs_signature.label,CP_classification_data$unique_motifs_signature.label,
	EC_classification_data$unique_motifs_signature.label,HP_classification_data$unique_motifs_signature.label,
	MH_classification_data$unique_motifs_signature.label,NG_classification_data$unique_motifs_signature.label)
	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data.basic_group <- generate.offsetted.dataset.array(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset.array(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset.array(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset.array.electric(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset.array.electric(processed_classification_data.wga, length_vector, selected_starts, background_type)
		processed_classification_data.label <- generate.offsetted.dataset.array.label(processed_classification_data.label, length_vector, selected_starts, background_type)
		#processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data.basic_group <- generate.offsetted.dataset.array(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset.array(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset.array(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset.array.electric(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset.array.electric(processed_classification_data.wga, length_vector, selected_starts, background_type)
		processed_classification_data.label <- generate.offsetted.dataset.array.label(processed_classification_data.label, length_vector, selected_starts, background_type)
	}

	motifs_signature.label_data_frame=as.data.frame(processed_classification_data.label[1,])
	colnames(motifs_signature.label_data_frame) <- c("label")
	sample_num<- data.frame(label=NA,train_idx=NA)
	max_sample_size <- min(table(motifs_signature.label_data_frame[,1])) # Find smallest category

	unique_label=unique(motifs_signature.label_data_frame[,1])
	print(unique_label)
	# filter 4mc 5mc 6ma to same number
	result<-foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data_label <- subset(motifs_signature.label_data_frame, label==label_type)
		rowname_subtype_training_data_label<-rownames(subtype_training_data_label)
		train_idx <- sample(rowname_subtype_training_data_label, size=max_sample_size, replace=FALSE) # Subsample every type to smallest category
        new_row <- data.frame(label=label_type,train_idx=train_idx)
        sample_num <- rbind(sample_num, new_row)
	}
	sample_num<-sample_num[-1,]

	sample_num$train_idx <- as.integer(sample_num$train_idx) 
	processed_classification_data.label <- processed_classification_data.label[sample_num$train_idx]
	processed_classification_data <- processed_classification_data[,sample_num$train_idx]
	processed_classification_data.basic_group <- processed_classification_data.basic_group[,sample_num$train_idx]
	processed_classification_data.sd_diff<-processed_classification_data.sd_diff[,sample_num$train_idx]
	processed_classification_data.nat<-processed_classification_data.nat[,,sample_num$train_idx]
	processed_classification_data.wga<-processed_classification_data.wga[,,sample_num$train_idx]
	processed_classification_data.label<-array(processed_classification_data.label,c(1,length(processed_classification_data.label)))


	processed_classification_data.basic_group[processed_classification_data.basic_group=='A']=0
	processed_classification_data.basic_group[processed_classification_data.basic_group=='C']=1
	processed_classification_data.basic_group[processed_classification_data.basic_group=='T']=2
	processed_classification_data.basic_group[processed_classification_data.basic_group=='G']=3
	len=dim(processed_classification_data.basic_group)
	processed_classification_data.basic_group=array(as.integer(processed_classification_data.basic_group),c(len[1],len[2]))

	processed_classification_data.label[processed_classification_data.label=='5mC_-12']=0 
	processed_classification_data.label[processed_classification_data.label=='5mC_-11']=1
	processed_classification_data.label[processed_classification_data.label=='5mC_-10']=2
	processed_classification_data.label[processed_classification_data.label=='5mC_-9']=3
	processed_classification_data.label[processed_classification_data.label=='5mC_-8']=4
	processed_classification_data.label[processed_classification_data.label=='5mC_-7']=5
	processed_classification_data.label[processed_classification_data.label=='5mC_-6']=6   

	processed_classification_data.label[processed_classification_data.label=='4mC_-12']=7 
	processed_classification_data.label[processed_classification_data.label=='4mC_-11']=8
	processed_classification_data.label[processed_classification_data.label=='4mC_-10']=9
	processed_classification_data.label[processed_classification_data.label=='4mC_-9']=10
	processed_classification_data.label[processed_classification_data.label=='4mC_-8']=11
	processed_classification_data.label[processed_classification_data.label=='4mC_-7']=12
	processed_classification_data.label[processed_classification_data.label=='4mC_-6']=13  

	processed_classification_data.label[processed_classification_data.label=='6mA_-12']=14
	processed_classification_data.label[processed_classification_data.label=='6mA_-11']=15
	processed_classification_data.label[processed_classification_data.label=='6mA_-10']=16
	processed_classification_data.label[processed_classification_data.label=='6mA_-9']=17
	processed_classification_data.label[processed_classification_data.label=='6mA_-8']=18
	processed_classification_data.label[processed_classification_data.label=='6mA_-7']=19
	processed_classification_data.label[processed_classification_data.label=='6mA_-6']=20 
	processed_classification_data.label=array(as.integer(processed_classification_data.label))
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	saveRDS(processed_classification_data.label, file = "my_train_models/all_46_basic_group/train_data_balance_label_1_28.rds")
	saveRDS(processed_classification_data, file = "my_train_models/all_46_basic_group/train_data_balance_1_28.rds")
	saveRDS(processed_classification_data.sd_diff, file = "my_train_models/all_46_basic_group/train_data_balance_sd_diff_1_28.rds")
	saveRDS(processed_classification_data.basic_group, file = "my_train_models/all_46_basic_group/train_data_balance_basic_group_1_28.rds")
	saveRDS(processed_classification_data.nat, file = "my_train_models/all_46_basic_group/train_data_balance_nat_1_28.rds")
	saveRDS(processed_classification_data.wga, file = "my_train_models/all_46_basic_group/train_data_balance_wga_1_28.rds")
	return(excluded_motifs)
}

evaluate.classifiers.performance.basic_group.single.array.test <- function(BA_classification_data,CP_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)
	#BA_processed_classification_data <- BA_classification_data$unique_motifs_signature.mean_diff
	processed_classification_data<-cbind(BA_classification_data$unique_motifs_signature.mean_diff,CP_classification_data$unique_motifs_signature.mean_diff)

	#BA_processed_classification_data.sd_diff <- BA_classification_data$unique_motifs_signature.sd_diff
	processed_classification_data.sd_diff<-cbind(BA_classification_data$unique_motifs_signature.sd_diff,CP_classification_data$unique_motifs_signature.sd_diff)

	#BA_processed_classification_data.nat <- BA_classification_data$unique_motifs_signature.nat
	dim_BA=dim(BA_classification_data$unique_motifs_signature.nat)

	dim_CP=dim(CP_classification_data$unique_motifs_signature.nat)

	dim_All=dim_BA[3]+dim_CP[3]
	processed_classification_data.nat<-array(c(BA_classification_data$unique_motifs_signature.nat,CP_classification_data$unique_motifs_signature.nat),c(128,22,dim_All))

	#BA_processed_classification_data.wga <- BA_classification_data$unique_motifs_signature.wga
	processed_classification_data.wga<-array(c(BA_classification_data$unique_motifs_signature.wga,CP_classification_data$unique_motifs_signature.wga),c(128,22,dim_All))

	#BA_processed_classification_data.basic_group <- BA_classification_data$unique_motifs_signature.basic_group
	processed_classification_data.basic_group<-cbind(BA_classification_data$unique_motifs_signature.basic_group,CP_classification_data$unique_motifs_signature.basic_group)

	#BA_processed_classification_data.label <- BA_classification_data$unique_motifs_signature.label
	processed_classification_data.label<-cbind(BA_classification_data$unique_motifs_signature.label,CP_classification_data$unique_motifs_signature.label)

	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data.basic_group <- generate.offsetted.dataset.array(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset.array(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset.array(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset.array.electric(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset.array.electric(processed_classification_data.wga, length_vector, selected_starts, background_type)
		processed_classification_data.label <- generate.offsetted.dataset.array.label(processed_classification_data.label, length_vector, selected_starts, background_type)
		#processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data.basic_group <- generate.offsetted.dataset.array(processed_classification_data.basic_group, length_vector, selected_starts, background_type)
		processed_classification_data <- generate.offsetted.dataset.array(processed_classification_data, length_vector, selected_starts, background_type)
		processed_classification_data.sd_diff <- generate.offsetted.dataset.array(processed_classification_data.sd_diff, length_vector, selected_starts, background_type)
		processed_classification_data.nat <- generate.offsetted.dataset.array.electric(processed_classification_data.nat, length_vector, selected_starts, background_type)
		processed_classification_data.wga <- generate.offsetted.dataset.array.electric(processed_classification_data.wga, length_vector, selected_starts, background_type)
		processed_classification_data.label <- generate.offsetted.dataset.array.label(processed_classification_data.label, length_vector, selected_starts, background_type)
	}	

	motifs_signature.label_data_frame=as.data.frame(processed_classification_data.label[1,])
	colnames(motifs_signature.label_data_frame) <- c("label")
	sample_num<- data.frame(label=NA,train_idx=NA)
	max_sample_size <- min(table(motifs_signature.label_data_frame[,1])) # Find smallest category
	unique_label=unique(motifs_signature.label_data_frame[,1])
	print(unique_label)
	# filter 4mc 5mc 6ma to same number
	result<-foreach(label_type=unique_label, .combine=rbind) %do% {
		subtype_training_data_label <- subset(motifs_signature.label_data_frame, label==label_type)
		rowname_subtype_training_data_label<-rownames(subtype_training_data_label)
		train_idx <- sample(rowname_subtype_training_data_label, size=max_sample_size, replace=FALSE) # Subsample every type to smallest category
        new_row <- data.frame(label=label_type,train_idx=train_idx)
        sample_num <- rbind(sample_num, new_row)
	}
	sample_num<-sample_num[-1,]

	sample_num$train_idx <- as.integer(sample_num$train_idx) 
	processed_classification_data.label <- processed_classification_data.label[sample_num$train_idx]
	processed_classification_data <- processed_classification_data[,sample_num$train_idx]
	processed_classification_data.basic_group <- processed_classification_data.basic_group[,sample_num$train_idx]
	processed_classification_data.sd_diff<-processed_classification_data.sd_diff[,sample_num$train_idx]
	processed_classification_data.nat<-processed_classification_data.nat[,,sample_num$train_idx]
	processed_classification_data.wga<-processed_classification_data.wga[,,sample_num$train_idx]
	processed_classification_data.label<-array(processed_classification_data.label,c(1,length(processed_classification_data.label)))
	

	processed_classification_data.basic_group[processed_classification_data.basic_group=='A']=0
	processed_classification_data.basic_group[processed_classification_data.basic_group=='C']=1
	processed_classification_data.basic_group[processed_classification_data.basic_group=='T']=2
	processed_classification_data.basic_group[processed_classification_data.basic_group=='G']=3
	len=dim(processed_classification_data.basic_group)
	processed_classification_data.basic_group=array(as.integer(processed_classification_data.basic_group),c(len[1],len[2]))

	processed_classification_data.label[processed_classification_data.label=='5mC_-12']=0 
	processed_classification_data.label[processed_classification_data.label=='5mC_-11']=1
	processed_classification_data.label[processed_classification_data.label=='5mC_-10']=2
	processed_classification_data.label[processed_classification_data.label=='5mC_-9']=3
	processed_classification_data.label[processed_classification_data.label=='5mC_-8']=4
	processed_classification_data.label[processed_classification_data.label=='5mC_-7']=5
	processed_classification_data.label[processed_classification_data.label=='5mC_-6']=6   

	processed_classification_data.label[processed_classification_data.label=='4mC_-12']=7 
	processed_classification_data.label[processed_classification_data.label=='4mC_-11']=8
	processed_classification_data.label[processed_classification_data.label=='4mC_-10']=9
	processed_classification_data.label[processed_classification_data.label=='4mC_-9']=10
	processed_classification_data.label[processed_classification_data.label=='4mC_-8']=11
	processed_classification_data.label[processed_classification_data.label=='4mC_-7']=12
	processed_classification_data.label[processed_classification_data.label=='4mC_-6']=13  

	processed_classification_data.label[processed_classification_data.label=='6mA_-12']=14
	processed_classification_data.label[processed_classification_data.label=='6mA_-11']=15
	processed_classification_data.label[processed_classification_data.label=='6mA_-10']=16
	processed_classification_data.label[processed_classification_data.label=='6mA_-9']=17
	processed_classification_data.label[processed_classification_data.label=='6mA_-8']=18
	processed_classification_data.label[processed_classification_data.label=='6mA_-7']=19
	processed_classification_data.label[processed_classification_data.label=='6mA_-6']=20 
	processed_classification_data.label=array(as.integer(processed_classification_data.label))
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	saveRDS(processed_classification_data.label, file = "my_train_models/all_46_basic_group/train_data_balance_label_1_255.rds")
	saveRDS(processed_classification_data, file = "my_train_models/all_46_basic_group/train_data_balance_1_25.rds")
	saveRDS(processed_classification_data.sd_diff, file = "my_train_models/all_46_basic_group/train_data_balance_sd_diff_1_25.rds")
	saveRDS(processed_classification_data.basic_group, file = "my_train_models/all_46_basic_group/train_data_balance_basic_group_1_25.rds")
	saveRDS(processed_classification_data.nat, file = "my_train_models/all_46_basic_group/train_data_balance_nat_1_25.rds")
	saveRDS(processed_classification_data.wga, file = "my_train_models/all_46_basic_group/train_data_balance_wga_1_25.rds")
	return(excluded_motifs)
}

evaluate.classifiers.performance.ALL <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)

	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data, attr(BA_classification_data, annotation_type))
	BF_processed_classification_data <- prepare.classifier.input(BF_classification_data, attr(BF_classification_data, annotation_type))
	CP_processed_classification_data <- prepare.classifier.input(CP_classification_data, attr(CP_classification_data, annotation_type))
	EC_processed_classification_data <- prepare.classifier.input(EC_classification_data, attr(EC_classification_data, annotation_type))
	HP_processed_classification_data <- prepare.classifier.input(HP_classification_data, attr(HP_classification_data, annotation_type))
	MH_processed_classification_data <- prepare.classifier.input(MH_classification_data, attr(MH_classification_data, annotation_type))
	NG_processed_classification_data <- prepare.classifier.input(NG_classification_data, attr(NG_classification_data, annotation_type))

	processed_classification_data<-rbind(BA_processed_classification_data,BF_processed_classification_data,CP_processed_classification_data,EC_processed_classification_data,
	HP_processed_classification_data,MH_processed_classification_data,NG_processed_classification_data)
	motifs_possible_mod_pos<-c(2,2,3,0,1,4,4,4,5,1,2,1,2,0,5,2,3,2,4,2,2,1,1,2,3,2,4,2,4,2,3,5,3,0,2,4,6,3,2,1,2,2,4,4,2,2)
	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
		#processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
	}	
	# Generate exclusion list
	#下面的源代码没注释
	#motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	#classification_data_label<-rbind(attr(BA_classification_data, "annotation_motif")$label)
	#motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences

	#train ALL
	classifier_datasets <- prepare.classifier.dataset.notest(processed_classification_data, excluded_motifs, 1, balancing, motif_center_summary)
	#saveRDS(classifier_datasets, file = "my_train_models/all_46/All46_train_data_balance.rds")
	meta_classifier_data <- evaluate.classifier(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,NA)
	return(excluded_motifs)
}

evaluate.classifiers.performance.LOOCV.draw <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)

	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data, attr(BA_classification_data, annotation_type))
	BF_processed_classification_data <- prepare.classifier.input(BF_classification_data, attr(BF_classification_data, annotation_type))
	CP_processed_classification_data <- prepare.classifier.input(CP_classification_data, attr(CP_classification_data, annotation_type))
	EC_processed_classification_data <- prepare.classifier.input(EC_classification_data, attr(EC_classification_data, annotation_type))
	HP_processed_classification_data <- prepare.classifier.input(HP_classification_data, attr(HP_classification_data, annotation_type))
	MH_processed_classification_data <- prepare.classifier.input(MH_classification_data, attr(MH_classification_data, annotation_type))
	NG_processed_classification_data <- prepare.classifier.input(NG_classification_data, attr(NG_classification_data, annotation_type))

	processed_classification_data<-rbind(BA_processed_classification_data,BF_processed_classification_data,CP_processed_classification_data,EC_processed_classification_data,
	HP_processed_classification_data,MH_processed_classification_data,NG_processed_classification_data)
	motifs_possible_mod_pos<-c(2,2,3,0,1,4,4,4,5,1,2,1,2,0,5,2,3,2,4,2,2,1,1,2,3,2,4,2,4,2,3,5,3,0,2,4,6,3,2,1,2,2,4,4,2,2)
	# Construct offsetted dataset
	if(length(selected_starts)==1){ # If only one window tested
		# Do not specify offset in label
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
		#processed_classification_data <- subset(processed_classification_data, select=c("label","motif",as.character(seq(selected_starts,selected_starts + (length_vector - 1)))))
	}else{
		processed_classification_data <- generate.offsetted.dataset(processed_classification_data, length_vector, selected_starts, background_type)
	}	
	# Generate exclusion list
	#下面的源代码没注释
	#motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	#classification_data_label<-rbind(attr(BA_classification_data, "annotation_motif")$label)
	#motif_summary <- motif_summary[motif_summary$explicit_motif %in% unique(attr(classification_data, "annotation_motif")$label),] # Remove untested motifs
	excluded_motifs <- t(combn(motif_summary$explicit_motif, nb_motifs))
	#excluded_motifs <- remove.impossible.exclusion(excluded_motifs, motif_summary) # Delete combination removing all modification type occurrences

	#train ALL
	#classifier_datasets <- prepare.classifier.dataset.notest(processed_classification_data, excluded_motifs, 1, balancing, motif_center_summary)
	#meta_classifier_data <- evaluate.classifier.predict(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,NA,NA)

	#draw ALL
	registerDoMC(nbCPU)
	meta_classifier_data <- foreach(idx=seq(1,25), .options.multicore=list(preschedule=TRUE), .combine=rbind) %dopar% { #nrow(excluded_motifs)
		# Split datsets in training/testing
		classifier_datasets <- prepare.classifier.dataset(processed_classification_data, excluded_motifs, idx, balancing, motif_center_summary)
		# Train and test classifier
		classifier_results <- evaluate.classifier.predict(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,excluded_motifs[idx,],motifs_possible_mod_pos[idx])
		return(classifier_results)
	}
	registerDoSEQ()
	draw.classification.results(meta_classifier_data,"predict_on_predict_mod_balance", "my_predict/loocv/")
	sink("train/LOOCV.result1.txt")
	print(meta_classifier_data)
	sink()
	registerDoMC(nbCPU)
	meta_classifier_data <- foreach(idx=seq(26,nrow(excluded_motifs)), .options.multicore=list(preschedule=TRUE), .combine=rbind) %dopar% { #nrow(excluded_motifs)
		# Split datsets in training/testing
		classifier_datasets <- prepare.classifier.dataset(processed_classification_data, excluded_motifs, idx, balancing, motif_center_summary)
		# Train and test classifier
		classifier_results <- evaluate.classifier.predict(classifier_datasets$training_data, classifier_datasets$testing_data, list_classifier,excluded_motifs[idx,],motifs_possible_mod_pos[idx])
		return(classifier_results)
	}
	registerDoSEQ()
	draw.classification.results(meta_classifier_data,"predict_on_predict_mod_balance_2", "my_predict/loocv/")
	sink("train/LOOCV.result2.txt")
	print(meta_classifier_data)
	sink()
	return(meta_classifier_data)
}

evaluate.classifiers.performance.LOOCV.tsne <- function(BA_classification_data,BF_classification_data,CP_classification_data,EC_classification_data,HP_classification_data,MH_classification_data,NG_classification_data
, annotation_type, motif_summary, motif_center_summary, length_vector, selected_starts, nb_motifs, nbCPU, balancing, list_classifier, background_type="multi"){
	library(VGAM) # TODO dup for speed-up
	library(MASS)
	library(caret)
	library(mda)
	library(klaR)
	library(nnet)
	library(kernlab)
	library(e1071) # Discard?
	library(randomForest)

	BA_processed_classification_data <- prepare.classifier.input(BA_classification_data, attr(BA_classification_data, annotation_type))
	BF_processed_classification_data <- prepare.classifier.input(BF_classification_data, attr(BF_classification_data, annotation_type))
	CP_processed_classification_data <- prepare.classifier.input(CP_classification_data, attr(CP_classification_data, annotation_type))
	EC_processed_classification_data <- prepare.classifier.input(EC_classification_data, attr(EC_classification_data, annotation_type))
	HP_processed_classification_data <- prepare.classifier.input(HP_classification_data, attr(HP_classification_data, annotation_type))
	MH_processed_classification_data <- prepare.classifier.input(MH_classification_data, attr(MH_classification_data, annotation_type))
	NG_processed_classification_data <- prepare.classifier.input(NG_classification_data, attr(NG_classification_data, annotation_type))

	processed_classification_data<-rbind(BA_processed_classification_data,BF_processed_classification_data,CP_processed_classification_data,EC_processed_classification_data,
	HP_processed_classification_data,MH_processed_classification_data,NG_processed_classification_data)
	tmp_data <- subset(processed_classification_data, select=c(as.character(seq(-14,7))))
	tmp_label <-  subset(processed_classification_data, select=c("motif"))
	tmp_label <- tmp_label %>% 
  	mutate(ID=row_number())
	tsne_out = Rtsne(
	as.matrix(tmp_data),
	dims = 2,
	pca = T,
	max_iter = 1000,
	theta = 0.4,
	perplexity = 20,
	verbose = F
	) # 进行t-SNE降维分析

	tsne_result = as.data.frame(tsne_out$Y)
	colnames(tsne_result) = c("tSNE1","tSNE2")
	tsne_result <- tsne_result %>%
    mutate(ID=row_number())
	tsne_result <- tsne_result %>%
    inner_join(tmp_label, by="ID")
	saveRDS(tsne_result,file = 'TSNE/tsne_result.rds')

	draw_result<-ggplot(tsne_result,aes(tSNE1,tSNE2,color=motif)) +
	geom_point()+
	facet_grid(~motif)
	
	ggsave(filename = "TSNE/tsne_motif.pdf", plot = draw_result, device = "pdf")
	return(processed_classification_data)
}

merge.classification.results <- function(list_performance_classifiers){
	classification_results <- NULL

	# Merge summaries
	classification_results$summary <- foreach(classification_result=list_performance_classifiers, .combine=rbind) %do% {
		return(classification_result$summary)
	}

	# Merge details
	classification_results$details <- foreach(classification_result=list_performance_classifiers, .combine=rbind) %do% {
		return(classification_result$details)
	}

	return(classification_results)
}

draw.classifier.performance.realistic <- function(performance_classifier, graph_suffix, motif_center_summary){
	meta_classifier_summary <- performance_classifier$summary
	meta_classifier_details <- performance_classifier$details

	nb_exclude_motifs <- length(levels(meta_classifier_details$excluded_motifs))
	nb_classifiers <- length(levels(meta_classifier_details$classifier))
	nb_classes <- length(levels(meta_classifier_details$Prediction))
	nb_position <- length(unique(str_split(meta_classifier_details$Prediction,"_",simplify=TRUE)[,2]))
	nb_mod <- length(unique(str_split(meta_classifier_details$Prediction,"_",simplify=TRUE)[,1]))

	gp <- ggplot(meta_classifier_summary) +
		geom_point(aes(x=excluded_motifs,y=accuracy)) +
		facet_wrap(~classifier) +
		labs(title="Classifier Performance", subtitle="By classifier") +
		labs(x="Tested motif (exclude from training)", y="Accuracy") +
		theme(axis.text.x=element_text(angle=-30, hjust=0))
	pdf(paste0("Classifier_statistics_byMethod_acc_",graph_suffix,"_v1.pdf"), width=4 + ceiling(sqrt(nb_classifiers)) * nb_exclude_motifs * 0.15, height=4 + floor(sqrt(nb_classifiers)) * nb_exclude_motifs * 0.15)
	print(gp)
	dev.off()

	gp <- ggplot(meta_classifier_summary) +
		geom_point(aes(x=classifier,y=accuracy)) +
		facet_wrap(~excluded_motifs) +
		labs(title="Classifier Performance", subtitle="By motif") +
		labs(x="Classifier", y="Accuracy") +
		theme(axis.text.x=element_text(angle=-30, hjust=0))
	pdf(paste0("Classifier_statistics_byMotif_acc_",graph_suffix,"_v1.pdf"), width=4 + ceiling(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.35, height=4 + floor(sqrt(nb_exclude_motifs)) * 2)
	print(gp)
	dev.off()

	# Tidy prediction & summarize for motif
	summary_meta_classifier_details <- meta_classifier_details %>%
		group_by(excluded_motifs, classifier, Reference) %>%
		filter(sum(Freq)>0) %>%
		mutate(explicit_motif=as.character(excluded_motifs)) %>%
		inner_join(subset(motif_center_summary, select=c(mod_pos, mod_type, explicit_motif, possible_mod_pos)), by=c("explicit_motif"="explicit_motif")) %>%
		mutate(corresponding_mod_type=str_split(Prediction,"_",simplify=TRUE)[,1]) %>%
		mutate(corresponding_mod_pos=(possible_mod_pos - (as.integer(str_split(Prediction,"_",simplify=TRUE)[,2]) + 9)) + 1) %>% # possible_mod_pos is 0-based and mod_pos 1-based
		mutate(score=(Freq/sum(Freq))*100) %>%
		mutate(prediction=ifelse(max(score)==score,TRUE,FALSE)) %>%
		mutate(expectation=ifelse(corresponding_mod_pos==mod_pos & corresponding_mod_type==mod_type,ifelse(prediction,"bold.italic","italic"),ifelse(prediction,"bold","plain"))) %>%
		mutate(TP=ifelse(corresponding_mod_pos==mod_pos & corresponding_mod_type==mod_type,TRUE,FALSE)) %>%
		mutate(tmp=ifelse(prediction & TP,TRUE,FALSE))

	gp <- ggplot(summary_meta_classifier_details, aes(x=corresponding_mod_pos, y=classifier)) +
		geom_tile(aes(fill=score), colour="white") +
		geom_text(aes(label=round(score), fontface=expectation)) +
		facet_grid(corresponding_mod_type~excluded_motifs, scales="free_x") +
		scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
		scale_x_continuous(breaks=seq(min(summary_meta_classifier_details$mod_pos), max(summary_meta_classifier_details$mod_pos), by=1)) +
		labs(title="Motif Characterization", subtitle="By excluded motif") +
		labs(x="Modification Position", y="Modification Type", fill="Percentage Prediction")
	pdf(paste0("Classifier_performance_DNAmod_Position_",graph_suffix,"_v1.pdf"), width=4 + nb_exclude_motifs * nb_position * 0.30, height=4 + nb_classifiers * nb_mod * 0.15)
	print(gp)
	dev.off()

	summary2_meta_classifier_details <- summary_meta_classifier_details %>%
		group_by(excluded_motifs, classifier, TP) %>%
		summarize(global_score=sum(score), .groups="drop_last")

	gp <- ggplot(summary2_meta_classifier_details, aes(TP, classifier)) +
		geom_tile(aes(fill=global_score), colour="white") +
		geom_text(aes(label=round(global_score))) +
		facet_wrap(~excluded_motifs) +
		scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
		labs(title="Prediction frequency", subtitle="By excluded motif") +
		labs(x="TP/FP Prediction", y="Classifier", fill="Percentage") +
		theme(axis.text.x=element_text(angle=-30, hjust=0))
	pdf(paste0("Classifier_performance_",graph_suffix,"_summary_v1.pdf"), width=4 + ceiling(sqrt(nb_exclude_motifs)) * 2 * 0.7, height=4 + floor(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.2)
	print(gp)
	dev.off()
	
	gp <- ggplot(subset(summary2_meta_classifier_details, TP==TRUE)) +
		stat_summary(aes(x=classifier, y=global_score), fun.data="mean_cl_boot") +
		geom_quasirandom(aes(x=classifier, y=global_score, col=excluded_motifs), alpha=0.9) +
		scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
		labs(title="Prediction frequency", subtitle="By excluded motif") +
		labs(x="Classifier", y="TPR") +
		theme(axis.text.x=element_text(angle=-30, hjust=0))
	pdf(paste0("Classifier_performance_",graph_suffix,"_summary_v2.pdf"), width=4 + nb_classifiers + nb_exclude_motifs * 0.05, height=5 + nb_exclude_motifs * 0.05)
	print(gp)
	dev.off()

	final_stat <- summary_meta_classifier_details %>%
		group_by(classifier, excluded_motifs) %>%
		summarize(res=ifelse(any(tmp),"Match","Error"), .groups="drop_last") %>%
		group_by(classifier, res) %>%
		summarize(res2=n(), .groups="drop_last") %>%
		spread(res,res2) %>%
		mutate(Accuracy=Match*100/(Match+ifelse(is.na(Error),0,Error)))
	print(final_stat)
}

draw.classifier.performance <- function(meta_classifier_summary, meta_classifier_details, graph_suffix){
	nb_exclude_motifs <- length(levels(meta_classifier_details$excluded_motifs))
	nb_classifiers <- length(levels(meta_classifier_details$classifier))
	nb_classes <- length(levels(meta_classifier_details$Prediction))

	gp <- ggplot(meta_classifier_summary) +
		geom_point(aes(x=excluded_motifs,y=accuracy)) +
		facet_wrap(~classifier) +
		labs(title="Classifier Performance", subtitle="By classifier") +
		labs(x="Tested motif (exclude from training)", y="Accuracy") +
		theme(axis.text.x=element_text(angle=-30, hjust=0))
	pdf(paste0("Classifier_statistics_byMethod_acc_",graph_suffix,"_v1.pdf"), width=4 + ceiling(sqrt(nb_classifiers)) * nb_exclude_motifs * 0.15, height=4 + floor(sqrt(nb_classifiers)) * nb_exclude_motifs * 0.15)
	print(gp)
	dev.off()

	gp <- ggplot(meta_classifier_summary) +
		geom_point(aes(x=classifier,y=accuracy)) +
		facet_wrap(~excluded_motifs) +
		labs(title="Classifier Performance", subtitle="By motif") +
		labs(x="Classifier", y="Accuracy") +
		theme(axis.text.x=element_text(angle=-30, hjust=0))
	pdf(paste0("Classifier_statistics_byMotif_acc_",graph_suffix,"_v1.pdf"), width=4 + ceiling(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.35, height=4 + floor(sqrt(nb_exclude_motifs)) * 1)
	print(gp)
	dev.off()

	# Mark good predictions
	meta_classifier_details$TP <- ifelse(as.character(meta_classifier_details$Prediction)==as.character(meta_classifier_details$Reference),TRUE,FALSE)

	if(length(strsplit(as.character(meta_classifier_details$Reference[1]),"_")[[1]])==2){ # If type and pos; Not great
		# Add relative position information
		meta_classifier_details$Ref_Pos <- do.call(rbind,strsplit(as.character(meta_classifier_details$Reference),"_"))[,2]
		meta_classifier_details$Ref_Pos <- as.factor(meta_classifier_details$Ref_Pos) # Bad order
		meta_classifier_details$Ref_Pos <- factor(meta_classifier_details$Ref_Pos, levels=mixedsort(levels(meta_classifier_details$Ref_Pos)))
		meta_classifier_details$Ref_Mod <- do.call(rbind,strsplit(as.character(meta_classifier_details$Reference),"_"))[,1]
		meta_classifier_details$Ref_Mod <- as.factor(meta_classifier_details$Ref_Mod) # Bad order
		meta_classifier_details$Ref_Mod <- factor(meta_classifier_details$Ref_Mod, levels=mixedsort(levels(meta_classifier_details$Ref_Mod)))
		nb_position <- length(unique(meta_classifier_details$Ref_Pos))
		nb_mod <- length(unique(meta_classifier_details$Ref_Mod))
		# Add relative position information
		meta_classifier_details$Pred_Pos <- do.call(rbind,strsplit(as.character(meta_classifier_details$Prediction),"_"))[,2]
		meta_classifier_details$Pred_Pos <- as.factor(meta_classifier_details$Pred_Pos) # Bad order
		meta_classifier_details$Pred_Pos <- factor(meta_classifier_details$Pred_Pos, levels=mixedsort(levels(meta_classifier_details$Pred_Pos)))
		meta_classifier_details$Pred_Mod <- do.call(rbind,strsplit(as.character(meta_classifier_details$Prediction),"_"))[,1]
		meta_classifier_details$Pred_Mod <- as.factor(meta_classifier_details$Pred_Mod) # Bad order
		meta_classifier_details$Pred_Mod <- factor(meta_classifier_details$Pred_Mod, levels=mixedsort(levels(meta_classifier_details$Pred_Mod)))

		# Transform Freq to real percentage
		meta_classifier_details <- meta_classifier_details %>%
			group_by(Ref_Pos, excluded_motifs, classifier) %>%
			mutate(Freq=(Freq*100)/sum(Freq, na.rm=TRUE))

		# By DNA modification type and offset
		summary_meta_classifier_details <- meta_classifier_details %>%
			group_by(Ref_Mod, Ref_Pos, classifier, excluded_motifs, TP) %>%
			summarize(sum_Freq=sum(Freq, na.rm=TRUE), .groups="drop_last") %>%
			mutate(group=paste0(TP,".",Ref_Mod))
		summary_meta_classifier_details$group <- as.factor(summary_meta_classifier_details$group) # Bad order
		summary_meta_classifier_details$group <- factor(summary_meta_classifier_details$group, levels=mixedsort(levels(summary_meta_classifier_details$group)))

		gp <- ggplot(summary_meta_classifier_details, aes(group, classifier)) +
			geom_tile(aes(fill=sum_Freq), colour="white") +
			geom_text(aes(label=round(sum_Freq))) +
			facet_grid(Ref_Pos~excluded_motifs) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="TP/FP Position and DNA mod. type Prediction", y="Classifier", fill="Percentage") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_DNAmod_Position_",graph_suffix,"_v1.pdf"), width=4 + nb_exclude_motifs * (nb_mod * 2) * 0.25, height=4 + nb_classifiers * nb_position * 0.15)
		print(gp)
		dev.off()

		# Global
		summary_meta_classifier_details <- meta_classifier_details %>%
			group_by(classifier, excluded_motifs, TP) %>%
			summarize(sum_Freq=sum(Freq, na.rm=TRUE)/nb_position, .groups="drop_last")

		gp <- ggplot(summary_meta_classifier_details, aes(TP, classifier)) +
			geom_tile(aes(fill=sum_Freq), colour="white") +
			geom_text(aes(label=round(sum_Freq))) +
			facet_wrap(~excluded_motifs) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="TP/FP Prediction", y="Classifier", fill="Percentage") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_",graph_suffix,"_summary_v1.pdf"), width=4 + ceiling(sqrt(nb_exclude_motifs)) * 2 * 0.7, height=4 + floor(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.2)
		print(gp)
		dev.off()
		
		gp <- ggplot(subset(summary_meta_classifier_details, TP==TRUE)) +
			stat_summary(aes(x=classifier, y=sum_Freq), fun.data="mean_cl_boot") +
			geom_quasirandom(aes(x=classifier, y=sum_Freq, col=excluded_motifs), alpha=0.9) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="Classifier", y="TPR") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_",graph_suffix,"_summary_v2.pdf"), width=4 + nb_classifiers + nb_exclude_motifs * 0.05, height=5 + nb_exclude_motifs * 0.05)
		print(gp)
		dev.off()
	
		# By offset
		# Mark good predictions
		meta_classifier_details$TP <- ifelse(as.character(meta_classifier_details$Ref_Pos)==as.character(meta_classifier_details$Pred_Pos),TRUE,FALSE)

		summary_meta_classifier_details <- meta_classifier_details %>%
			group_by(Ref_Pos, classifier, excluded_motifs, TP) %>%
			summarize(sum_Freq=sum(Freq, na.rm=TRUE), .groups="drop_last") %>%
			mutate(group=paste0(TP,".",Ref_Pos))
		summary_meta_classifier_details$group <- as.factor(summary_meta_classifier_details$group) # Bad order
		summary_meta_classifier_details$group <- factor(summary_meta_classifier_details$group, levels=mixedsort(levels(summary_meta_classifier_details$group)))
		
		gp <- ggplot(summary_meta_classifier_details, aes(group, classifier)) +
			geom_tile(aes(fill=sum_Freq), colour="white") +
			geom_text(aes(label=round(sum_Freq))) +
			facet_wrap(~excluded_motifs) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="TP/FP Position Prediction", y="Classifier", fill="Percentage") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_Position_",graph_suffix,"_v1.pdf"), width=4 + ceiling(sqrt(nb_exclude_motifs)) * nb_classes * 0.2, height=4 + floor(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.2)
		print(gp)
		dev.off()

		# By DNA modification type
		# Mark good predictions
		meta_classifier_details$TP <- ifelse(as.character(meta_classifier_details$Ref_Mod)==as.character(meta_classifier_details$Pred_Mod),TRUE,FALSE)

		summary_meta_classifier_details <- meta_classifier_details %>%
			group_by(Ref_Mod, classifier, excluded_motifs, TP) %>%
			summarize(sum_Freq=sum(Freq, na.rm=TRUE)/nb_position, .groups="drop_last") %>%
			mutate(group=paste0(TP,".",Ref_Mod))
		summary_meta_classifier_details$group <- as.factor(summary_meta_classifier_details$group) # Good order

		gp <- ggplot(summary_meta_classifier_details, aes(group, classifier)) +
			geom_tile(aes(fill=sum_Freq), colour="white") +
			geom_text(aes(label=round(sum_Freq))) +
			facet_wrap(~excluded_motifs) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="TP/FP DNA mod. type Prediction", y="Classifier", fill="Percentage") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_DNAmod_",graph_suffix,"_v1.pdf"), width=4 + ceiling(sqrt(nb_exclude_motifs)) * nb_classes * 0.1 , height=4 + floor(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.2)
		print(gp)
		dev.off()
	}else{
		# Transform Freq to real percentage
		meta_classifier_details <- meta_classifier_details %>%
			group_by(excluded_motifs, classifier) %>%
			mutate(Freq=(Freq*100)/sum(Freq, na.rm=TRUE))

		# By DNA modification type
		summary_meta_classifier_details <- meta_classifier_details %>%
			group_by(classifier, excluded_motifs, TP, Reference) %>%
			summarize(sum_Freq=sum(Freq, na.rm=TRUE), .groups="drop_last") %>%
			mutate(group=paste0(TP,".",Reference))
		summary_meta_classifier_details$group <- as.factor(summary_meta_classifier_details$group) # Good order

		gp <- ggplot(summary_meta_classifier_details, aes(x=group, y=classifier)) +
			geom_tile(aes(fill=sum_Freq), colour="white") +
			geom_text(aes(label=round(sum_Freq))) +
			facet_wrap(~excluded_motifs) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="Prediction/Reference", y="Classifier", fill="Percentage") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_DNAmod_",graph_suffix,"_v1.pdf"), width=ceiling(sqrt(nb_exclude_motifs)) * nb_classes^2 * 0.5, height=floor(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.25)
		print(gp)
		dev.off()

		# Global
		summary_meta_classifier_details <- meta_classifier_details %>%
			group_by(classifier, excluded_motifs, TP) %>%
			summarizse(sum_Freq=sum(Freq, na.rm=TRUE))

		gp <- ggplot(summary_meta_classifier_details, aes(TP, classifier)) +
			geom_tile(aes(fill=sum_Freq), colour="white") +
			geom_text(aes(label=round(sum_Freq))) +
			facet_wrap(~excluded_motifs) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="TP/FP Prediction", y="Classifier", fill="Percentage") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_",graph_suffix,"_summary_v1.pdf"), width=ceiling(sqrt(nb_exclude_motifs)) * 2 * 1, height=floor(sqrt(nb_exclude_motifs)) * nb_classifiers * 0.25)
		print(gp)
		dev.off()

		gp <- ggplot(subset(summary_meta_classifier_details, TP==TRUE)) +
			stat_summary(aes(x=classifier, y=sum_Freq), fun.data="mean_cl_boot") +
			geom_quasirandom(aes(x=classifier, y=sum_Freq, col=excluded_motifs), alpha=0.9) +
			scale_fill_gradientn(colours=rev(brewer.pal(11,'Spectral'))) +
			labs(title="Prediction frequency", subtitle="By excluded motif") +
			labs(x="Classifier", y="TPR") +
			theme(axis.text.x=element_text(angle=-30, hjust=0))
		pdf(paste0("Classifier_performance_",graph_suffix,"_summary_v2.pdf"), width=nb_classifiers, height=5)
		print(gp)
		dev.off()
	}
}

draw.classifier.results <- function(performance_classifier, motif_summary, base_name){
	meta_classifier_details <- performance_classifier$details
	default_offset <- -9 # TODO change if position signal corrected

	# Split categories
	meta_classifier_details$Ref_Pos <- do.call(rbind,strsplit(as.character(meta_classifier_details$Reference),"_"))[,2]
	meta_classifier_details$Ref_Pos <- as.numeric(meta_classifier_details$Ref_Pos) # Bad order
	meta_classifier_details$Ref_Mod <- do.call(rbind,strsplit(as.character(meta_classifier_details$Reference),"_"))[,1]
	meta_classifier_details$Ref_Mod <- as.factor(meta_classifier_details$Ref_Mod) # Bad order
	meta_classifier_details$Ref_Mod <- factor(meta_classifier_details$Ref_Mod, levels=mixedsort(levels(meta_classifier_details$Ref_Mod)))
	meta_classifier_details$Pred_Pos <- do.call(rbind,strsplit(as.character(meta_classifier_details$Prediction),"_"))[,2]
	meta_classifier_details$Pred_Pos <- as.numeric(meta_classifier_details$Pred_Pos) # Bad order
	meta_classifier_details$Pred_Mod <- do.call(rbind,strsplit(as.character(meta_classifier_details$Prediction),"_"))[,1]
	meta_classifier_details$Pred_Mod <- as.factor(meta_classifier_details$Pred_Mod) # Bad order
	meta_classifier_details$Pred_Mod <- factor(meta_classifier_details$Pred_Mod, levels=mixedsort(levels(meta_classifier_details$Pred_Mod)))

	# Analyze classifier results by motif
	motif_summary$explicit_motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
	classifier_results <- foreach(idx_motif=seq_along(motif_summary$explicit_motif)) %do% {
		motif <- motif_summary$motif[idx_motif]
		explicit_motif <- motif_summary$explicit_motif[idx_motif]
		mod_pos <- motif_summary$mod_pos[idx_motif]
		mod_type <- motif_summary$mod_type[idx_motif]

		# Limit positions to evaluate
		# TODO refine position from motif signatures

		# list_DNA_mod_pos <- gregexpr("[CA]", motif, perl=TRUE)[[1]] # Only test C and A because only expecting 4mC, 5mC or 6mA
		list_DNA_mod_pos <- seq(1:nchar(motif)) # Seems better
		offset <- default_offset - (mod_pos - 1)
		list_DNA_mod_pos <- list_DNA_mod_pos + offset

		# subset_meta_classifier_details <- subset(meta_classifier_details, excluded_motifs==motif)
		subset_meta_classifier_details <- subset(meta_classifier_details, excluded_motifs==explicit_motif & Ref_Pos %in% list_DNA_mod_pos) # Similar

		predicted_mod_type <- subset_meta_classifier_details %>%
			group_by(excluded_motifs, classifier, Ref_Pos) %>% # Ref_Pos==Mod_Type * Offset
			mutate(Perc=(Freq*100)/sum(Freq)) %>%
			group_by(excluded_motifs, classifier, Ref_Pos, Pred_Mod) %>%
			summarize(tmp_score=sum(Perc), .groups="drop_last") %>% # Per offseted subset
			group_by(excluded_motifs, classifier, Pred_Mod) %>%
			summarize(score_mod=sum(tmp_score)/n(), .groups="drop_last") %>%
			dplyr::rename(Mod_Type=Pred_Mod)

		# Duplicate but needed for mod pos
		results_mod_type <- predicted_mod_type %>% # Summarize across offseted subset
			filter(score_mod==max(score_mod)) # Keep best prediction

		subset_meta_classifier_details <- merge(subset_meta_classifier_details, results_mod_type)

		# TODO refine mod. pos range from motif signature
		
		predicted_mod_pos <- subset_meta_classifier_details %>%
			group_by(excluded_motifs, classifier, Ref_Pos) %>%
			mutate(Perc=(Freq*100)/sum(Freq)) %>%
			mutate(relative_Ref_Pos=Ref_Pos + 9, relative_Pred_Pos=Pred_Pos + 9) %>% # Convert to position relative to DNA mod.
			mutate(Motif_Pos=relative_Ref_Pos + mod_pos) %>% # Convert offseted data to motif position, i.e. what we know in real usage 
			mutate(Mod_Pos=Motif_Pos - relative_Pred_Pos) %>%
			filter(Ref_Mod==Mod_Type) %>%
			group_by(excluded_motifs, classifier, Mod_Pos) %>%
			summarize(score_pos=sum(Perc), .groups="drop_last")

		return(list(type=as.data.frame(predicted_mod_type), pos=as.data.frame(predicted_mod_pos)))
	}

	# Tidy results
	tmp_pos <- foreach(res=classifier_results, .combine=rbind) %do% {
		return(res$pos)
	}

	tmp_type <- foreach(res=classifier_results, .combine=rbind) %do% {
		return(res$type)
	}
	classifier_results <- NULL
	classifier_results$pos <- tmp_pos
	classifier_results$type <- tmp_type

	# Modification Type
	classifier_results$type <- merge(classifier_results$type, motif_summary, by.x="excluded_motifs", by.y="explicit_motif") # Add mod_pos
	best_mod_prediction <- classifier_results$type %>%
		group_by(excluded_motifs, classifier) %>%
		filter(score_mod>0) %>%
		filter(score_mod==max(score_mod))
		
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	gp <- ggplot(classifier_results$type) +
		geom_tile(aes(y=classifier, x=Mod_Type, fill=score_mod)) +
		geom_point(aes(y=classifier, x=mod_type, colour="black"), size=0.5) +
		geom_point(data=best_mod_prediction, aes(y=classifier, x=Mod_Type, colour="red"), size=0.5) +
		facet_wrap(~excluded_motifs) +
		scale_fill_gradientn(name="DNA Mod.\nScore", colours=myPalette(100)) +
		scale_colour_manual(name="DNA Mod.\nType", values=c("black"="black","red"="red"), labels=c("Reality","Prediction")) +
		labs(title="DNA Modification Type Prediction per Classifier") +
		labs(x="DNA Modification Type", y="Classifier") +
		theme_bw()
	pdf(paste0("Classifier_results_DNAmod_",base_name,"_v1.pdf"), width=15, height=10)
	print(gp)
	dev.off()

	# Modification Position
	classifier_results$pos <- merge(classifier_results$pos, motif_summary, by.x="excluded_motifs", by.y="explicit_motif") # Add mod_pos
	best_pos_prediction <- classifier_results$pos %>%
		group_by(excluded_motifs, classifier) %>%
		filter(score_pos>0) %>%
		filter(score_pos==max(score_pos))
		
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	gp <- ggplot(classifier_results$pos) +
		geom_tile(aes(y=classifier, x=Mod_Pos, fill=score_pos)) +
		geom_point(aes(y=classifier, x=mod_pos, colour="black"), size=0.5) +
		geom_point(data=best_pos_prediction, aes(y=classifier, x=Mod_Pos, colour="red"), size=0.5) +
		facet_wrap(~excluded_motifs) +
		scale_fill_gradientn(name="Position\nScore", colours=myPalette(100)) +
		scale_colour_manual(name="DNA Mod.\nPosition", values=c("black"="black","red"="red"), labels=c("Expectation","Prediction")) +
		scale_x_continuous(breaks=seq(min(classifier_results$pos$Mod_Pos), max(classifier_results$pos$Mod_Pos))) +
		labs(title="DNA Modification Position Prediction per Classifier") +
		labs(x="Position in Motif", y="Classifier") +
		theme_bw()
	pdf(paste0("Classifier_results_Position_",base_name,"_v1.pdf"), width=15, height=10)
	print(gp)
	dev.off()

	summary_best_pos_prediction <- best_pos_prediction %>%
		group_by(excluded_motifs, classifier) %>%
		mutate(error=Mod_Pos - mod_pos)

	gp <- ggplot(summary_best_pos_prediction) +
		geom_jitter(aes(x=classifier, y=error, colour=col_motif), height=0) +
		scale_y_continuous(breaks=seq(min(summary_best_pos_prediction$error), max(summary_best_pos_prediction$error))) +
		scale_colour_identity(name="Motifs", breaks=motif_summary$col_motif, labels=motif_summary$explicit_motif, guide="legend") +
		labs(title="Summary error in DNA Modification Position Prediction per Classifier") +
		labs(x="Classifier", y="Distance to expected position") +
		theme(axis.text.x=element_text(angle=-30, hjust=0), legend.position="bottom")
	pdf(paste0("Classifier_results_summary_Position_",base_name,"_v1.pdf"), width=8)
	print(gp)
	dev.off()
}

# For unknown motifs 
select.motif <- function(classification_data, motif){
	motif_annotation <- attr(classification_data, "annotation_motif")
	row_to_keep <- grepl(motif, as.character(motif_annotation$id))

	subset_classification_data <- classification_data[row_to_keep,]
	attr(subset_classification_data, "annotation_motif") <- attr(classification_data, "annotation_motif")[row_to_keep,]
	# attr(classification_data, "annotation_mod") <- attr(classification_data, "annotation_mod")[row_to_keep,]
	# attr(classification_data, "annotation_dir") <- attr(classification_data, "annotation_dir")[row_to_keep,]
	# if("annotation_strain" %in% names(attributes(classification_data))){
	# 	attr(classification_data, "annotation_strain") <- attr(classification_data, "annotation_strain")[row_to_keep,]
	# }

	return(subset_classification_data)
}
prepare.classification.data.old <- function(methylation_signal, strain_id, motif_summary, genome, min_cov, filter_iso, iupac_nc, nbCPU, data_type="train"){
	left_signal <- -1 # Remove nothing
	right_signal <- -1 # Remove nothing
	error_margin <- -1 # Remove nothing
	expected_signal_left <- -6
	expected_signal_right <- -1
	signal_margin <- 8 # Huge margin needed for offseted windows

	motifs_signature <- extract.motifs.signature.old(methylation_signal, genome, motif_summary, iupac_nc, left_signal, right_signal, error_margin, expected_signal_left, expected_signal_right, signal_margin, filter_iso, min_cov, nbCPU)
	if(data_type=="train"){
		motifs_signature <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif_summary$explicit_motif[match(motif,motif_summary$motif)],"_",motif_summary$mod_type[match(motif,motif_summary$motif)],"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,mean_diff,distance)) %>%
			spread(distance, mean_diff)
	}else if(data_type=="classify"){
		motifs_signature <- motifs_signature %>%
			mutate(id=paste0(gsub("_","",strain_id),"_",motif,"_",gsub("_","",contig),"_",pos_motif,"_",dir,"_",strand)) %>%
			dplyr::select(c(id,mean_diff,distance)) %>%
			spread(distance, mean_diff)
	}
	rownames(motifs_signature) <- motifs_signature$id
	motifs_signature$id <- NULL
	motifs_signature <- motifs_signature[complete.cases(motifs_signature),] # TODO check if can be handle better

	return(motifs_signature)
}
prepare.annotation.data.old <- function(motifs_signature, ann_type, motif_summary=NA){
	if(ann_type=="motif"){
		motif_summary$old_motif <- motif_summary$motif # TODO modify motif_summary
		motif_summary$motif <- paste0(substr(motif_summary$motif,1,motif_summary$mod_pos-1), as.character(motif_summary$mod_type), substr(motif_summary$motif,motif_summary$mod_pos+1,nchar(motif_summary$motif)))
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,2]) # Motif
		annotation <- merge(annotation, subset(motif_summary, select=c(motif, col_motif)), by.x=c("label"), by.y=c("motif"))
		names(annotation)[names(annotation)=='col_motif'] <- 'col_label'
		annotation$col_label <- as.factor(annotation$col_label)
	}else if(ann_type=="mod_type"){
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,3]) # Modification types
		annotation <- merge(annotation, data.frame(mod_type=c("4mC","5mC","6mA"), col_label=c("#36D2A0","#444CF0","#FF5733")), by.x=c("label"), by.y=c("mod_type"))
		annotation$col_label <- as.factor(annotation$col_label)
	}else if(ann_type=="dir"){
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,6]) # Strand types
		annotation <- merge(annotation, data.frame(dir=c("fwd","rev"), col_label=c("#444CF0","#FF5733")), by.x=c("label"), by.y=c("dir"))
		annotation$col_label <- as.factor(annotation$col_label)
	}else if(ann_type=="classify"){
		annotation <- data.frame(id=rownames(motifs_signature), label=do.call(rbind,strsplit(rownames(motifs_signature),"_"))[,2]) # Motif
		annotation <- merge(annotation, subset(motif_summary, select=c(motif, col_motif)), by.x=c("label"), by.y=c("motif"))
		names(annotation)[names(annotation)=='col_motif'] <- 'col_label'
		annotation$col_label <- as.factor(annotation$col_label)
	}
	annotation <- annotation[match(rownames(motifs_signature), annotation$id),]

	return(annotation)
}
classify.detected.motifs.time <- function(methylation_signal, strain_id, motif_center_summary, model, genome, min_cov, keepIsolated, iupac_nc, nbCPU){
	## Details:
	# We train the model on 7 offseted 12 bp windows -2:+2
	# We approximate mod. pos. by looking at signature center but often +/-1 offset
	# Using approx. as input we fall in training window and output will return improve position + type

	# Default model preprocessing
	length_vector <- 12 # Depend on training parameter; TODO Change if changed for training
	vector_start <- -9 # Considere mod. pos. approximation as input

	# Prepare classification data
	motif_center_summary$mod_pos <- motif_center_summary$possible_mod_pos + 1 # Convert to 1-based; mod_pos used as signature center downstream of prepare.classification.data
	#classification_data <- prepare.classification.data.old(methylation_signal, strain_id, motif_center_summary, genome, min_cov, keepIsolated, iupac_nc, nbCPU, "classify")
	
	motifs_signature <- prepare.classification.data.basic_group.time(methylation_signal,basecall, strain_id, motif_summary, genome,genome, min_cov, keepIsolated, iupac_nc, nbCPU)
	unique_motifs_signature.mean_diff<-motifs_signature$motifs_signature.mean_diff
	unique_motifs_signature.basic_group<-motifs_signature$motifs_signature.basic_group

	#attr(classification_data, "annotation_motif") <- prepare.annotation.data.old(classification_data, "classify", motif_center_summary)

	attr(unique_motifs_signature.mean_diff, "annotation_motif") <- prepare.annotation.data.old(unique_motifs_signature.mean_diff, "classify", motif_center_summary)
	attr(unique_motifs_signature.basic_group, "annotation_motif") <- prepare.annotation.data.old(unique_motifs_signature.basic_group, "classify", motif_center_summary)
	# Predict modification type & position for each motifs
	print_message("  Classify motif(s)")
	classification_results <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.mean_diff, motif)
		subset_classification_data.basic <- select.motif(unique_motifs_signature.basic_group, motif)
		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		features_classification_data.basic <- subset(subset_classification_data.basic, select=c(signature_features))

		colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns

		colnames(features_classification_data.basic) <- c(paste0("pos",seq(1,ncol(features_classification_data.basic)))) # Rename position/features
		features_classification_data.basic$label <- attr(subset_classification_data.basic, "annotation_motif")$label
		features_classification_data.basic$label <- as.factor(paste0(features_classification_data.basic$label,"_",vector_start))
		features_classification_data.basic <- features_classification_data.basic %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data.basic)-1)))

		# Characterize each motif occurrences
		print_message(paste0("    Classifing ",motif_center_summary$motif[idx_motif]))
		if(!is.null(model$modelInfo$library)){
			if(model$modelInfo$library=="randomForest"){
				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}else if(model$modelInfo$library=="klaR"){
				library(klaR) # Not used

				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}else if(model$modelInfo$library=="nnet"){
				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}
		}else if(!is.null(model$modelInfo$label)){
			if(model$modelInfo$label=="k-Nearest Neighbors"){
				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}
		}else if(any(grepl("nnet",class(model)))){			
			print(features_classification_data$label[0])
			prediction_results <- predict(model, subset(features_classification_data, select=-c(label)), type="raw")
			saveRDS(prediction_results, file =  paste0("my_train_models/all_46_basic_group/",motif_center_summary$motif[idx_motif],".rds"))
			saveRDS(features_classification_data.basic, file =  paste0("my_train_models/all_46_basic_group/",motif_center_summary$motif[idx_motif],"_basic.rds"))
			#print(class(prediction_results))
			#print(dim(subset(features_classification_data, select=c(label))))
			#print(dim(prediction_results))
		}else if(any(grepl("randomForest",class(model)))){
			prediction_results <- predict(model,subset(features_classification_data, select=-c(label)))
		}else if(any(grepl("knn3",class(model)))){
			prediction_results <- predict(model, subset(features_classification_data, select=-c(label)), type="class")
		}

		# Tidy prediction & summarize for motif
		prediction_summary <- data.frame(table(prediction_results))
		colnames(prediction_summary) <- c("label","count")

		# Generate complete data.frame
		base_df <- expand.grid(mod=c("4mC","5mC","6mA"), pos=seq(-12,-6)) %>%
			mutate(label=paste0(mod,"_",pos)) %>%
			arrange(label) %>%
			dplyr::select(label)
		
		prediction_summary <- merge(prediction_summary, base_df, by=c("label"), all=TRUE) %>% # Add combination potentially missing
			mutate(count=ifelse(is.na(count),0,count)) %>% # Put 0 when missing
			mutate(motif=motif_center_summary$motif[idx_motif]) %>%
			mutate(mod_type=str_split(label,"_",simplify=TRUE)[,1]) %>%
			mutate(mod_pos=possible_mod_pos - (as.integer(str_split(label,"_",simplify=TRUE)[,2]) - vector_start)) %>% # Leave it 0-based
			mutate(score=(count/sum(count))*100)

		return(prediction_summary)
	}

	return(classification_results)
}

classify.detected.motifs <- function(methylation_signal, strain_id, motif_center_summary, model, genome, min_cov, keepIsolated, iupac_nc, nbCPU){
	## Details:
	# We train the model on 7 offseted 12 bp windows -2:+2
	# We approximate mod. pos. by looking at signature center but often +/-1 offset
	# Using approx. as input we fall in training window and output will return improve position + type

	# Default model preprocessing
	length_vector <- 12 # Depend on training parameter; TODO Change if changed for training
	vector_start <- -9 # Considere mod. pos. approximation as input

	# Prepare classification data
	motif_center_summary$mod_pos <- motif_center_summary$possible_mod_pos + 1 # Convert to 1-based; mod_pos used as signature center downstream of prepare.classification.data
	#classification_data <- prepare.classification.data.old(methylation_signal, strain_id, motif_center_summary, genome, min_cov, keepIsolated, iupac_nc, nbCPU, "classify")
	
	motifs_signature <- prepare.classification.data.basic_group.time(methylation_signal,basecall, strain_id, motif_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU)
	unique_motifs_signature.mean_diff<-motifs_signature$motifs_signature.mean_diff
	unique_motifs_signature.basic_group<-motifs_signature$motifs_signature.basic_group

	#attr(classification_data, "annotation_motif") <- prepare.annotation.data.old(classification_data, "classify", motif_center_summary)

	attr(unique_motifs_signature.mean_diff, "annotation_motif") <- prepare.annotation.data.old(unique_motifs_signature.mean_diff, "classify", motif_center_summary)
	attr(unique_motifs_signature.basic_group, "annotation_motif") <- prepare.annotation.data.old(unique_motifs_signature.basic_group, "classify", motif_center_summary)
	# Predict modification type & position for each motifs
	print_message("  Classify motif(s)")
	classification_results <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(classification_data, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns

		# Characterize each motif occurrences
		print_message(paste0("    Classifing ",motif_center_summary$motif[idx_motif]))
		if(!is.null(model$modelInfo$library)){
			if(model$modelInfo$library=="randomForest"){
				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}else if(model$modelInfo$library=="klaR"){
				library(klaR) # Not used

				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}else if(model$modelInfo$library=="nnet"){
				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}
		}else if(!is.null(model$modelInfo$label)){
			if(model$modelInfo$label=="k-Nearest Neighbors"){
				prediction_results <- predict(model, subset(features_classification_data, select=-c(label)))
			}
		}else if(any(grepl("nnet",class(model)))){			
			print(features_classification_data$label[0])
			prediction_results <- predict(model, subset(features_classification_data, select=-c(label)), type="raw")
			saveRDS(prediction_results, file =  paste0("my_train_models/all_46_basic_group/",motif_center_summary$motif[idx_motif],".rds"))
			#print(class(prediction_results))
			#print(dim(subset(features_classification_data, select=c(label))))
			#print(dim(prediction_results))
		}else if(any(grepl("randomForest",class(model)))){
			prediction_results <- predict(model,subset(features_classification_data, select=-c(label)))
		}else if(any(grepl("knn3",class(model)))){
			prediction_results <- predict(model, subset(features_classification_data, select=-c(label)), type="class")
		}

		# Tidy prediction & summarize for motif
		prediction_summary <- data.frame(table(prediction_results))
		colnames(prediction_summary) <- c("label","count")

		# Generate complete data.frame
		base_df <- expand.grid(mod=c("4mC","5mC","6mA"), pos=seq(-12,-6)) %>%
			mutate(label=paste0(mod,"_",pos)) %>%
			arrange(label) %>%
			dplyr::select(label)
		
		prediction_summary <- merge(prediction_summary, base_df, by=c("label"), all=TRUE) %>% # Add combination potentially missing
			mutate(count=ifelse(is.na(count),0,count)) %>% # Put 0 when missing
			mutate(motif=motif_center_summary$motif[idx_motif]) %>%
			mutate(mod_type=str_split(label,"_",simplify=TRUE)[,1]) %>%
			mutate(mod_pos=possible_mod_pos - (as.integer(str_split(label,"_",simplify=TRUE)[,2]) - vector_start)) %>% # Leave it 0-based
			mutate(score=(count/sum(count))*100)

		return(prediction_summary)
	}

	return(classification_results)
}

classify.detected.motifs.get_test_data <- function(methylation_signal, strain_id, motif_center_summary, model, genome, min_cov, keepIsolated, iupac_nc, nbCPU){
	## Details:
	# We train the model on 7 offseted 12 bp windows -2:+2
	# We approximate mod. pos. by looking at signature center but often +/-1 offset
	# Using approx. as input we fall in training window and output will return improve position + type

	# Default model preprocessing
	length_vector <- 12 # Depend on training parameter; TODO Change if changed for training
	vector_start <- -9 # Considere mod. pos. approximation as input

	# Prepare classification data
	motif_center_summary$mod_pos <- motif_center_summary$possible_mod_pos + 1 # Convert to 1-based; mod_pos used as signature center downstream of prepare.classification.data
	classification_data <- prepare.classification.data(methylation_signal, strain_id, motif_center_summary, genome, min_cov, keepIsolated, iupac_nc, nbCPU, "classify")
	attr(classification_data, "annotation_motif") <- prepare.annotation.data(classification_data, "classify", motif_center_summary)
	# Predict modification type & position for each motifs
	print_message("  Classify motif(s)")
	classification_results <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(classification_data, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}
	saveRDS(classification_results, file = "my_train_models/all_46/evalute_data_NO.rds")

	return(classification_results)
}


classify.detected.motifs.get_test_data.basic_group <- function(methylation_signal,basecall, strain_id, motif_center_summary, model, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU){
	## Details:
	# We train the model on 7 offseted 12 bp windows -2:+2
	# We approximate mod. pos. by looking at signature center but often +/-1 offset
	# Using approx. as input we fall in training window and output will return improve position + type

	# Default model preprocessing
	length_vector <- 22 # Depend on training parameter; TODO Change if changed for training
	vector_start <- -14 # Considere mod. pos. approximation as input
	#length_vector <- 22 # Depend on training parameter; TODO Change if changed for training
	#vector_start <- -14 # Considere mod. pos. approximation as input
	# Prepare classification data
	motif_center_summary$mod_pos <- motif_center_summary$possible_mod_pos + 1 # Convert to 1-based; mod_pos used as signature center downstream of prepare.classification.data
	classification_data <- prepare.classification.data.basic_group.time(methylation_signal,basecall, strain_id, motif_center_summary, genome,genome_rev, min_cov, keepIsolated, iupac_nc, nbCPU, "classify")

	unique_motifs_signature.mean_diff<-classification_data$motifs_signature.mean_diff
	unique_motifs_signature.std_diff<-classification_data$motifs_signature.std_diff
	unique_motifs_signature.time_diff<-classification_data$motifs_signature.time_diff
	unique_motifs_signature.basic_group<-classification_data$motifs_signature.basic_group
	unique_motifs_signature.t_test_pval<-classification_data$motifs_signature.t_test_pval
	unique_motifs_signature.u_test_pval<-classification_data$motifs_signature.u_test_pval
	unique_motifs_signature.avg_basequality<-classification_data$motifs_signature.avg_basequality
	unique_motifs_signature.base_num_fraction<-classification_data$motifs_signature.base_num_fraction
	unique_motifs_signature.g_num<-classification_data$motifs_signature.g_num
	unique_motifs_signature.g_value<-classification_data$motifs_signature.g_value
	attr(unique_motifs_signature.mean_diff, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.basic_group, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.std_diff, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.time_diff, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.t_test_pval, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.u_test_pval, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.avg_basequality, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.base_num_fraction, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.g_num, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	attr(unique_motifs_signature.g_value, "annotation_motif") <- prepare.annotation.data(classification_data$motifs_signature.mean_diff,"classify", motif_center_summary) # ok if not unique_motifs_signature
	#attr(classification_data, "annotation_motif") <- prepare.annotation.data(classification_data, "classify", motif_center_summary)
	# Predict modification type & position for each motifs
	print_message("  Classify motif(s)")
	classification_results <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.mean_diff, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.basic_group <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.basic_group, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.std_diff <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.std_diff, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.time_diff <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.time_diff, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.t_test_pval <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.t_test_pval, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.u_test_pval <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.u_test_pval, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.avg_basequality <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.avg_basequality, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.base_num_fraction <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.base_num_fraction, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.g_num <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.g_num, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}

	classification_results.g_value <- foreach(idx_motif=seq(1, nrow(motif_center_summary)), .combine=rbind) %do% {
		# Select motif & approximate modification position
		motif <- paste0(gsub("_","",strain_id),"_",motif_center_summary$motif[idx_motif])
		possible_mod_pos <- motif_center_summary$possible_mod_pos[idx_motif]

		# Select corresponding broad signature
		subset_classification_data <- select.motif(unique_motifs_signature.g_value, motif)

		# Select feature for prediction
		signature_features <- as.character(seq(vector_start, vector_start + (length_vector - 1)))
		features_classification_data <- subset(subset_classification_data, select=c(signature_features))
		#下面要解开注释
		#colnames(features_classification_data) <- c(paste0("pos",seq(1,ncol(features_classification_data)))) # Rename position/features
		colnames(features_classification_data) <- seq(1,ncol(features_classification_data))
		features_classification_data$label <- attr(subset_classification_data, "annotation_motif")$label
		features_classification_data$label <- as.factor(paste0(features_classification_data$label,"_",vector_start))
		#下面要解开注释
		#features_classification_data <- features_classification_data %>% dplyr::select("label",paste0("pos",seq(1,ncol(features_classification_data)-1))) # Reorder columns
		features_classification_data <- features_classification_data %>% dplyr::select("label",seq(1,ncol(features_classification_data)-1))
		features_classification_data$possible_mod_pos <- possible_mod_pos
		return(features_classification_data)
	}
	baterial<-"NG"
	saveRDS(classification_results, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_12_code_filter_iso.rds"))
	saveRDS(classification_results.basic_group, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_basic_group_12_code_filter_iso.rds"))
	saveRDS(classification_results.std_diff, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_std_diff_12_code_filter_iso.rds"))
	saveRDS(classification_results.time_diff, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_time_diff_12_code_filter_iso.rds"))
	saveRDS(classification_results.t_test_pval, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_t_test_pval_12_code_filter_iso.rds"))
	saveRDS(classification_results.u_test_pval, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_u_test_pval_12_code_filter_iso.rds"))
	saveRDS(classification_results.avg_basequality, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_avg_basequality_12_code_filter_iso.rds"))
	saveRDS(classification_results.base_num_fraction, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_base_num_fraction_12_code_filter_iso.rds"))
	saveRDS(classification_results.g_num, file = paste0("my_train_models/all_46_basic_group/evalute_data_", baterial, "_g_num_12_code_filter_iso.rds"))
	#saveRDS(classification_results.g_value, file = "my_train_models/all_46_basic_group/evalute_data_BA_g_value_12_code_filter_iso.rds")
	return(classification_results)
}

draw.classification.results <- function(classification_results, base_name, path_output=NA){
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	gp <- ggplot(classification_results, aes(x=as.factor(mod_pos), y=mod_type)) +
		geom_tile(aes(fill=score)) +
		geom_text(aes(label=round(score))) +
		facet_wrap(~motif) +
		scale_fill_gradientn(colours=myPalette(100), limits=c(0,100)) +
		coord_cartesian(expand=FALSE) +
		labs(title=paste0(base_name," motifs classification results")) +
		labs(x="Motif position (0-based)", y="Modification Type", fill="Prediction\nScore (%)") +
		guides(fill=guide_colorbar(nbin=100, draw.ulim=FALSE, draw.llim=FALSE))

	nb_row <- round(sqrt(length(unique(classification_results$motif))))

	plot_name <- paste0("Motifs_classification_",base_name,".pdf")
	if(!is.na(path_output)){
		plot_name <- paste0(path_output,plot_name)
	}
	pdf(plot_name, width=9 + (nb_row+1)*0.5, height=3 + nb_row*0.5)
	print(gp)
	dev.off()
}

write.best.classification.results <- function(classification_results, base_name, path_output){
	best_prediction <- classification_results %>%
		group_by(motif) %>%
		mutate(is_within=ifelse(mod_pos >= 0 & mod_pos <= nchar(motif) - 1,TRUE,FALSE)) %>%
		mutate(replacing=str_split(mod_type,"",simplify=TRUE)[,3]) %>%
		filter(mod_pos>=0 & mod_pos<=3) %>%
		rowwise() %>%
		mutate(safe_mod_pos=ifelse(is_within, mod_pos+2, 1)) %>%
		mutate(safe_base=cbind("Z",str_split(motif,"",simplify=TRUE))[,safe_mod_pos]) %>%
		mutate(is_consistent=ifelse(is_within & replacing==safe_base,TRUE,FALSE)) %>%
		ungroup() %>%
		filter(is_consistent) %>%
		# filter(score>20)
		group_by(motif) %>%
		filter(score==max(score)) %>%
		mutate(clean_motif=paste0(substr(motif, 1,mod_pos),mod_type,substr(motif, mod_pos+2,nchar(motif)))) %>%
		mutate(clear_mod_pos=mod_pos+1) %>%
		mutate(clear_score=round(score, 2)) %>%
		dplyr::select(clean_motif, motif, clear_mod_pos, mod_type, clear_score)
	colnames(best_prediction) <- c("Characterized_motif","Motif","Predicted_position","Predicted_type","Prediction_score")

	file_name <- paste0("Motifs_classification_",base_name,".tsv")
	output_file_name <- paste0(path_output,file_name)

	# con <- file(output_file_name, open="wt")
	# writeLines(paste("# This only shows the prediction with highest scores."), con)
	# writeLines(paste("# Please consult the associated .pdf for the full prediction results."), con)
	write.table(best_prediction, output_file_name, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	# close(con)
}

######################
## Visualization IGV
######################

find.motifs <- function(genome, motif_summary, iupac_nc){
	g_seq <- readDNAStringSet(genome) # Not design for multiple contig

	genome_annotation <- foreach(direction=c("fwd","rev"), .combine=rbind) %do% { # Process both strand
		# List all modified bases and keep duplicate
		sub_genome_annotation <- foreach(idx_motif=seq(1,nrow(motif_summary)), .combine=rbind) %do% {
			original_motif <- motif_summary$motif[idx_motif]
			mod_pos <- motif_summary$mod_pos[idx_motif]
			len_motif <- nchar(original_motif)

			if(direction=="rev"){
				motif <- reverseComplement(DNAString(original_motif)) # Double checked
				mod_pos <- (len_motif - mod_pos) + 1
			}else{
				motif <- original_motif
			}

			pos_motif <- vmatchPattern(motif, g_seq, fixed=FALSE)[[1]]@start + (mod_pos - 1)

			return(data.frame(position=pos_motif, motif=as.character(original_motif)))
		}
		sub_genome_annotation$dir <- direction

		return(sub_genome_annotation)
	}
	genome_annotation$dir <- as.factor(genome_annotation$dir)

	return(genome_annotation)
}

generate.motif.annotation.gff <- function(genome, motif_summary, iupac_nc, output_name){
	name_chr <- names(readDNAStringSet(genome))

	genome_annotation <- find.motifs(genome, motif_summary, iupac_nc)
	motif_summary <- motif_summary %>%
		mutate(motif_len=nchar(motif))
	genome_annotation <- merge(genome_annotation, motif_summary)

	genome_annotation <- genome_annotation %>% group_by(motif, dir) %>% mutate(id=row_number())

	stifle <- foreach(direction=c("fwd","rev")) %do% {
		dir_genome_annotation <- subset(genome_annotation, dir==direction)
		if(direction=="fwd"){
			start_motifs <- dir_genome_annotation$position - (dir_genome_annotation$mod_pos - 1)
			end_motifs <- start_motifs + (dir_genome_annotation$motif_len - 1)
		}else{
			start_motifs <- dir_genome_annotation$position - (((dir_genome_annotation$motif_len - dir_genome_annotation$mod_pos) + 1) - 1)
			end_motifs <- start_motifs + (dir_genome_annotation$motif_len - 1)
		}
		strand_motifs <- ifelse(dir_genome_annotation$dir=="fwd","+","-")

		id_motifs <- paste0("ID=",dir_genome_annotation$motif,"_",dir_genome_annotation$id)
		name_motifs <- paste0("Name=",dir_genome_annotation$motif)
		color_motifs <- paste0("color=",dir_genome_annotation$col_motif)
		attributes_motifs <- paste0(id_motifs,";",name_motifs,";",color_motifs)
		
		gff_annotation <- paste(name_chr,".","motif",start_motifs,end_motifs,".",strand_motifs,".",attributes_motifs,sep="\t")

		path_output_file <- paste0(output_name,"_",direction,".gff")
		cat(paste0("##gff-version 3.2.1\n"), file=path_output_file)
		write.table(gff_annotation, file=path_output_file, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
}

save.stat.wig <- function(stat_data, base_output_name, stat_val, genome){
	name_chr <- names(readDNAStringSet(genome))
	
	if(length(grep("_test_",stat_val))==1){ # Need log10 transformation  # NEW
		stat_data <- stat_data %>%
			mutate_(score=paste0("-log10(",stat_val,")"))

		stat_data$score[stat_data$score>200] <- 200 # Capping
		stat_data$score <- round(stat_data$score) # Rounding
		param_wig <- paste0("track type=wiggle_0 windowingFunction=none featureVisibilityWindow=-1 viewLimits=0:200\n")
	}else if(length(grep("p",stat_val))==1){ # Need log10 transformation  # NEW
		stat_data <- stat_data %>%
			mutate_(score=paste0("-log10(",stat_val,")"))

		stat_data$score[stat_data$score>200] <- 200 # Capping
		stat_data$score <- round(stat_data$score) # Rounding
		param_wig <- paste0("track type=wiggle_0 windowingFunction=none featureVisibilityWindow=-1 viewLimits=0:200\n")
	}else{
		stat_data <- stat_data %>%
			mutate_(score=paste0(stat_val))

		stat_data$score <- round(stat_data$score, 2) # Rounding
		param_wig <- paste0("track type=wiggle_0 windowingFunction=none featureVisibilityWindow=-1 viewLimits=-9:15\n") #viewLimits=0:200
	}

	stat_data <- subset(stat_data, !is.na(score)) # Remove NA

	stifle <- foreach(direction=c("fwd","rev")) %do% {
		path_output_file <- paste0(base_output_name,"_",direction,".wig")
		cat(paste0(param_wig), file=path_output_file)
		cat(paste0("variableStep chrom=",name_chr,"\n"), file=path_output_file, append=TRUE) # TODO change chr if multiple contig
		write.table(subset(stat_data, select=c("position","score"), dir==direction), file=path_output_file, append=TRUE, row.names=FALSE, col.names=FALSE)
	}
}

#  ______ _____ _____    ___              _           _     
#  | ___ \  _  /  __ \  / _ \            | |         (_)    
#  | |_/ / | | | /  \/ / /_\ \_ __   __ _| |_   _ ___ _ ___ 
#  |    /| | | | |     |  _  | '_ \ / _` | | | | / __| / __|
#  | |\ \\ \_/ / \__/\ | | | | | | | (_| | | |_| \__ \ \__ \
#  \_| \_|\___/ \____/ \_| |_/_| |_|\__,_|_|\__, |___/_|___/
#                                            __/ |          
#                                           |___/           

evaluate.feature.ROC <- function(main_subset_stat_data, start, end, threshold, ROC_classes){
 # main_subset_stat_data <- subset(stat_data, dir=="fwd" & strand=="t")

 	processed_dir <- unique(main_subset_stat_data$dir)
	main_subset_ROC_classes <- subset(ROC_classes, dir==ifelse(processed_dir=="fwd","+","-")) # Extract corresponding dir of ROC_classes

	ROC_classes_motifs <- subset(main_subset_ROC_classes, grepl("^motif_NAT",main_subset_ROC_classes$feature)) # motif_NAT and motif_FAUX same ranges so only need one.
	ranges_motifs_features <- GRanges(
		seqnames="genome", # TODO change genome_annotation, add contig
		ranges=IRanges(start=ROC_classes_motifs$start, end=ROC_classes_motifs$end),
		strand=ROC_classes_motifs$dir # Is strand specific
	)
	ROC_classes_ctrl <- subset(main_subset_ROC_classes, grepl("^ctrl$",main_subset_ROC_classes$feature))
	ranges_ctrl_features <- GRanges(
		seqnames="genome", # TODO change genome_annotation, add contig
		ranges=IRanges(start=ROC_classes_ctrl$start, end=ROC_classes_ctrl$end),
		strand=ROC_classes_ctrl$dir # Is strand specific
	)

	subset_stat_data <- subset(main_subset_stat_data, select=c(position, dir), data_type=="NAT" & -log10(p)>=threshold)
	ranges_nat_passing_threshold <- GRanges(
		seqnames=rep("genome",nrow(subset_stat_data)), # TODO change main_subset_stat_data, add contig
		ranges=IRanges(start=subset_stat_data$position, end=subset_stat_data$position),
		strand=ifelse(subset_stat_data$dir=="fwd","+","-") # Is strand specific
	)
	subset_stat_data <- subset(main_subset_stat_data, select=c(position, dir), data_type=="FAUX" & -log10(p)>=threshold)
	ranges_faux_passing_threshold <- GRanges(
		seqnames=rep("genome",nrow(subset_stat_data)), # TODO change main_subset_stat_data, add contig
		ranges=IRanges(start=subset_stat_data$position, end=subset_stat_data$position),
		strand=ifelse(subset_stat_data$dir=="fwd","+","-") # Is strand specific
	)

	idx <- unique(findOverlaps(ranges_motifs_features,ranges_nat_passing_threshold)@from) # peaks are width 1, so type="within" not needed
	main_subset_ROC_classes$detected[grepl("^motif_NAT",main_subset_ROC_classes$feature)][idx] <- rep(TRUE,length(idx))

	idx <- unique(findOverlaps(ranges_motifs_features,ranges_faux_passing_threshold)@from)
	main_subset_ROC_classes$detected[grepl("^motif_FAUX",main_subset_ROC_classes$feature)][idx] <- rep(TRUE,length(idx))

	idx <- unique(findOverlaps(ranges_ctrl_features,ranges_nat_passing_threshold)@from)
	main_subset_ROC_classes$detected[grepl("^ctrl$",main_subset_ROC_classes$feature)][idx] <- rep(TRUE,length(idx))

	return(main_subset_ROC_classes)
}

summarize.feature.ROC <- function(detailed_parameters, stat_data, ROC_parameters, ROC_classes){
	threshold <- detailed_parameters[["threshold"]]
	# min_cov <- detailed_parameters[["min_cov"]]

	ROC_res <- ddply(stat_data, .(dir,strand), evaluate.feature.ROC, start=ROC_parameters$primary$start, end=ROC_parameters$primary$end, threshold=threshold, ROC_classes=ROC_classes, .parallel=FALSE)

	list_motif_feature <- levels(ROC_res$feature)[grepl("^motif_",levels(ROC_res$feature))]
	list_motif_feature <- unique(gsub("motif_NAT_|motif_FAUX_","",list_motif_feature))
	final_ROC_res <- foreach(motif=c("all",list_motif_feature), .combine=rbind) %do% {
		if(motif=="all"){
			subset_motif_ROC_res <- ROC_res[grepl("^motif_",ROC_res$feature),] # Keep all motifs; remove ctrl
		}else{
			subset_motif_ROC_res <- ROC_res[grepl(paste0("_",motif,"$"),ROC_res$feature),] # Keep one motif; remove ctrl
		}
		subset_ctrl_ROC_res <- ROC_res[grepl("^ctrl$",ROC_res$feature),]

		if(nrow(subset_motif_ROC_res)>0){
			NAT_res <- table(subset_motif_ROC_res[grepl("_NAT_",subset_motif_ROC_res$feature),]$detected)
			FAUX_res <- table(subset_motif_ROC_res[grepl("_FAUX_",subset_motif_ROC_res$feature),]$detected)
			ctrl_res <- table(subset_ctrl_ROC_res$detected)

			TP <- ifelse("TRUE" %in% names(NAT_res), NAT_res["TRUE"], 0)
			FN <- ifelse("FALSE" %in% names(NAT_res), NAT_res["FALSE"], 0)
			TN <- ifelse("FALSE" %in% names(FAUX_res), FAUX_res["FALSE"], 0)
			FP <- ifelse("TRUE" %in% names(FAUX_res), FAUX_res["TRUE"], 0)
			ctrl_TN <- ifelse("FALSE" %in% names(ctrl_res), ctrl_res["FALSE"], 0)
			ctrl_FP <- ifelse("TRUE" %in% names(ctrl_res), ctrl_res["TRUE"], 0)
			FP2 <- FP+ctrl_FP
			TN2 <- TN+ctrl_TN

			results <- data.frame(
				motif=motif, TP=TP, FN=FN, TN=TN, FP=FP, FP2=FP2, TN2=TN2
			)

			return(results)
		}else{
			results <- data.frame(
				motif=motif, TP=0, FN=0, TN=0, FP=0, FP2=0, TN2=0
			)

			return(results)
		}
	}

	return(final_ROC_res)
}

new.summarize.feature.ROC <- function(detailed_parameters, stat_data, ROC_parameters, ROC_classes, list_exception_motifs){
	threshold <- detailed_parameters[["threshold"]]
	# min_cov <- detailed_parameters[["min_cov"]]

	# Strand not used; TODO will need to integrate different strand if used (i.e. will return ROC_classes twice).
	# ROC_res2 <- ddply(stat_data, .(dir,strand), evaluate.feature.ROC, start=ROC_parameters$primary$start, end=ROC_parameters$primary$end, threshold=threshold, ROC_classes=ROC_classes, .parallel=FALSE)
	ROC_res <- stat_data %>% group_by(dir,strand) %>% do(evaluate.feature.ROC(.,start=ROC_parameters$primary$start, end=ROC_parameters$primary$end, threshold=threshold, ROC_classes=ROC_classes))

	list_motif_feature <- levels(ROC_res$feature)[grepl("^motif_",levels(ROC_res$feature))]
	list_motif_feature <- unique(gsub("motif_NAT_|motif_FAUX_","",list_motif_feature))
	final_ROC_res <- foreach(motif=c("all","trusted",list_motif_feature), .combine=rbind) %do% {
		if(motif=="all"){
			subset_motif_ROC_res <- ROC_res[grepl("^motif_",ROC_res$feature),] # Keep all motifs; remove ctrl
		}else if(motif=="trusted" & !any(is.na(list_exception_motifs))){
			subset_motif_ROC_res <- ROC_res[!gsub("motif_NAT_|motif_FAUX_","",ROC_res$feature) %in% list_exception_motifs,] # Remove unconfident motifs
		}else{
			subset_motif_ROC_res <- ROC_res[grepl(paste0("_",motif,"$"),ROC_res$feature),] # Keep one motif; remove ctrl
		}
		subset_ctrl_ROC_res <- ROC_res[grepl("^ctrl$",ROC_res$feature),]

		if(nrow(subset_motif_ROC_res)>0){
			NAT_query <- subset_motif_ROC_res[grepl("_NAT_",subset_motif_ROC_res$feature),]
			FAUX_query <- subset_motif_ROC_res[grepl("_FAUX_",subset_motif_ROC_res$feature),]

			# Matching control
			NAT_res <- table(NAT_query$detected)
			FAUX_res <- table(FAUX_query$detected)

			TP <- ifelse("TRUE" %in% names(NAT_res), NAT_res["TRUE"], 0)
			FN <- ifelse("FALSE" %in% names(NAT_res), NAT_res["FALSE"], 0)
			TN <- ifelse("FALSE" %in% names(FAUX_res), FAUX_res["FALSE"], 0)
			FP <- ifelse("TRUE" %in% names(FAUX_res), FAUX_res["TRUE"], 0)

			# Realistic control
			ctrl_res <- table(subset_ctrl_ROC_res$detected)
			ctrl_TN <- ifelse("FALSE" %in% names(ctrl_res), ctrl_res["FALSE"], 0)
			ctrl_FP <- ifelse("TRUE" %in% names(ctrl_res), ctrl_res["TRUE"], 0)
			TN2 <- TN+ctrl_TN
			FP2 <- FP+ctrl_FP

			# Realistic control with fixed ratio
			target_ratio <- 1/3
			nb_subsampling <- 100

			nb_pos <- nrow(NAT_query)
			nb_neg <- nrow(FAUX_query)
			nb_ctrl <- nrow(subset_ctrl_ROC_res)
			min_ratio <- nb_pos/(nb_pos + nb_neg + nb_ctrl)
			if(min_ratio<target_ratio){ # Random sampling of negative
				max_nb_neg <- nb_pos * 2
				combined_FAUXctrl_query <- rbind(FAUX_query, subset_ctrl_ROC_res)
				subsampling_res <- foreach(idx=seq(1, nb_subsampling), .combine=rbind) %do% {
					subsample_FAUXctrl_res <- table(combined_FAUXctrl_query$detected[sample(nrow(combined_FAUXctrl_query), max_nb_neg, replace=FALSE)])
					tmp_TN <- ifelse("FALSE" %in% names(subsample_FAUXctrl_res), subsample_FAUXctrl_res["FALSE"], 0)
					tmp_FP <- ifelse("TRUE" %in% names(subsample_FAUXctrl_res), subsample_FAUXctrl_res["TRUE"], 0)

					return(data.frame(tmp_TN=tmp_TN, tmp_FP=tmp_FP))
				} %>%
					summarize(mean_TN=mean(tmp_TN), sd_TN=sd(tmp_TN), mean_FP=mean(tmp_FP), sd_FP=sd(tmp_FP), .groups="drop_last")

				TP3 <- TP
				FN3 <- FN
				TN3 <- floor(subsampling_res$mean_TN)
				FP3 <- floor(subsampling_res$mean_FP)
				TP3_sd <- 0
				FN3_sd <- 0
				TN3_sd <- round(subsampling_res$sd_TN, 2)
				FP3_sd <- round(subsampling_res$sd_FP, 2)
			}else if(min_ratio>=target_ratio){ # Random sampling of positive
				max_nb_pos <- (nb_neg + nb_ctrl)/2
				subsampling_res <- foreach(idx=seq(1, nb_subsampling), .combine=rbind) %do% {
					subsample_NAT_res <- table(NAT_query$detected[sample(nrow(NAT_query), max_nb_pos, replace=FALSE)])
					tmp_TP <- ifelse("TRUE" %in% names(subsample_NAT_res), subsample_NAT_res["TRUE"], 0)
					tmp_FN <- ifelse("FALSE" %in% names(subsample_NAT_res), subsample_NAT_res["FALSE"], 0)

					return(data.frame(tmp_TP=tmp_TP, tmp_FN=tmp_FN))
				} %>%
					summarize(mean_TP=mean(tmp_TP), sd_TP=sd(tmp_TP), mean_FN=mean(tmp_FN), sd_FN=sd(tmp_FN), .groups="drop_last")

				TP3 <- floor(subsampling_res$mean_TP)
				FN3 <- floor(subsampling_res$mean_FN)
				TN3 <- TN2
				FP3 <- FP2
				TP3_sd <- round(subsampling_res$sd_TP, 2)
				FN3_sd <- round(subsampling_res$sd_FN, 2)
				TN3_sd <- 0
				FP3_sd <- 0
			}

			results <- data.frame(
				motif=motif, TP=TP, FN=FN, TN=TN, FP=FP, FP2=FP2, TN2=TN2, TP3=TP3, FN3=FN3, TN3=TN3, FP3=FP3, TP3_sd=TP3_sd, FN3_sd=FN3_sd, TN3_sd=TN3_sd, FP3_sd=FP3_sd
			)

			return(results)
		}else{
			results <- data.frame(
				motif=motif, TP=0, FN=0, TN=0, FP=0, FP2=0, TN2=0, TP3=0, FN3=0, TN3=0, FP3=0, TP3_sd=0, FN3_sd=0, TN3_sd=0, FP3_sd=0
			)

			return(results)
		}
	}

	return(final_ROC_res)
}

submeta.ROC <- function(parameters, stat_data, ROC_parameters, ROC_classes, list_exception_motifs){
	smoothing_win_size <- as.numeric(unique(parameters["smoothing_win_size"]))
	nbCPU <- ROC_parameters$primary$nbCPU

	stat_data <- stat_data[!is.na(stat_data$mean_diff),] # TODO what's the effect of smoothing?
	print(paste0("   Smoothing window: ",smoothing_win_size))
	registerDoMC(min(4,nbCPU))
	stat_data <- ddply(stat_data, .(dir, strand, data_type), rollingFunction, win_size=smoothing_win_size, stat_val=ROC_parameters$primary$stat_val, smooth_func=ROC_parameters$primary$smooth_func, .parallel=TRUE) # NEW
	registerDoSEQ()

	# Perform peak detection if TRUE
	if(ROC_parameters$primary$half_win_size>0){
		print(paste0("   Peaks detection: ",ROC_parameters$primary$half_win_size))
		registerDoMC(min(4,nbCPU))
		stat_data <- ddply(stat_data, .(dir, strand, data_type), mark.peaks, half_win_size=ROC_parameters$primary$half_win_size, .parallel=TRUE)
		registerDoSEQ()
	}

	# Choose threshold values | scale variate depending on smoothing
	nthreshold <- unique(parameters$threshold)
	threshold_range <- range(-log10(stat_data$p), na.rm=TRUE)
	step_size <- (threshold_range[2] - threshold_range[1])/nthreshold
	tmp_parameters <- as.list(parameters)
	detailed_parameters <- list()
	for(i in names(tmp_parameters)){
		detailed_parameters[[i]] <- unique(tmp_parameters[[i]])
	}
	detailed_parameters$threshold <- seq(threshold_range[1], threshold_range[2], by=step_size) # stat>=threshold so no need of max + step_size
	detailed_parameters <- expand.grid(detailed_parameters)

	if("peak" %in% colnames(stat_data)){ # If peak detection enable
		stat_data <- subset(stat_data, peak==TRUE) # Remove non-peak from dataset
	}

	print(paste0("      Summarize results per threshold."))
	registerDoMC(nbCPU)
	# res <- ddply(detailed_parameters, .(threshold), summarize.feature.ROC, stat_data=stat_data, ROC_parameters=ROC_parameters, ROC_classes=ROC_classes, .parallel=TRUE)
	res <- ddply(detailed_parameters, .(threshold), new.summarize.feature.ROC, stat_data=stat_data, ROC_parameters=ROC_parameters, ROC_classes=ROC_classes, list_exception_motifs=list_exception_motifs, .parallel=TRUE)
	registerDoSEQ()

	ROC_data <- merge(detailed_parameters, res) # Useless?

	return(ROC_data)
}

genome.annotation.ROC <- function(genome, motif_summary, iupac_nc, seq_params, start, end){
	# Find motifs ranges
	genome_annotation <- find.motifs(genome, motif_summary, iupac_nc) # TODO how duplicate are handled?
	genome_annotation <- subset(genome_annotation, position >= start & position < end)
	colnames(genome_annotation)[colnames(genome_annotation)=="motif"] <- "feature"
	genome_annotation$llen <- ifelse(genome_annotation$dir=="fwd",seq_params$fwd_llen,seq_params$rev_llen) # TODO Remove if change signal alignment
	genome_annotation$rlen <- ifelse(genome_annotation$dir=="fwd",seq_params$fwd_rlen,seq_params$rev_rlen)
	genome_annotation$start <- genome_annotation$position - genome_annotation$llen
	genome_annotation$end <- genome_annotation$position + genome_annotation$rlen
	genome_annotation <- subset(genome_annotation, select=-c(position,llen,rlen))
	genome_annotation$dir <- revalue(genome_annotation$dir, c("fwd"="+", "rev"="-"))

	# Prepare ROC classes for NAT (TP, FN) and FAUX (FP, TN)
	tmp <- genome_annotation
	tmp$feature <- paste0("motif_FAUX_",tmp$feature)
	genome_annotation$feature <- paste0("motif_NAT_",genome_annotation$feature)
	genome_annotation <- rbind(genome_annotation, tmp)
	genome_annotation$feature <- as.factor(genome_annotation$feature)
	genome_annotation$detected <- FALSE

	# Reorder columns
	genome_annotation <- subset(genome_annotation, select=c(start,end,dir,feature,detected))

	return(genome_annotation)
}

genome.ctrl.ROC <- function(genome_annotation, genome, seq_params, start, end){
	motifs_GR <- GRanges(seqnames="genome", # TODO change stat_data, add contig
						ranges=IRanges(start=genome_annotation$start, end=genome_annotation$end),
						strand=genome_annotation$dir) # Is strand specific
	empty_GR <- gaps(motifs_GR, start=1, end=width(readDNAStringSet(genome)))
	empty_GR <- empty_GR[strand(empty_GR)!="*"] # Remove gaps no strand entry
	genome_ctrl <- as.data.frame(empty_GR)
	genome_ctrl <- droplevels(genome_ctrl)
	genome_ctrl$start[genome_ctrl$start < start] <- start
	genome_ctrl$end[genome_ctrl$end > end] <- end
	genome_ctrl[genome_ctrl$end==end | genome_ctrl$start==start,] <- genome_ctrl[genome_ctrl$end==end | genome_ctrl$start==start,] %>%
		mutate(width=(end+1) - start)

	len_motif_range <- seq_params$fwd_llen + seq_params$fwd_rlen + 1
	tmp <- genome_ctrl %>%
		filter(width>=len_motif_range) %>%
		mutate(nbChunk=floor(width/len_motif_range), len_margin=width %% len_motif_range)

	set.seed(101)
	tmp$start_margin <- by(tmp, 1:nrow(tmp), function(x) sample(seq_len(x$len_margin+1),1) - 1)
	res <- by(tmp, 1:nrow(tmp), function(x) {
			vec <- seq(x$start+x$start_margin, x$end, by=len_motif_range-1)
			return(data.frame(start=vec[seq_along(vec[-length(vec)])], end=vec[seq_along(vec[-length(vec)])+1], dir=x$strand))
		}
	)
	genome_ctrl <- do.call(rbind.fill,res)
	genome_ctrl$feature <- as.factor("ctrl")
	genome_ctrl$detected <- FALSE

	# Reorder columns
	genome_ctrl <- subset(genome_ctrl, select=c(start,end,dir,feature,detected))

	return(genome_ctrl)
}

filter.not.queryable <- function(genome_combined, stat_data){
 # genome_combined <- ROC_classes

	query_gr <- GRanges(
		seqnames="genome", # TODO change stat_data, add contig
		ranges=IRanges(start=genome_combined$start, end=genome_combined$end),
		strand=genome_combined$dir # Is strand specific
	)

	stat_pos_data <- subset(stat_data, data_type=="NAT")
	stat_pos_gr <- GRanges(
		seqnames="genome", # TODO change stat_data, add contig
		ranges=IRanges(start=stat_pos_data$position, end=stat_pos_data$position),
		strand=ifelse(stat_pos_data$dir=="fwd","+","-") # Is strand specific
	)
	stat_neg_data <- subset(stat_data, data_type=="FAUX")
	stat_neg_gr <- GRanges(
		seqnames="genome", # TODO change stat_data, add contig
		ranges=IRanges(start=stat_neg_data$position, end=stat_neg_data$position),
		strand=ifelse(stat_neg_data$dir=="fwd","+","-") # Is strand specific
	)

	# tmp <- findOverlaps(query_gr, stat_pos_gr, type="any", select="all")
	genome_combined$nb_pos <- countOverlaps(query_gr, stat_pos_gr, type="any")
	genome_combined$nb_neg <- countOverlaps(query_gr, stat_neg_gr, type="any")
	genome_combined <- genome_combined %>%
		filter(nb_pos>5 & nb_neg>5) %>%
		dplyr::select(-c(nb_pos, nb_neg))

	return(genome_combined)
}

define.class.ROC <- function(genome, motif_summary, iupac_nc, ROC_parameters){
	seq_params <- ROC_parameters$seq_params
	start <- ROC_parameters$primary$start
	end <- ROC_parameters$primary$end

	genome_annotation <- genome.annotation.ROC(genome, motif_summary, iupac_nc, seq_params, start, end)

	genome_ctrl <- genome.ctrl.ROC(genome_annotation, genome, seq_params, start, end)

	genome_combined <- rbind(genome_annotation, genome_ctrl)

	return(genome_combined)
}

meta.ROC <- function(stat_pos_data, stat_neg_data, genome, motif_summary, iupac_nc, ROC_parameters, list_exception_motifs=NA){
	start <- ROC_parameters$primary$start
	end <- ROC_parameters$primary$end

	print(paste0("Select genomic subset from ",start," to ",end,"."))
	stat_pos_data$data_type <- "NAT"
	stat_neg_data$data_type <- "FAUX"
	stat_data <- rbind(stat_pos_data,stat_neg_data)
	stat_data <- subset(stat_data, position >= start & position < end)
	stat_data$data_type <- as.factor(stat_data$data_type)

	ROC_classes <- define.class.ROC(genome, motif_summary, iupac_nc, ROC_parameters)
	ROC_classes <- filter.not.queryable(ROC_classes, stat_data) # Remove region without signal from ROC_classes

	print("Processing sample with various smoothing and filtering.")
	ROC_data <- ddply(ROC_parameters$secondary, .(smoothing_win_size), submeta.ROC, stat_data=stat_data, ROC_parameters=ROC_parameters, ROC_classes=ROC_classes, list_exception_motifs=list_exception_motifs, .parallel=FALSE)

	# Add ploting information for motifs and total
	special_groups <- data.frame(motif=c("all","trusted"), mod_type=c(NA,NA), mod_pos=c(NA,NA), col_motif=c("#000000","#999999"), stringsAsFactors=FALSE)
	tmp_motif_summary <- rbind(special_groups, subset(motif_summary, select=c(motif, mod_type, mod_pos, col_motif)))
	ROC_data <- merge(tmp_motif_summary, ROC_data, by.x="motif", by.y="motif")

	print("Done.")
	return(ROC_data)
}

prepare.parameters.ROC <- function(smoothing_win_sizes, peak_filtering_radii, peak_widths, win_sizes, offsets, start_points, coverages=NA, thresholds=NA){
	if(is.na(min(coverages))){coverages <- 4} # Min produce NA if in list
	if(is.na(min(thresholds))){thresholds <- seq(0,32,0.1)} #0.1 32
	parameters <- expand.grid(list(smoothing_win_size=smoothing_win_sizes, peak_filtering_radius=peak_filtering_radii, 
		threshold=thresholds, peak_width=peak_widths, win_size=win_sizes, start_point=start_points, offset=offsets, min_cov=coverages))
}

compute.ROC.stats <- function(ROC_data){
	# Edge case:
		# PPV https://stats.stackexchange.com/questions/1773/what-are-correct-values-for-precision-and-recall-in-edge-cases
		# MCC https://lettier.github.io/posts/2016-08-05-matthews-correlation-coefficient.html
	# TPR: sensitivity
	# TNR: specificity

	scored_ROC_data <- ROC_data %>%
		group_by(motif) %>%
		mutate(FPR=FP/(FP + TN)) %>%
		mutate(TPR=TP/(TP + FN)) %>%
		mutate(TNR=TN/(TN + FP)) %>%
		mutate(FNR=1 - TPR) %>%
		mutate(PPV=ifelse(TP+FP==0,1,TP/(TP + FP))) %>% # With edge cases
		mutate(NPV=ifelse(TN+FN==0,1,TN/(TN + FN))) %>% # With edge cases
		mutate(FDR=1 - PPV) %>%
		mutate(FOR=1 - NPV) %>%
		mutate(ACC=(TP + TN)/(TP + TN + FP + FN)) %>%
		mutate(F1=((1 + 1^2)*TP)/((1 + 1^2)*TP + 1^2*(FP + FN))) %>%
		mutate(F2=((1 + 2^2)*TP)/((1 + 2^2)*TP + 2^2*(FP + FN))) %>%
		mutate(F05=((1 + 0.5^2)*TP)/((1 + 0.5^2)*TP + 0.5^2*(FP + FN))) %>%
		mutate(MCC=ifelse(TP+FP==0 | TN+FN==0,0,((TP*TN) - (FP*FN))/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)))) %>% # With edge cases
		mutate(BM=TPR + TNR - 1) %>%
		mutate(MK=PPV + NPV - 1) %>%
		mutate(PLR=TPR/FPR) %>%
		mutate(NLR=FNR/TNR) %>%
		mutate(DOR=PLR/NLR) %>%
		# With realistic control
		mutate(FPR2=FP2/(FP2 + TN2)) %>%
		mutate(TNR2=TN2/(TN2 + FP2)) %>%
		mutate(PPV2=ifelse(TP==0 & FP2==0,1,TP/(TP + FP2))) %>% # With edge cases
		mutate(NPV2=ifelse(TN2==0 & FN==0,1,TN2/(TN2 + FN))) %>% # With edge cases
		mutate(FDR2=1 - PPV2) %>%
		mutate(FOR2=1 - NPV2) %>%
		mutate(ACC2=(TP + TN2)/(TP + TN2 + FP2 + FN)) %>%
		mutate(F1_2=((1 + 1^2)*TP)/((1 + 1^2)*TP + 1^2*(FP2 + FN))) %>%
		mutate(F2_2=((1 + 2^2)*TP)/((1 + 2^2)*TP + 2^2*(FP2 + FN))) %>%
		mutate(F05_2=((1 + 0.5^2)*TP)/((1 + 0.5^2)*TP + 0.5^2*(FP2 + FN))) %>%
		mutate(MCC2=ifelse(TP+FP2==0 | TN2+FN==0,0,((TP*TN2) - (FP2*FN))/sqrt((TP + FP2)*(TP + FN)*(TN2 + FP2)*(TN2 + FN)))) %>% # With edge cases
		mutate(BM2=TPR + TNR2 - 1) %>%
		mutate(MK2=PPV2 + NPV2 - 1) %>%
		mutate(PLR2=TPR/FPR2) %>%
		mutate(NLR2=FNR/TNR2) %>%
		mutate(DOR2=PLR2/NLR2) %>%
		# With realistic control and frequency normalization
		mutate(FPR3=FP3/(FP3 + TN3)) %>%
		mutate(TPR3=TP3/(TP3 + FN3)) %>%
		mutate(TNR3=TN3/(TN3 + FP3)) %>%
		mutate(FNR3=1 - TPR3) %>%
		mutate(PPV3=ifelse(TP3==0 & FP3==0,1,TP3/(TP3 + FP3))) %>% # With edge cases
		mutate(NPV3=ifelse(TN3==0 & FN3==0,1,TN3/(TN3 + FN3))) %>% # With edge cases
		mutate(FDR3=1 - PPV3) %>%
		mutate(FOR3=1 - NPV3) %>%
		mutate(ACC3=(TP3 + TN3)/(TP3 + TN3 + FP3 + FN3)) %>%
		mutate(F1_3=((1 + 1^2)*TP3)/((1 + 1^2)*TP3 + 1^2*(FP3 + FN3))) %>%
		mutate(F2_3=((1 + 2^2)*TP3)/((1 + 2^2)*TP3 + 2^2*(FP3 + FN3))) %>%
		mutate(F05_3=((1 + 0.5^2)*TP3)/((1 + 0.5^2)*TP3 + 0.5^2*(FP3 + FN3))) %>%
		mutate(MCC3=ifelse(TP3+FP3==0 | TN3+FN3==0,0,((TP3*TN3) - (FP3*FN3))/sqrt((TP3 + FP3)*(TP3 + FN3)*(TN3 + FP3)*(TN3 + FN3)))) %>% # With edge cases
		mutate(BM3=TPR3 + TNR3 - 1) %>%
		mutate(MK3=PPV3 + NPV3 - 1) %>%
		mutate(PLR3=TPR3/FPR3) %>%
		mutate(NLR3=FNR3/TNR3) %>%
		mutate(DOR3=PLR3/NLR3)

	return(scored_ROC_data)
}

compute.AUC <- function(ROC_data){
	myAUC <- function(sample_ROC_data){
		sample_ROC_data <- sample_ROC_data %>%
			arrange(desc(threshold))
		
		# Matching control
		FPR <- sample_ROC_data$FPR
		TPR <- sample_ROC_data$TPR
		dFPR <- c(diff(FPR), 0)
		dTPR <- c(diff(TPR), 0)
		auroc <- sum(TPR * dFPR) + sum(dTPR * dFPR)/2

		PPV <- sample_ROC_data$PPV
		dPPV <- c(diff(PPV), 0)
		aupr <- sum(PPV * dTPR) + sum(dTPR * dPPV)/2

		# Realistic control
		FPR <- sample_ROC_data$FPR2
		TPR <- sample_ROC_data$TPR
		dFPR <- c(diff(FPR), 0)
		dTPR <- c(diff(TPR), 0)
		auroc2 <- sum(TPR * dFPR) + sum(dTPR * dFPR)/2

		PPV <- sample_ROC_data$PPV2
		dPPV <- c(diff(PPV), 0)
		aupr2 <- sum(PPV * dTPR) + sum(dTPR * dPPV)/2

		# Realistic control with fixed ratio
		FPR <- sample_ROC_data$FPR3
		TPR <- sample_ROC_data$TPR3
		dFPR <- c(diff(FPR), 0)
		dTPR <- c(diff(TPR), 0)
		auroc3 <- sum(TPR * dFPR) + sum(dTPR * dFPR)/2

		PPV <- sample_ROC_data$PPV3
		dPPV <- c(diff(PPV), 0)
		aupr3 <- sum(PPV * dTPR) + sum(dTPR * dPPV)/2

		return(data.frame(auroc=auroc, aupr=aupr, auroc2=auroc2, aupr2=aupr2, auroc3=auroc3, aupr3=aupr3))
	}

	if("dataset_name" %in% colnames(ROC_data)){
		# res <- ddply(ROC_data, .(dataset_name,min_cov,smoothing_win_size,motif), myAUC, .parallel=FALSE)
		res <- ROC_data %>% group_by(dataset_name,min_cov,smoothing_win_size,motif) %>% do(myAUC(.))
	}else{
		# res <- ddply(ROC_data, .(min_cov,smoothing_win_size,motif), myAUC, .parallel=FALSE)
		res <- ROC_data %>% group_by(min_cov,smoothing_win_size,motif) %>% do(myAUC(.))
	}
	# res <- res %>% spread(motif,auroc)

	return(res)
}

plot.ROC <- function(ROC_data, base_output_name, range_zoom){
	ROC_data <- ROC_data %>% arrange(threshold)

	threshold_palette <- colorRampPalette(brewer.pal(12, "Paired"))

	gp_th <- ggplot(subset(ROC_data, motif=="all")) +
		geom_point(aes(x=FPR, y=TPR, colour=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		scale_colour_gradientn(colours=threshold_palette(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		geom_abline(intercept=0, slope=1) + 
		theme(legend.justification=c(1,0), legend.position=c(1,0)) +
		labs(title="ROC curves", colour="Threshold", shape="Smoothing\nWin. Size") +
		labs(x="FPR", y="TPR")
	pdf(paste0("ROC_v1_",base_output_name,".pdf"))
	print(gp_th)
	dev.off()

	gp_th2 <- ggplot(subset(ROC_data, motif=="all")) +
		geom_point(aes(x=FPR2, y=TPR, colour=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		scale_colour_gradientn(colours=threshold_palette(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		geom_abline(intercept=0, slope=1) + 
		theme(legend.justification=c(1,0), legend.position=c(1,0)) +
		labs(title="ROC curves", colour="Threshold", shape="Smoothing\nWin. Size") +
		labs(x="FPR", y="TPR")
	pdf(paste0("ROC_v2_",base_output_name,".pdf"))
	print(gp_th2)
	dev.off()

	gp_zth <- ggplot(subset(ROC_data, motif=="all")) +
		geom_point(aes(x=FPR, y=TPR, colour=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		scale_colour_gradientn(colours=threshold_palette(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		geom_abline(intercept=0, slope=1) + 
		theme(legend.justification=c(1,0), legend.position=c(1,0)) +
		labs(title="ROC curves", subtitle=paste0("Zoomed "), colour="Threshold", shape="Smoothing\nWin. Size") + #TODO precise zoom
		labs(x="FPR", y="TPR") +
		coord_cartesian(xlim=range_zoom$x,ylim=range_zoom$y)
	pdf(paste0("ROC_zoom_",base_output_name,".pdf"))
	print(gp_zth)
	dev.off()

	labels <- unique(subset(ROC_data, select=c("motif","col_motif")))
	labels[] <- lapply(labels, as.character) # convert to character

	gp_motif <- ggplot(subset(ROC_data, smoothing_win_size==0)) +
		geom_point(aes(x=FPR, y=TPR, colour=motif), size=1) +
		geom_abline(intercept=0, slope=1) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		theme(legend.justification=c(1,0), legend.position=c(1,0)) +
		labs(title="ROC curves by motifs", subtitle="No smoothing") +
		labs(x="FPR", y="TPR")
	pdf(paste0("ROC_detail_v1_",base_output_name,".pdf"))
	print(gp_motif)
	dev.off()

	gp_motif2 <- ggplot(subset(ROC_data, smoothing_win_size==0)) +
		geom_point(aes(x=FPR2, y=TPR, colour=motif), size=1) +
		geom_abline(intercept=0, slope=1) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		theme(legend.justification=c(1,0), legend.position=c(1,0)) +
		labs(title="ROC curves by motifs", subtitle="No smoothing") +
		labs(x="FPR", y="TPR")
	pdf(paste0("ROC_detail_v2_",base_output_name,".pdf"))
	print(gp_motif2)
	dev.off()

	gp_motif2 <- ggplot(ROC_data) +
		geom_point(aes(x=FPR2, y=TPR, colour=motif), size=1) +
		geom_abline(intercept=0, slope=1) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		facet_wrap(~ smoothing_win_size) +
		labs(title="ROC curves by motifs") +
		labs(x="FPR", y="TPR")
	pdf(paste0("ROC_detail_v3_",base_output_name,".pdf"), width=15)
	print(gp_motif2)
	dev.off()

	gp_motif2 <- ggplot(ROC_data) +
		geom_point(aes(x=FPR2, y=TPR, colour=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		scale_colour_gradientn(colours=threshold_palette(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		geom_abline(intercept=0, slope=1) +
		facet_wrap(~ motif) +
		labs(title="ROC curves by motifs") +
		labs(x="FPR", y="TPR", colour="Threshold", shape="Smoothing\nWin. Size")
	pdf(paste0("ROC_detail_v4_",base_output_name,".pdf"), height=15, width=16)
	print(gp_motif2)
	dev.off()

	gp_zmotif <- ggplot(subset(ROC_data, smoothing_win_size==0)) +
		geom_point(aes(x=FPR, y=TPR, colour=motif), size=1) +
		geom_abline(intercept=0, slope=1) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		theme(legend.justification=c(1,0), legend.position=c(1,0)) +
		labs(title="ROC curves by motifs", subtitle=paste0("Zoomed & No smoothing")) + #TODO precise zoom
		labs(x="FPR", y="TPR") +
		coord_cartesian(xlim=range_zoom$x,ylim=range_zoom$y)
	pdf(paste0("ROC_detail_zoom_",base_output_name,".pdf"))
	print(gp_zmotif)
	dev.off()

	gp_PR <- ggplot(ROC_data) +
		geom_path(aes(x=TPR, y=PPV2, col=as.factor(smoothing_win_size)), size=0.5, alpha=0.7) +
		facet_wrap(~motif) +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
		labs(title="PR curves", subtitle="By Smoothing Win. Size", colour="Smoothing\nWin. Size") +
		labs(x="TPR", y="PPV2")
	pdf(paste0("ROC_PRxSmooth_curves_",base_output_name,".pdf"), width=10, height=10)
	print(gp_PR)
	dev.off()

	gp_PR <- ggplot(ROC_data) +
		geom_point(aes(x=TPR, y=PPV2, col=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		facet_wrap(~motif) +
		scale_colour_gradientn(colours=colorRampPalette(brewer.pal(12, "Paired"))(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
		labs(title="PR curves", subtitle="By Threshold", colour="Threshold", shape="Smoothing\nWin. Size") +
		labs(x="TPR", y="PPV2")
	pdf(paste0("ROC_PRxThreshold_curves_",base_output_name,".pdf"), width=10, height=10)
	print(gp_PR)
	dev.off()

	gp_motif2 <- ggplot(ROC_data) +
		geom_point(aes(x=FPR2, y=BM2, colour=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		scale_colour_gradientn(colours=threshold_palette(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		facet_wrap(~ motif) +
		labs(colour="Threshold", shape="Smoothing\nWin. Sizee") +
		labs(x="FPR2", y="BM2")
	pdf(paste0("ROC_BMxFPR_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_motif2)
	dev.off()

	gp_motif3 <- ggplot(ROC_data) +
		geom_point(aes(x=FPR2, y=MCC2, colour=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		scale_colour_gradientn(colours=threshold_palette(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		facet_wrap(~ motif) +
		labs(colour="Threshold", shape="Smoothing\nWin. Size") +
		labs(x="FPR2", y="MCC2")
	pdf(paste0("ROC_MCCxFPR_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_motif3)
	dev.off()

	gp_motif2 <- ggplot(ROC_data) +
		geom_path(aes(x=threshold, y=BM2, colour=as.factor(smoothing_win_size)), size=0.5, alpha=0.7) +
		facet_wrap(~ motif) +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		labs(colour="Smoothing\nWin. Size") +
		labs(x="Threshold", y="BM2")
	pdf(paste0("ROC_BMxThreshold_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_motif2)
	dev.off()

	gp_motif3 <- ggplot(ROC_data) +
		geom_path(aes(x=threshold, y=MCC2, colour=as.factor(smoothing_win_size)), size=0.5, alpha=0.7) +
		facet_wrap(~ motif) +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		labs(colour="Smoothing\nWin. Size") +
		labs(x="Threshold", y="MCC2")
	pdf(paste0("ROC_MCCxThreshold_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_motif3)
	dev.off()

	res <- compute.AUC(ROC_data)
	gp_auc <- ggplot(res) +
		geom_path(aes(x=smoothing_win_size, y=auroc, col=motif)) +
		geom_point(aes(x=smoothing_win_size, y=auroc, col=motif), size=1) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		labs(title="Area under ROC curves by motifs", subtitle=paste0("AUROC by Smoothing Win. Size (0 for no smoothing)")) +
		labs(x="Smoothing Win. Size", y="AUROC")
	pdf(paste0("ROC_AUROC_",base_output_name,".pdf"), height=10, width=10)
	print(gp_auc)
	dev.off()
	gp_auc <- ggplot(res) +
		geom_path(aes(x=smoothing_win_size, y=aupr, col=motif)) +
		geom_point(aes(x=smoothing_win_size, y=aupr, col=motif), size=1) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		labs(title="Area under PR curves by motifs", subtitle=paste0("AUPR by Smoothing Win. Size (0 for no smoothing)")) +
		labs(x="Smoothing Win. Size", y="AUPR")
	pdf(paste0("ROC_AUPR_",base_output_name,".pdf"), height=10, width=10)
	print(gp_auc)
	dev.off()

	tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$BM2),])
	res <- tmp_ROC_data[tmp_ROC_data[, .I[BM2==max(BM2)], by=list(motif,smoothing_win_size)]$V1] # Can have multiple row occurrences, TODO min?
	gp_bm <- ggplot(res) +
		geom_path(aes(x=smoothing_win_size, y=BM2, col=motif)) +
		geom_point(aes(x=smoothing_win_size, y=BM2, col=motif)) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		labs(title="Informedness by motifs", subtitle=paste0("BM by Smoothing Win. Size (0 for no smoothing)")) +
		labs(x="Smoothing Win. Size", y="BM")
	pdf(paste0("ROC_BMxSmooth_",base_output_name,".pdf"), height=10, width=10)
	print(gp_bm)
	dev.off()
	gp_bm <- ggplot(res) +
		geom_point(aes(x=threshold, y=BM2, colour=as.factor(smoothing_win_size))) +
		facet_grid(motif~., scales="free_y") +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		theme(strip.text.y=element_text(angle=0)) +
		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
		labs(title="Informedness by motifs", subtitle=paste0("BM by Best Threshold (0 for no smoothing)")) +
		labs(x="Best Threshold", y="BM", colour="Smoothing\nWin. Size")
	pdf(paste0("ROC_BMxThreshold_",base_output_name,".pdf"), height=15, width=7)
	print(gp_bm)
	dev.off()

	tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$MCC2),])
	res <- tmp_ROC_data[tmp_ROC_data[, .I[MCC2==max(MCC2)], by=list(motif,smoothing_win_size)]$V1] # Can have multiple row occurrences, TODO min?
	gp_mcc <- ggplot(res) +
		geom_path(aes(x=smoothing_win_size, y=MCC2, col=motif)) +
		geom_point(aes(x=smoothing_win_size, y=MCC2, col=motif)) +
		scale_colour_manual(name="Colors:", breaks=labels$motif, values=labels$col_motif, labels=labels$motif) +
		labs(title="Matthews correlation coefficient by motifs", subtitle=paste0("MCC by Smoothing Win. Size (0 for no smoothing)")) +
		labs(x="Smoothing Win. Size", y="MCC")
	pdf(paste0("ROC_MCCxSmooth_",base_output_name,".pdf"), height=10, width=10)
	print(gp_mcc)
	dev.off()
	gp_mcc <- ggplot(res) +
		geom_point(aes(x=threshold, y=MCC2, colour=as.factor(smoothing_win_size))) +
		facet_grid(motif~., scales="free_y") +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		theme(strip.text.y=element_text(angle=0)) +
		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
		labs(title="Matthews correlation coefficient by motifs", subtitle=paste0("MCC by Best Threshold (0 for no smoothing)")) +
		labs(x="Best Threshold", y="MCC", colour="Smoothing\nWin. Size")
	pdf(paste0("ROC_MCCxThreshold_",base_output_name,".pdf"), height=15, width=7)
	print(gp_mcc)
	dev.off()

	gp_motif3 <- ggplot(ROC_data) +
		geom_path(aes(x=threshold, y=FDR2, colour=as.factor(smoothing_win_size)), size=0.5, alpha=0.7) +
		facet_wrap(~ motif) +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		labs(title="FDR curves",colour="Smoothing\nWin. Size") +
		labs(x="Threshold", y="FDR2")
	pdf(paste0("ROC_FDRxThreshold_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_motif3)
	dev.off()

	gp_Fscores <- ggplot(ROC_data) +
		geom_path(aes(x=threshold, y=F05_2, col=as.factor(smoothing_win_size)), size=0.5, alpha=0.7) +
		facet_wrap(~motif) +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
		labs(title=expression(F[0.5]~scores), subtitle="By Smoothing Win. Size", colour="Smoothing\nWin. Size") +
		labs(x="Threshold", y=expression(F[0.5]))
	pdf(paste0("ROC_F05xThreshold_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_Fscores)
	dev.off()
	gp_Fscores <- ggplot(ROC_data) +
		geom_path(aes(x=threshold, y=F1_2, col=as.factor(smoothing_win_size)), size=0.5, alpha=0.7) +
		facet_wrap(~motif) +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
		labs(title=expression(F[1]~scores), subtitle="By Smoothing Win. Size", colour="Smoothing\nWin. Size") +
		labs(x="Threshold", y=expression(F[1]))
	pdf(paste0("ROC_F1xThreshold_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_Fscores)
	dev.off()

	gp_Fscores <- ggplot(ROC_data) +
		geom_point(aes(x=TP, y=F05_2, col=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		facet_wrap(~motif, scales="free_x") +
		scale_colour_gradientn(colours=colorRampPalette(brewer.pal(12, "Paired"))(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
		labs(title=expression(F[0.5]~scores), subtitle="By Threshold", colour="Threshold", shape="Smoothing\nWin. Size") +
		labs(x="TP", y=expression(F[0.5]))
	pdf(paste0("ROC_F05xTP_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_Fscores)
	dev.off()
	gp_Fscores <- ggplot(ROC_data) +
		geom_point(aes(x=TP, y=F1_2, col=threshold, shape=as.factor(smoothing_win_size)), size=1, alpha=0.5) +
		facet_wrap(~motif, scales="free_x") +
		scale_colour_gradientn(colours=colorRampPalette(brewer.pal(12, "Paired"))(100)) +
		scale_shape_manual(values=seq_along(levels(as.factor(ROC_data$smoothing_win_size)))) +
		theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
		labs(title=expression(F[1]~scores), subtitle="By Threshold", colour="Threshold", shape="Smoothing\nWin. Size") +
		labs(x="TP", y=expression(F[1]))
	pdf(paste0("ROC_F1xTP_curves_",base_output_name,".pdf"), height=10, width=10)
	print(gp_Fscores)
	dev.off()

	tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$F05_2),])
	res <- tmp_ROC_data[tmp_ROC_data[, .I[F05_2==max(F05_2)], by=list(motif,smoothing_win_size)]$V1] # Can have multiple row occurrences, TODO min?
	gp_mcc <- ggplot(res) +
		geom_point(aes(x=threshold, y=F05_2, colour=as.factor(smoothing_win_size))) +
		facet_grid(motif~., scales="free_y") +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		theme(strip.text.y=element_text(angle=0)) +
		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
		labs(title=expression(F[0.5]~scores~by~motifs), subtitle=expression(F[0.5]*" by Best Threshold (0 for no smoothing)")) +
		labs(x="Best Threshold", y=expression(F[0.5]), colour="Smoothing\nWin. Size")
	pdf(paste0("ROC_F05xThreshold_",base_output_name,".pdf"), height=15, width=7)
	print(gp_mcc)
	dev.off()
	tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$F1_2),])
	res <- tmp_ROC_data[tmp_ROC_data[, .I[F1_2==max(F1_2)], by=list(motif,smoothing_win_size)]$V1] # Can have multiple row occurrences, TODO min?
	gp_mcc <- ggplot(res) +
		geom_point(aes(x=threshold, y=F1_2, colour=as.factor(smoothing_win_size))) +
		facet_grid(motif~., scales="free_y") +
		scale_colour_brewer(palette="Spectral", direction=-1) +
		theme(strip.text.y=element_text(angle=0)) +
		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
		labs(title=expression(F[1]~scores~by~motifs), subtitle=expression(F[1]*" by Best Threshold (0 for no smoothing)")) +
		labs(x="Best Threshold", y=expression(F[1]), colour="Smoothing\nWin. Size")
	pdf(paste0("ROC_F1xThreshold_",base_output_name,".pdf"), height=15, width=7)
	print(gp_mcc)
	dev.off()
}

merge.ROC <- function(list_ROC_data, list_ROC_name){
	merged_ROC_data <- foreach(ROC_data_idx=seq(1,length(list_ROC_data)), .combine=rbind) %do% {
		ROC_data <- list_ROC_data[[ROC_data_idx]]
		ROC_data$dataset_name <- list_ROC_name[ROC_data_idx]
		ROC_data <- as.data.frame(ROC_data) # Avoid unnecessary data.table information

		return(ROC_data)
	}
	merged_ROC_data$dataset_name <- as.factor(merged_ROC_data$dataset_name) # After rbind to avoid conflicting factor levels
	merged_ROC_data$dataset_name <- ordered(merged_ROC_data$dataset_name, levels=list_ROC_name)

	return(merged_ROC_data)
}

compare.ROC <- function(base_output_name, list_ROC_data, list_ROC_name, list_exception_motifs=NA){
	ROC_data <- merge.ROC(list_ROC_data, list_ROC_name)

	# AUCs
	res <- compute.AUC(ROC_data) %>% arrange(smoothing_win_size)

	# AUCs Matching control
	{
		gp_auc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=auroc, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=auroc, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Area under ROC curves by motifs", subtitle=paste0("AUROC by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="AUROC", colour="Type dataset:") +
			theme(legend.position="bottom", legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_AUROC_",base_output_name,".pdf"), height=10, width=13)
		print(gp_auc)
		dev.off()
		gp_auc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=aupr, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=aupr, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Area under PR curves by motifs", subtitle=paste0("AUPR by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="AUPR", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_AUPR_",base_output_name,".pdf"), height=10, width=13)
		print(gp_auc)
		dev.off()
	}

	# AUCs Realistic control
	{
		gp_auc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=auroc2, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=auroc2, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Area under ROC curves by motifs", subtitle=paste0("AUROC by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="AUROC2", colour="Type dataset:") +
			theme(legend.position="bottom", legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_AUROC2_",base_output_name,".pdf"), height=10, width=13)
		print(gp_auc)
		dev.off()
		gp_auc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=aupr2, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=aupr2, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Area under PR curves by motifs", subtitle=paste0("AUPR by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="AUPR2", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_AUPR2_",base_output_name,".pdf"), height=10, width=13)
		print(gp_auc)
		dev.off()
	}

	# AUCs Realistic control with fixed ratio
	{
		gp_auc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=auroc3, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=auroc3, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Area under ROC curves by motifs", subtitle=paste0("AUROC by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="AUROC3", colour="Type dataset:") +
			theme(legend.position="bottom", legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_AUROC3_",base_output_name,".pdf"), height=10, width=13)
		print(gp_auc)
		dev.off()
		gp_auc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=aupr3, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=aupr3, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Area under PR curves by motifs", subtitle=paste0("AUPR by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="AUPR3", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_AUPR3_",base_output_name,".pdf"), height=10, width=13)
		print(gp_auc)
		dev.off()
	}

	# Matching control
	{
		tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$BM2),])
		res <- tmp_ROC_data[tmp_ROC_data[, .I[BM2==max(BM2)], by=list(dataset_name,motif,smoothing_win_size)]$V1] %>%
			arrange(smoothing_win_size)
		gp_bm <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=BM2, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=BM2, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Informedness by motifs", subtitle=paste0("BM by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="BM", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_BM_",base_output_name,".pdf"), height=10, width=13)
		print(gp_bm)
		dev.off()

		tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$MCC2),])
		res <- tmp_ROC_data[tmp_ROC_data[, .I[MCC2==max(MCC2)], by=list(dataset_name,motif,smoothing_win_size)]$V1] %>%
			arrange(smoothing_win_size)
		gp_mcc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=MCC2, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=MCC2, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Matthews correlation coefficient by motifs", subtitle=paste0("MCC by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="MCC", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_MCC_",base_output_name,".pdf"), height=10, width=13)
		print(gp_mcc)
		dev.off()

		res <- ROC_data %>%
			filter(smoothing_win_size==5) %>%
			arrange(threshold)
		gp_pr <- ggplot(res) +
			geom_path(aes(x=TPR, y=PPV, col=as.factor(dataset_name)), size=0.5, alpha=0.7) +
			facet_wrap(~motif) +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
			labs(title="PR curves", subtitle="By Data Types", colour="Data\nTypes") +
			labs(x="TPR", y="PPV")
		pdf(paste0("ROC_summary_PR_curves_",base_output_name,".pdf"), width=10, height=10)
		print(gp_pr)
		dev.off()
	}
	
	# Realistic control
	{
		tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$BM2),])
		res <- tmp_ROC_data[tmp_ROC_data[, .I[BM2==max(BM2)], by=list(dataset_name,motif,smoothing_win_size)]$V1] %>%
			arrange(smoothing_win_size)
		gp_bm <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=BM2, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=BM2, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Informedness by motifs", subtitle=paste0("BM by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="BM2", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_BM2_",base_output_name,".pdf"), height=10, width=13)
		print(gp_bm)
		dev.off()

		tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$MCC2),])
		res <- tmp_ROC_data[tmp_ROC_data[, .I[MCC2==max(MCC2)], by=list(dataset_name,motif,smoothing_win_size)]$V1] %>%
			arrange(smoothing_win_size)
		gp_mcc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=MCC2, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=MCC2, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Matthews correlation coefficient by motifs", subtitle=paste0("MCC by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="MCC2", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_MCC2_",base_output_name,".pdf"), height=10, width=13)
		print(gp_mcc)
		dev.off()

		res <- ROC_data %>%
			filter(smoothing_win_size==5) %>%
			arrange(threshold)
		gp_pr <- ggplot(res) +
			geom_path(aes(x=TPR, y=PPV2, col=as.factor(dataset_name)), size=0.5, alpha=0.7) +
			facet_wrap(~motif) +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
			labs(title="PR curves", subtitle="By Data Types", colour="Data\nTypes") +
			labs(x="TPR", y="PPV2")
		pdf(paste0("ROC_summary_PR2_curves_",base_output_name,".pdf"), width=10, height=10)
		print(gp_pr)
		dev.off()
	}

	# Realistic control with fixed ratio
	{
		tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$BM3),])
		res <- tmp_ROC_data[tmp_ROC_data[, .I[BM3==max(BM3)], by=list(dataset_name,motif,smoothing_win_size)]$V1] %>%
			arrange(smoothing_win_size)
		gp_bm <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=BM3, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=BM3, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Informedness by motifs", subtitle=paste0("BM by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="BM3", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_BM3_",base_output_name,".pdf"), height=10, width=13)
		print(gp_bm)
		dev.off()

		tmp_ROC_data <- as.data.table(ROC_data[!is.na(ROC_data$MCC3),])
		res <- tmp_ROC_data[tmp_ROC_data[, .I[MCC3==max(MCC3)], by=list(dataset_name,motif,smoothing_win_size)]$V1] %>%
			arrange(smoothing_win_size)
		gp_mcc <- ggplot(res) +
			geom_path(aes(x=smoothing_win_size, y=MCC3, col=dataset_name)) +
			geom_point(aes(x=smoothing_win_size, y=MCC3, col=dataset_name), size=1) +
			facet_wrap(~motif, scales="free_y") +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			scale_x_continuous(breaks=c(0,seq(3,21,4)), minor_breaks=c(0,seq(3,21,2))) +
			labs(title="Matthews correlation coefficient by motifs", subtitle=paste0("MCC by Smoothing Win. Size (0 for no smoothing)")) +
			labs(x="Smoothing Win. Size", y="MCC3", colour="Type dataset:") +
			theme(legend.position="bottom",legend.direction="horizontal") +
			# theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
			guides(colour=guide_legend(nrow=1))
		pdf(paste0("ROC_summary_MCC3_",base_output_name,".pdf"), height=10, width=13)
		print(gp_mcc)
		dev.off()

		res <- ROC_data %>%
			filter(smoothing_win_size==5) %>%
			arrange(threshold)
		gp_pr <- ggplot(res) +
			geom_path(aes(x=TPR3, y=PPV3, col=as.factor(dataset_name)), size=0.5, alpha=0.7) +
			facet_wrap(~motif) +
			scale_colour_brewer(palette="Spectral", direction=-1) +
			theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
			labs(title="PR curves", subtitle="By Data Types", colour="Data\nTypes") +
			labs(x="TPR3", y="PPV3")
		pdf(paste0("ROC_summary_PR3_curves_",base_output_name,".pdf"), width=10, height=10)
		print(gp_pr)
		dev.off()

		myPalette <- subset(res, select=c(motif,col_motif)) %>% unique()
		named_color_vector <- myPalette$col_motif
		names(named_color_vector) <- myPalette$motif
		gp_pr <- ggplot(res) +
			geom_path(aes(x=TPR3, y=PPV3, col=motif), size=0.5, alpha=0.7) +
			facet_wrap(~dataset_name) +
			scale_colour_manual(values=named_color_vector) +
			theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
			labs(title="PR curves", subtitle="By Motif", colour="Motifs") +
			labs(x="TPR3", y="PPV3") +
			guides(colour=guide_legend(ncol=3))
		pdf(paste0("ROC_summary_PR3_curves_byMotif_",base_output_name,".pdf"), width=10, height=10)
		print(gp_pr)
		dev.off()

		if(!is.na(list_exception_motifs[1])){
			gp_pr <- ggplot(subset(res, ! motif %in% list_exception_motifs)) +
				geom_path(aes(x=TPR3, y=PPV3, col=motif), size=0.5, alpha=0.7) +
				facet_wrap(~dataset_name) +
				scale_colour_manual(values=named_color_vector) +
				theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal") +
				labs(title="PR curves", subtitle="By Motif", colour="Motifs") +
				labs(x="TPR3", y="PPV3") +
				guides(colour=guide_legend(ncol=3))
			pdf(paste0("ROC_summary_PR3_curves_byMotif2_",base_output_name,".pdf"), width=10, height=10)
			print(gp_pr)
			dev.off()
		}
	}
}







