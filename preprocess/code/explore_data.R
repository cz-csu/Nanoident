load.libraries.extract <- function(){
	library("rhdf5")
	library("foreach")
	library("doMC")
	# library("doFuture")
	library("optparse")
	# library("progressr")
	library("progress")
	library("Biostrings")
	library("stringr")
}
suppressMessages(load.libraries.extract())
h5ls("dataset/unzip_data/MinION_BA_NAT/0/DESKTOP_T1OHN23_20171108_FAH33570_MN20624_mux_scan_Pool_Mock9b_NAT_43075_read_4_ch_404_strand.fast5")