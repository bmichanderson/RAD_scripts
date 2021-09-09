
#############
# Author: B. Anderson based partly on O'Leary et al. 2018 DOI: 10.1111/mec.14792
# Date: Sep 2021
# Description: write a file of individuals (one per line) that do not pass a loci missing filter for a VCF
#############


## define a helper function for calls
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to return a file listing individuals that do not pass a filter\n")
		cat("Usage: Rscript indv_miss.R -l limit -p prefix\n")
		cat("Option:\n")
		cat("\t-l\tThe max missing limit to filter on\n")
		cat("\t-p\tThe file name prefix used when running vcftools\n")
        cat("\t\tNOTE: this script expects the following file to be present in the working directory:\n")
        cat("\t\t*.imiss\n")
	} else {
		cat(help_message)
	}
}


## parse command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
    prefix_present <- FALSE
    limit_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-l") {
			limit <- args[index + 1]
            limit_present <- TRUE
		} else if (args[index] == "-p") {
            prefix <- args[index + 1]
            prefix_present <- TRUE
        } else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
    }
}
if (!all(prefix_present, limit_present)) {
    stop(help(), call. = FALSE)
}


## load file and filter
imiss <- read.table(paste0(prefix, ".imiss"), header = TRUE, stringsAsFactors = FALSE)
to_remove <- imiss[imiss$F_MISS > limit, ]$INDV
write.table(to_remove, paste0(prefix, ".to_remove"), col.names = FALSE, row.names = FALSE, quote = FALSE)
