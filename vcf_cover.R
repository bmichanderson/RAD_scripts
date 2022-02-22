
##########
# Author: Ben Anderson
# Date: Feb 2022
# Description: report the number of loci and SNPs present in a VCF (arg1) in >= percent (arg2) samples
# Note: this needs a VCF from ipyrad, with naming convention for loci as "locus_position"
##########


# parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	stop("Missing arguments for VCF file and percent cover", call. = FALSE)
} else {
	vcf_file <- args[1]
	percent <- as.numeric(args[2])
}


# load the library, analyse the VCF, and report
suppressMessages(library(vcfR))
myvcf <- read.vcfR(vcf_file, verbose = FALSE)
gt <- extract.gt(myvcf, element = "GT")
greater <- apply(gt, MARGIN = 1, function(x) { 100 * sum(!is.na(x)) / ncol(gt) >= percent })
nsnps <- sum(greater)
nloci <- length(unique(sub("_.*", "", (rownames(as.data.frame(greater[greater]))))))
cat(paste0("The VCF contains ", nloci, " loci and ", nsnps, " SNPs in >= ", percent, "% of samples\n"))
