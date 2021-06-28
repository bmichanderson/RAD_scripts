
##########
# Authors: Rachel Binks and Ben Anderson, based on scripts from Mastretta-Yanes et al. 2015
# Date: June 2021
# Description: run Mastretta-Yanes et al. error estimation for technical replicates
##########


###############
############### NOTE: We need to figure out how to get a locus/allele error rate file set up
###############




# set parameters needed for output

#param1 <- X		# param1 description


# a helper function for when the script is called without arguments
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to run Mastretta-Yanes et al. error estimation for technical",
			"replicates using a vcf file\n")
		cat("Usage: Rscript MY_replicate_error.R -o output -s sample_file -v vcf_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name [default output.txt]\n")
		cat("\t-s\tThe sample names file if vcf labels differ, with the format:\n",
			"\t\tvcf_label\t*tab*\tsample_name\n",
			"\t\tNote: the name should have an \"_R\" at the end for a replicate\n")
		cat("\t-v\tThe vcf file to be analysed\n")
	} else {
		cat(help_message)
	}
}


# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1

	outfile <- "output.txt"
	samples_present <- FALSE
	vcf_present <- FALSE

	for (index in 1:length(args)) {
		if (args[index] == "-o") {
			outfile <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present = TRUE
			vcf_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! samples_present) {
	cat("No argument for sample names, so we assume duplicates are named with an \"_R\"\n")
}

if (! vcf_present) {
	stop(help("Missing argument for vcf file!\n"), call. = FALSE)
}


# load required libraries
library(vcfR)
library(adegenet)


# read in the input files
if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
}

vcf <- read.vcfR(vcf_file, verbose = FALSE)
gl <- vcfR2genlight(vcf)		# convert to genlight object


# assign correct names to the genlight object, if needed
if (samples_present) {
	indNames(gl) <- sample_table$V2[match(indNames(gl), sample_table$V1)]
}


# NOTE: I don't think we should do any filtering for missing data, as part of the error
# test is to see whether loci are called in one replicate but not the other


# determine names of replicate pairs and combine into a dataframe
reps <- grep("_R$", indNames(gl), value = TRUE)
samps <- indNames(gl)[match(sub("_R$", "", reps), indNames(gl))]
pairs <- cbind(reps, samps)
pairs <- pairs[rowSums(is.na(pairs)) < 1, ]	# remove rows with missing pairs
npairs <- nrow(pairs)


######## Loci and allele error rates VS. SNP error rates
# Locus error rate: if a locus is in one replicate but not the other,
# allele error rate: if the locus is in both but differs by at least one SNP, and
# SNP error rate: if a SNP is found in both and differs (there may be multiple SNPs per locus)

# If we use VCF input, then we can only get locus error rate if there is one SNP per locus, otherwise
# we are measuring SNP 'missingness'; also
# allele error rate == SNP error rate if only one SNP has been output per locus

