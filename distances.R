
##########
# Author: Ben Anderson
# Date: Nov 2021
# Modified: Oct 2023 (removed plotting and simplified)
# Description: calculate and output distances between samples from a VCF file
##########


# load required libraries
suppressMessages(library(adegenet))
suppressMessages(library(ape))
suppressMessages(library(pofadinr))
suppressMessages(library(vcfR))


# Define a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate distances from an input VCF file, then output the matrix\n")
		cat("Usage: Rscript distances.R -d dist_method -o output -v vcf_file [-s samples]\n")
		cat("Options:\n")
		cat("\t-d\tThe distance method: \"G\"ENPOFAD [default], \"E\"uclidean, \"M\"ATCHSTATES\n")
		cat("\t-o\tThe output file name prefix [default \"output\"]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-s\tA file with sampleID and sample name (tab separated, one per line) [optional]\n")
	} else {
		cat(help_message)
	}
}


# Parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	dist_method <- 1
	output <- "output"
	vcf_present <- FALSE
	samples_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-d") {
			dist_method <- args[index + 1]
		} else if (args[index] == "-o") {
			output <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present <- TRUE
			vcf_file <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! vcf_present) {
	stop(help("Missing argument for vcf file!\n"), call. = FALSE)
}


# read in the input
vcf <- read.vcfR(vcf_file, verbose = FALSE)
cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
	length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")

if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
}


# calculate distances between samples
if (dist_method == "G") {
	dnabin <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE)
	# we need to change the default behaviour of converting missing into "n"!!!!!
	temp <- as.character(dnabin)
	temp[temp == "n"] <- "?"
	dnabin <- as.DNAbin(temp)
	individuals <- rownames(dnabin)
	## GENPOFAD (allows better comparison between hets and homo; ambiguity codes)
	distance <- dist.snp(dnabin, model = "GENPOFAD")
	dist_suffix <- "_distGENPOFAD.nex"
	method <- paste0("GENPOFAD distances from ", dim(dnabin)[2], " SNPs")
} else if (dist_method == "E") {
	genl <- vcfR2genlight(vcf)
	individuals <- indNames(genl)
	## Euclidean
	distance <- dist(as.matrix(genl))
	dist_suffix <- "_distEuclidean.nex"
	method <- paste0("Euclidean distances from ", nLoc(genl), " SNPs")
} else if (dist_method == "M") {
	dnabin <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE)
	# we need to change the default behaviour of converting missing into "n"!!!!!
	temp <- as.character(dnabin)
	temp[temp == "n"] <- "?"
	dnabin <- as.DNAbin(temp)
	individuals <- rownames(dnabin)
	## MATCHSTATES (another way to use ambiguity codes)
	distance <- dist.snp(dnabin, model = "MATCHSTATES")
	dist_suffix <- "_distMATCHSTATES.nex"
	method <- paste0("MATCHSTATES distances from ", dim(dnabin)[2], " SNPs")
} else {
	stop(help("Distance method incorrectly specified! Quitting...\n"), call. = FALSE)
}


# replace sample names if provided
if (samples_present) {
	for (index in seq_len(length(individuals))) {
		id <- individuals[index]
		if (id %in% sample_table$V1) {
			individuals[index] <- sample_table$V2[match(id, sample_table$V1)]
		} else {
			cat(paste0("Sample ", id, " not found in sample table\n"))
		}
	}

	temp_mat <- as.matrix(distance)
	rownames(temp_mat) <- individuals
	distance <- as.dist(temp_mat)
}


# Output the distance matrix
taxa <- individuals
taxa_block <- paste0("BEGIN TAXA;\n\tDIMENSIONS NTAX=", length(taxa), ";\n\t",
	"TAXLABELS ", paste(taxa, collapse = " "), ";\nEND;\n")
dist_block <- paste0("BEGIN DISTANCES;\n\tFORMAT\n\t\tTRIANGLE=BOTH\n\t\tDIAGONAL\n\t\t",
	"LABELS=LEFT\n\t;\n\tMATRIX\n")
outfile <- file(paste0(output, dist_suffix), open = "w")
writeLines("#NEXUS", con = outfile)
writeLines(taxa_block, con = outfile)
writeLines(dist_block, con = outfile)
write.table(as.matrix(distance), file = outfile, col.names = FALSE,
	append = TRUE, quote = FALSE)
writeLines("\t;\nEND;\n", con = outfile)
close(outfile)
