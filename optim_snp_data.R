##########
# Author: B. Anderson
# Date: Nov 2023
# Description: evaluate multiple SNP datasets to optimise clustering threshold
# Note: the minimum arguments are multiple VCF files (one per parameter combination)
#	It is helpful to also indicate the association between the files and the parameter values (-p),
#	specify which sample IDs are replicates (-r), a samples and populations file to limit which
#	to analyse (-s), as well as	a custom stats file from ipyrad (if you used that) for additional output (-i)
##########


# load libraries
suppressMessages(library(adegenet))
suppressMessages(library(vcfR))


##########
# define functions

## a help function for when the script is called incorrectly or without arguments
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to evaluate SNP datasets called under different parameter settings\n\n")
		cat("Usage: Rscript optim_snp_data.R <optional args> VCF_file1 VCF_file2 ...\n")
		cat("Arguments:\n")
		cat("\t...\tPaths to VCF files in the order for plotting/comparing (e.g. sequentially increasing)\n")
		cat("\t-i\tA custom-formatted ipyrad stats file with a header and columns containing specific values:\n")
		cat("\t\tcolumn1 = sample ID; column2 = clustering threshold; column3 = autosomal heterozygosity\n")
		cat("\t-p\tA text file with parameter values for each input VCF in the order provided,",
			"one per line (e.g. \"85\" or \"M2\")\n")
		cat("\t\tIf not provided, the parameter values will be displayed arbitrarily\n")
		cat("\t-r\tA text file with sample IDs of replicate pairs (tab separated, one per line)\n")
		cat("\t-s\tA text file with sample IDs and populations (tab separated, one per line)\n")
		cat("\t\tNote: you can optionally include decimal latitude and longitude as the third and fourth columns\n")
		cat("\t\tIf a samples file is provided, only the specified samples in the VCF files will be analysed\n")
	} else {
		cat(help_message)
	}
}


## a function to read in a VCF file and convert to a genlight
## returns a genlight
vcftogenl <- function(vcf_file) {
	vcf <- read.vcfR(vcf_file, verbose = FALSE)
	suppressWarnings(vcfR2genlight(vcf))
}


## a function to compute the number of loci shared by >=80% of samples
## -- based on the metric used in Paris et al. (2017) --
## returns a number
paris <- function(genl) {
	mat <- as.matrix(genl)
	keep_snps <- colSums(! is.na(mat)) / nrow(mat) >= 0.80
	sub_genl <- genl[, keep_snps]
	length(unique(sub_genl@chromosome))
}


## functions to run error rate assessments
## -- based on metrics used in Mastretta-Yanes et al. (2015) --
## returns a dataframe with columns corresponding to locus, allele and snp error rates
mastretta_error <- function(genl, reps) {
	error_df <- as.data.frame(matrix(nrow = nrow(reps), ncol = 3))
	colnames(error_df) <- c("locus", "allele", "snp")
	for (pair_row in seq_len(nrow(reps))) {
		pair_genl <- genl[match(c(reps[pair_row, ]), genl@ind.names)]
		# locus
		# loci found in one rep but not the other divided by total potential shared
		na_pos <- NA.posi(pair_genl)				# NA positions in the genlight
		na_rep <- na_pos[[1]]						# NA positions in the rep
		rep_genl <- pair_genl[1, -na_rep]			# remove NA
		rep_loci <- unique(rep_genl@chromosome)		# find unique loci names
		na_samp <- na_pos[[2]]						# NA positions in the samp
		samp_genl <- pair_genl[2, -na_samp]			# remove NA
		samp_loci <- unique(samp_genl@chromosome)	# find unique loci names
		diffs <- length(setdiff(rep_loci, samp_loci)) + length(setdiff(samp_loci, rep_loci))
		in_both <- intersect(rep_loci, samp_loci)
		locus_error <- 100 * diffs / (diffs + length(in_both))	# diffs divided by possible comparisons

		# allele
		# differences between loci as a whole
		chrom_list <- rep(NA, length = length(in_both))
		allele_genl <- pair_genl[, pair_genl@chromosome %in% in_both]
		comp_mat <- as.matrix(allele_genl)
		comp_mat[is.na(comp_mat)] <- 3			# change NA to a number that won't match
		index <- 1
		allele_errors <- 0
		for (snp in seq_len(ncol(comp_mat))) {
			if (as.character(allele_genl@chromosome[snp]) %in% chrom_list) {	# if we have seen this locus already
				if (! hit) {	# if we haven't found an error yet
					if (comp_mat[1, snp] != comp_mat[2, snp]) {
						allele_errors <- allele_errors + 1
						hit <- TRUE
					}
				}
			} else {	# new locus
				chrom_list[index] <- as.character(allele_genl@chromosome[snp])	# add it to the list
				index <- index + 1
				if (comp_mat[1, snp] != comp_mat[2, snp]) {
					allele_errors <- allele_errors + 1
					hit <- TRUE
				} else {
					hit <- FALSE
				}
			}
		}
		allele_error <- 100 * allele_errors / (length(in_both))

		# SNP
		# differences between every shared SNP
		comp_mat <- comp_mat[, colSums(comp_mat == 3) == 0]		# completely remove any NA comparisons
		snp_errors <- 0
		for (snp in seq_len(ncol(comp_mat))) {
			if (comp_mat[1, snp] != comp_mat[2, snp]) {
				snp_errors <- snp_errors + 1
			}
		}
		snp_error <- 100 * snp_errors / ncol(comp_mat)

		# capture the values
		error_df[pair_row, ] <- c(locus_error, allele_error, snp_error)
	}
	error_df
}


## functions to run intrapopulation distances
## -- based on a metric used in Mastretta-Yanes et al. (2015) --
## returns a dataframe with columns of number of individuals and average distances
mastretta_dist <- function(genl, pop_table) {
	dist_obj <- dist(as.matrix(genl))		# choose a better metric?
	# normalize the distances
	dist_mat <- as.matrix(dist_obj / max(dist_obj))
	# extract distances between individuals of the same population
	pops <- unique(pop_table[, 2])
	dists_df <- as.data.frame(matrix(nrow = length(pops), ncol = 2))
	colnames(dists_df) <- c("number_ind", "mean_intra_dist")
	for (popnum in seq_len(length(pops))) {
		pop <- pops[popnum]
		pos <- pop_table[, 2] %in% pop			# find the positions which match that pop
		dmat <- dist_mat[pos, pos]			# extract a smaller matrix of samples from that pop
		dists <- dmat[lower.tri(dmat)]			# keep the lower triangle values
		avg_dist <- mean(dists)				# calculate an average within population distance
		dists_df[popnum, ] <- c(nrow(dmat), avg_dist)
	}
	dists_df
}




###########
# process the call
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	extra <- 1
	catch <- TRUE
	ipyrad_present <- FALSE
	ipyrad_file <- ""
	params_present <- FALSE
	params_file <- ""
	reps_present <- FALSE
	reps_file <- ""
	samples_present <- FALSE
	samples_file <- ""
	for (index in seq_len(length(args))) {
		if (args[index] == "-i") {
			ipyrad_present <- TRUE
			ipyrad_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-p")  {
			params_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-r")  {
			reps_present <- TRUE
			reps_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-s")  {
			samples_present <- TRUE
			samples_file <- args[index + 1]
			catch <- FALSE
		} else {
			if (catch) {
				catch_args[extra] <- args[index]
				extra <- extra + 1
			} else {
				catch <- TRUE
			}
		}
	}
}
if (length(catch_args) < 2) {
	stop(help("Missing VCF files!\n"), call. = FALSE)
}


# read in the input
if (ipyrad_present) {
	het_stats <- read.table(ipyrad_file, sep = "\t", header = TRUE)
}

if (params_present) {
	params <- as.character(read.table(params_file, header = FALSE)[, 1])
} else {
	params <- unlist(lapply(seq_len(length(catch_args)), function(x) c("P", x)))
}

if (reps_present) {
	reps_table <- read.table(reps_file, sep = "\t", header = FALSE)
}

if (samples_present) {
	samples_table <- read.table(samples_file, sep = "\t", header = FALSE)
	sampleids <- as.character(samples_table[, 1])
	pops <- unique(as.character(samples_table[, 2]))
	if (ncol(samples_table) == 4) {
		locations_present <- TRUE
	} else {
		locations_present <- FALSE
	}
}





# execute the analyses





# plot

if (ipyrad_present) {
	boxplot(het_stats[, 3] * 100 ~ het_stats[, 2],
		main = "Autosomal heterozygosity\n(ambiguous bases in consensus sequences within samples)",
		ylab = "Heterozygous sites (%)", xlab = "Clustering threshold (%)")
}
