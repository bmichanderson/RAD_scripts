
##########
# Author: B.M. Anderson
# Date: April 2022; May 2025 (cleaned up and adjusted calculations)
# Description: assess tree height and population divergences for setting priors for SNAPP/ER
# Note: the input Nexus file can be produced with `vcf_to_nexus.py`
##########


# load required library (for reading in Nexus)
suppressMessages(library(ape))


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to assess divergences in a SNAPP/ER input Nexus file\n")
		cat("Usage: Rscript SNAPP_divergence.R -s {SNAPP/ER_infile.nex} -p {pops_file}\n")
		cat("Options:\n")
		cat("\t-s\tSNAPP/ER input Nexus file\n")
		cat("\t-p\tPopulations file, with sampleID and population tab-separated, one per line\n")
	} else {
		cat(help_message)
	}
}


# A function to calculate pairwise sequence divergence from 012 data (biallelic)
# Note: this is calculated only on variable sites, SNPs (if using that), so values are high
# It compares each row to every other row, and sums absolute distances
# Genotypes are: 0 = 00, 1 = 01, 2 = 11
# Distances are between phased alleles, so there are 2 comparisons per site
# For example, if sample1 has "1", then sample2 having "0", "1", and "2" = distances of 0+1, 0+0, and 1+0
# Missing data ("?") is ignored (sites dropped for that comparison)
# It then calculates substitutions per site by dividing by number of comparisons (2 per site)
# The argument is a matrix with samples as rows and SNPs as columns
# Returns a (lower tri) matrix with substitutions per site between samples (rows and columns)
sequence_diversity <- function(snp_mat) {
	missing <- vector("numeric")
	out_mat <- matrix(nrow = nrow(snp_mat), ncol = nrow(snp_mat))
	rownames(out_mat) <- rownames(snp_mat)
	colnames(out_mat) <- rownames(snp_mat)
	index <- 2
	for (row in seq_len(nrow(snp_mat) - 1)) {
		for (comp_row in seq(index, nrow(snp_mat))) {
			not_miss <- snp_mat[row, ] != "?" & snp_mat[comp_row, ] != "?"
			missing <- append(missing, (ncol(snp_mat) - sum(not_miss)) / ncol(snp_mat))
			comps <- snp_mat[c(row, comp_row), not_miss]
			distance <- sum(abs(as.numeric(comps[1, ]) - as.numeric(comps[2, ])))
			out_mat[comp_row, row] <-  distance / (2 * ncol(comps))
		}
		index <- index + 1
	}
	cat(paste0("Average missing comparisons proportion: ", round(mean(missing, na.rm = TRUE), 3), "\n"))
	out_mat
}


# a function to calculate lambda (used for SNAPP/ER)
# Lambda is a parameter for the Yule Prior for the species tree and branch lengths
# This equation is from "A Rough Guide to SNAPP"
# lambda = 1/h * (sum(k=1 to n-1) k/(n*(n - k)))
# h = expected(prior) tree height, n = number of species
get_lambda <- function(height, num) {
	values <- vector("numeric")
	for (kval in seq_len(num - 1)) {
		values <- append(values, kval / (num * (num - kval)))
	}
	sum(values, na.rm = TRUE) / height
}


# Read in the data

# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	snapp_present <- FALSE
	pops_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-p") {
			pops_present <- TRUE
			pops_file <- args[index + 1]
		} else if (args[index] == "-s") {
			snapp_present <- TRUE
			snapp_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! snapp_present || ! pops_present) {
	stop(help("Missing argument for SNAPP/ER file and/or populations file!\n"), call. = FALSE)
}


# read in the NEXUS SNP matrix (automatically converted to a list indexed by sample)
mynex <- read.nexus.data(snapp_file)
mymat <- matrix(unlist(mynex), nrow = length(mynex), ncol = length(mynex[[1]]), byrow = TRUE)
rownames(mymat) <- names(mynex)


# read in the populations file as a table
pops_table <- read.table(pops_file, header = FALSE, sep = "\t")
colnames(pops_table) <- c("sample", "pop")


# calculate pairwise sequence divergence as substitutions per site
# this calculation takes the longest (but not too bad)
subs_mat <- sequence_diversity(mymat)


# report possible lambda values for different numbers of species (2–15)
# uses maximum divergence divided by two (height of the tree)
height <- round(max(subs_mat, na.rm = TRUE) / 2, 4)
cat(paste0("\nExpected tree height: ", height, " substitutions per site\n"))
cat(paste0("Average divergence between all samples: ",
	round(mean(subs_mat, na.rm = TRUE), 4)), "substitutions per site\n")
cat("\nLambda (parameter for Yule Prior for species tree) by number of species:\n")
cat("Species\tLambda\n")
for (num in seq(2, 15)) {
	lambda <- round(get_lambda(height, num), 1)
	cat(paste0(num, "\t", lambda, "\n"))
}


# report average divergences within populations
cat("\nAverage substitutions per site between members of a population:\n")
pops_to_subset <- unique(pops_table$pop[pops_table$sample %in% rownames(mymat)])

result_df <- as.data.frame(matrix(nrow = length(pops_to_subset), ncol = 1))
rownames(result_df) <- pops_to_subset
colnames(result_df) <- c("Avg_intra_subs")

for (pop in pops_to_subset) {
	samples <- pops_table$sample[pops_table$pop == pop]
	comp_mat <- subs_mat[samples, samples]
	result_df[pop, ] <- round(mean(comp_mat, na.rm = TRUE), 3)
}
result_df
cat(paste0("Min:\t", min(result_df, na.rm = TRUE), "\n"))
cat(paste0("Avg:\t", round(mean(result_df[, 1], na.rm = TRUE), 3), "\n"))
cat(paste0("Max:\t", max(result_df, na.rm = TRUE), "\n"))


# report a distance matrix between populations and average between pop distances
cat("\nAverage substitutions per site between populations:\n")
result_mat <- matrix(NA, nrow = length(pops_to_subset), ncol = length(pops_to_subset))
rownames(result_mat) <- pops_to_subset
colnames(result_mat) <- pops_to_subset
index <- 2
for (pop_num in seq_len(length(pops_to_subset) - 1)) {
	for (pop2_num in seq(index, length(pops_to_subset))) {
		samples1 <- pops_table$sample[pops_table$pop == pops_to_subset[pop_num]]
		samples2 <- pops_table$sample[pops_table$pop == pops_to_subset[pop2_num]]
		comp_mat1 <- subs_mat[samples1, samples2]
		comp_mat2 <- subs_mat[samples2, samples1]
		result_mat[pops_to_subset[pop2_num], pops_to_subset[pop_num]] <- round(mean(c(comp_mat1, comp_mat2), na.rm = TRUE), 3)
	}
	index <- index + 1
}
result_mat
cat(paste0("Min:\t", min(result_mat, na.rm = TRUE), "\n"))
cat(paste0("Avg:\t", round(mean(result_mat, na.rm = TRUE), 3), "\n"))
cat(paste0("Max:\t", max(result_mat, na.rm = TRUE), "\n"))
