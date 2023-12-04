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
		cat("\t\tIf a samples file is provided, only the specified samples in the VCF files will be analysed\n")
	} else {
		cat(help_message)
	}
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


## a function (and subfunctions) to run error rate assessments
## -- based on metrics used in Mastretta-Yanes et al. (2015) --
## arguments are a genlight and a table with the replicate pair sampleIDs, one pair per row (2 columns)
## returns a dataframe with columns corresponding to locus, allele and snp error rates
mastretta_error <- function(genl, reps) {
	### a function to grab a list of loci present (at least one SNP called) for a sample index in a genlight
	getloci <- function(genl, sampind) {
		na_positions <- NA.posi(genl)[[sampind]]
		unique(genl[sampind, -na_positions]@chromosome)
	}

	### a function to find differences between two samples for every locus
	### the genlight should only have two samples
	### returns the number of loci that differ
	compare_alleles <- function(genl) {
		allele_errors <- 0
		for (chrom in unique(genl@chromosome)) {
			comp_mat <- as.matrix(genl[, genl@chromosome == chrom])
			comp_mat[is.na(comp_mat)] <- 3		# change to a number that won't match
			for (col in seq_len(ncol(comp_mat))) {
				if (comp_mat[1, col] != comp_mat[2, col]) {
					allele_errors <- allele_errors + 1
					break
				}
			}
		}
		allele_errors
	}

	### a function to find differences between two samples for every SNP
	### the genlight should only have two samples
	### returns the number of called SNPs that differ, and the number called in both
	compare_snps <- function(genl) {
		snp_errors <- 0
		comp_mat <- as.matrix(genl)
		comp_mat <- comp_mat[, colSums(is.na(comp_mat)) == 0]	# remove SNPs with at least one uncalled
		for (col in seq_len(ncol(comp_mat))) {
			if (comp_mat[1, col] != comp_mat[2, col]) {
				snp_errors <- snp_errors + 1
			}
		}
		c(snp_errors, ncol(comp_mat))
	}

	### make a dataframe for output
	error_df <- as.data.frame(matrix(nrow = nrow(reps), ncol = 3))
	colnames(error_df) <- c("locus", "allele", "snp")

	### cycle through the replicate pairs
	for (pair_row in seq_len(nrow(reps))) {
		#### subset the genl for that pair
		pair_genl <- genl[match(c(reps[pair_row, ]), genl@ind.names)]

		#### locus error: loci found in one rep but not the other divided by total potential shared
		rep_loci <- getloci(pair_genl, 1)
		samp_loci <- getloci(pair_genl, 2)
		num_diffs <- length(setdiff(rep_loci, samp_loci)) + length(setdiff(samp_loci, rep_loci))
		in_both <- intersect(rep_loci, samp_loci)
		locus_error <- 100 * num_diffs / (num_diffs + length(in_both))

		#### allele error: whether shared loci differ in sequence (any SNP, including no call in one)
		share_genl <- pair_genl[, pair_genl@chromosome %in% in_both]
		allele_errors <- compare_alleles(share_genl)
		allele_error <- 100 * allele_errors / (length(in_both))

		#### snp error: differences between any shared SNP divided by total shared
		snp_output <- compare_snps(share_genl)
		snp_errors <- snp_output[1]
		snps_shared <- snp_output[2]
		snp_error <- 100 * snp_errors / snps_shared

		#### capture the values
		error_df[pair_row, ] <- c(locus_error, allele_error, snp_error)
	}
	error_df
}


## a function to evaluate within-population distances
## -- based on the suggestion in Mastretta-Yanes et al. (2015) --
## returns a dataframe with columns of number of individuals and average distances
pop_dist <- function(vcfr, sampleids, pops) {
	### load libraries
	suppressMessages(library(ape))
	suppressMessages(library(pofadinr))

	### calculate distances for the VCF samples
	dnabin <- suppressMessages(vcfR2DNAbin(vcfr, consensus = TRUE, extract.haps = FALSE))
	temp <- as.character(dnabin)
	temp[temp == "n"] <- "?"	# replace missing with a character recognised by pofadinr
	dnabin <- as.DNAbin(temp)
	distances <- dist.snp(dnabin, model = "GENPOFAD")
	dist_mat <- as.matrix(distances)

	### extract distances between individuals of the same population
	pop_names <- unique(pops)
	dists_df <- as.data.frame(matrix(nrow = length(pop_names), ncol = 2))
	colnames(dists_df) <- c("number_ind", "mean_intra_dist")
	for (popnum in seq_len(length(pop_names))) {
		pop <- pop_names[popnum]
		samples <- sampleids[pops == pop]
		sub_mat <- dist_mat[rownames(dist_mat) %in% samples, colnames(dist_mat) %in% samples]
		avg_dist <- mean(sub_mat)
		dists_df[popnum, ] <- c(nrow(sub_mat), avg_dist)
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
	reps_table[] <- lapply(reps_table, as.character)
}

if (samples_present) {
	samples_table <- read.table(samples_file, sep = "\t", header = FALSE)
	sampleids <- as.character(samples_table[, 1])
	pops <- as.character(samples_table[, 2])		# corresponding sequence to samples (not unique)
}


##########
# execute the analyses

## create the dataframes for holding the results


## cycle through the VCF files
for (vcf_file in catch_args) {
	vcfr <- read.vcfR(vcf_file, verbose = FALSE)
	genl <- suppressMessages(vcfR2genlight(vcfr))

	## if there is a samples file, subset the genl
	## also, calculate within-population distances
	if (samples_present) {
		genl <- genl[match(sampleids, genl@ind.names)]
		pop_dists <- pop_dist(vcfr, sampleids, pops)
	}

	## calculate Paris
	paris_num <- paris(genl)

	## if there is a replicates file, calculate error rates
	if (reps_present) {
		merror_df <- mastretta_error(genl, reps_table)
	}
}




# If there is a samples file, plot within-population distances
if (samples_present) {
######
	boxplot(100 * mydist_table$V3 ~ mydist_table$V1,
		main = "Within Population Distances",
		ylab = "Average GENPOFAD distance"
	)
}



##############################

if (merror_present) {
	# plot Mastretta-Yanes errors
	pdf(file = paste0(outpre, "_merror.pdf"))
	myerror_table <- read.table(merror_file, sep = "\t", header = FALSE)
	boxplot(myerror_table$V2 ~ myerror_table$V1,
		main = paste0("Locus error rates\n",
			"(loci present in one rep but not the other, relative to total in either/both)"),
		ylab = "Locus error rate (%)", xlab = "Clustering threshold (%)")
	boxplot(myerror_table$V3 ~ myerror_table$V1,
		main = "Allele error rates\n(loci present in both that differ, relative to total in both)",
		ylab = "Allele error rate (%)", xlab = "Clustering threshold (%)")
	boxplot(myerror_table$V4 ~ myerror_table$V1,
		main = "SNP error rates\n(SNPs called in both that differ, relative to total called in both)",
		ylab = "SNP error rate (%)", xlab = "Clustering threshold (%)")
	invisible(dev.off)
}

if (mparis_present) {
	# plot Paris et al counts of loci and SNPs in >= 80% of individuals
	pdf(file = paste0(outpre, "_mparis.pdf"))
	count_table <- read.table(mparis_file, sep = "\t", header = FALSE)
	plot(count_table$V1, count_table$V2, main = "Total number of loci recovered in >= 80% of individuals",
		ylab = "Number of loci", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(count_table$V1), max(count_table$V1), by = 1))
	plot(count_table$V1, count_table$V3, main = "Total number of SNPs recovered in >= 80% of individuals",
		ylab = "Number of SNPs", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(count_table$V1), max(count_table$V1), by = 1))
	invisible(dev.off())
}


############################


# plot

if (ipyrad_present) {
	boxplot(het_stats[, 3] * 100 ~ het_stats[, 2],
		main = "Autosomal heterozygosity\n(ambiguous bases in consensus sequences within samples)",
		ylab = "Heterozygous sites (%)", xlab = "Clustering threshold (%)")
}
