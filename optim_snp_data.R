##########
# Author: B. Anderson
# Date: Nov-Dec 2023
# Description: evaluate multiple SNP datasets to optimise assembly parameters
# Note: the minimum arguments are multiple VCF files (one per parameter combination)
#	It is helpful to also indicate the association between the files and the parameter values (-p),
#	specify which sample IDs are replicate pairs (-r), a samples and populations file to limit which
#	to analyse (-s), as well as a custom stats file from ipyrad (if you used that) for autosomal output (-i)
# Note2: this script relies on genlight objects, which only uses biallelic SNPs (may underestimate heterozygosity)
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
		cat("\t-o\tOutput prefix [default: \"out\"]\n")
		cat("\t-p\tA text file with parameter values for each input VCF in the order provided,",
			"one per line (e.g. \"85\" or \"M2\")\n")
		cat("\t\tIf not provided, the parameter values will be displayed arbitrarily\n")
		cat("\t-r\tA text file with sample IDs of replicate pairs (tab separated, one per line)\n")
		cat("\t-s\tA text file with sample IDs and populations (tab separated, one per line)\n")
		cat("\t\tIf a samples file is provided, only those specified samples will be analysed\n")
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


## a function to run error rate assessments
## -- based on metrics used in Mastretta-Yanes et al. (2015) --
## arguments are a genlight and a table with the replicate pair sampleIDs, one pair per row (2 columns)
## returns a matrix with columns corresponding to locus, allele and snp error rates
mastretta_error <- function(genl, reps) {
	errors <- matrix(nrow = nrow(reps), ncol = 3)
	cat(paste0("Comparing replicate pair (of ", nrow(reps), "): "))
	for (pair_row in seq_len(nrow(reps))) {
		cat(paste0(pair_row, " "))
		### subset the genl for that pair and make a matrix
		pair_genl <- genl[match(c(reps[pair_row, ]), genl@ind.names)]
		comp_mat <- as.matrix(pair_genl)
		na1 <- is.na(comp_mat[1, ])
		na2 <- is.na(comp_mat[2, ])
		na_both <- na1 & na2
		pair_genl <- pair_genl[, !na_both]	# remove NA in both
		comp_mat <- as.matrix(pair_genl)
		na1 <- is.na(comp_mat[1, ])
		na2 <- is.na(comp_mat[2, ])

		### locus error: loci found in one rep but not the other, divided by total potential shared
		loci1 <- unique(pair_genl[1, !na1]@chromosome)
		loci2 <- unique(pair_genl[2, !na2]@chromosome)
		diffs <- length(setdiff(loci1, loci2)) + length(setdiff(loci2, loci1))
		in_both <- intersect(loci1, loci2)
		locus_error <- 100 * diffs / (diffs + length(in_both))

		### allele error: whether shared loci differ in sequence (any SNP, including no call in one)
		pair_genl <- pair_genl[, pair_genl@chromosome %in% in_both]
		comp_mat <- as.matrix(pair_genl)
		comp_mat[is.na(comp_mat)] <- 3		# so NA can be compared
		diff_pos <- comp_mat[1, ] != comp_mat[2, ]
		allele_errors <- length(unique(pair_genl[, diff_pos]@chromosome))		# so that multiple in an allele aren't counted
		allele_error <- 100 * allele_errors / (length(in_both))

		### snp error: differences between any shared SNP divided by total shared
		comp_mat <- comp_mat[, colSums(comp_mat == 3) == 0]		# remove any missing
		diffs <- sum(comp_mat[1, ] != comp_mat[2, ])
		snp_error <- 100 * diffs / ncol(comp_mat)

		### capture the values
		errors[pair_row, ] <- c(locus_error, allele_error, snp_error)
	}
	cat("\n")
	errors
}


## a function to calculate within-population distances
## -- based on the suggestion in Mastretta-Yanes et al. (2015) --
## returns a vector of average distances
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
	dists <- vector()
	index <- 1
	for (pop in pop_names) {
		samples <- sampleids[pops == pop]
		if (length(samples) < 2) {		# can't compute
			next
		}
		sub_mat <- dist_mat[rownames(dist_mat) %in% samples, colnames(dist_mat) %in% samples]
		sub_mat <- sub_mat[lower.tri(sub_mat)]		# keep lower triangle
		avg_dist <- mean(sub_mat)
		dists[index] <- avg_dist
		index <- index + 1
	}
	dists
}


## a function to calculate heterozygosity
## returns a dataframe with the observed heterozygosity for all samples
het_est <- function(genl) {
	comp_mat <- as.matrix(genl)
	as.data.frame(rowSums(comp_mat == 1, na.rm = TRUE) / (ncol(comp_mat) - rowSums(is.na(comp_mat))))
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
	out_pref <- "out"
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
		} else if (args[index] == "-o")  {
			out_pref <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-p")  {
			params_present <- TRUE
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
} else {
	cat("No ipyrad file detected, so calculating heterozygosity only from SNPs\n")
}

if (params_present) {
	params <- as.character(read.table(params_file, header = FALSE)[, 1])
} else {
	params <- unlist(lapply(seq_len(length(catch_args)), function(x) paste0("P", x)))
}

if (reps_present) {
	reps_table <- read.table(reps_file, sep = "\t", header = FALSE)
	reps_table[] <- lapply(reps_table, as.character)
} else {
	cat("No replicates designated, so skipping error rate calculations\n")
}

if (samples_present) {
	samples_table <- read.table(samples_file, sep = "\t", header = FALSE)
	sampleids <- as.character(samples_table[, 1])
	pops <- as.character(samples_table[, 2])		# corresponding sequence to samples (not unique)
} else {
	cat("No samples file detected, so skipping population metrics and using all samples\n")
}


##########
# execute the analyses

## initiate the matrices for holding the results
size_mat <- matrix(ncol = 3)
colnames(size_mat) <- c("param", "num_loci", "num_snps")
paris_mat <- matrix(ncol = 2)
colnames(paris_mat) <- c("param", "num_loci")
error_mat <- matrix(ncol = 4)
colnames(error_mat) <- c("param", "locus_err", "allele_err", "snp_err")
dists_pop <- matrix(ncol = 2)
colnames(dists_pop) <- c("param", "dist")
het_obs <- matrix(ncol = 3)
colnames(het_obs) <- c("sample", "param", "het")

## cycle through the VCF files
index <- 1
for (vcf_file in catch_args) {
	cat(paste0("\nProcessing file ", as.character(index), " of ", as.character(length(catch_args)), "...\n"))
	vcfr <- read.vcfR(vcf_file, verbose = FALSE)
	genl <- suppressWarnings(vcfR2genlight(vcfr))
	cat(paste0("Read in and converted the VCF to a genlight with ", nInd(genl), " samples, ",
		length(unique(genl@chromosome)), " loci, and ",
		nLoc(genl), " SNPs\n"))

	## if there is a samples file, subset the genl
	if (samples_present) {
		keep_samples <- sampleids
		keep_pops <- pops
		genl <- genl[match(keep_samples, genl@ind.names)]
		cat(paste0("Subsetted the genlight to ", nInd(genl), " samples\n"))
	}

	## if there is a replicates file, calculate error rates
	## and remove some samples afterward
	if (reps_present) {
		cat("Running error rate assessments using replicates\n")
		# check that the rep pairs are present in the genlight
		rep_mat <- matrix(ncol = 2)
		for (pair_row in seq_len(nrow(reps_table))) {
			sample1 <- reps_table[pair_row, 1]
			sample2 <- reps_table[pair_row, 2]
			if (all(sample1 %in% genl@ind.names, sample2 %in% genl@ind.names)) {
				if (pair_row == 1) {
					rep_mat[1, ] <- c(sample1, sample2)
				} else {
					rep_mat <- rbind(rep_mat, c(sample1, sample2))
				}
			} else {
				c(paste0("No sample(s) found for replicate pair: ", sample1, " ", sample2, "\n"))
			}
		}
		merror <- mastretta_error(genl, rep_mat)
		out_mat <- as.matrix(cbind(rep(params[index], nrow(merror)),
			merror))
		if (index == 1) {
			error_mat <- out_mat
		} else {
			error_mat <- rbind(error_mat, out_mat)
		}

		## remove one of the reps from the genlight
		genl <- genl[-match(rep_mat[, 2], genl@ind.names)]
		cat(paste0("Subsetted the genlight to ", nInd(genl), " samples\n"))
		## if there is a samples file, remove the reps (and corresponding pops) from those
		if (samples_present) {
			drop_indices <- keep_samples %in% rep_mat[, 2]
			keep_samples <- keep_samples[!drop_indices]
			keep_pops <- keep_pops[!drop_indices]
		}
	}

	## summarise number of loci and SNPs
	if (index == 1) {
		size_mat[1, ] <- c(params[index], length(unique(genl@chromosome)), nLoc(genl))
	} else {
		size_mat <- rbind(size_mat, c(params[index], length(unique(genl@chromosome)), nLoc(genl)))
	}

	## calculate Paris
	cat("Assessing how many loci are shared by at least 80% of samples\n")
	paris_num <- paris(genl)
	if (index == 1) {
		paris_mat[1, ] <- c(params[index], paris_num)
	} else {
		paris_mat <- rbind(paris_mat, c(params[index], paris_num))
	}

	## calculate heterozygosity
	cat("Calculating SNP heterozygosity\n")
	het_df <- het_est(genl)
	out_mat <- as.matrix(cbind(rownames(het_df),
		rep(params[index], nrow(het_df)),
		het_df[, 1]))
	if (index == 1) {
		het_obs <- out_mat
	} else {
		het_obs <- rbind(het_obs, out_mat)
	}

	## if there is a samples file, calculate within population distances
	if (samples_present) {
		suppressMessages(library(SNPfiltR))
		cat("Calculating within population genetic (GENPOFAD) distances\n")
		### make the vcfr smaller by removing samples
		### note that vcfr objects have an extra row for format
		### additionally, drop SNPs that are now monomorphic
		indices <- match(keep_samples, genl@ind.names)
		indices <- c(0, indices)
		indices <- indices + 1
		vcfr <- vcfr[is.biallelic(vcfr), indices]
		vcfr <- min_mac(vcfr, min.mac = 1)
		dists <- pop_dist(vcfr, keep_samples, keep_pops)
		out_mat <- as.matrix(cbind(rep(params[index], length(dists)),
			dists))
		if (index == 1) {
			dists_pop <- out_mat
		} else {
			dists_pop <- rbind(dists_pop, out_mat)
		}
	}

	## increment
	index <- index + 1
}


##########
# save the results to text files
write.table(size_mat, file = paste0(out_pref, "_size.tab"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(paris_mat, file = paste0(out_pref, "_paris.tab"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(error_mat, file = paste0(out_pref, "_error.tab"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(het_obs, file = paste0(out_pref, "_snphet.tab"), quote = FALSE, row.names = FALSE, sep = "\t")
if (samples_present) {
	write.table(dists_pop, file = paste0(out_pref, "_distspop.tab"), quote = FALSE, row.names = FALSE, sep = "\t")
}


##########
# plot
cat("\nPlotting...\n")

## start the pdf
pdf(paste0(out_pref, "_optim.pdf"), width = 10, height = 10)

## size
plot(NULL, xlim = c(1, nrow(size_mat)), ylim = c(0, max(as.numeric(size_mat[, 2])) * 1.05),
	main = "Total number of loci recovered",
	xaxt = "n", ylab = "Number of loci", xlab = "")
points(seq_len(nrow(size_mat)), as.numeric(size_mat[, 2]), pch = 16)
axis(1, at = seq_len(nrow(size_mat)), labels = size_mat[, 1])
plot(NULL, xlim = c(1, nrow(size_mat)), ylim = c(0, max(as.numeric(size_mat[, 3])) * 1.05),
	main = "Total number of biallelic SNPs recovered",
	xaxt = "n", ylab = "Number of biallelic SNPs", xlab = "")
points(seq_len(nrow(size_mat)), as.numeric(size_mat[, 3]), pch = 16)
axis(1, at = seq_len(nrow(size_mat)), labels = size_mat[, 1])

## Paris
plot(NULL, xlim = c(1, nrow(paris_mat)), ylim = c(0, max(as.numeric(paris_mat[, 2])) * 1.05),
	main = "Total number of loci recovered in >= 80% of individuals",
	xaxt = "n", ylab = "Number of loci", xlab = "")
points(seq_len(nrow(paris_mat)), as.numeric(paris_mat[, 2]), pch = 16)
axis(1, at = seq_len(nrow(paris_mat)), labels = paris_mat[, 1])

## heterozygosity
if (ipyrad_present) {
	boxplot(as.numeric(het_stats[, 3]) * 100 ~ het_stats[, 2],
		ylim = c(0, max(as.numeric(het_stats[, 3]) * 100) * 1.05),
		main = "Autosomal heterozygosity\n(ambiguous bases in consensus sequences within samples)",
		ylab = "Heterozygous sites (%)", xlab = "Clustering threshold (%)")
	boxplot(as.numeric(het_obs[, 3]) * 100 ~ het_obs[, 2],
		ylim = c(0, max(as.numeric(het_obs[, 3]) * 100) * 1.05),
		main = "Observed SNP heterozygosity",
		ylab = "Heterozygosity (%)", xlab = "")
} else {
	boxplot(as.numeric(het_obs[, 3]) * 100 ~ het_obs[, 2],
		ylim = c(0, max(as.numeric(het_obs[, 3]) * 100) * 1.05),
		main = "Observed SNP heterozygosity",
		ylab = "Heterozygosity (%)", xlab = "")
}

## error
if (reps_present) {
	boxplot(as.numeric(error_mat[, 2]) ~ error_mat[, 1],
		main = paste0("Locus error rates\n",
			"(loci present in one rep but not the other, relative to total in either/both)"),
		ylim = c(0, max(as.numeric(error_mat[, 2])) * 1.05),
		ylab = "Locus error rate (%)", xlab = "")
	boxplot(as.numeric(error_mat[, 3]) ~ error_mat[, 1],
		main = "Allele error rates\n(loci present in both that differ in sequence, relative to total in both)",
		ylim = c(0, max(as.numeric(error_mat[, 3])) * 1.05),
		ylab = "Allele error rate (%)", xlab = "")
	boxplot(as.numeric(error_mat[, 4]) ~ error_mat[, 1],
		main = "SNP error rates\n(SNPs called in both that differ, relative to total called in both)",
		ylim = c(0, max(as.numeric(error_mat[, 4])) * 1.05),
		ylab = "SNP error rate (%)", xlab = "")
}

## within population distances
if (samples_present) {
	boxplot(as.numeric(dists_pop[, 2]) ~ dists_pop[, 1],
		main = "Within Population Genetic Distances",
		ylim = c(0, max(as.numeric(dists_pop[, 2]), na.rm = TRUE) * 1.05),
		ylab = "Average GENPOFAD distance", xlab = "")
}

## stop the pdf
invisible(dev.off)
