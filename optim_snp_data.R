##########
# Author: B.M. Anderson
# Date: Nov–Dec 2023
# Modified: Mar 2025 (added fasta input; updated Paris for polymorphic loci and per population; adjusted plotting)
#	Apr 2025 (adjusted population distance calculations and reporting, and Paris reporting; adjusted het site calc)
#	Feb 2026 (added saving intermediate results and memory cleanup to try to avoid/mitigate crashes)
#	Mar 2026 (adjusted intermediate saves and changed execution to optionally run on those text files)
# Description: evaluate multiple SNP datasets to optimise assembly parameters
# The minimum arguments are multiple VCF files (one per parameter combination) OR a flag to use text files (-e)
#	Optionally, provide text files to:
#	indicate the association between the files and the parameter values (-p), one label per line;
#	specify which sampleIDs are replicate pairs (-r), tab-separated and one per line (second column will be dropped after);
#	enable population calculations and limit analysis to specific sampleIDs (-s), tab-separated with population ID;
#	plot a custom stats file from ipyrad (if you used that) for autosomal heterozygosity output (-i)
#	calculate rough autosomal heterozygosity from alignments, one fasta filename per line (-f)
# Note: this script relies on genlight objects, which only use biallelic SNPs (so may underestimate SNP heterozygosity)
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
		cat("\t-e\tFlag to execute using previously-generated text files with the output prefix (no VCF files)\n")
		cat("\t\tNote: this requires that there be a parameters file (-p) set based on the output\n")
		cat("\t-f\tA file with names of fasta files corresponding to the VCF files, one per line\n")
		cat("\t-i\tA custom-formatted ipyrad stats file with a header and columns containing specific values:\n")
		cat("\t\t\"clust\" = clustering threshold; \"sample\" = sample ID; \"heterozygosity\" = autosomal heterozygosity\n")
		cat("\t-o\tOutput prefix [default: \"out\"]\n")
		cat("\t-p\tA text file with parameter values for each input VCF in the order provided,",
			"one per line (e.g. \"85\" or \"M2\")\n")
		cat("\t\tIf not provided, the parameter values will be arbitrary\n")
		cat("\t-r\tA text file with sampleIDs of replicate pairs (tab separated, one per line)\n")
		cat("\t\tThe second column will be dropped after error rate analysis\n")
		cat("\t-s\tA text file with sampleIDs and populations (tab separated, one per line)\n")
		cat("\t\tIf a samples file is provided, only those specified samples will be analysed\n")
	} else {
		cat(help_message)
	}
}


## a function to compute the number of polymorphic loci shared by >= 80% of samples
## -- based on the metric used in Paris et al. (2017) --
## argument is a genlight object
## returns a number
paris <- function(genl) {
	len_unique <- function(x) length(unique(x[! is.na(x)]))
	mat <- as.matrix(genl)		# rows = individuals, cols = loci
	keep_cols <- apply(mat, 2, len_unique) > 1		# keep loci with multiple genotypes
	mat <- mat[, keep_cols]
	if (!is.null(dim(mat))) {
		keep_snps <- colSums(!is.na(mat)) / nrow(mat) >= 0.80
		length(unique(genl[, keep_snps]@chromosome))
	} else {
		0
	}
}


## a function to run error rate assessments
## -- based on metrics used in Mastretta-Yanes et al. (2015) --
## arguments are a genlight object and a table with the replicate pair sampleIDs, one pair per row (2 columns)
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
		pair_genl <- pair_genl[, !na_both]	# remove loci with NA in both
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
		comp_mat[is.na(comp_mat)] <- -9		# so NA can be compared
		diff_pos <- comp_mat[1, ] != comp_mat[2, ]
		allele_errors <- length(unique(pair_genl[, diff_pos]@chromosome))		# so that multiple in an allele aren't counted
		allele_error <- 100 * allele_errors / (length(in_both))

		### snp error: differences between any shared SNP divided by total shared
		comp_mat <- comp_mat[, colSums(comp_mat == -9) == 0]		# remove any missing
		diffs <- sum(comp_mat[1, ] != comp_mat[2, ])
		snp_error <- 100 * diffs / ncol(comp_mat)

		### capture the values
		errors[pair_row, ] <- c(locus_error, allele_error, snp_error)
	}
	cat("\n")
	errors
}


## a function to calculate within-population genetic distances
## -- based on the suggestion in Mastretta-Yanes et al. (2015) --
## arguments are a vcfR object, a vector of sampleIDs and a corresponding vector of pop names
## returns a vector of average distances
pop_dist <- function(vcfr, sampleids, pops) {
	### load libraries
	suppressMessages(require(ape))
	suppressMessages(require(pofadinr))

	### cycle through the populations, calculating distances
	pop_names <- unique(pops)
	dists <- vector()
	index <- 1
	cat("Populations with fewer than three samples not used (designated with brackets)\n")
	cat(paste0("Population (of ", length(pop_names), "): "))
	for (pop in pop_names) {
		samples <- sampleids[pops == pop]
		if (length(samples) < 3) {		# average distance not meaningful
			cat(paste0("(", index, ") "))
			index <- index + 1
			next
		}
		cat(paste0(index, " "))
		sub_vcfr <- vcfr[sample = samples]
		dnabin <- suppressMessages(vcfR2DNAbin(sub_vcfr, consensus = TRUE, extract.haps = FALSE))
		dnabin <- as.character(dnabin)
		dnabin[dnabin == "n"] <- "?"	# replace missing with a character recognised by pofadinr
		dnabin <- as.DNAbin(dnabin)
		distances <- dist.snp(dnabin, model = "GENPOFAD")
		dist_mat <- as.matrix(distances)
		sub_mat <- dist_mat[lower.tri(dist_mat)]		# keep lower triangle
		avg_dist <- mean(sub_mat)
		dists[index] <- avg_dist
		index <- index + 1
	}
	cat("\n")
	### return the vector
	dists
}


## a function to calculate SNP heterozygosity
## argument is a genlight object
## returns a dataframe with the observed heterozygosity for all samples
het_est <- function(genl) {
	comp_mat <- as.matrix(genl)
	as.data.frame(rowSums(comp_mat == 1, na.rm = TRUE) / (ncol(comp_mat) - rowSums(is.na(comp_mat))))
}


## a function to calculate observed heterozygosity across sites
## this is the proportion of sites that are hets (ambiguities)
het_measure <- function(x) {
	y <- x[!is.na(x)]		# only account for non NA
	mylen <- length(y)
	if (mylen == 0) {
		NA
	} else {
		sum(y %in% c("k", "m", "r", "s", "w", "y", "b", "d", "h", "v")) / mylen
	}
}


###########
# process the call
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
}

catch_args <- vector("list")
extra <- 1
catch <- TRUE
skip_VCF <- FALSE
ipyrad_present <- FALSE
ipyrad_file <- ""
out_pref <- "out"
params_present <- FALSE
params_file <- ""
reps_present <- FALSE
reps_file <- ""
samples_present <- FALSE
samples_file <- ""
fasta_present <- FALSE
fasta_file <- ""
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
	} else if (args[index] == "-f") {
		fasta_present <- TRUE
		fasta_file <- args[index + 1]
		catch <- FALSE
	} else if (args[index] == "-e") {
		skip_VCF <- TRUE
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

if (length(catch_args) < 1) {
	if (skip_VCF) {
		cat("Using pre-generated output for execution instead of VCF files\n")
	} else {
		stop(help("Missing VCF file(s)!\n"), call. = FALSE)
	}
}


# read in the input
if (ipyrad_present) {
	het_stats <- read.table(ipyrad_file, sep = "\t", header = TRUE)
} else {
	cat("No ipyrad file detected, so calculating heterozygosity only from SNPs and/or fasta files\n")
}

if (params_present) {
	params <- as.character(read.table(params_file, header = FALSE)[, 1])
} else if (skip_VCF) {
	stop(help("Missing parameters file for execution without VCFs!\n"), call. = FALSE)
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
	pops <- as.character(samples_table[, 2])		# in same order as samples (not unique)
} else {
	cat("No samples file detected, so skipping population metrics and using all samples\n")
}

if (fasta_present) {
	fastas <- as.character(read.table(fasta_file, header = FALSE)[, 1])
} else {
	cat("No fasta files designated, so skipping fasta heterozygosity calculations\n")
}


##########
# execute the analyses

## initiate the matrices for holding the results and designate output files
error_mat <- matrix(ncol = 4)
colnames(error_mat) <- c("param", "locus_err", "allele_err", "snp_err")
error_file <- paste0(out_pref, "_error.txt")

size_mat <- matrix(ncol = 3)
colnames(size_mat) <- c("param", "num_loci", "num_snps")
size_file <- paste0(out_pref, "_size.txt")

paris_mat <- matrix(ncol = 2)
colnames(paris_mat) <- c("param", "num_loci")
paris_file <- paste0(out_pref, "_paris.txt")

paris_pop_mat <- matrix(ncol = 2)
colnames(paris_pop_mat) <- c("param", "num_loci")
paris_pop_file <- paste0(out_pref, "_paris_pop.txt")

het_mat <- matrix(ncol = 3)
colnames(het_mat) <- c("param", "sample", "het")
het_file <- paste0(out_pref, "_het.txt")

dists_mat <- matrix(ncol = 2)
colnames(dists_mat) <- c("param", "dist")
dists_file <- paste0(out_pref, "_dists.txt")

fastahet_mat <- matrix(ncol = 3)
colnames(fastahet_mat) <- c("param", "sample", "het")
fastahet_file <- paste0(out_pref, "_fastahet.txt")

## cycle through the VCF files
index <- 1
for (vcf_file in catch_args) {
	cat(paste0("\nProcessing VCF file ", as.character(index), " of ", as.character(length(catch_args)), "...\n"))
	vcfr <- read.vcfR(vcf_file, verbose = FALSE)

	## if there is a samples file, subset the vcfR
	if (samples_present) {
		vcfr_samples <- colnames(vcfr@gt)[2: length(colnames(vcfr@gt))]		# the first column is FORMAT
		keep_samples <- sampleids[sampleids %in% vcfr_samples]
		keep_pops <- pops[sampleids %in% vcfr_samples]
		vcfr <- vcfr[sample = keep_samples]
		cat(paste0("Subsetted the VCF to ", ncol(vcfr@gt) - 1, " samples\n"))
		not_found <- sampleids[! sampleids %in% vcfr_samples]
		if (length(not_found) > 0) {
			cat(paste0("Could not find ", length(not_found), " samples from the samples table in the VCF\n"))
			cat("Samples: ")
			cat(not_found, sep = " ")
			cat("\n")
		}
	}

	## convert to genlight
	genl <- suppressWarnings(vcfR2genlight(vcfr))
	cat(paste0("Converted the VCF to a genlight with ", nInd(genl), " samples, ",
		length(unique(genl@chromosome)), " loci, and ", nLoc(genl), " SNPs\n"))

	## if there is a replicates file, calculate error rates
	## and remove some samples afterward
	if (reps_present) {
		cat("Running error rate assessments using replicates\n")
		## check that the rep pairs are present in the genlight
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
		## calculate error and record it
		merror <- mastretta_error(genl, rep_mat)
		out_mat <- as.matrix(cbind(rep(params[index], nrow(merror)), merror))
		write.table(out_mat, file = error_file, append = TRUE, sep = "\t",
			row.names = FALSE, col.names = FALSE, quote = FALSE)		# output to file
		error_mat <- rbind(error_mat, out_mat)
		error_mat <- error_mat[!rowSums(is.na(error_mat)) == ncol(error_mat), ]		# remove rows with all NA (the first one)

		## remove one of the reps from the genlight
		reps_to_drop <- unique(rep_mat[, 2])
		genl <- genl[-match(reps_to_drop, genl@ind.names)]
		cat(paste0("Subsetted the genlight to ", nInd(genl), " samples\n"))
		## if there is a samples file, remove the reps (and corresponding pops) from those and from the vcfr
		if (samples_present) {
			drop_indices <- keep_samples %in% reps_to_drop
			keep_samples <- keep_samples[!drop_indices]
			keep_pops <- keep_pops[!drop_indices]
			vcfr <- vcfr[sample = keep_samples]
		}
	}

	## summarise number of loci and SNPs
	size_row <- c(params[index], length(unique(genl@chromosome)), nLoc(genl))
	cat(size_row, file = size_file, append = TRUE, sep = "\t")		# output to file
	cat("\n", file = size_file, append = TRUE)
	size_mat <- rbind(size_mat, size_row)
	size_mat <- size_mat[!rowSums(is.na(size_mat)) == ncol(size_mat), ]		# remove rows with all NA (the first one)

	## calculate Paris
	cat("Assessing how many polymorphic loci are shared by at least 80% of samples\n")
	paris_num <- paris(genl)
	paris_row <- c(params[index], paris_num)
	cat(paris_row, file = paris_file, append = TRUE, sep = "\t")		# output to file
	cat("\n", file = paris_file, append = TRUE)
	paris_mat <- rbind(paris_mat, paris_row)
	paris_mat <- paris_mat[!rowSums(is.na(paris_mat)) == ncol(paris_mat), ]		# remove rows with all NA (the first one)

	if (samples_present) {
		cat("Also assessing polymorphic loci shared by at least 80% of samples per population\n")
		paris_pop_counts <- vector()
		pop_names <- sort(unique(keep_pops))
		unused <- vector()
		pop_index <- 1
		cat(paste0("Population (of ", length(pop_names), "): "))
		for (pop in pop_names) {
			samples <- keep_samples[keep_pops == pop]
			if (length(samples) > 4) {
				cat(paste0(pop_index, " "))
				subgenl <- genl[match(samples, genl@ind.names)]
				paris_num_pop <- paris(subgenl)
				paris_pop_counts[pop_index] <- paris_num_pop
			} else {
				cat(paste0("(", pop_index, ") "))
				unused <- append(unused, pop_names[pop_index])
			}
			pop_index <- pop_index + 1
		}
		cat("\n")
		if (length(unused) > 0) {
			cat(paste0("The following populations had fewer than 5 individuals and were not used: ",
				paste0(unused, collapse = " "), "\n"))
		}
		out_mat <- as.matrix(cbind(rep(params[index], length(paris_pop_counts)), paris_pop_counts))
		write.table(out_mat, file = paris_pop_file, append = TRUE, sep = "\t",
			row.names = FALSE, col.names = FALSE, quote = FALSE)		# output to file
		paris_pop_mat <- rbind(paris_pop_mat, out_mat)
		paris_pop_mat <- paris_pop_mat[!rowSums(is.na(paris_pop_mat)) == ncol(paris_pop_mat), ]		# remove rows with all NA (the first one)
	}

	## calculate heterozygosity
	cat("Calculating SNP heterozygosity\n")
	het_df <- het_est(genl)
	out_mat <- as.matrix(cbind(rep(params[index], nrow(het_df)), rownames(het_df), het_df[, 1]))
	write.table(out_mat, file = het_file, append = TRUE, sep = "\t",
		row.names = FALSE, col.names = FALSE, quote = FALSE)		# output to file
	het_mat <- rbind(het_mat, out_mat)
	het_mat <- het_mat[!rowSums(is.na(het_mat)) == ncol(het_mat), ]		# remove rows with all NA (the first one)

	## if there is a samples file, calculate within population distances
	if (samples_present) {
		cat("Calculating within population genetic (GENPOFAD) distances\n")
		dists <- pop_dist(vcfr, keep_samples, keep_pops)
		out_mat <- as.matrix(cbind(rep(params[index], length(dists)), dists))
		write.table(out_mat, file = dists_file, append = TRUE, sep = "\t",
			row.names = FALSE, col.names = FALSE, quote = FALSE)		# output to file
		dists_mat <- rbind(dists_mat, out_mat)
		dists_mat <- dists_mat[!rowSums(is.na(dists_mat)) == ncol(dists_mat), ]		# remove rows with all NA (the first one)
	}

	## remove large variables and cleanup memory
	rm(vcfr, genl)
	gc()

	## increment
	index <- index + 1
}

## if fasta files are present, cycle through them
if (fasta_present) {
	suppressMessages(require(ape))
	index <- 1
	for (fasta_file in fastas) {
		cat(paste0("\nProcessing fasta file ", as.character(index), " of ", as.character(length(fastas)), "...\n"))
		fasta <- read.dna(fasta_file, format = "fasta")

		## if samples and/or replicates were included, drop any that are remaining in the fasta file
		keep_samples <- rownames(fasta)
		input <- length(keep_samples)
		if (samples_present) {
			keep_samples <- sampleids[sampleids %in% keep_samples]
		}
		if (reps_present) {
			reps_to_drop <- unique(reps_table[, 2])
			drop_indices <- keep_samples %in% reps_to_drop
			keep_samples <- keep_samples[!drop_indices]
		}
		output <- length(keep_samples)
		if (output < input) {
			fasta <- fasta[keep_samples, ]
			cat(paste0("Subsetted the fasta DNAbin to ", length(rownames(fasta)), " samples\n"))
		}

		## convert to matrix and correct the DNA so that missing data is coded as NA
		mymat <- as.matrix(as.character(fasta))
		mymat[mymat == "?"] <- NA
		mymat[mymat == "n"] <- NA
		mymat[mymat == "-"] <- NA

		## calculate observed autosomal heterozygosity
		hets <- apply(mymat, 1, het_measure)
		het_df <- as.data.frame(hets)
		out_mat <- as.matrix(cbind(rep(params[index], nrow(het_df)), rownames(het_df), het_df[, 1]))
		write.table(out_mat, file = fastahet_file, append = TRUE, sep = "\t",
			row.names = FALSE, col.names = FALSE, quote = FALSE)		# output to file
		fastahet_mat <- rbind(fastahet_mat, out_mat)
		fastahet_mat <- fastahet_mat[!rowSums(is.na(fastahet_mat)) == ncol(fastahet_mat), ]		# remove rows with all NA (the first one)

		## remove large variables and cleanup memory
		rm(fasta, mymat)
		gc()

		## increment
		index <- index + 1
	}
}


##########
# plot
cat("\nPlotting ")

# determine if using variables in this run, or reading in from previous
if (skip_VCF) {
	cat("using data from previous runs in output text files\n")

	if (file.exists(error_file)) {
		error_mat <- as.matrix(read.table(error_file, sep = "\t", header = FALSE))
		colnames(error_mat) <- c("param", "locus_err", "allele_err", "snp_err")
	}

	size_mat <- as.matrix(read.table(size_file, sep = "\t", header = FALSE))
	colnames(size_mat) <- c("param", "num_loci", "num_snps")

	paris_mat <- as.matrix(read.table(paris_file, sep = "\t", header = FALSE))
	colnames(paris_mat) <- c("param", "num_loci")

	if (file.exists(paris_pop_file)) {
		paris_pop_mat <- as.matrix(read.table(paris_pop_file, sep = "\t", header = FALSE))
		colnames(paris_pop_mat) <- c("param", "num_loci")
	}

	het_mat <- as.matrix(read.table(het_file, sep = "\t", header = FALSE))
	colnames(het_mat) <- c("param", "sample", "het")

	if (file.exists(dists_file)) {
		dists_mat <- as.matrix(read.table(dists_file, sep = "\t", header = FALSE))
		colnames(dists_mat) <- c("param", "dist")
	}

	if (file.exists(fastahet_file)) {
		fastahet_mat <- as.matrix(read.table(fastahet_file, sep = "\t", header = FALSE))
		colnames(fastahet_mat) <- c("param", "sample", "het")
	}
} else {
	cat("using only data from the VCFs/fasta files of this run\n")
}

# adjust margins to fit custom labels on x-axis
par(mar = c(7, 4, 4, 1) + 0.1)

# start the pdf
pdf(paste0(out_pref, "_plots.pdf"), width = 10, height = 10)

## size
size_mat <- size_mat[match(params, size_mat[, 1]), ]
plot(NULL, xlim = c(1, nrow(size_mat)), ylim = c(0, max(as.numeric(size_mat[, 2])) * 1.05),
	main = "Total number of loci recovered",
	xaxt = "n", ylab = "Number of loci", xlab = "")
points(seq_len(nrow(size_mat)), as.numeric(size_mat[, 2]), pch = 16)
axis(1, at = seq_len(nrow(size_mat)), labels = FALSE)
labels <- size_mat[, 1]
yrange <- max(as.numeric(size_mat[, 2])) * 1.05
text(seq_len(nrow(size_mat)), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
	labels = labels, xpd = TRUE)

plot(NULL, xlim = c(1, nrow(size_mat)), ylim = c(0, max(as.numeric(size_mat[, 3])) * 1.05),
	main = "Total number of biallelic SNPs recovered",
	xaxt = "n", ylab = "Number of biallelic SNPs", xlab = "")
points(seq_len(nrow(size_mat)), as.numeric(size_mat[, 3]), pch = 16)
axis(1, at = seq_len(nrow(size_mat)), labels = FALSE)
yrange <- max(as.numeric(size_mat[, 3])) * 1.05
text(seq_len(nrow(size_mat)), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
	labels = labels, xpd = TRUE)

## Paris
paris_mat <- paris_mat[match(params, paris_mat[, 1]), ]
plot(NULL, xlim = c(1, nrow(paris_mat)), ylim = c(0, max(as.numeric(paris_mat[, 2])) * 1.05),
	main = "Number of polymorphic loci recovered in >= 80% of all samples",
	xaxt = "n", ylab = "Number of loci", xlab = "")
points(seq_len(nrow(paris_mat)), as.numeric(paris_mat[, 2]), pch = 16)
labels <- paris_mat[, 1]
yrange <- max(as.numeric(paris_mat[, 2])) * 1.05
axis(1, at = seq_len(nrow(paris_mat)), labels = FALSE)
text(seq_len(nrow(paris_mat)), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
	labels = labels, xpd = TRUE)

if (nrow(paris_pop_mat) > 1) {
	paris_pop_df <- as.data.frame(paris_pop_mat)
	mylevels <- paris_pop_df[, 1][match(params, paris_pop_df[, 1])]
	paris_pop_df[, 1] <- factor(paris_pop_df[, 1], levels = mylevels)
	paris_pop_df[, 2] <- as.numeric(paris_pop_df[, 2])

	boxplot(paris_pop_df[, 2] ~ paris_pop_df[, 1],
		main = "Number of polymorphic loci recovered in >= 80% of samples within each population",
		ylim = c(0, max(paris_pop_df[, 2], na.rm = TRUE) * 1.05),
		ylab = "Number of loci", xlab = "", xaxt = "n")
	yrange <- max(paris_pop_df[, 2], na.rm = TRUE) * 1.05
	axis(1, at = seq_len(length(mylevels)), labels = FALSE)
	text(seq_len(length(mylevels)), par("usr")[3] - (0.03 * yrange),
		srt = 45, adj = 1, labels = mylevels, xpd = TRUE)
}

## heterozygosity
het_df <- as.data.frame(het_mat)
mylevels <- het_df[, 1][match(params, het_df[, 1])]
het_df[, 1] <- factor(het_df[, 1], levels = mylevels)
het_df[, 3] <- as.numeric(het_df[, 3])

boxplot(het_df[, 3] * 100 ~ het_df[, 1],
	ylim = c(0, max(het_df[, 3] * 100) * 1.05),
	main = "Observed SNP heterozygosity",
	ylab = "Heterozygosity (%)", xlab = "", xaxt = "n")
yrange <- max(het_df[, 3] * 100) * 1.05
axis(1, at = seq_len(length(mylevels)), labels = FALSE)
text(seq_len(length(mylevels)), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
	labels = mylevels, xpd = TRUE)

if (ipyrad_present) {
	mylevels <- unique(het_stats[, "clust"])
	het_stats[, "clust"] <- factor(het_stats[, "clust"], levels = mylevels)
	het_stats[, "heterozygosity"] <- as.numeric(het_stats[, "heterozygosity"])

	boxplot(het_stats[, "heterozygosity"] * 100 ~ het_stats[, "clust"],
		ylim = c(0, max(het_stats[, "heterozygosity"] * 100) * 1.05),
		main = "Autosomal heterozygosity (ipyrad stats)\n(ambiguous bases in consensus sequences within samples)",
		ylab = "Heterozygous sites (%)", xlab = "", xaxt = "n")
	yrange <- max(het_stats[, "heterozygosity"] * 100) * 1.05
	axis(1, at = seq_len(length(mylevels)), labels = FALSE)
	text(seq_len(length(mylevels)), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
		labels = mylevels, xpd = TRUE)
}

if (nrow(fastahet_mat) > 1) {
	fastahet_df <- as.data.frame(fastahet_mat)
	mylevels <- fastahet_df[, 1][match(params, fastahet_df[, 1])]
	fastahet_df[, 1] <- factor(fastahet_df[, 1], levels = mylevels)
	fastahet_df[, 3] <- as.numeric(fastahet_df[, 3])

	boxplot(fastahet_df[, 3] * 100 ~ fastahet_df[, 1],
		ylim = c(0, max(fastahet_df[, 3] * 100) * 1.05),
		main = "Observed autosomal heterozygosity (fasta file)",
		ylab = "Heterozygosity (%)", xlab = "", xaxt = "n")
	yrange <- max(fastahet_df[, 3] * 100) * 1.05
	axis(1, at = seq_len(length(mylevels)), labels = FALSE)
	text(seq_len(length(mylevels)), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
		labels = mylevels, xpd = TRUE)
}

## error
if (nrow(error_mat) > 1) {
	error_df <- as.data.frame(error_mat)
	mylevels <- error_df[, 1][match(params, error_df[, 1])]
	error_df[, 1] <- factor(error_df[, 1], levels = mylevels)
	error_df[, 2] <- as.numeric(error_df[, 2])
	error_df[, 3] <- as.numeric(error_df[, 3])
	error_df[, 4] <- as.numeric(error_df[, 4])

	boxplot(error_df[, 2] ~ error_df[, 1],
		main = paste0("Locus error rates\n",
			"(loci present in one rep but not the other, relative to total in either/both)"),
		ylim = c(0, max(error_df[, 2]) * 1.05),
		ylab = "Locus error rate (%)", xlab = "", xaxt = "n")
	yrange <- max(error_df[, 2]) * 1.05
	axis(1, at = seq_len(length(mylevels)), labels = FALSE)
	text(seq_len(length(mylevels)), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
		labels = mylevels, xpd = TRUE)

	boxplot(error_df[, 3] ~ error_df[, 1],
		main = "Allele error rates\n(loci present in both that differ in sequence, relative to total in both)",
		ylim = c(0, max(error_df[, 3]) * 1.05),
		ylab = "Allele error rate (%)", xlab = "", xaxt = "n")
	axis(1, at = seq_len(length(mylevels)), labels = FALSE)
	yrange <- max(error_df[, 3]) * 1.05
	text(seq_len(length(unique((error_df[, 1])))), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
		labels = mylevels, xpd = TRUE)

	boxplot(error_df[, 4] ~ error_df[, 1],
		main = "SNP error rates\n(SNPs called in both that differ, relative to total called in both)",
		ylim = c(0, max(error_df[, 4]) * 1.05),
		ylab = "SNP error rate (%)", xlab = "", xaxt = "n")
	axis(1, at = seq_len(length(mylevels)), labels = FALSE)
	yrange <- max(error_df[, 4]) * 1.05
	text(seq_len(length(unique((error_df[, 1])))), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
		labels = mylevels, xpd = TRUE)
}

## within population distances
if (nrow(dists_mat) > 1) {
	dists_df <- as.data.frame(dists_mat)
	mylevels <- dists_df[, 1][match(params, dists_df[, 1])]
	dists_df[, 1] <- factor(dists_df[, 1], levels = mylevels)
	dists_df[, 2] <- as.numeric(dists_df[, 2])

	boxplot(dists_df[, 2] ~ dists_df[, 1],
		main = "Within Population Genetic Distances",
		ylim = c(0, max(dists_df[, 2], na.rm = TRUE) * 1.05),
		ylab = "Average GENPOFAD distance", xlab = "", xaxt = "n")
	yrange <- max(dists_df[, 2], na.rm = TRUE) * 1.05
	axis(1, at = seq_len(length(mylevels)), labels = FALSE)
	text(seq_len(length(unique((dists_df[, 1])))), par("usr")[3] - (0.03 * yrange), srt = 45, adj = 1,
		labels = mylevels, xpd = TRUE)
}

# stop the pdf
invisible(dev.off())
