
##########
# Authors: Rachel Binks and Ben Anderson, based on scripts from Mastretta-Yanes et al. 2015
# Date: June-July-Aug 2021
# Description: run Mastretta-Yanes et al. error estimation for technical replicates
##########


# Define any parameters that may need to be changed by users

min_present_rate <- 0.04	# the minimum proportion of samples a SNP must be in (for distances, not error)


# a helper function for when the script is called without arguments
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to run Mastretta-Yanes et al. error estimation for technical",
			"replicates and/or intrapopulation distances using a vcf file\n")
		cat("Usage: Rscript MY_replicate_error.R -o output_pre -r reps_include -s sample_file -v vcf_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default \"output\"]\n")
		cat("\t-r\tWhether the data includes technical replicates (yes [default] or no)\n")
		cat("\t-s\tThe sample names file and pops, with the format:\n",
			"\t\tvcf_label\tsample_name\tpopulation\n",
			"\t\tNote: replicate names should have an \"_R\" after the same text for the other rep\n",
			"\t\tNote: if sample names already match vcf, then you only need one column + pop\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t\tNote: only biallelic SNPs will be kept\n")
	} else {
		cat(help_message)
	}
}



###########
# Read in and format the data
###########


# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1

	outpre <- "output"
	reps_present <- TRUE
	samples_present <- FALSE
	vcf_present <- FALSE

	for (index in 1: length(args)) {
		if (args[index] == "-o") {
			outpre <- args[index + 1]
		} else if (all(args[index] == "-r", args[index + 1] == "no")) {
			reps_present <- FALSE
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
	stop(help("Missing argument for sample/pops file!\n"), call. = FALSE)
}

if (! vcf_present) {
	stop(help("Missing argument for vcf file!\n"), call. = FALSE)
}


# load required libraries
suppressMessages(library(vcfR))
suppressMessages(library(adegenet))


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
myvcf <- read.vcfR(vcf_file, verbose = FALSE)
genl <- suppressWarnings(vcfR2genlight(myvcf))		# convert to genlight object
#genl			# print a summary of the object


# assign correct sample names and population labels to the genlight object
# and filter the genlight to only retain samples present in the sample table
# format1 = id  sample_name     pop		(if vcf ids need to be renamed)
# format2 = id  pop				(if vcf ids = sample names already)
len_original <- length(genl@ind.names)
original_names <- genl@ind.names
if (ncol(sample_table) != 2) {
	genl@ind.names <- sample_table$V2[match(genl@ind.names, sample_table$V1)]
	populations <- sample_table$V3[match(genl@ind.names, sample_table$V2)]
} else {
	genl@ind.names <- sample_table$V1[match(genl@ind.names, sample_table$V1)]
	populations <- sample_table$V2[match(genl@ind.names, sample_table$V1)]
}
if (length(genl@ind.names[is.na(genl@ind.names)]) > 0) {		# if there are any NAs (non-matches)
	original_names <- original_names[!is.na(genl@ind.names)]	# to make sure vcf data matches
	populations <- populations[!is.na(genl@ind.names)]
	genl <- genl[!is.na(genl@ind.names), ]
}
len_new <- length(genl@ind.names)
cat("Retained", length(genl@loc.names), "SNPs from the VCF file\n")
cat("Removed", len_original - len_new, "samples not present in the input sample table\n")


######## Locus, allele and SNP error rates
# Locus error rate: if a locus is in one replicate but not the other relative to total present
## NOTE: previously, and with MY, this is relative to the total number of loci in the dataset
## I think it might make more sense to be relative to the total possible comparisons between the two
# allele error rate: if the locus is in both but differs (by at least one SNP) relative to total present, and
# SNP error rate: if a SNP is found in both and differs (there may be multiple SNPs per locus) relative to total
# If we use VCF input, then we get locus error rate from looking at "chromosome" field
# If we only output one SNP per locus, allele error rate == SNP error rate

# determine names of replicate pairs if reps are present, then run error rate calculations
if (reps_present) {
	cat("Running error estimation\n")
	reps <- grep("_R", genl@ind.names, value = TRUE)
	nonreps <- grep("_R", genl@ind.names, value = TRUE, invert = TRUE)
	samps <- nonreps[pmatch(sub("_R.*", "", reps), nonreps)]
	pairs <- cbind(reps, samps)
	pairs <- pairs[rowSums(is.na(pairs)) < 1, ]		# remove rows with missing pairs

	# calculate error rates for each pair
	error_df <- as.data.frame(matrix(nrow = nrow(pairs), ncol = 3))
	colnames(error_df) <- c("Locus_error", "Allele_error", "SNP_error")
	for (pair_row in 1: nrow(pairs)) {
		cat("Processing replicate pair", pair_row, "of", nrow(pairs), "\n")
		pair_genl <- genl[match(c(pairs[pair_row,]), genl@ind.names)]	# subset genlight for the pair

		# locus
		# use the "chromosome" tag to determine if a locus is present
		NApos <- NA.posi(pair_genl)			# NA positions in the genlight
		NA_rep <- NApos[[1]]				# NA positions in the rep
		rep_genl <- pair_genl[1, -NA_rep]		# remove NA
		rep_loci <- unique(rep_genl@chromosome)		# find unique loci names
		NA_samp <- NApos[[2]]				# NA positions in the samp
		samp_genl <- pair_genl[2, -NA_samp]		# remove NA
		samp_loci <- unique(samp_genl@chromosome)	# find unique loci names
		diffs <- length(setdiff(rep_loci, samp_loci)) + length(setdiff(samp_loci, rep_loci)) 	# differences
		in_both <- intersect(rep_loci, samp_loci)
		locus_error <- 100 * diffs / (diffs + length(in_both))	# diffs divided by possible comparisons

		# allele
		# loci in both
		allele_genl <- pair_genl[, pair_genl@chromosome %in% in_both]	# subset genlight for loci in both
		chrom_list <- rep(NA, length = length(in_both))
		comp_mat <- as.matrix(allele_genl)
		comp_mat[is.na(comp_mat)] <- 3			# change NA to a number that won't match
		index <- 1
		allele_errors <- 0
		for (snp in 1: ncol(comp_mat)) {
			if (as.character(allele_genl@chromosome[snp]) %in% chrom_list) {	# if we have seen this locus already
				if (hit == 0) {				# if we haven't found an error yet for it
					if (comp_mat[1, snp] != comp_mat[2, snp]) {
						allele_errors <- allele_errors + 1
						hit <- 1
					}
				}
			} else {								# new locus
				chrom_list[index] <- as.character(allele_genl@chromosome[snp])	# add it to the list
				index <- index + 1
				if (comp_mat[1, snp] != comp_mat[2, snp]) {
					allele_errors <- allele_errors + 1
					hit <- 1
				} else {
					hit <- 0
				}
			}
		}
		allele_error <- 100 * allele_errors / (length(in_both))

		# SNP
		# all SNPs comparable, but no NA as differences
		comp_mat <- as.matrix(allele_genl)
		comp_mat <- comp_mat[, colSums(is.na(comp_mat)) == 0]	# completely remove any NA comparisons
		snp_errors <- 0
		for (snp in 1: ncol(comp_mat)) {
			if (comp_mat[1, snp] != comp_mat[2, snp]) {
				snp_errors <- snp_errors + 1
			}
		}
		snp_error <- 100 * snp_errors / ncol(comp_mat)

		# Now capture these values for that pair
		error_df[pair_row,] <- c(locus_error, allele_error, snp_error)

		# report
		cat("\tProcessed", length(union(rep_loci, samp_loci)), "total loci, of which", 
			diffs, "were present in one but not the other\n")
		cat("\tProcessed", length(in_both), "loci present in both, of which", 
			allele_errors, "had at least one SNP difference\n")
		cat("\tProcessed", ncol(comp_mat), "SNPs present in both, of which",
			snp_errors, "differed\n")
	}
	# Write the resulting error rates to a file for plotting after
	write.table(error_df, file = paste0(outpre, "_error_table.tab"), sep = "\t", 
			quote = FALSE, row.names = FALSE)
}


# This next step is not yet implemented, nor is it clear whether it should be (independent plotting)
##########
# create PCoA
##########

# create a PCoA to display how well populations/replicates are clustering
#cat("Running a PCoA based on the genlight\n")

# Remove monomorphic SNPs (?)
#gl_matrix <- as.matrix(genl)
#loc_list <- array(NA, length(genl@loc.names))
#for (index in 1: length(genl@loc.names)) {
#	row <- gl_matrix[, index]
#	if (all(row == 0, na.rm = TRUE) | all(row == 2, na.rm = TRUE)
#	| all(row == 1, na.rm = TRUE) | all(is.na(row))) {
#		loc_list[index] <- genl@loc.names[index]
#	}
#}
#loc_list <- loc_list[!is.na(loc_list)]
#if (length(loc_list) > 0) {
#	cat("Dropping", length(loc_list), "monomorphic or all NA loci\n")
#	pca_genl <- genl[, is.na(match(genl@loc.names, loc_list))]		# the ones that don't match
#} else {
#	cat("There are no monomorphic loci\n")
#}
#cat("Retained", length(pca_genl@loc.names), "SNPs\n")

# run the PCA
#pca <- glPca(pca_genl, nf = 4, loadings = FALSE, parallel = TRUE, n.cores = 4)

# note the percentage variation from the first 8 eigenvalues
#perc_eig <- 100 * pca$eig / sum(pca$eig)
#sumvar <- sum(perc_eig[1:8])

# plot the PCoA for visual check of population/replicate clustering

###########
###########


# Remove replicates if they are present
if (reps_present) {
	cat("Dropping", length(reps), "replicates\n")
	original_names <- original_names[! genl@ind.names %in% reps]
	populations <- populations[! genl@ind.names %in% reps]
	genl <- genl[! genl@ind.names %in% reps, ]
}


# Calculate how many loci and SNPs are present in >= 80% of individuals
# this is to emulate some of the ideas in Paris et al. 2017
check_mat <- as.matrix(genl)
keep_snps <- colSums(! is.na(check_mat)) / nrow(check_mat) >= 0.80
check_genl <- genl[, keep_snps]
nsnps <- sum(keep_snps)
nloci <- length(unique(check_genl@chromosome))
cat("There were", nloci, "loci and", nsnps,"SNPs present in more than 80% of individuals\n")
# write to a table for plotting
count_df <- as.data.frame(matrix(ncol = 2, nrow = 1))
colnames(count_df) <- c("loci", "snps")
count_df[1, ] <- c(nloci, nsnps)
write.table(count_df, file = paste0(outpre, "_count80.tab"), sep = "\t", quote = FALSE, row.names = FALSE)


# Remove monomorphic SNPs
gl_matrix <- as.matrix(genl)
loc_list <- array(NA, length(genl@loc.names))
for (index in 1: length(genl@loc.names)) {
	row <- gl_matrix[, index]
	if (all(row == 0, na.rm = TRUE) | all(row == 2, na.rm = TRUE)
	| all(row == 1, na.rm = TRUE) | all(is.na(row))) {
		loc_list[index] <- genl@loc.names[index]
	}
}
loc_list <- loc_list[!is.na(loc_list)]
if (length(loc_list) > 0) {
	cat("Dropping", length(loc_list), "monomorphic or all NA loci\n")
	genl <- genl[, is.na(match(genl@loc.names, loc_list))]		# the ones that don't match
} else {
	cat("There are no monomorphic loci\n")
}
cat("Retained", length(genl@loc.names), "SNPs\n")



##########
# Genetic distance and populations
##########


# If possible, assess genetic distances between individuals in the same populations
cat("Assessing intrapopulation Euclidean distances between individuals\n")

# Filter on missing data
check_mat <- as.matrix(genl)
keep_snps <- colSums(! is.na(check_mat)) / nrow(check_mat) >= min_present_rate
genl <- genl[, keep_snps]
cat("Kept", sum(keep_snps),"SNPs after filtering for a minimum sample coverage of", 
	min_present_rate, "\n")


# Calculate distances between individuals
# Use dist and Euclidean distance, which *should* omit calculations with missing values
dist_obj <- dist(as.matrix(genl))

# Normalize the distances
dist_mat <- as.matrix(dist_obj / max(dist_obj))

# Check that labels match (should)
if (! all(original_names == rownames(dist_mat))) {
	cat("Whoops, names aren't matching as expected\n")
	cat(populations)
	populations <- populations[match(rownames(dist_mat), original_names)]
	cat(populations)
}

# Extract distances between individuals of the same population
pops <- unique(populations)
dists_df <- as.data.frame(matrix(nrow = length(pops), ncol = 2))
colnames(dists_df) <- c("number_ind", "mean_intra_dist")
for (popnum in 1: length(pops)) {
	pop <- pops[popnum]
	pos <- populations %in% pop			# find the positions which match that pop
	dmat <- dist_mat[pos, pos]			# extract a smaller matrix of samples from that pop
	dists <- dmat[lower.tri(dmat)]			# keep the lower triangle values
	avg_dist <- mean(dists)				# calculate an average within population distance
	dists_df[popnum,] <- c(nrow(dmat), avg_dist)
}

# Write the resulting file for plotting after
write.table(dists_df, file = paste0(outpre, "_dist_table.tab"), sep = "\t", quote = FALSE, row.names = FALSE)


# Report that the script has finished
if (reps_present) {
	cat("Finished running MY replicate error assessment and intrapopulation Euclidean distances\n\n")
} else {
	cat("Finished running intrapopulation Euclidean distances\n\n")
}

