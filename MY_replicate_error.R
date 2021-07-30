
##########
# Authors: Rachel Binks and Ben Anderson, based on scripts from Mastretta-Yanes et al. 2015
# Date: June-July 2021
# Description: run Mastretta-Yanes et al. error estimation for technical replicates
##########


# Define any parameters that may need to be changed by users

min_present_rate <- 0.7		# the minimum proportion of samples a SNP must be in (for distances, not error)



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
		} else if (args[index] == "-r") {
			if (args[index + 1] == "no") {
				reps_present <- FALSE
			}
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
genl <- vcfR2genlight(myvcf)		# convert to genlight object
genl			# print a summary of the object

# assign correct sample names and population labels to the genlight object
# and filter the genlight to only retain samples present in the sample table
# format1 = id  sample_name     pop		(if vcf ids need to be renamed)
# format2 = id  pop				(if vcf ids = sample names already)
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


# NOTE: I don't think we should do any filtering for missing data, as part of the error
# test is to see whether loci are called in one replicate but not the other


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
		pair_genl <- genl[match(c(pairs[pair_row,]), genl@ind.names)]		# subset the genlight for the pair

		# locus
		# use the "chromosome" tag to determine if a locus is present
		NApos <- NA.posi(pair_genl)			# NA positions in the genlight
		NA_rep <- NApos[[1]]				# NA positions in the rep
		rep_genl <- pair_genl[1, -NA_rep]		# remove NA
		rep_loci <- unique(rep_genl@chromosome)		# find unique loci names
		NA_samp <- NApos[[2]]				# NA positions in the samp
		samp_genl <- pair_genl[2, -NA_samp]		# remove NA
		samp_loci <- unique(samp_genl@chromosome)	# find unique loci names
		errors <- length(setdiff(rep_loci, samp_loci)) + length(setdiff(samp_loci, rep_loci)) 	# differences
		in_both <- intersect(rep_loci, samp_loci)
		loc_error <- 100 * errors / (errors + length(in_both))	# errors divided by total possible comparisons

		# allele
		# loci in both
		alle_genl <- pair_genl[, pair_genl@chromosome %in% in_both]	# subset the genlight for loci in both
		chrom_list <- rep(NA, length = length(in_both))
		comp_mat <- as.matrix(alle_genl)
		comp_mat[is.na(comp_mat)] <- 3			# change NA to a number that won't match
		index <- 1
		errors <- 0
		for (snp in 1: ncol(comp_mat)) {
			if (as.character(alle_genl@chromosome[snp]) %in% chrom_list) {	# if we have seen this locus already
				if (hit == 0) {				# if we haven't found an error yet for it
					if (comp_mat[1, snp] != comp_mat[2, snp]) {
						errors <- errors + 1
						hit <- 1
					}
				}
			} else {								# new locus
				chrom_list[index] <- as.character(alle_genl@chromosome[snp])	# add it to the list
				index <- index + 1
				if (comp_mat[1, snp] != comp_mat[2, snp]) {
					errors <- errors + 1
					hit <- 1
				} else {
					hit <- 0
				}
			}
		}
		alle_error <- 100 * errors / (length(in_both))

		# SNP
		# all SNPs comparable, but no NA as differences
		comp_mat <- as.matrix(alle_genl)
		comp_mat <- comp_mat[, colSums(is.na(comp_mat)) == 0]		# completely remove any NA comparisons
		errors <- 0
		for (snp in 1: ncol(comp_mat)) {
			if (comp_mat[1, snp] != comp_mat[2, snp]) {
				errors <- errors + 1
			}
		}
		snp_error <- 100 * errors / ncol(comp_mat)

		# Now capture these values for that pair
		error_df[pair_row,] <- c(loc_error, alle_error, snp_error)
	}

	# Write the resulting error rates to a file for plotting after
	write.table(error_df, file = paste0(outpre, "_error_table.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
}


##########
# Genetic distance and populations
##########

# If possible, assess genetic distances between individuals in the same populations

# First, remove replicates if they are present
if (reps_present) {
	cat("Dropping", length(reps), "replicates\n")
	original_names <- original_names[! genl@ind.names %in% reps]
	populations <- populations[! genl@ind.names %in% reps]
	genl <- genl[! genl@ind.names %in% reps, ]
}


# Filter on missing data
cat("Filtering for a minimum sample coverage of", min_present_rate, "to keep a SNP\n")
check_mat <- as.matrix(genl)
keep_snps <- colSums(! is.na(check_mat)) / nrow(check_mat) >= min_present_rate
genl <- genl[, keep_snps]


# Calculate distances between individuals
# Use dist and Euclidean distance, which should omit calculations with missing values
cat("Calculating Euclidean distance between individuals\n")
dist_obj <- dist(as.matrix(genl))

# Normalize the distances
dist_mat <- as.matrix(dist_obj / max(dist_obj))

# Check that labels match (should)
if (! all(original_names == rownames(dist_mat))) {
	new_names <- original_names[match(rownames(dist_mat), original_names)]
	populations <- populations[match(rownames(dist_mat), original_names)]
}

# Extract distances between individuals of the same population
pops <- unique(populations)
dists_df <- as.data.frame(matrix(nrow = length(pops), ncol = 2))
colnames(dists_df) <- c("number_ind", "mean_intra_dist")
for (popnum in 1: length(pops)) {
	pop <- pops[popnum]
	pos <- populations == pops[popnum]		# find the positions which match that pop
	dmat <- dist_mat[pos, pos]			# extract a smaller matrix of samples from that pop
	dists <- dmat[lower.tri(dmat)]			# keep the lower triangle values
	avg_dist <- mean(dists)				# calculate an average within population distance
	dists_df[popnum,] <- c(nrow(dmat), avg_dist)
}

# Write the resulting file for plotting after
write.table(dists_df, file = paste0(outpre, "_dist_table.tab"), sep = "\t", quote = FALSE, row.names = FALSE)



##########
# create PCoA
##########

# create a PCoA to display how well populations are clustering and proportion of variation explained
cat("Running a PCoA based on the genlight\n")
pca_df <- as.data.frame(matrix(nrow = 1, ncol = 4))
colnames(pca_df) <- c("PC1", "PC2", "PC3", "PC4")
pca <- glPca(genl, nf = 4, loadings = FALSE, parallel = TRUE, n.cores = 4)

# percentage variation from eigenvalues
perc_eig <- 100 * pca$eig / sum(pca$eig)

# capture the first four PCs % contribution
pca_df[1, ] <- perc_eig[1:4]

# Write the resulting file and the PCoA scores with populations for plotting after
write.table(pca_df, file = paste0(outpre, "_pca_table.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pca$scores, file = paste0(outpre, "_pca_scores.tab"), sep = "\t", quote = FALSE, 
		row.names = TRUE, col.names = TRUE)
write.table(populations, file = paste0(outpre, "_pca_pops.tab"), sep = "\t", quote = FALSE, 
		row.names = FALSE, col.names = FALSE)



# Report that the script has finished
if (reps_present) {
	cat("Finished running MY replicate error assessment and intrapopulation Euclidean distances and PCoA\n")
} else {
	cat("Finished running MY intrapopulation Euclidean distances and PCoA\n")
}

