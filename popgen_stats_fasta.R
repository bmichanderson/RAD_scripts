
##########
# Author: Ben Anderson
# Date: Jan 2022
# Modified: May 2022
# Description: calculate various popgen stats for a fasta alignment
##########


# load required libraries
suppressMessages(library(ape))


# Create an ambiguity code dataframe
amb_df <- data.frame(matrix(nrow = 10, ncol = 2))
rownames(amb_df) <- c("a", "c", "g", "t", "k", "m", "r", "s", "w", "y")
amb_df[, 1] <- c("a", "c", "g", "t", "g", "a", "a", "c", "a", "c")
amb_df[, 2] <- c("a", "c", "g", "t", "t", "c", "g", "g", "t", "t")


# Define functions

## a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate various popgen stats from an input fasta alignment\n")
		cat("Usage: Rscript popgen_stats_fasta.R -o out_pref -f fasta_file -s samples_file -t [run Fst calculations]\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-f\tThe fasta file to be analysed (names matching samples file exactly)\n")
		cat("\t-s\tSamples and populations as tab-delimited sample IDs and pops, one per line\n")
		cat("\t-t\tRun pairwise Fst calculations between populations (default: do not)\n")
		cat("\t\tNOTE: Fst calculations should probably not be done on full alignments due to linkage\n")
	} else {
		cat(help_message)
	}
}

## a function to calculate standard error
## see: https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
## and: https://www.rdocumentation.org/packages/plotrix/versions/3.8-2/topics/std.error
stde <- function(x) {
	sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))
}

## a function to plot a barplot with error bars
## see: https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
mybarplot <- function(data_val, data_se, names, main) {
	myplot <- barplot(data_val, names = names, las = 2,
		main = main,
		ylim = c(min(c(0, 1.2 * (data_val - data_se * 2))), 1.2 * max(data_val + data_se * 2)))
	arrows(myplot, y0 = data_val + data_se * 2, y1 = data_val - data_se * 2,
		angle = 90, code = 3, length = 0.1)
}

## a function to calculate n choose 2
## this is the number of ways to choose two elements from a set
## see e.g.: https://www.reddit.com/r/learnmath/comments/3bg21t/what_exactly_is_n_choose_2_in_probability/
nchoose2 <- function(n) {
	n * (n - 1) / 2
}

## a function to calculate observed heterozygosity
## this is the proportion of genotypes/sites that are hets (ambiguities)
het_measure <- function(x) {
	y <- x[!is.na(x)]		# only account for non NA
	mylen <- length(y)
	if (mylen == 0) {
		NA
	} else {
		sum(! y %in% c("a", "c", "g", "t")) / mylen
	}
}

## a function to calculate gene diversity
## this is the average number of differences between sequences
## or the probability of choosing an allele that is different
## or the expected heterozygosity
## we can calculate this with either:
## 1) Nei & Roychoudhury 1974:
##		pi = (n / n - 1) * (1 - sum[ (ni / n)^2 ])
## 2) Hohenlohe et al. 2010:
##		pi = 1 - sum[ (ni choose 2) / (n choose 2) ]
## where ni = count of allele i in the sample; n = sum of all allele counts
## we'll use Hohenlohe
gd_measure <- function(x) {
	y <- x[!is.na(x)]		# only account for non NA
	mylen <- length(y)
	if (mylen == 0) {
		NA
	} else {
		# count alleles
		alleles <- vector("character", mylen * 2)
		index <- 1
		for (genotype in y) {
			alleles[index] <- amb_df[genotype, 1]
			alleles[index + 1] <- amb_df[genotype, 2]
			index <- index + 2
		}
		counts <- vector("numeric", length(unique(alleles)))
		index <- 1
		for (allele in unique(alleles)) {
			counts[index] <- sum(alleles == allele)
			index <- index + 1
		}
		# return the gene diversity
		1 - sum(nchoose2(counts)) / nchoose2(sum(counts))
	}
}

## a function to count alleles sampled at each locus for later Fst
countallele <- function(x) {
	y <- x[!is.na(x)]
	mylen <- length(y)
	if (length(y) == 0) {
		NA
	} else {
		# return the number of alleles sampled
		mylen * 2
	}
}


# Parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	out_pref <- "output"
	fasta_present <- FALSE
	samples_present <- FALSE
	run_fst <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			out_pref <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else if (args[index] == "-f") {
			fasta_present <- TRUE
			fasta_file <- args[index + 1]
		} else if (args[index] == "-t") {
			run_fst <- TRUE
			cat("Will run Fst calculations\n")
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}
if (any(c(! fasta_present, ! samples_present))) {
	stop(help("Missing argument for fasta file and/or samples file!\n"), call. = FALSE)
}


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
fasta <- read.dna(fasta_file, format = "fasta")


# convert to matrix and correct the DNA so that missing data is coded as NA
mymat <- as.matrix(as.character(fasta))
mymat[mymat == "?"] <- NA
mymat[mymat == "n"] <- NA
mymat[mymat == "-"] <- NA


# set the populations to the order of the rows
# this requires that the fasta file have the same sample names as in the table
for (indiv in rownames(mymat)) {
	if (! indiv %in% sample_table$V1) {
		stop(help("Sample missing from samples file! Stopping\n"), call. = FALSE)
	}
}
populations <- sample_table$V2[match(rownames(mymat), sample_table$V1)]


# set up a summary dataframe to store the resulting calculation outputs
summary <- as.data.frame(matrix(ncol = 8, nrow = length(unique(populations))))
colnames(summary) <- c("Samples", "PercentMissing", "Ho", "Ho_SE", "He",
	"He_SE", "Fis", "Fis_SE")
rownames(summary) <- unique(populations)


# set up dataframes to store the gene diversity and n choose 2 for each locus for each pop
if (run_fst) {
	gened_df <- as.data.frame(matrix(ncol = ncol(mymat), nrow = length(unique(populations))))
	rownames(gened_df) <- unique(populations)
	nchoose2_df <- as.data.frame(matrix(ncol = ncol(mymat), nrow = length(unique(populations))))
	rownames(nchoose2_df) <- unique(populations)
}

# cycle through populations and calculate stats
index <- 1
for (pop in unique(populations)) {
	cat("Analyzing population", index, "of", length(unique(populations)), "\n")
	submat <- mymat[populations == pop, ]

	# calculate population sample size and amount of missing data
	summary[pop, "Samples"] <- nrow(submat)
	summary[pop, "PercentMissing"] <- 100 * sum(is.na(submat)) / length(submat)

	# calculate observed heterozygosity
	hets <- apply(submat, 2, het_measure)
	summary[pop, "Ho"] <- mean(hets, na.rm = TRUE)
	summary[pop, "Ho_SE"] <- stde(hets)

	# calculate gene/nucleotide diversity (expected heterozygosity)
	gened <- apply(submat, 2, gd_measure)
	summary[pop, "He"] <- mean(gened, na.rm = TRUE)
	summary[pop, "He_SE"] <- stde(gened)

	if (run_fst) {
		# store the gene diversity
		gened_df[pop, ] <- gened
		# calculate n choose 2 for alleles per locus
		nchoose2_df[pop, ] <- nchoose2(apply(submat, 2, countallele))
	}

	# calculate the fixation index inbreeding coefficient Fis
	# from Nei 1977:
	#	Fis = (He - Ho) / He = 1 - Ho / He
	# note: this only uses sites with non-zero gene diversity (He), not all sites
	indices <- which(gened > 0)
	fis <- 1 - hets[indices] / gened[indices]
	summary[pop, "Fis"] <- mean(fis, na.rm = TRUE)
	summary[pop, "Fis_SE"] <- stde(fis)

	# increment the pop index
	index <- index + 1
}


# export the table to file
write.table(summary, file = paste0(out_pref, "_summary.txt"),
	quote = FALSE, row.names = TRUE)


# in the rare event there are NaN values in the summary, change them to zero for plotting
summary[is.na(summary)] <- 0


# Create plots
# start creating a pdf
pdf(paste0(out_pref, "_summary.pdf"), width = 10, height = 10)

# graph the sample sizes by pop
barplot(summary$Samples, las = 2, main = "Sample size",
	names = rownames(summary))

# graph the amount of missing data
barplot(summary$PercentMissing, las = 2, main = "Percentage Missing",
	names = rownames(summary))

# graph the values by population
# also, use a confidence interval of +/- 2 * SE
mybarplot(summary$Ho, summary$Ho_SE, names = rownames(summary),
	main = "Ho, observed heterozygosity")
mybarplot(summary$He, summary$He_SE, names = rownames(summary),
	main = "He, estimated gene diversity\n(expected heterozygosity)")
mybarplot(summary$Fis, summary$Fis_SE, names = rownames(summary),
	main = "Inbreeding coefficient Fis\n(1 - Ho / He)")

# stop creating the pdf
invisible(dev.off())


# Run Fst calculations if requested
# Hohenlohe et al. 2010:
#	Fst = 1 - sumj [ (nj choose 2) * gdj ] / [ gdp * sumj (nj choose 2) ]
# where nj is the number of alleles sampled in pop j,
# gdj is the gene diversity in pop j,
# and gdp is the gene diversity of the pooled pops
if (run_fst) {
	pops <- unique(populations)
	numpops <- length(pops)
	fst_mat <- matrix(ncol = numpops, nrow = numpops)
	rownames(fst_mat) <- pops
	colnames(fst_mat) <- pops
	nchoose2_mat <- as.matrix(nchoose2_df)
	gened_mat <- as.matrix(gened_df)
	for (index in seq_len(numpops - 1)) {
		cat("Calculating Fst for pop", index, "versus ")
		pop <- pops[index]
		for (index2 in (index + 1): numpops) {
			cat(index2, "")
			pop2 <- pops[index2]
			# determine total gene diversity
			submat <- mymat[populations %in% c(pop, pop2), ]
			total_he <- apply(submat, 2, gd_measure)
			# only use sites with gene diversity
			indices <- which(total_he > 0)
			num_pop <- nchoose2_mat[pop, indices] * gened_mat[pop, indices]
			num_pop2 <- nchoose2_mat[pop2, indices] * gened_mat[pop2, indices]
			denom <- nchoose2_mat[pop, indices] + nchoose2_mat[pop2, indices]
			fst <- 1 - (num_pop + num_pop2) / (total_he[indices] * denom)
			# NOTE: this can produce negative Fst values for individual sites
			fst_mat[index, index2] <- mean(fst, na.rm = TRUE)
		}
		cat("\n")
	}

	# make the resulting pairwise matrix a lower triangle matrix
	myfsts <- fst_mat
	myfsts[lower.tri(myfsts)] <- t(myfsts)[lower.tri(myfsts)]
	myfsts[upper.tri(myfsts)] <- NA
	diag(myfsts) <- NA

	# export the pairwise Fst values to file
	write.table(myfsts, file = paste0(out_pref, "_Fst.txt"),
		quote = FALSE, row.names = TRUE)
}
