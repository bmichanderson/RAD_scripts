
##########
# Author: Ben Anderson, with ideas from Rachel Binks
# Date: Dec 2021
# Modified: May 2022, Oct 2022 (combined with fasta)
# Description: calculate various popgen stats for a VCF file or a fasta alignment
# Note:	at least one of the VCF or fasta files must be present
##########


# load required libraries
suppressMessages(library(vcfR))
suppressMessages(library(ape))
suppressMessages(library(adegenet))
suppressMessages(library(hierfstat))
suppressMessages(library(StAMPP))


# Create an ambiguity code dataframe
amb_df <- data.frame(matrix(nrow = 10, ncol = 2))
rownames(amb_df) <-	c("a", "c", "g", "t", "k", "m", "r", "s", "w", "y")
amb_df[, 1] <- 		c("a", "c", "g", "t", "g", "a", "a", "c", "a", "c")
amb_df[, 2] <- 		c("a", "c", "g", "t", "t", "c", "g", "g", "t", "t")


#####################
# Functions
#####################


## a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate various popgen stats from an input VCF file or fasta alignment\n")
		cat("Usage: Rscript popgen_stats.R -s samples_file [-o out_pref] [-v vcf_file] [-f fasta_file]",
			"[-t run Fst calculations]\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-f\tThe fasta alignment to be analysed\n")
		cat("\t-s\tSamples and populations as a tab-delimited file of sample IDs and pops, one per line\n")
		cat("\t-t\tRun pairwise Fst calculations between populations (only for VCF) [default: do not]\n")
	} else {
		cat(help_message)
	}
}


## a function to calculate standard error
### see: https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
### and: https://www.rdocumentation.org/packages/plotrix/versions/3.8-2/topics/std.error
stde <- function(x) {
	sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))
}


## a function to plot a barplot with error bars of 2 * SE
### see: https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
mybarplot <- function(data_val, data_se, names, main) {
	myplot <- barplot(data_val, names = names, las = 2,
		main = main,
		ylim = c(min(c(0, 1.2 * (data_val - data_se * 2))), 1.2 * max(data_val + data_se * 2)))
	arrows(myplot, y0 = data_val + data_se * 2, y1 = data_val - data_se * 2,
		angle = 90, code = 3, length = 0.1)
}


## a function to calculate n choose 2
### this is the number of ways to choose two elements from a set
### see e.g.: https://www.reddit.com/r/learnmath/comments/3bg21t/what_exactly_is_n_choose_2_in_probability/
nchoose2 <- function(n) {
	n * (n - 1) / 2
}


## a function to calculate observed heterozygosity
### this is the proportion of genotypes/sites that are hets (ambiguities)
het_measure <- function(x) {
	y <- x[!is.na(x)]		# only account for non NA
	mylen <- length(y)
	if (mylen == 0) {
		NA
	} else {
		sum(! y %in% c("a", "c", "g", "t")) / mylen
	}
}


## a function to calculate gene diversity for a fasta alignment
### this is the average number of differences between sequences
### or the probability of choosing an allele that is different
### or the expected heterozygosity
### we can calculate this with either:
### 1) Nei & Roychoudhury 1974:
###		pi = (n / n - 1) * (1 - sum[ (ni / n)^2 ])
### 2) Hohenlohe et al. 2010:
###		pi = 1 - sum[ (ni choose 2) / (n choose 2) ]
### where ni = count of allele i in the sample; n = sum of all allele counts
### we'll use Hohenlohe
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


#####################
# Execution
#####################


# Parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	out_pref <- "output"
	vcf_present <- FALSE
	fasta_present <- FALSE
	samples_present <- FALSE
	run_fst <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			out_pref <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present <- TRUE
			vcf_file <- args[index + 1]
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

if (! samples_present) {
	stop(help("Missing argument for samples file!\n"), call. = FALSE)
}
if (all(c(! vcf_present, ! fasta_present))) {
	stop(help("Missing argument for VCF file or fasta file!\n"), call. = FALSE)
}


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
if (vcf_present) {
	vcf <- read.vcfR(vcf_file, verbose = FALSE)
	cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
		length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")
}
if (fasta_present) {
	fasta <- read.dna(fasta_file, format = "fasta")
}


#####################
# VCF
#####################


# process the VCF file, if present, and calculate stats
if (vcf_present) {
	## check that the samples in the VCF are in the table
	for (indiv in colnames(vcf@gt)[2: ncol(vcf@gt)]) {
		if (! indiv %in% sample_table$V1) {
			stop(help("VCF sample missing from samples file! Stopping\n"), call. = FALSE)
		}
	}


	## calculate amount of missing data, and average per pop
	gt <- extract.gt(vcf, element = "GT")
	missing <- apply(gt, MARGIN = 2, function(x) { sum(is.na(x)) })
	missing <- 100 * missing / nrow(vcf)
	misspops <- sample_table$V2[match(names(missing), sample_table$V1)]
	names(missing) <- misspops
	meanmiss <- tapply(missing, names(missing), mean)
	sample_size <- table(names(missing))


	## convert the VCF into a hierfstat dataframe with pop and genotype
	geni <- vcfR2genind(vcf, return.alleles = TRUE)
	populations <- sample_table$V2[match(rownames(geni@tab), sample_table$V1)]
	geni@pop <- as.factor(populations)
	my_hfst <- genind2hierfstat(geni)


	## calculate basic stats and create a summary dataframe
	cat("Calculating and plotting basic popgen stats for the VCF\n")
	bstats <- basic.stats(my_hfst)
	### if there was only one pop, hierfst adds a dummy pop,
	### which breaks my code, so only keep the first column
	if (length(unique(populations)) == 1) {
		bstats$Ho <- data.frame(bstats$Ho[, 1])
		bstats$Hs <- data.frame(bstats$Hs[, 1])
		bstats$Fis <- data.frame(bstats$Fis[, 1])
		colnames(bstats$Ho) <- populations[1]
	}
	summary <- as.data.frame(matrix(ncol = 8, nrow = length(unique(populations))))
	colnames(summary) <- c("Samples", "PercentMissing", "Ho", "Ho_SE", "Hs",
		"Hs_SE", "Fis", "Fis_SE")
	rownames(summary) <- colnames(bstats$Ho)


	## populate the dataframe with values
	summary[, "Samples"] <- sample_size
	summary[, "PercentMissing"] <- meanmiss
	summary[, "Ho"] <- apply(bstats$Ho, 2, function(x) mean(x, na.rm = TRUE))
	summary[, "Ho_SE"] <- apply(bstats$Ho, 2, stde)
	summary[, "Hs"] <- apply(bstats$Hs, 2, function(x) mean(x, na.rm = TRUE))
	summary[, "Hs_SE"] <- apply(bstats$Hs, 2, stde)
	summary[, "Fis"] <- apply(bstats$Fis, 2, function(x) mean(x, na.rm = TRUE))
	summary[, "Fis_SE"] <- apply(bstats$Fis, 2, stde)


	## export the table to file
	write.table(summary, file = paste0(out_pref, "_summary_VCF.txt"),
		quote = FALSE, row.names = TRUE)


	## in the rare event there are NaN values in the summary, change them to zero for plotting
	summary[is.na(summary)] <- 0


	## Create plots
	### start creating a pdf
	pdf(paste0(out_pref, "_summary_VCF.pdf"), width = 10, height = 10)

	### graph the sample sizes by pop
	barplot(summary$Samples, las = 2, main = "Sample size",
		names = rownames(summary))

	### graph the amount of missing data
	barplot(summary$PercentMissing, las = 2, main = "Mean percentage missing",
		names = rownames(summary))

	### graph the values by population
	mybarplot(summary$Ho, summary$Ho_SE, names = rownames(summary),
		main = "Ho, observed heterozygosity")
	mybarplot(summary$Hs, summary$Hs_SE, names = rownames(summary),
		main = "Hs, estimated gene diversity\n(expected heterozygosity)")
	mybarplot(summary$Fis, summary$Fis_SE, names = rownames(summary),
		main = "Inbreeding coefficient Fis\n(1 - Ho / Hs)")

	### stop creating the pdf
	invisible(dev.off())


	## Run Fst calculations if requested
	if (run_fst) {
		### calculate pairwise Fst
		cat("Calculating pairwise Fst with StAMPP\n")

		### The faster way is via a genlight and StAMPP
		### convert the vcf to a genlight, then add pop
		### NOTE: this will remove sites that are not biallelic
		genl <- vcfR2genlight(vcf)
		pop(genl) <- populations
		cat(paste0("The genlight has ", nLoc(genl), " loci\n"))

		### now use the genlight in StAMPP
		### this calculates Fst following Weir and Cockerham 1984
		### if wanting to get confidence intervals, one could use e.g.
		### nboots = 100, percent = 95, nclusters = 4
		fsts <- stamppFst(genl, nboots = 1)

		### export the pairwise Fst values to file
		write.table(fsts, file = paste0(out_pref, "_StAMPP_Fst.txt"),
			quote = FALSE, row.names = TRUE)

		### A slower way is using the hierfstat package
		### this calculates Fst following Weir and Cockerham 1984
		### since this doesn't convert to genlight, it keeps all sites
		#cat("And more slowly with hierfstat...\n")
		#h_fsts <- pairwise.WCfst(my_hfst)

		### make the resulting square matrix a lower triangle matrix
		#h_fsts[upper.tri(h_fsts)] <- NA

		### export the pairwise Fst values to file
		#write.table(h_fsts, file = paste0(out_pref, "_hierfstat_Fst.txt"),
		#	quote = FALSE, row.names = TRUE)
	}
}


#####################
# fasta
#####################


# process the fasta file, if present, and calculate stats
if (fasta_present) {
	## convert to matrix and correct the DNA so that missing data is coded as NA
	mymat <- as.matrix(as.character(fasta))
	mymat[mymat == "?"] <- NA
	mymat[mymat == "n"] <- NA
	mymat[mymat == "-"] <- NA


	## set the populations to the order of the rows
	## this requires that the fasta file have the same sample names as in the table
	for (indiv in rownames(mymat)) {
		if (! indiv %in% sample_table$V1) {
			stop(help("Fasta sample missing from samples file! Stopping\n"), call. = FALSE)
		}
	}
	populations <- sample_table$V2[match(rownames(mymat), sample_table$V1)]


	## set up a summary dataframe to store the resulting calculation outputs
	summary <- as.data.frame(matrix(ncol = 8, nrow = length(unique(populations))))
	colnames(summary) <- c("Samples", "PercentMissing", "Ho", "Ho_SE", "He",
		"He_SE", "Fis", "Fis_SE")
	rownames(summary) <- unique(populations)


	## cycle through populations and calculate stats
	cat("Calculating and plotting basic popgen stats for the fasta\n")
	index <- 1
	for (pop in unique(populations)) {
		cat("Analyzing population", index, "of", length(unique(populations)), "\n")
		submat <- mymat[populations == pop, ]

		### calculate population sample size and amount of missing data
		summary[pop, "Samples"] <- nrow(submat)
		summary[pop, "PercentMissing"] <- 100 * sum(is.na(submat)) / length(submat)

		### calculate observed heterozygosity
		hets <- apply(submat, 2, het_measure)
		summary[pop, "Ho"] <- mean(hets, na.rm = TRUE)
		summary[pop, "Ho_SE"] <- stde(hets)

		### calculate gene/nucleotide diversity (expected heterozygosity)
		gened <- apply(submat, 2, gd_measure)
		summary[pop, "He"] <- mean(gened, na.rm = TRUE)
		summary[pop, "He_SE"] <- stde(gened)

		### calculate the fixation index inbreeding coefficient Fis
		### from Nei 1977:
		###	Fis = (He - Ho) / He = 1 - Ho / He
		### note: this only uses sites with non-zero gene diversity (He), not all sites
		indices <- which(gened > 0)
		fis <- 1 - hets[indices] / gened[indices]
		summary[pop, "Fis"] <- mean(fis, na.rm = TRUE)
		summary[pop, "Fis_SE"] <- stde(fis)

		### increment the pop index
		index <- index + 1
	}


	## export the table to file
	write.table(summary, file = paste0(out_pref, "_summary_fasta.txt"),
		quote = FALSE, row.names = TRUE)


	## in the rare event there are NaN values in the summary, change them to zero for plotting
	summary[is.na(summary)] <- 0


	## Create plots
	### start creating a pdf
	pdf(paste0(out_pref, "_summary_fasta.pdf"), width = 10, height = 10)

	### graph the sample sizes by pop
	barplot(summary$Samples, las = 2, main = "Sample size",
		names = rownames(summary))

	### graph the amount of missing data
	barplot(summary$PercentMissing, las = 2, main = "Percentage Missing",
		names = rownames(summary))

	### graph the values by population
	mybarplot(summary$Ho, summary$Ho_SE, names = rownames(summary),
		main = "Ho, observed heterozygosity")
	mybarplot(summary$He, summary$He_SE, names = rownames(summary),
		main = "He, estimated gene diversity\n(expected heterozygosity)")
	mybarplot(summary$Fis, summary$Fis_SE, names = rownames(summary),
		main = "Inbreeding coefficient Fis\n(1 - Ho / He)")

	### stop creating the pdf
	invisible(dev.off())
}
