
##########
# Author: B.M. Anderson, with ideas from Rachel Binks
# Date: Dec 2021
# Modified: May 2022, Oct 2022 (combined with fasta), Oct 2023 (added Reich Fst),
#	Mar 2025 (adjust Reich calcs and made a separate function; added argument to select Fst type; corrected boxplot)
#	Apr 2025 (changed fasta behaviour to take a list of multiple files and average across; assume one population)
#	Apr 2025 (added support for individual sample measures and adjusted output and calcs for more potential hets)
# Description: calculate various popgen stats for a VCF file or a set of fasta alignments specified in a text file
# Note:	either the VCF or fasta files must be present
# Note: if specifying fasta input, it is assumed the alignments come from a single population
#	(one set of stats calculated); the samples file can be used to designate which samples are of interest;
#	other samples (if present) will be dropped from the fasta alignments
##########


# load required libraries
suppressMessages(library(vcfR))
suppressMessages(library(ape))
suppressMessages(library(adegenet))
suppressMessages(library(hierfstat))
suppressMessages(library(StAMPP))


#####################
# Functions
#####################


## a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate various popgen stats from an input VCF file or list of fasta alignments\n")
		cat("Usage: Rscript popgen_stats.R -s samples_file [-o out_pref] [-v vcf_file] [-f fasta_list]",
			"[-t flag to run Fst calculations] [-r Fst type to run]\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-f\tThe text file list of fasta alignments to be analysed (from a single population)\n")
		cat("\t-s\tSamples and populations as a tab-delimited file of sample IDs and pops, one per line\n")
		cat("\t-t\tRun pairwise Fst calculations between populations (only for VCF) [default: do not]\n")
		cat("\t-r\tType of Fst calculation: \"s\" for StAMPP Weir and Cockerham 1984 or",
			"\"r\" for Reich et al. 2009 [default: both]\n")
	} else {
		cat(help_message)
	}
}


## a function to calculate standard error
### see: https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
stde <- function(x) {
	sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))
}


## a function to plot a barplot with error bars of 2 * SE
### see: https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
mybarplot <- function(data_val, data_se, names, main) {
	myplot <- barplot(data_val, names = names, las = 2,
		main = main,
		ylim = c(min(c(0, 1.2 * min(data_val - data_se * 2))),
			max(c(0, 1.2 * max(data_val + data_se * 2)))))
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
		sum(y %in% c("k", "m", "r", "s", "w", "y", "b", "d", "h", "v")) / mylen
	}
}


## a function to calculate gene/nucleotide diversity (pi) for a population
### this represents the average number of differences between sequences
### or the probability of choosing an allele that is different
### or the expected heterozygosity
### Two formulas for calculating this use counts of alleles at a locus in a population
### where ni = count of allele i, and n = sum of all allele counts
### (assuming HWE)
### 1) Nei & Roychoudhury 1974:	pi = n * (1 - sum[(ni / n)^2]) / (n - 1)
### and
### 2) Hohenlohe et al. 2010: pi = 1 - sum[(ni choose 2) / (n choose 2)]
### The metric can be averaged across sites/loci, but separate estimates should be unlinked!
### If there is complete inbreeding, the sample number should be individuals rather than loci

### Create an ambiguity code dataframe
amb_df <- data.frame(matrix(nrow = 14, ncol = 3))
rownames(amb_df) <-	c("a", "c", "g", "t", "k", "m", "r", "s", "w", "y", "b", "d", "h", "v")
amb_df[, 1] <- 		c("a", "c", "g", "t", "g", "a", "a", "c", "a", "c", "c", "a", "a", "a")
amb_df[, 2] <- 		c("a", "c", "g", "t", "t", "c", "g", "g", "t", "t", "g", "g", "c", "c")
amb_df[, 3] <- 		c("", "", "", "", "", "", "", "", "", "", "t", "t", "t", "g")

### Use (2) Hohenlohe et al. 2010 for the function:
gd_measure <- function(x) {
	y <- x[!is.na(x)]		# remove missing
	mylen <- length(y)
	if (mylen == 0) {
		NA
	} else {
		# count alleles
		alleles <- vector("character")
		index <- 1
		for (genotype in y) {
			alleles[index] <- amb_df[genotype, 1]
			alleles[index + 1] <- amb_df[genotype, 2]
			if (genotype %in% c("b", "d", "h", "v")) {
				alleles[index + 2] <- amb_df[genotype, 3]
				index <- index + 3
			} else {
				index <- index + 2
			}
		}
		counts <- vector("numeric", length(unique(alleles)))
		index <- 1
		for (allele in unique(alleles)) {
			counts[index] <- sum(alleles == allele)
			index <- index + 1
		}
		# return pi
		1 - sum(nchoose2(counts)) / nchoose2(sum(counts))
	}
}


## a function to calculate Reich et al. 2009 Fst estimates
### Another way to calculate Fst that accounts for smaller/different
### sample sizes was put forward by Reich et al. 2009 and shown to
### be less biased by sample size in Willing et al. 2012
### This measure was implemented for a genlight object here:
### https://github.com/jessicarick/reich-fst/blob/master/reich_fst.R
### A modified version of that calculation is included here to avoid
### the need to use other R packages
### Note: due to precision in R, error can accumulate when there is
### addition after division, so reduce this where possible and
### round the result to 4 decimal places
reich_fst <- function(genl, populations) {
	pops <- unique(populations)
	fsts_reich <- matrix(nrow = length(pops),
		ncol = length(pops),
		dimnames = list(pops, pops))
	index <- 1
	cat(paste0("Calculating pairwise Fst for populations",
		" (", length(pops), "):"))
	for (pop in pops) {
		cat(paste0(" ", index))
		pop1mat <- as.matrix(genl[genl@pop == pop])
		a1 <- apply(pop1mat, 2, function(x) sum(x, na.rm = TRUE))
		n1 <- apply(pop1mat, 2, function(x) 2 * sum(!is.na(x)))
		h1 <- (a1 * (n1 - a1)) / (n1 * (n1 - 1))
		for (pop2 in pops[-1: -index]) {
			pop2mat <- as.matrix(genl[genl@pop == pop2])
			a2 <- apply(pop2mat, 2, function(x) sum(x, na.rm = TRUE))
			n2 <- apply(pop2mat, 2, function(x) 2 * sum(!is.na(x)))
			h2 <- (a2 * (n2 - a2)) / (n2 * (n2 - 1))
			bign <- (((a1 * n2) - (a2 * n1)) / (n1 * n2))^2 -
				(a1 * (n1 - a1)) / (n1 * n1 * (n1 - 1)) - (a2 * (n2 - a2)) / (n2 * n2 * (n2 - 1))
			bigd <- bign + h1 + h2
			fst_r <- sum(bign, na.rm = TRUE) / sum(bigd, na.rm = TRUE)
			fsts_reich[pop2, pop] <- round(fst_r, digits = 4)
		}
		index <- index + 1
	}
	cat("\n")
	# return the matrix
	fsts_reich
}


#####################
# Execution
#####################


# Parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) { # nolint
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	extra <- 1
	catch <- TRUE
	out_pref <- "output"
	vcf_present <- FALSE
	fasta_present <- FALSE
	samples_present <- FALSE
	run_fst <- FALSE
	type_fst <- "both"
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			out_pref <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-v") {
			vcf_present <- TRUE
			vcf_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-f") {
			fasta_present <- TRUE
			fasta_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-t") {
			run_fst <- TRUE
			cat("Will run Fst calculations\n")
		} else if (args[index] == "-r") {
			type_fst <- as.character(args[index + 1])
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

if (! samples_present) {
	stop(help("Missing argument for samples file!\n"), call. = FALSE)
}
if (all(c(! vcf_present, ! fasta_present))) {
	stop(help("Missing argument for VCF file or fasta file list!\n"), call. = FALSE)
}


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
sample_table[2] <- lapply(sample_table[2], as.character)		# convert if populations are numbers
if (vcf_present) {
	vcf <- read.vcfR(vcf_file, verbose = FALSE)
	cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
		length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")
}
if (fasta_present) {
	# read in the list and store (a single column of strings)
	fasta_list <- read.table(fasta_file, header = FALSE)[, 1]
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
	missing <- apply(gt, MARGIN = 2, function(x) {
			sum(is.na(x))
		}
	)
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
	summary <- as.data.frame(matrix(ncol = 9, nrow = length(unique(populations))))
	colnames(summary) <- c("Samples", "Loci", "PercentMissing", "Ho", "Ho_SE", "Hs",
		"Hs_SE", "Fis", "Fis_SE")
	rownames(summary) <- colnames(bstats$Ho)

	## populate the dataframe with values
	summary[, "Samples"] <- sample_size
	summary[, "SNPs"] <- nrow(bstats$Ho)
	summary[, "Mean_percent_miss"] <- meanmiss
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
	barplot(summary$Samples, las = 2, main = "Samples",
		names = rownames(summary))

	### graph the amount of missing data
	barplot(summary$Mean_percent_miss, las = 2, main = "Average missing (%)",
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
		### convert to genlight
		genl <- vcfR2genlight(vcf)
		pop(genl) <- populations
		cat(paste0("The genlight has ", nLoc(genl), " loci\n"))

		if (type_fst == "r") {
			### calculate Reich Fst
			cat("Calculating pairwise Fst using an adaptation of Reich et al. 2009\n")
			fsts_reich <- reich_fst(genl, populations)
			### export the pairwise Fst values to file
			write.table(fsts_reich, file = paste0(out_pref, "_Reich_Fst.txt"),
				quote = FALSE, row.names = TRUE)
		} else if (type_fst == "s") {
			### calculate Weir and Cockerham Fst with StAMPP
			cat("Calculating pairwise Fst using StAMPP\n")
			fsts <- stamppFst(genl, nboots = 1)
			### export the pairwise Fst values to file
			write.table(fsts, file = paste0(out_pref, "_StAMPP_Fst.txt"),
				quote = FALSE, row.names = TRUE)
		} else {
			### calculate Weir and Cockerham Fst with StAMPP
			cat("Calculating pairwise Fst using StAMPP\n")
			fsts <- stamppFst(genl, nboots = 1)
			### export the pairwise Fst values to file
			write.table(fsts, file = paste0(out_pref, "_StAMPP_Fst.txt"),
				quote = FALSE, row.names = TRUE)
			### calculate Reich Fst
			cat("Calculating pairwise Fst using an adaptation of Reich et al. 2009\n")
			fsts_reich <- reich_fst(genl, populations)
			### export the pairwise Fst values to file
			write.table(fsts_reich, file = paste0(out_pref, "_Reich_Fst.txt"),
				quote = FALSE, row.names = TRUE)
		}
	}
}


#####################
# fasta
#####################


# process the fasta files (individual loci), if present, and calculate stats
if (fasta_present) {
	cat("Calculating basic popgen stats for the fasta files (loci)\n")
	allpops <- vector()

	## make a matrix to hold the results per locus (rows) for the population
	## Columns:
	## 1: number of samples for that locus
	## 2: number of sites for that locus
	## 3: number of polymorphic alignment positions for that locus
	## 4: average percentage missing sites across samples for that locus
	## 5: average observed heterozygosity for genotypes in that locus
	## 6: average nucleotide diversity for sites in that locus
	## 7: average inbreeding coefficient for sites in the locus
	results <- matrix(0, ncol = 7, nrow = length(fasta_list))
	colnames(results) <- c("Samples", "Sites", "Polymorphic", "Missing", "Ho", "pi", "Fis")
	index <- 1		# track each locus

	## make a second matrix to hold the measures for each sample
	## Columns:
	## 1: count of loci containing that sample
	## 2: count of sites across loci containing that sample
	## 3: count of polymorphic/heterozygous sites across loci containing that sample
	## 4: count of missing sites across loci containing that sample
	## 5: sum across loci of average observed heterozygosity across sites per locus
	results_samp <- matrix(0, ncol = 5, nrow = nrow(sample_table))
	colnames(results_samp) <- c("Loci", "Sites", "Polymorphic", "Missing", "Ho")

	## cycle through the list of files
	incongruence <- FALSE
	for (fasta in fasta_list) {
		### if there are many files, report progress
		if (index %% 1000 == 0) {
			cat("Processed", index, "files of", as.character(length(fasta_list)), "\n")
		}

		### read in the alignment and convert to matrix
		fasta_align <- read.dna(fasta, format = "fasta")
		mymat <- as.matrix(as.character(fasta_align))

		### check and filter the alignment to samples in the sample table
		drop_indices <- vector()
		drop_index <- 1
		for (indiv in rownames(mymat)) {
			if (! indiv %in% sample_table$V1) {
				drop_indices <- append(drop_indices, drop_index)
			}
			drop_index <- drop_index + 1
		}
		if (length(drop_indices) > 0) {
			incongruence <- TRUE
			mymat <- mymat[-drop_indices, ]
		}

		### check there is one population left
		populations <- sample_table$V2[match(rownames(mymat), sample_table$V1)]
		if (length(unique(populations)) > 1) {
				cat("WARNING: more than one population present in", fasta, "; stats will be wrong!\n")
		}
		allpops <- unique(c(allpops, populations))

		### correct the DNA so that missing data is coded as NA
		mymat[mymat == "?"] <- NA
		mymat[mymat == "n"] <- NA
		mymat[mymat == "-"] <- NA

		### check all samples are present; if there are missing, make NA rows
		samples_present <- 0
		for (sample_index in seq_len(nrow(sample_table))) {
			sample <- sample_table$V1[sample_index]
			if (! sample %in% rownames(mymat)) {
				curr_names <- rownames(mymat)
				mymat <- rbind(mymat, rep(NA, ncol(mymat)))
				rownames(mymat) <- c(curr_names, sample)
			} else {
				samples_present <- samples_present + 1
				# increment the locus and site counts for present samples
				results_samp[sample_index, "Loci"] <- results_samp[[sample_index, "Loci"]] + 1
				results_samp[sample_index, "Sites"] <- results_samp[[sample_index, "Sites"]] + ncol(mymat)
			}
		}

		### order the matrix by the sample table
		mymat <- mymat[as.character(sample_table$V1), ]

		### calculate stats per sample
		polymorphic <- apply(mymat, 1, function(x) sum(x %in% c("k", "m", "r", "s", "w", "y", "b", "d", "h", "v")))
		results_samp[, "Polymorphic"] <- results_samp[, "Polymorphic"] + polymorphic
		missing <- rowSums(is.na(mymat))
		missing[which(missing == ncol(mymat))] <- 0		# only count missing when not all NA (loci count)
		results_samp[, "Missing"] <- results_samp[, "Missing"] + missing
		hets <- polymorphic / (ncol(mymat) - missing)
		results_samp[, "Ho"] <- results_samp[, "Ho"] + hets

		### calculate stats across the population samples
		results[index, "Samples"] <- samples_present
		results[index, "Sites"] <- ncol(mymat)
		genotypes <- apply(mymat, 2, function(x) length(unique(x[!is.na(x)])))
		results[index, "Polymorphic"] <- sum(genotypes > 1)
		results[index, "Missing"] <- 100 * sum(is.na(mymat)) / length(mymat)
		hets <- apply(mymat, 2, het_measure)		# observed heterozygosity across sites
		results[index, "Ho"] <- mean(hets, na.rm = TRUE)
		gened <- apply(mymat, 2, gd_measure)		# expected heterozygosity across sites
		results[index, "pi"] <- mean(gened, na.rm = TRUE)

		### calculate the fixation index inbreeding coefficient Fis
		### from Nei 1977: Fis = (pi - Ho) / pi
		### Note: this only uses sites with non-zero gene diversity (pi), not all sites
		indices <- which(gened > 0)
		fis <- (gened[indices] - hets[indices]) / gened[indices]
		results[index, "Fis"] <- mean(fis, na.rm = TRUE)

		### increment the index and start on the next locus
		index <- index + 1
	}

	## check populations across files
	if (length(unique(allpops)) > 1) {
		cat("WARNING: more than one population present across fasta files; stats will be wrong!\n")
	}

	## report if there were extra samples not included in calculations
	if (incongruence) {
		cat("Samples were dropped from input alignments if missing from the samples table\n")
	}

	## set up a summary dataframe to store the overall results
	summary <- as.data.frame(matrix(ncol = 17, nrow = 1))
	colnames(summary) <- c("Samples", "Loci", "Sites", "Polymorphic_pop", "Missing_percent",
		"Ho", "Ho_se", "pi", "pi_se", "Fis", "Fis_se", "Samples_loci", "Samples_sites", "Samples_polymorphic",
		"Samples_missing", "Samples_Ho", "Samples_Ho_se")
	rownames(summary) <- allpops[1]

	## calculate overall results
	summary[1, "Samples"] <- mean(results[, "Samples"])
	summary[1, "Loci"] <- length(fasta_list)
	summary[1, "Sites"] <- sum(results[, "Sites"])
	summary[1, "Polymorphic_pop"] <- sum(results[, "Polymorphic"])
	summary[1, "Missing_percent"] <- mean(results[, "Missing"], na.rm = TRUE)
	summary[1, "Ho"] <- mean(results[, "Ho"], na.rm = TRUE)			# observed heterozygosity
	summary[1, "Ho_se"] <- stde(results[, "Ho"])			# Ho standard error
	summary[1, "pi"] <- mean(results[, "pi"], na.rm = TRUE)			# gene diversity or expected heterozygosity
	summary[1, "pi_se"] <- stde(results[, "pi"])			# pi standard error
	summary[1, "Fis"] <- mean(results[, "Fis"], na.rm = TRUE)			# inbreeding coefficient
	summary[1, "Fis_se"] <- stde(results[, "Fis"])			# Fis standard error

	summary[1, "Samples_loci"] <- mean(results_samp[, "Loci"])
	summary[1, "Samples_sites"] <- mean(results_samp[, "Sites"])
	summary[1, "Samples_polymorphic"] <- mean(results_samp[, "Polymorphic"])
	summary[1, "Samples_missing"] <- mean(results_samp[, "Missing"])
	heterozygosity <- results_samp[, "Ho"] / results_samp[, "Loci"]
	summary[1, "Samples_Ho"] <- mean(heterozygosity, na.rm = TRUE)
	summary[1, "Samples_Ho_se"] <- stde(heterozygosity)

	## export the table to file
	write.table(summary, file = paste0(out_pref, "_summary_fasta.txt"),
		quote = FALSE, row.names = TRUE)
}
