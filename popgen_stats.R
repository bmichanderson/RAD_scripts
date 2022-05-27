
##########
# Author: Ben Anderson with ideas from Rachel Binks
# Date: Dec 2021
# Modified: May 2022
# Description: calculate various popgen stats for a VCF file
##########


# load required libraries
suppressMessages(library(vcfR))
suppressMessages(library(adegenet))
suppressMessages(library(hierfstat))
suppressMessages(library(StAMPP))


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate various popgen stats from an input VCF file\n")
		cat("Usage: Rscript popgen_stats.R -o out_pref -v vcf_file -s samples_file -t [run Fst calculations]\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-s\tSamples and populations as tab-delimited sample IDs and pops, one per line\n")
		cat("\t-t\tRun pairwise Fst calculations between populations (default: do not)\n")
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
		} else if (args[index] == "-t") {
			run_fst <- TRUE
			cat("Will run Fst calculations\n")
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}
if (any(c(! vcf_present, ! samples_present))) {
	stop(help("Missing argument for vcf file and/or samples file!\n"), call. = FALSE)
}


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
vcf <- read.vcfR(vcf_file, verbose = FALSE)
cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
	length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")


# check that the samples in the VCF are in the table
for (indiv in colnames(vcf@gt)[2: ncol(vcf@gt)]) {
	if (! indiv %in% sample_table$V1) {
		stop(help("Sample missing from samples file! Stopping\n"), call. = FALSE)
	}
}

# calculate amount of missing data, and average per pop
gt <- extract.gt(vcf, element = "GT")
missing <- apply(gt, MARGIN = 2, function(x) { sum(is.na(x)) })
missing <- 100 * missing / nrow(vcf)
misspops <- sample_table$V2[match(names(missing), sample_table$V1)]
names(missing) <- misspops
meanmiss <- tapply(missing, names(missing), mean)
sample_size <- table(names(missing))


# convert the VCF into a hierfstat dataframe with pop and genotype
geni <- vcfR2genind(vcf, return.alleles = TRUE)
populations <- sample_table$V2[match(rownames(geni@tab), sample_table$V1)]
geni@pop <- as.factor(populations)
my_hfst <- genind2hierfstat(geni)


# calculate basic stats and create a summary dataframe
cat("Calculating and plotting basic popgen stats\n")
bstats <- basic.stats(my_hfst)
# if there was only one pop, hierfst adds a dummy pop,
# which breaks my code, so only keep the first column
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


# populate the dataframe with values
summary[, "Samples"] <- sample_size
summary[, "PercentMissing"] <- meanmiss
summary[, "Ho"] <- apply(bstats$Ho, 2, function(x) mean(x, na.rm = TRUE))
summary[, "Ho_SE"] <- apply(bstats$Ho, 2, stde)
summary[, "Hs"] <- apply(bstats$Hs, 2, function(x) mean(x, na.rm = TRUE))
summary[, "Hs_SE"] <- apply(bstats$Hs, 2, stde)
summary[, "Fis"] <- apply(bstats$Fis, 2, function(x) mean(x, na.rm = TRUE))
summary[, "Fis_SE"] <- apply(bstats$Fis, 2, stde)


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
barplot(summary$PercentMissing, las = 2, main = "Mean percentage missing",
	names = rownames(summary))

# graph the values by population
mybarplot(summary$Ho, summary$Ho_SE, names = rownames(summary),
	main = "Ho, observed heterozygosity")
mybarplot(summary$Hs, summary$Hs_SE, names = rownames(summary),
	main = "Hs, estimated gene diversity\n(expected heterozygosity)")
mybarplot(summary$Fis, summary$Fis_SE, names = rownames(summary),
	main = "Inbreeding coefficient Fis\n(1 - Ho / Hs)")

# stop creating the pdf
invisible(dev.off())


# Run Fst calculations if requested
if (run_fst) {
	# calculate pairwise Fst
	cat("Calculating pairwise Fst\n")

	# The faster way is via a genlight and StAMPP
	# convert the vcf to a genlight, then add pop
	# NOTE: this will remove sites that are not biallelic
	genl <- vcfR2genlight(vcf)
	pop(genl) <- populations

	# now use the genlight in StAMPP
	# this calculates Fst following Weir and Cockerham 1984
	# if wanting to get confidence intervals, one could use e.g.
	# nboots = 100, percent = 95, nclusters = 4
	fsts <- stamppFst(genl, nboots = 1)

	# export the pairwise Fst values to file
	write.table(fsts, file = paste0(out_pref, "_StAMPP_Fst.txt"),
		quote = FALSE, row.names = TRUE)


	# A slower way is using the hierfstat package
	# this calculates Fst following Weir and Cockerham 1984
	# since this doesn't convert to genlight, it keeps all sites
	cat("And more slowly with hierfstat...\n")
	h_fsts <- pairwise.WCfst(my_hfst)

	# make the resulting square matrix a lower triangle matrix
	h_fsts[upper.tri(h_fsts)] <- NA

	# export the pairwise Fst values to file
	write.table(h_fsts, file = paste0(out_pref, "_hierfstat_Fst.txt"),
		quote = FALSE, row.names = TRUE)
}
