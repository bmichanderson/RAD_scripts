
##########
# Author: Ben Anderson with ideas from Rachel Binks
# Date: Dec 2021
# Description: calculate various popgen stats for a VCF file
##########


# load required libraries
suppressMessages(library(adegenet))
suppressMessages(library(vcfR))
suppressMessages(library(StAMPP))
suppressMessages(library(hierfstat))


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate various popgen stats from an input VCF file\n")
		cat("Usage: Rscript popgen_stats.R -o out_pref -v vcf_file -s samples_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-s\tSamples and populations as tab-delimited sample IDs and pops, one per line\n")
	} else {
		cat(help_message)
	}
}


# a function to calculate standard error
# see: https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
# and: https://www.rdocumentation.org/packages/plotrix/versions/3.8-2/topics/std.error
se <- function(x) {
	sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))
}


# a function to add an error bar
# from https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
errorbar <- function(x, y, upper, lower = upper, length = 0.1, ...) {
	arrows(x, y + upper, x, y - lower, angle = 90, code = 3, length = length, ...)
}


# a function to plot a barplot with error bars
mybarplot <- function(vec_val, vec_se, names, ...) {
	myplot <- barplot(vec_val, names = names, las = 2, ...)
	errorbar(myplot, vec_val, vec_se)
}


# Read in and format the data

# parse the command line
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
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			out_pref <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present <- TRUE
			vcf_file <- args[index + 1]
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
if (vcf_present) {
	vcf <- read.vcfR(vcf_file, verbose = FALSE)
	cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
		length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")
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
geni <- vcfR2genind(vcf)
populations <- sample_table$V2[match(rownames(geni@tab), sample_table$V1)]
geni@pop <- as.factor(populations)
my_hfst <- genind2hierfstat(geni)


# calculate basic stats and create a summary dataframe
print("Calculating and plotting basic popgen stats")
bstats <- basic.stats(my_hfst)
summary <- as.data.frame(matrix(ncol = 6, nrow = length(unique(populations))))
colnames(summary) <- c("ho", "hose", "hs", "hsse", "fis", "fisse")
rownames(summary) <- colnames(bstats$Ho)


# populate the dataframe with values
summary$ho <- apply(bstats$Ho, 2, function(x) mean(x, na.rm = TRUE))
summary$hose <- apply(bstats$Ho, 2, se)
summary$hs <- apply(bstats$Hs, 2, function(x) mean(x, na.rm = TRUE))
summary$hsse <- apply(bstats$Hs, 2, se)
summary$fis <- apply(bstats$Fis, 2, function(x) mean(x, na.rm = TRUE))
summary$fisse <- apply(bstats$Fis, 2, se)


# export the table to file
write.table(summary, file = paste0(out_pref, "_summary.txt"),
			quote = FALSE, row.names = TRUE)


# start creating a pdf of plots
pdf(paste0(out_pref, "_summary.pdf"), width = 10, height = 10)


# graph the amount of missing data
barplot(meanmiss, las = 2, main = "Mean missing data")


# graph the sample sizes by pop
barplot(sample_size, las = 2, main = "Sample size")


# graph the values by population
# also, use a confidence interval of +/- 2 * SE
mybarplot(summary$ho, summary$hose * 2, names = rownames(summary),
			main = "Ho, observed heterozygosity",
			ylim = c(0, 1.2 * max(summary$ho)))
mybarplot(summary$hs, summary$hsse * 2, names = rownames(summary),
			main = "Hs, observed gene diversity",
			ylim = c(0, 1.2 * max(summary$hs)))
mybarplot(summary$fis, summary$fisse * 2, names = rownames(summary),
			main = "Fis = 1 - Ho/Hs",
			ylim = c(1.2 * min(summary$fis), 1.2 * max(summary$fis)))


# stop creating the pdf
invisible(dev.off())


# calculate pairwise Fst
print("Calculating and plotting pariwise Fst")

# first, convert the vcf to a genlight, then add pop
genl <- vcfR2genlight(vcf)
pop(genl) <- populations

# now use the genlight in StAMPP
fsts <- stamppFst(genl, nboots = 100, percent = 95, nclusters = 4)


# order the resulting pairwise matrix by row
# NOTE: could make it the order of the input pop file?
myfsts <- fsts$Fsts
ord <- sort(rownames(myfsts))
myfsts[upper.tri(myfsts)] <- t(myfsts)[upper.tri(myfsts)]
myfsts <- myfsts[ord, ord]
myfsts[upper.tri(myfsts)] <- NA


# export the pairwise Fst values to file
write.table(myfsts, file = paste0(out_pref, "_Fst.txt"),
			quote = FALSE, row.names = FALSE)


# start creating a pdf of plots
pdf(paste0(out_pref, "_fst.pdf"), width = 10, height = 10)


# plot a heatmap
par(mar = c(5, 4, 4, 5) + 0.1)
image(myfsts, col = hcl.colors(n = 100, palette = "greens", rev = TRUE),
		axes = FALSE, main = "Pairwise Fst")
axis(1, at = seq(0, 1, length.out = ncol(myfsts)),
	labels = colnames(myfsts), las = 2, lwd.ticks = 0)
axis(4, at = seq(0, 1, length.out = ncol(myfsts)),
	labels = colnames(myfsts), las = 2, lwd.ticks = 0)


# stop creating the pdf
invisible(dev.off())
