
##########
# Author: Ben Anderson
# Date: Jan 2022
# Description: calculate various popgen stats for a fasta alignment
##########


# load required libraries
suppressMessages(library(ape))


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate various popgen stats from an input fasta alignment\n")
		cat("Usage: Rscript popgen_stats_fasta.R -o out_pref -f fasta_file -s samples_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-f\tThe fasta file to be analysed (names matching samples file exactly)\n")
		cat("\t-s\tSamples and populations as tab-delimited sample IDs and pops, one per line\n")
		cat("\t-t\tRun Fst calculations (default: do not)\n")
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


# a function to calculate observed heterozygosity
# this is the proportion of samples at each site that are hets
het_measure <- function(x) {
		y <- x[!is.na(x)]		# only account for non NA
		if (length(y) == 0) {
			NA
		} else {
			sum(!(y %in% c("a", "c", "g", "t"))) / length(y)
		}
	}


# a function to calculate gene diversity
# this is the average number of differences between sequences
# or the probability of choosing an allele that is different
# we can calculate this with either:
# 1) Nei & Roychoudhury 1974: (n/n-1) * (1 - sum[(ni/n)^2])
# 2) Hohenlohe et al. 2010: 1 - sum[(ni-1)*ni]/(n-1)*n
# where ni = count of allele i; n = sum of allele counts
# we'll use 2
gd_measure <- function(x) {
	y <- x[!is.na(x)]		# only account for non NA
	if (length(y) == 0) {
		NA
	} else {
		# first, turn the genotypes into two alleles
		z <- vector("character", length(y) * 2)
		index <- 1
		for (gt in y) {
			if (gt %in% c("a", "c", "g", "t")) {		# homo
				z[index] <- gt
				z[index + 1] <- gt
				index <- index + 2
			} else {		# het
				if (gt == "k") {
					b1 <- "g"
					b2 <- "t"
				} else if (gt == "m") {
					b1 <- "a"
					b2 <- "c"
				} else if (gt == "r") {
					b1 <- "a"
					b2 <- "g"
				} else if (gt == "s") {
					b1 <- "c"
					b2 <- "g"
				} else if (gt == "w") {
					b1 <- "a"
					b2 <- "t"
				} else if (gt == "y") {
					b1 <- "c"
					b2 <- "t"
				} else {
					stop(help("Unrecognized base!\n"), call. = FALSE)
				}
				z[index] <- b1
				z[index + 1] <- b2
				index <- index + 2
			}
		}
		# now, count alleles and multiply
		n <- length(z)
		calcs <- 0
		for (allele in unique(z)) {
			count <- sum(z == allele)
			calc <- (count - 1) * count
			calcs <- calcs + calc
		}
		# return the gene diversity measure
		1 - (calcs / ((n - 1) * n))
	}
}


# a function to count alleles at each locus for later Fst
# it returns the result of 2 * (n choose 2) = (n - 1) * n
n_measure <- function(x) {
	y <- x[!is.na(x)]		# only account for non NA
	if (length(y) == 0) {
		NA
	} else {
		n <- length(y) * 2
		# return n choose 2 (without the 1/2 because cancelled later)
		(n - 1) * n
	}
}


# a function to calculate the hierfstat version of gene diversity
# hierfstat uses an additional correction from Nei 1987
# it also bases the stat on genotypic number rather than allelic
# hs = (n/n-1) * (1 - sum[(ni/n)^2] - ho/2n)
# I'm now no longer using this
hs_measure <- function(x) {
	y <- x[!is.na(x)]		# only account for non NA
	if (length(y) == 0) {
		NA
	} else {
		# first, calculate observed heterozygosity
		ho <- sum(!(y %in% c("a", "c", "g", "t"))) / length(y)
		# second, turn the genotypes into alleles for frequency calcs
		z <- vector("character", length(y) * 2)
		index <- 1
		for (gt in y) {
			if (gt %in% c("a", "c", "g", "t")) {		# homo
				z[index] <- gt
				z[index + 1] <- gt
				index <- index + 2
			} else {		# het
				if (gt == "k") {
					b1 <- "g"
					b2 <- "t"
				} else if (gt == "m") {
					b1 <- "a"
					b2 <- "c"
				} else if (gt == "r") {
					b1 <- "a"
					b2 <- "g"
				} else if (gt == "s") {
					b1 <- "c"
					b2 <- "g"
				} else if (gt == "w") {
					b1 <- "a"
					b2 <- "t"
				} else if (gt == "y") {
					b1 <- "c"
					b2 <- "t"
				} else {
					stop(help("Unrecognized base!\n"), call. = FALSE)
				}
				z[index] <- b1
				z[index + 1] <- b2
				index <- index + 2
			}
		}
		# now, count alleles and multiply
		n <- length(z)
		calcs <- 0
		for (allele in unique(z)) {
			count <- sum(z == allele)
			calc <- (count / n) ^ 2
			calcs <- calcs + calc
		}
		# return the gene diversity measure
		# hs = (n/n-1) * (1 - sum[(ni/n)^2] - ho/2n)
		# this time, we only use sample size based on genotype
		ng <- length(y)
		(ng / (ng - 1)) * (1 - calcs - ho / (2 * ng))
	}
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


# correct the DNA so that missing data is coded as NA
mymat <- as.matrix(as.character(fasta))
mymat[mymat == "?"] <- NA
mymat[mymat == "n"] <- NA
mymat[mymat == "-"] <- NA
mydnabin <- as.DNAbin(mymat)


# set the rownames of the DNAbin object to populations
# this requires that the fasta file have the same sample names as in the table
populations <- sample_table$V2[match(rownames(mydnabin), sample_table$V1)]
rownames(mydnabin) <- populations


# set up a summary dataframe to store the resulting calculation outputs
summary <- as.data.frame(matrix(ncol = 8, nrow = length(unique(populations))))
colnames(summary) <- c("n", "miss", "ho", "hose", "he",
					"hese", "fis", "fisse")
rownames(summary) <- unique(populations)


# set up convenience dataframes to store the n counts/choose for each locus for each pop (for Fst)
# and to store the gene diversities for each locus for each pop
if (run_fst) {
	n_choose <- as.data.frame(matrix(ncol = ncol(mydnabin), nrow = length(unique(populations))))
	rownames(n_choose) <- unique(populations)
	gd_df <- as.data.frame(matrix(ncol = ncol(mydnabin), nrow = length(unique(populations))))
	rownames(gd_df) <- unique(populations)
}

# cycle through populations and calculate stats
index <- 1
for (pop in unique(populations)) {
	cat("Analyzing population", index, "of", length(unique(populations)), "\n")
	index <- index + 1
	thisdnabin <- mydnabin[rownames(mydnabin) == pop, ]
	pop_sequences <- as.character(thisdnabin)

	# calculate sample size and amount of missing data
	samples <- length(rownames(thisdnabin))
	summary[pop, "n"] <- samples
	num_miss <- sum(is.na(pop_sequences))
	perc_miss <- 100 * num_miss / length(pop_sequences)
	summary[pop, "miss"] <- perc_miss

	# calculate observed heterozygosity
	hets <- matrix(ncol = length(thisdnabin[1, ]), nrow = 1)
	rownames(hets) <- c("ho")
	hets["ho", ] <- apply(pop_sequences, 2, het_measure)
	summary[pop, "ho"] <- mean(hets["ho", ], na.rm = TRUE)
	summary[pop, "hose"] <- se(hets["ho", ])

	# calculate gene/nucleotide diversity (expected heterozygosity) and store
	gd <- matrix(ncol = length(thisdnabin[1, ]), nrow = 2)
	rownames(gd) <- c("he", "hs")
	gd["he", ] <- apply(pop_sequences, 2, gd_measure)
	summary[pop, "he"] <- mean(gd["he", ], na.rm = TRUE)
	summary[pop, "hese"] <- se(gd["he", ])

	if (run_fst) {
		# store the diversity
		gd_df[pop, ] <- gd["he", ]
		# calculate the allele counts per locus and choose 2, and store
		n_choose[pop, ] <- apply(pop_sequences, 2, n_measure)
	}

	# calculate the fixation index inbreeding coefficient Fis
	# from Nei 1977, this is: (hs - ho)/hs
	# or 1 - ho/hs
	# note: this only uses sites with gene diversity (not all sites)
	fis <- matrix(ncol = length(thisdnabin[1, ]), nrow = 2)
	rownames(fis) <- c("fis", "fiss")
	for (col in seq_len(ncol(gd))) {
		if (is.na(gd["he", col])) {
			fis["fis", col] <- NA
		} else {
			if (gd["he", col] == 0) {		# can't divide by zero or use site
				fis["fis", col] <- NA
			} else {
				fis["fis", col] <- 1 - hets["ho", col] / gd["he", col]
			}
		}
	}
	summary[pop, "fis"] <- mean(fis["fis", ], na.rm = TRUE)
	summary[pop, "fisse"] <- se(fis["fis", ])
}


# export the table to file
write.table(summary, file = paste0(out_pref, "_summary.txt"),
			quote = FALSE, row.names = TRUE)


# start creating a pdf of plots
pdf(paste0(out_pref, "_summary.pdf"), width = 10, height = 10)


# graph the sample sizes by pop
barplot(summary$n, las = 2, main = "Sample size", names = rownames(summary))


# graph the amount of missing data
barplot(summary$miss, las = 2, main = "Missing data", names = rownames(summary))


# graph the values by population
# also, use a confidence interval of +/- 2 * SE
mybarplot(summary$ho, summary$hose * 2, names = rownames(summary),
			main = "Ho, observed heterozygosity",
			ylim = c(0, 1.2 * max(summary$ho + summary$hose * 2)))
mybarplot(summary$he, summary$hese * 2, names = rownames(summary),
			main = "He, estimated gene diversity\n(expected heterozygosity)",
			ylim = c(0, 1.2 * max(summary$he + summary$hese * 2)))
mybarplot(summary$fis, summary$fisse * 2, names = rownames(summary),
			main = "Inbreeding coefficicient Fis\n(1 - Ho/He)",
			ylim = c(1.2 * min(c(0, summary$fis - summary$fisse * 2)),
					1.2 * max(c(0, summary$fis + summary$fisse * 2))))


# stop creating the pdf
invisible(dev.off())


if (run_fst) {
	# find a way to calculate pairwise Fst for sequence data
	# complicated --> figure out how to apply Weir and Cockerham 1984...
	# another way would be Fst = 1 - Hs/Ht (??)
	# this doesn't correct for small/different sample sizes
	# we could try to implement the approach taken by Hohenlohe et al. 2010
	# Fst = 1 - [(nj - 1) * nj * hej + (nk - 1) * nk * hek] / heT * [(nj - 1) * nj + (nk - 1) * nk]
	# where heT is gd_measure of the combined two populations
	pops <- unique(populations)
	fst_mat <- matrix(ncol = length(pops), nrow = length(pops))
	rownames(fst_mat) <- pops
	colnames(fst_mat) <- pops
	for (index in seq_len(length(pops) - 1)) {
		cat("Calculating Fst for pop", index, "versus ")
		pop <- pops[index]
		for (index2 in (index + 1): length(pops)) {
			cat(index2, "")
			pop2 <- pops[index2]
			# determine total diversity
			thisdnabin <- mydnabin[rownames(mydnabin) %in% c(pop, pop2), ]
			pop_sequences <- as.character(thisdnabin)
			total_he <- apply(pop_sequences, 2, gd_measure)
			total_he[total_he == 0] <- NA				# set to NA to avoid 0 in denominator
			# only want to look at sites with non-NA
			indices <- which(!is.na(total_he))
			het <- total_he[indices]
			he1 <- gd_df[pop, ][indices]
			n1 <- n_choose[pop, ][indices]
			he2 <- gd_df[pop2, ][indices]
			n2 <- n_choose[pop2, ][indices]
			# calculate pairwise Fst for each locus
			fsts <- 1 - ((n1 * he1 + n2 * he2) / (het * (n1 + n2)))
			# record the mean
			fst_mat[index, index2] <- mean(as.matrix(fsts), na.rm = TRUE)
		}
		cat("\n")
	}


	# order the resulting pairwise matrix by row
	myfsts <- fst_mat
	ord <- sort(rownames(myfsts))
	myfsts[lower.tri(myfsts)] <- t(myfsts)[lower.tri(myfsts)]
	myfsts <- myfsts[ord, ord]
	myfsts[upper.tri(myfsts)] <- NA
	diag(myfsts) <- 0


	# export the pairwise Fst values to file
	write.table(myfsts, file = paste0(out_pref, "_Fst.txt"),
				quote = FALSE, row.names = TRUE)


	# start creating a pdf of plots
	pdf(paste0(out_pref, "_fst.pdf"), width = 10, height = 10)


	# plot a heatmap
	# for making legend, see https://stackoverflow.com/a/13355440 and https://stackoverflow.com/a/70522655
	par(mar = c(5, 4, 4, 5) + 0.1)
	# heatmap
	image(myfsts, col = c("#FFFFFF", hcl.colors(n = 100, palette = "greens", rev = TRUE)),
			axes = FALSE, main = expression("Pairwise F"[ST]), useRaster = TRUE)
	axis(1, at = seq(0, 1, length.out = ncol(myfsts)),
		labels = colnames(myfsts), las = 2, lwd.ticks = 0)
	axis(4, at = seq(0, 1, length.out = ncol(myfsts)),
		labels = colnames(myfsts), las = 2, lwd.ticks = 0)
	# legend
	subx <- grconvertX(c(0, 0.2), from = "user", to = "ndc")
	suby <- grconvertY(c(0.5, 1), from = "user", to = "ndc")
	op <- par(fig = c(subx, suby),
			mar = c(1, 1, 1, 0),
			new = TRUE)
	legend_colours <- as.raster(c(hcl.colors(n = 100, palette = "greens"), "#FFFFFF"))
	legend_seq <- seq(0, max(myfsts, na.rm = TRUE), length = 5)
	legend_labels <- format(round(legend_seq, 2), nsmall = 2)
	plot(x = c(0, 2), y = c(0, 1), type = "n",
		axes = FALSE, xlab = "", ylab = "", main = "")
	axis(side = 4, at = seq(0, 1, length = 5), pos = 1, labels = FALSE,
		col = 0, col.ticks = 1)
	mtext(legend_labels, side = 4, line = -0.5, at = seq(0, 1, length = 5), las = 2)
	rasterImage(legend_colours, xleft = 0, ybottom = 0, xright = 1, ytop = 1)
	par(op)


	# stop creating the pdf
	invisible(dev.off())
}
