
##########
# Author: B.M. Anderson
# Date: Feb 2022
# Modified: Nov 2023 (account for pop names as numbers), Apr 2025 (reduced screen output)
# Description: compare a distance matrix (e.g. Fst) to geographic distance
##########


# load required libraries
suppressMessages(library(geosphere))
suppressMessages(library(vegan))


# Define functions

## a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to compare distances to geographic distances\n")
		cat("Usage: Rscript isolation_by_distance.R -o out_pref -d dist_file -l loc\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-d\tThe distance matrix to be compared\n")
		cat("\t-l\tThe localities of each unit in the matrix, tab delimited as sample_ID (identically labelled) lat lon\n")
		cat("\t-m\tFlag to run a Mantel test [default do not]\n")
		cat("\t-r\tFlag to run a linear regression on transformed data [default do not]\n")
		cat("\t-s\tA file with samples to include (unique names, one per line) [optional]")
		cat("\t-z\tConvert negative genetic distances to 0 [default: do not]")
	} else {
		cat(help_message)
	}
}


# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	out_pref <- "output"
	dist_present <- FALSE
	loc_present <- FALSE
	samples_present <- FALSE
	run_mantel <- FALSE
	run_regression <- FALSE
	convert_zero <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			out_pref <- args[index + 1]
		} else if (args[index] == "-d") {
			dist_present <- TRUE
			dist_file <- args[index + 1]
		} else if (args[index] == "-l") {
			loc_present <- TRUE
			loc_file <- args[index + 1]
		} else if (args[index] == "-m") {
			run_mantel <- TRUE
		} else if (args[index] == "-r") {
			run_regression <- TRUE
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			sample_file <- args[index + 1]
		} else if (args[index] == "-z") {
			convert_zero <- TRUE
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}
if (any(c(! dist_present, ! loc_present))) {
	stop(help("Missing argument for distance matrix and/or locations file!\n"), call. = FALSE)
}


# read in the data; filter the distance matrix to only keep samples if specified
dist_mat <- read.csv(dist_file, sep = " ", header = TRUE, check.names = FALSE)	# in case populations are numbers
loc_table <- read.table(loc_file, sep = "\t", header = FALSE)
if (samples_present) {
	sample_table <- read.table(sample_file, header = FALSE)
	sample_table[1] <- lapply(sample_table[1], as.character)		# convert if populations are numbers
	dist_mat <- dist_mat[rownames(dist_mat) %in% sample_table[, 1], colnames(dist_mat) %in% sample_table[, 1]]
}


# check that the items in the dist_mat have entries in the loc_table
if (any(! rownames(dist_mat) %in% loc_table$V1)) {
	stop(help("Entries in the distance matrix do not have location coordinates!\n"), call. = FALSE)
}


# convert the distance matrix to square (account for multiple input orientations)
if (is.na(dist_mat[1, 2])) {	# lower triangle matrix
	dist_mat[upper.tri(dist_mat)] <- t(dist_mat)[upper.tri(dist_mat)]
} else if (is.na(dist_mat[2, 1])) {		# upper triangle matrix
	dist_mat[lower.tri(dist_mat)] <- t(dist_mat)[lower.tri(dist_mat)]
} else {
	# assume a square distance matrix input
}
diag(dist_mat) <- NA


# convert negative distances to zeros, if requested
if (convert_zero) {
	dist_mat[dist_mat < 0] <- 0
}


# calculate pairwise geographic distances and store as a matrix
geo_mat <- matrix(nrow = nrow(dist_mat), ncol = ncol(dist_mat))
rownames(geo_mat) <- rownames(dist_mat)
colnames(geo_mat) <- colnames(dist_mat)
for (index in seq_len(nrow(dist_mat) - 1)) {
	entry <- rownames(dist_mat)[index]
	lat1 <- loc_table[, 2][loc_table[, 1] == entry]
	lon1 <- loc_table[, 3][loc_table[, 1] == entry]
	for (index2 in (index + 1): nrow(dist_mat)) {
		entry2 <- rownames(dist_mat)[index2]
		lat2 <- loc_table[, 2][loc_table[, 1] == entry2]
		lon2 <- loc_table[, 3][loc_table[, 1] == entry2]
		# calculate distance
		mydist <- distm(c(lon1, lat1), c(lon2, lat2), fun = distGeo)
		# add to the matrix as km distances
		geo_mat[index, index2] <- mydist / 1000
		geo_mat[index2, index] <- mydist / 1000
	}
}


# store the geo_mat as log-transformed as well
# if there are zero values, then there will be a problem, so do a log(x+1) transformation
if (sum(geo_mat == 0, na.rm = TRUE) > 0) {
	log_geo_mat <- log(geo_mat + 1)
} else {
	log_geo_mat <- log(geo_mat)
}


# if the distances are Fst, then make a correction for fitting to the log geo
alter_dist_mat <- dist_mat / (1 - dist_mat)


# run a Mantel test if requested
if (run_mantel) {
	mymant <- mantel(dist_mat, geo_mat, method = "pearson",
		permutations = 999, na.rm = TRUE)
	cat("Mantel test statistic from", mymant$permutations, "permutations is",
		mymant$statistic, "with a significance of", mymant$signif, "\n")
}


# run a linear regression if requested
if (run_regression) {
	lmodel <- lm(dist_mat[lower.tri(dist_mat, diag = FALSE)] ~
				geo_mat[lower.tri(geo_mat, diag = FALSE)])
	rsquared <- summary(lmodel)$r.squared
	adjrsquared <- format(round(summary(lmodel)$adj.r.squared, 2), nsmall = 2)
	pvalue <- format(summary(lmodel)$coefficients[2, 4], digits = 2)

	# start making a pdf
	pdf(paste0(out_pref, "_IBD.pdf"), width = 7, height = 7)

	# plot the distances by geographic distance with the regression line
	par(mar = c(5.1, 6.1, 5.1, 2.1), cex.lab = 2)
	xvar <- geo_mat[lower.tri(geo_mat, diag = FALSE)]
	yvar <- dist_mat[lower.tri(dist_mat, diag = FALSE)]
	plot(x = xvar, y = yvar, xlab = "Geographic distance (km)",
		ylab = expression("F"[ST]), cex.axis = 1.5, pch = 19)
	abline(lmodel)
	ltext <- bquote(atop(R^2 ~ "=" ~ .(adjrsquared), "P =" ~ .(pvalue)))
	legend("topleft", legend = ltext, bty = "n", cex = 1.5, xjust = 0)

	# if running with Fst, then may want to do another regression with log geo dist
	lmodel <- lm(alter_dist_mat[lower.tri(dist_mat, diag = FALSE)] ~
				log_geo_mat[lower.tri(geo_mat, diag = FALSE)])
	rsquared <- summary(lmodel)$r.squared
	adjrsquared <- format(round(summary(lmodel)$adj.r.squared, 2), nsmall = 2)
	pvalue <- format(summary(lmodel)$coefficients[2, 4], digits = 2)

	# plot the altered distances by log geographic distance with the regression line
	xvar <- log_geo_mat[lower.tri(geo_mat, diag = FALSE)]
	yvar <- alter_dist_mat[lower.tri(dist_mat, diag = FALSE)]
	plot(x = xvar, y = yvar, xlab = "Geographic distance (ln(km))",
		ylab = expression("F"[ST] * "/(1 - F"[ST] * ")"), cex.axis = 1.5, pch = 19)
	abline(lmodel)
	ltext <- bquote(atop(R^2 ~ "=" ~ .(adjrsquared), "P =" ~ .(pvalue)))
	legend("topleft", legend = ltext, bty = "n", cex = 1.5, xjust = 0)

	# stop the pdf
	invisible(dev.off())
}
