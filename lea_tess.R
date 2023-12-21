#######################
# Author: B. Anderson
# Date: Dec 2023
# Description: run population genetic structure analyses using LEA:snmf and optionally TESS3 with geographic data
#######################


# load libraries
suppressMessages(library("LEA"))
suppressMessages(library("tess3r"))


# a help function for when the script is called without arguments or incorrectly
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to run LEA:snmf with a Structure input file (from vcf_to_structure.py) and ")
		cat("optionally TESS3 with geographic data\n\n")
		cat("Usage: Rscript lea_tess.R <options> input.str\n")
		cat("Options:\n")
		cat("\t-k\tMaximum K to run to [default: 7]\n")
		cat("\t-r\tNumber of reps per K [default: 40]\n")
		cat("\t-s\tSamples file with tab-separated columns: sampleID, pop, latitude, longitude (no header) [optional]\n")
		cat("\t\tNOTE: the samples **must** be in the same order as in the Structure file\n")
		cat("\t-t\tFlag for whether to run TESS3 [optional]; if set, needs the samples file with geographic information\n")
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
	extra <- 1
	catch <- TRUE
	maxk <- 7
	reps <- 40
	samples_present <- FALSE
	samples_file <- ""
	run_tess <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-k") {
			maxk <- as.integer(args[index + 1])
			catch <- FALSE
		} else if (args[index] == "-r") {
			reps <- as.integer(args[index + 1])
			catch <- FALSE
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
			catch <- FALSE
		} else if (args[index] == "-t") {
			run_tess <- TRUE
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
if (length(catch_args) < 1) {
	stop(help("Missing Structure input file!\n"), call. = FALSE)
} else if (length(catch_args) > 1) {
	stop(help("Too many inputs specified! There should only be one Structure file\n"), call. = FALSE)
}


# read in and parse the files
if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
	if (ncol(sample_table) != 4) {
		stop(help("Samples table formatted improperly!\n"), call. = FALSE)
	}
}

str_file <- catch_args[[1]]
# the next command will write two new files (extensions *.geno and *.lfmm)
# the extra.col are sampleID and pop from the typical output from vcf_to_structure.py
struct2geno(str_file, ploidy = 2, FORMAT = 2, extra.row = 0, extra.col = 2)
geno_file <- paste0(str_file, ".geno")
lfmm_file <- paste0(str_file, ".lfmm")


if (run_tess) {
	if (! samples_present) {
		stop(help("Missing samples file with geographic information!\n"), call. = FALSE)
	}
	coords <- as.matrix(sample_table[, c(4, 3)])		# TESS3 wants lon lat
}


# run sparse non-negative matrix factorization (snmf) to estimate admixture
cat(paste0("Running snmf for ", reps,  " repetitions per K value from 1 to ", maxk, "...\n"))
invisible(capture.output(myproject <- snmf(geno_file, K = 1: maxk, project = "new", repetitions = reps,
	CPU = 8, entropy = TRUE, percentage = 0.2, ploidy = 2)))


# plot the cross-entropy for each K to estimate an optimal K
pdf("cross-entropy_lea.pdf", width = 7, height = 7)
plot(myproject, col = "#2c2a90", pch = 19)
invisible(dev.off())


# for each K from 2 to maxk, get the run with the lowest cross-entropy
for (kval in 2: maxk) {
	bestrun <- which.min(cross.entropy(myproject, K = kval))
	qmat <- Q(myproject, kval, bestrun)
	qmat <- format(qmat, scientific = FALSE)
	write.table(qmat, file = paste0("qfile_lea", kval, ".txt"), quote = FALSE,
		row.names = FALSE, col.names = FALSE)
}


# output the results for the top 10 runs per K as input to CLUMPAK
# to use CLUMPAK, we need to construct a format it recognises (mimic Structure)
if (samples_present) {
	sampleids <- sample_table[, 1]
	pops <- sample_table[, 2]
} else {
	sampleids <- seq_len(nrow(qmat))
	pops <- rep("1", nrow(qmat))
}
for (kval in 1: maxk) {
	entropy <- cross.entropy(myproject, K = kval)
	runs <- sort(entropy, index.return = TRUE)$ix[1:10]
	for (run in runs) {
		qmat <- Q(myproject, kval, run)
		qmat <- format(qmat, scientific = FALSE)
		outmat <- cbind(seq_len(nrow(qmat)), sampleids,
			rep(paste0("(", 5, ")"), nrow(qmat)), pops,
			rep(":", nrow(qmat)), qmat)
		outfile <- file(paste0("lea_", kval, "_", run, "_f"), open = "w")
		writeLines("# Fake Structure file\n", con = outfile)
		writeLines("Run parameters:", con = outfile)
		writeLines(paste0(" ", nrow(qmat), " individuals"), con = outfile)
		writeLines(paste0(" ", kval, " populations assumed\n"), con = outfile)
		writeLines("Inferred ancestry of individuals:", con = outfile)
		writeLines("    Label  (%Miss)  Pop:  Inferred clusters", con = outfile)
		write.table(outmat, file = outfile, append = TRUE,
			quote = FALSE, row.names = FALSE, col.names = FALSE)
		close(outfile)
	}
}


# use TESS3 to run another estimate using geographic information as well, if requested
if (run_tess) {
	# read the lfmm geno file as a matrix and convert missing data to a code tess3r likes
	input_mat <- as.matrix(read.table(lfmm_file))
	input_mat[input_mat == 9] <- NA

	# run
	# note: may fail if openMP is not installed on your system (in which case, delete that option)
	cat(paste0("Running tess3 gnmf for ", reps, " repetitions per K value from 1 to ", maxk, "...\n"))
	invisible(capture.output(tess_results <- tess3(input_mat, coord = coords, K = 1: maxk, ploidy = 2,
		rep = reps, openMP.core.num = 8, keep = "all", mask = 0.2)))

	# plot the cross-validation for each K
	pdf("cross-validation_tess.pdf", width = 7, height = 7)
	plot(tess_results, col = "#2c2a90", pch = 19,
		xlab = "K", ylab = "Cross-validation score",
		main = "Cross-validation including all runs (not just the best one)")
	invisible(dev.off())

	# output the top 10 runs per K (based on root mean squared error, following tess3r)
	for (kval in 1: maxk) {
		rootmse <- tess_results[[kval]]$rmse
		runs <- sort(rootmse, index.return = TRUE)$ix[1:10]
		for (run in runs) {
			qmat <- qmatrix(tess_results, kval, run)
			qmat <- format(qmat, scientific = FALSE)
			outmat <- cbind(seq_len(nrow(qmat)), sampleids,
				rep(paste0("(", 5, ")"), nrow(qmat)), pops,
				rep(":", nrow(qmat)), qmat)
			outfile <- file(paste0("tess_", kval, "_", run, "_f"), open = "w")
			writeLines("# Fake Structure file\n", con = outfile)
			writeLines("Run parameters:", con = outfile)
			writeLines(paste0(" ", nrow(qmat), " individuals"), con = outfile)
			writeLines(paste0(" ", kval, " populations assumed\n"), con = outfile)
			writeLines("Inferred ancestry of individuals:", con = outfile)
			writeLines("    Label  (%Miss)  Pop:  Inferred clusters", con = outfile)
			write.table(outmat, file = outfile, append = TRUE,
				quote = FALSE, row.names = FALSE, col.names = FALSE)
			close(outfile)
		}
	}

	# determine the best run per K, and map the interpolation
	suppressMessages(library("maps"))
	suppressMessages(library("rworldmap"))
	# note: tess_colours needs to be at least as long as maxk
	tess_colours <- c("tomato", "orange", "lightblue",
		"wheat", "olivedrab", "violet", "gold", "chartreuse",
		"blue", "red", "darkgreen", "pink")
	tess_palette <- CreatePalette(tess_colours, 9)
	for (kval in 2: maxk) {
		bestrun <- which.min(tess_results[[kval]]$rmse)
		qmat <- qmatrix(tess_results, kval, bestrun)
		pdf(paste0("map_tess_", kval, ".pdf"), width = 7, height = 7)
			invisible(capture.output(suppressWarnings(
				plot(qmat, coords, method = "map.max", interpol = FieldsKrigModel(10),
				main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude",
				resolution = c(500, 500), cex = .4,	col.palette = tess_palette)
			)))
		invisible(dev.off())
		qmat <- format(qmat, scientific = FALSE)
		write.table(qmat, file = paste0("qfile_tess", kval, ".txt"), quote = FALSE,
			row.names = FALSE, col.names = FALSE)
	}
}
