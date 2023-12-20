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
		cat("\t-s\tSamples file with tab-separated columns: sampleID, latitude, longitude (no header) [optional]\n")
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
	if (ncol(sample_table) != 3) {
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
	coords <- sample_table[, c(3, 2)]		# TESS3 wants lon lat
}


# run sparse non-negative matrix factorization (snmf) to estimate admixture
myproject <- snmf(geno_file, K = 1: maxk, project = "new", repetitions = reps,
	CPU = 8, entropy = TRUE, percentage = 0.2, ploidy = 2)


# plot the cross-entropy for each K to estimate an optimal K
pdf("cross-entropy.pdf", width = 7, height = 7)
plot(myproject, col = "#2c2a90", pch = 19)
invisible(dev.off())


# for each K from 2 to "maxk", get the run with the lowest cross-entropy
for (kval in 2: maxk) {
	bestrun <- which.min(cross.entropy(myproject, K = kval))
	qmat <- Q(myproject, kval, bestrun)
	write.table(qmat, file = paste0("qfile", kval, ".txt"), quote = FALSE,
		row.names = FALSE, col.names = FALSE)
}


# alternatively, output the results for the top 10 runs per K as input to CLUMPAK

##### NOTE: need pop labels and a table like in Structure output

for (kval in 1: maxk) {
	entropy <- cross.entropy(myproject, K = kval)
	runs <- sort(entropy, index.return = TRUE)$ix[1:10]
	for (run in runs) {
		qmat <- Q(myproject, kval, run)
		qmat <- format(qmat, scientific = FALSE)
		outmat <- cbind(seq_len(nrow(qmat)), rep(kval, nrow(qmat)),
			rep(paste0("(", 5, ")"), nrow(qmat)),
			rep(kval, nrow(qmat)), rep(":", nrow(qmat)), qmat)
		write.table(outmat, file = paste0("result_", kval, "_", run, "_f"),
			quote = FALSE, row.names = FALSE, col.names = FALSE)
	}
}


# use TESS

