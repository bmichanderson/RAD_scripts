
##########
# Author: Ben Anderson
# Date: April 2022
# Description: assess tree height and population divergences for setting priors for SNAPP/ER
# NOTE: the input Nexus file can be produced with `vcf_to_nexus.py`
##########


# load required library (for reading in Nexus)
suppressMessages(library(ape))


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to assess divergences in a SNAPP/ER input Nexus file\n")
		cat("Usage: Rscript SNAPP_divergence.R -s SNAPP/ER_infile.nex -p pops_file\n")
		cat("Options:\n")
		cat("\t-s\tSNAPP/ER input Nexus file\n")
		cat("\t-p\tPopulations file, with sampleID and population tab-separated, one per line\n")
	} else {
		cat(help_message)
	}
}


# a function to calculate pairwise sequence divergence from 012 data
# NOTE: this is calculated only on SNPs (if using that), so values are high
# it compares each row of a dataframe to every other row
# it sums absolute distances (a "1" vs "0", "1", "2" = distance of 1, 0, 1)
# missing data ("?") is ignored (sites dropped for that comparison)
# it then calculates substitutions per site by dividing by number of comparisons
# it returns a distance matrix
seq_div <- function(align_df) {
	out_mat <- matrix(nrow = nrow(align_df), ncol = nrow(align_df))
	rownames(out_mat) <- rownames(align_df)
	colnames(out_mat) <- rownames(align_df)
	index <- 2
	for (row in seq_len(nrow(align_df) - 1)) {
		for (comp_row in seq(index, nrow(align_df))) {
			not_miss <- align_df[row, ] != "?" & align_df[comp_row, ] != "?"
			comps <- align_df[c(row, comp_row), not_miss]
			dist <- sum(abs(as.numeric(comps[1, ]) - as.numeric(comps[2, ])))
			out_mat[comp_row, row] <-  dist / ncol(comps)
		}
		index <- index + 1
	}
	as.dist(out_mat)
}


# a function to calculate lambda
# this is from "A Rough Guide to SNAPP"
# lambda = 1/h * (sum(k=1 to n-1) k/(n*(n - k)))
# h = tree height, n = number of species
get_lambda <- function(height, num) {
	mysum <- 0
	for (kval in seq_len(num - 1)) {
		val <- kval / (num * (num - kval))
		mysum <- mysum + val
	}
	(1 / height) * mysum
}


# Read in the data

# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	snapp_present <- FALSE
	pops_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-p") {
			pops_present <- TRUE
			pops_file <- args[index + 1]
		} else if (args[index] == "-s") {
			snapp_present <- TRUE
			snapp_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! snapp_present || ! pops_present) {
	stop(help("Missing argument for SNAPP/ER file and/or populations file!\n"), call. = FALSE)
}


# read in the alignment as a list, and convert to a dataframe (samples as rows)
mynex <- read.nexus.data(snapp_file)
mydf <- t(as.data.frame(mynex))
rownames(mydf) <- names(mynex)		# for some reason it replaces "-" with "." otherwise


# read in the populations file as a table
pops_table <- read.table(pops_file, header = FALSE)
colnames(pops_table) <- c("sample", "pop")


# calculate pairwise sequence divergence as substitutions per site
# this calculation takes the longest (but not too bad)
subs_per_site <- seq_div(mydf)


# report possible lambda values for different numbers of species
# need maximum divergence divided by two (height of the tree)
height <- round(max(subs_per_site) / 2, 2)
cat(paste0("Tree height: ", height, "\n"))
cat("Lambda values for numbers of species:\n")
for (num in seq(2, 15)) {
	lambda <- round(get_lambda(height, num), 1)
	cat(paste0(num, "\t", lambda, "\n"))
}


# report average divergence within each population
cat("\nAverage substitutions per site between members of a population:\n")
pops_to_subset <- unique(pops_table$pop[pops_table$sample %in% rownames(mydf)])
subs_df <- as.data.frame(as.matrix(subs_per_site))
result_df <- as.data.frame(matrix(nrow = length(pops_to_subset), ncol = 1))
rownames(result_df) <- pops_to_subset
colnames(result_df) <- c("avg")
for (pop in pops_to_subset) {
	samples <- pops_table$sample[pops_table$pop == pop]
	sub_mat <- as.matrix(subs_df[rownames(subs_df) %in% samples, colnames(subs_df) %in% samples])
	result_df[pop, ] <- round(mean(as.dist(sub_mat)), 3)
}
result_df
summary(result_df)
