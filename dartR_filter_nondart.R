
##########
# Authors: Rachel Binks and Ben Anderson
# Date: June 2021
# Description: filter VCF fiile of SNPs from Stacks or ipyrad using the dartR package
##########


# set parameters needed for output

#param1 <- X		# param1 description


# a helper function for when the script is called without arguments
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to filter a VCF file of RADseq data\n")
		cat("Usage: Rscript dartR_filter_nondart.R -o output -s sample_file -v vcf_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-s\tThe sample names and pops file if vcf labels differ, with the format:\n",
			"\t\tvcf_label\t*tab*\tsample_name\t*tab*\tpopulation\n",
			"\t\tNote: the name should have an \"_R\" at the end for a replicate\n",
			"\t\tNote2: if sample names match vcf_labels, you only need one column then pop\n")
		cat("\t-v\tThe vcf file to be analysed\n")
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

	output <- "output"
	samples_present <- FALSE
	vcf_present <- FALSE

	for (index in 1:length(args)) {
		if (args[index] == "-o") {
			output <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present = TRUE
			vcf_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! samples_present) {
	stop(help("Missing argument for sample/pops file!\n"), call. = FALSE)
}

if (! vcf_present) {
	stop(help("Missing argument for vcf file!\n"), call. = FALSE)
}


# load required libraries
library(vcfR)
library(dartR)
library(ape)


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)

vcf <- read.vcfR(vcf_file, verbose = FALSE)
gl <- vcfR2genlight(vcf)		# convert to genlight object


# quick visualization
pdf(cat(output, "_initial.pdf"), paper = "A4")
plot(gl)
dev.off()


# assign correct sample names and population labels to the genlight object if needed
if (ncol(sample_table) != 2) {						# need to rename
	indNames(gl) <- sample_table$V2[match(indNames(gl), sample_table$V1)]
	pop(gl) <- sample_table$V3[match(indNames(gl), sample_table$V2)]
} else {
	pop(gl) <- sample_table$V2[match(indNames(gl), sample_table$V1)]
}


# determine names of replicate pairs
reps <- grep("_R$", indNames(gl), value = TRUE)
samps <- indNames(gl)[match(sub("_R$", "", reps), indNames(gl))]
#pairs <- cbind(reps, samps)
#pairs <- pairs[rowSums(is.na(pairs)) < 1, ]			# remove rows with missing pairs
#npairs <- nrow(pairs)


####
######### May want to implement the testing from Kate Farquaharson
####



# check for monomorphic loci (no variation) and make a list
gl_matrix <- as.matrix(gl)
loc_list <- array(NA, nLoc(gl))

for (index in 1: nLoc(gl)) {
	row <- gl_matrix[, index]

	if (all(row == 0, na.rm = TRUE) | all(row == 2, na.rm = TRUE) | all(row == 1, na.rm = TRUE) | all(is.na(row))) {
		loc_list[index] <- locNames(gl)[index]
	}
}

loc_list <- loc_list[!is.na(loc_list)]
if (length(loc_list) > 0) {
	cat("There are", length(loc_list), "monomorphic or all NA loci")
	print(loc_list)
}




######
########### NEED to find a way to add metadata that may be missing for filters
######



# remove the replicates from the dataset, filter monomorphic loci, and recalculate stats
### NOTE: this will remove more monomorphic loci than in the loc_list because dropping samples may increase 
### the number of monomorphic loci
gl <- gl.drop.ind(gl, ind.list = reps, mono.rm = TRUE, recalc = TRUE, verbose = 3)


# Check for dartR compliance
gl$other$loc.metrics <- as.data.frame(gl$other$loc.metrics)
gl <- gl.compliance.check(gl)


# Filter dataset
# print histograms for evaluating filter choices
pdf(cat(output, "_histograms_prefilter.pdf"), paper = "A4")

hist.callrate = hist(gl@other$loc.metrics$CallRate)	# proportion of samples called
hist.het = hist(gl@other$loc.metrics$FreqHets)		# frequence of heterozygotes
hist.homref = hist(gl@other$loc.metrics$FreqHomRef)	# homozygous for ref allele
hist.homsnp = hist(gl@other$loc.metrics$FreqHomSnp)	# homozygous for alt allele
hist.maf = hist(gl@other$loc.metrics$maf)		# minor allele frequency


######
######### Change plot prettiness?
######


plot(hist.callrate)
plot(hist.het)
plot(hist.homref)
plot(hist.homsnp)
plot(hist.maf)

dev.off()

















