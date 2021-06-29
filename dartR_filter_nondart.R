
##########
# Authors: Rachel Binks and Ben Anderson
# Date: June 2021
# Description: filter a VCF fiile of SNPs from Stacks or ipyrad using the dartR package
# NOTE: this only works for loci with two alleles; more diverse loci should be filtered in a different script
##########


# set parameters needed for output

remove_outgroupZ <- TRUE	# set to remove outgroup samples which have a "Z" in their names
remove_replicates <- TRUE	# set to remove replicate samples which have an "_R" at the end of their names
missing_loc_thresh <- 0.5	# set the call rate by locus threshold
missing_ind_thresh <- 0.2	# set the call rate by ind threshold
minor_allele_thresh <- 0.02	# set the minor allele threshold
lower_rdepth <- 5		# set lower read depth threshold
upper_rdepth <- 200		# set upper read depth threshold


# a helper function for when the script is called without arguments
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to filter a VCF file of RADseq data\n")
		cat("Usage: Rscript dartR_filter_nondart.R -o output -s sample_file -v vcf_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-s\tThe sample names and pops file if vcf labels differ, with the format:\n",
			"\t\tvcf_label\t*tab*\tsample_name\t*tab*\tpopulation\n",
			"\t\tNote: the name should have an \"_R\" in it after the same sequence in the other\n",
			"\t\tNote2: if sample names match vcf_labels, you only need one column then pop\n")
		cat("\t-v\tThe VCF file to be analysed\n",
			"\t\tNote: the VCF file will be filtered only for loci with two alleles\n") 
	} else {
		cat(help_message)
	}
}


# define a function to check for monomorphic loci (no variation) and return a list
# NOTE: this function differs in looking for all '1' as well, which isn't part of gl.filter.monomorphs
filt_mono <- function(my_genlight) {
	gl_matrix <- as.matrix(my_genlight)
	loc_list <- array(NA, nLoc(my_genlight))
	for (index in 1: nLoc(my_genlight)) {
		row <- gl_matrix[, index]
		if (all(row == 0, na.rm = TRUE) | all(row == 2, na.rm = TRUE) 
		| all(row == 1, na.rm = TRUE) | all(is.na(row))) {
			loc_list[index] <- locNames(my_genlight)[index]
		}
	}
	loc_list <- loc_list[!is.na(loc_list)]
	if (length(loc_list) > 0) {
		cat("There are", length(loc_list), "monomorphic or all NA loci\n")
		return(loc_list)
	} else {
		cat("There are no monomorphic loci\n")
		return(0)
	}
}


###############################
###############################


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

# NOTE: dartR importing of VCF via genind or genlight does NOT preserve reference SNP

#geni <- vcfR2genind(vcf)	# convert to genind
#genl <- gi2gl(geni)		# convert to genlight object
# OR
#genl <- gl.read.vcf(vcf)	# convert to genlight object after reading VCF
# OR
genl <- vcfR2genlight(vcf)	# convert to genlight
ploidy(genl) <- 2		# needed for dartR to recognize as SNP data


# assign correct sample names and population labels to the genlight object if needed
if (ncol(sample_table) != 2) {
	indNames(genl) <- sample_table$V2[match(indNames(genl), sample_table$V1)]
	pop(genl) <- sample_table$V3[match(indNames(genl), sample_table$V2)]
} else {
	pop(genl) <- sample_table$V2[match(indNames(genl), sample_table$V1)]
}


# run a compliance check to get proper headings and stats
genl <- gl.compliance.check(genl)


# quick visualization
pdf(paste(output, "_initial.pdf", sep = ""), paper = "A4")
gl.plot(genl)
dev.off()


## read depth
# first, let's calculate an average read depth for each locus
# NOTE: in ipyrad output (VCF 4.0), there is no "AD"=allele depth, so would have to parse out the data used by Kym 

depths <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
depths[depths==0] <- NA		# turn all zero depth values to NA for ignoring
avgdepths <- data.frame(round(rowMeans(depths, na.rm = TRUE), 1))
avgdepths <- subset(avgdepths, !(is.na(avgdepths)))

# ignore rows that are not in the genlight object
avgdepths <- subset(avgdepths, rownames(avgdepths) %in% locNames(genl))

# now create then fill the rdepth column with values
genl@other$loc.metrics$rdepth <- NA
genl@other$loc.metrics$rdepth <- avgdepths[match(locNames(genl), rownames(avgdepths)), ]

# now remove any loci that have an NA in read depth
rm_list <- locNames(genl)[is.na(genl@other$loc.metrics$rdepth)]
genl <- gl.drop.loc(genl, loc.list = rm_list)



# determine names of outgroup samples and remove them if requested,
if (remove_outgroupZ) {
	outgroups <- indNames(genl)[grep("Z", indNames(genl))]
	genl <- gl.drop.ind(genl, ind.list = outgroups, verbose = 3)
}


# determine names of replicate pairs
reps <- grep("_R", indNames(genl), value = TRUE)
fewer_names <- grep("_R", indNames(genl), value = TRUE, invert = TRUE)
samps <- fewer_names[pmatch(sub("_R.*", "", reps), fewer_names)]
#pairs <- cbind(reps, samps)
#pairs <- pairs[rowSums(is.na(pairs)) < 1, ]			# remove rows with missing pairs
#npairs <- nrow(pairs)


####
######### May want to implement the testing from Kate Farquaharson
####


# remove the replicates from the dataset, if requested
if (remove_replicates) {
	genl <- gl.drop.ind(genl, ind.list = reps, verbose = 3)
}


# drop monomorphic loci
loc_list <- filt_mono(genl)
genl <- gl.drop.loc(genl, loc.list = loc_list)


# run gl.filter.monomorphs to set flag, then re-calculate stats
genl <- gl.filter.monomorphs(genl, verbose = 3)
genl <- gl.recalc.metrics(genl, verbose = 3)


# Check for dartR compliance
#genl$other$loc.metrics <- as.data.frame(genl$other$loc.metrics)
genl <- gl.compliance.check(genl)
#names(genl@other$loc.metrics)


# Visualize prior to filtering
# print histograms for evaluating filter choices
pdf(paste(output, "_histograms_prefilter.pdf", sep = ""), paper = "A4")

hist.callrate = hist(genl@other$loc.metrics$CallRate)	# proportion of samples called
hist.het = hist(genl@other$loc.metrics$FreqHets)	# frequency of heterozygotes
hist.homref = hist(genl@other$loc.metrics$FreqHomRef)	# homozygous for ref allele
hist.homsnp = hist(genl@other$loc.metrics$FreqHomSnp)	# homozygous for alt allele
hist.maf = hist(genl@other$loc.metrics$maf)		# minor allele frequency

gl.report.rdepth(genl)

######
######### Change plot prettiness?
######

dev.off()


#### Apply filters

filter <- genl

# include a filter on reproducibility ? (see above)
#filter <- gl.filter.RepAvg(filter, threshold = 0.95, verbose = 3)

filter <- gl.filter.callrate(filter, method = "loc", threshold = missing_loc_thresh, verbose = 3, mono.rm=TRUE)

filter <- gl.filter.maf(filter, threshold = minor_allele_thresh, verbose = 3)

filter <- gl.filter.rdepth(filter, lower = lower_rdepth, upper = upper_rdepth, verbose = 3)
# This was not working because there were NA values in the read depth column (not sure why)

filter <- gl.filter.callrate(filter, method = "ind", threshold = missing_ind_thresh, verbose = 3)


# Visualize post filtering

data <- filter

pdf(paste(output, "_histograms_postfilter.pdf", sep = ""), paper = "A4")

hist.callrate = hist(data@other$loc.metrics$CallRate)	# proportion of samples called
hist.het = hist(data@other$loc.metrics$FreqHets)	# frequency of heterozygotes
hist.homref = hist(data@other$loc.metrics$FreqHomRef)	# homozygous for ref allele
hist.homsnp = hist(data@other$loc.metrics$FreqHomSnp)	# homozygous for alt allele
hist.maf = hist(data@other$loc.metrics$maf)		# minor allele frequency

gl.report.rdepth(data)

######
######### Change plot prettiness?
######

dev.off()


pdf(paste(output, "_pcoa_nj.pdf", sep = ""), paper = "A4")


# create a pcoa
pcoa <- gl.pcoa(data, nfactors = 3)
gl.pcoa.scree(pcoa)
gl.pcoa.plot(pcoa, data, labels = "pop", xaxis = 1, yaxis = 2)


# create NJ tree
njtree <- nj(dist(as.matrix(data)))
plot(njtree, typ = "unrooted", main = "nj tree", cex.main = 0.6, show.tip.label = FALSE, edge.width = 0)
tiplabels(data$pop, frame = "none", cex = 0.5)

dev.off()


# write pcoa scores
write.table(pcoa$scores, paste(output, "pcoascores.txt", quote = FALSE, sep = "\t")





