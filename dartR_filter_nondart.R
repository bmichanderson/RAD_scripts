
##########
# Authors: Rachel Binks and Ben Anderson
# Date: June - Sep 2021
# Description: filter a VCF file with the dartR package
# NOTE: this only works for biallelic loci
# NOTE: will downsample (randomly) to one SNP per locus
##########


## set output parameters (may adjust after a run)
remove_outgroups <- TRUE	# remove outgroup samples which have a "Z" in their names
remove_replicates <- FALSE	# remove replicate samples which have an "_R" at the end of the common text
repro_thresh <- 0.8			# set reproducibility threshold (NOTE: sensitive to number of reps included)
miss_loc_thresh <- 0.8		# set call rate by locus threshold (min proportion of samples with a locus)
miss_ind_thresh <- 0.2		# set call rate by ind threshold (min proportion of loci per individual)
minor_allele_thresh <- 0.01	# set minor allele threshold (min frequency of allele)
lower_rdepth <- 10			# set lower/upper read depth thresholds (min/max read depth for a SNP)
upper_rdepth <- 500
output_structure <- FALSE	# set what to output
output_faststructure <- FALSE
output_snapp <- FALSE
output_svdquartets <- FALSE
output_treemix <- FALSE


## load required libraries
suppressMessages(library(ape))
suppressMessages(library(adegenet))
suppressMessages(library(vcfR))
suppressMessages(library(dartR))


#######
# Define functions
#######

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to filter a VCF file of RADseq data\n")
		cat("Usage: Rscript dartR_filter_nondart.R -o output -s sample_file -v vcf_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-s\tThe tab-separated sample names and pops file if vcf labels differ, with the format:\n",
			"\t\tvcf_label\tsample_name\tpopulation\n",
			"\t\tNote: replicate names should have an \"_R\" in them after the same text as other rep\n",
			"\t\tNote: if sample names already match vcf_labels, you only need one column then pop\n")
		cat("\t-v\tThe VCF file to be analysed\n",
			"\t\tNote: only loci with two alleles will be kept for filtering\n")
	} else {
		cat(help_message)
	}
}


# a function to check for monomorphic loci (no variation) in a genlight object and return a list
# NOTE: this function differs in looking for all '1' as well, which isn't part of gl.filter.monomorphs
filt_mono <- function(my_genl) {
	gl_matrix <- as.matrix(my_genl)
	loc_list <- array(NA, length(my_genl@loc.names))
	for (index in seq_len(length(my_genl@loc.names))) {
		row <- gl_matrix[, index]
		if (any(all(row == 0, na.rm = TRUE), all(row == 2, na.rm = TRUE),
		all(row == 1, na.rm = TRUE), all(is.na(row)))) {
			loc_list[index] <- my_genl@loc.names[index]
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



#######
# Read in and format the data
#######


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

	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			output <- args[index + 1]
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

if (! samples_present) {
	stop(help("Missing argument for sample/pops file!\n"), call. = FALSE)
}

if (! vcf_present) {
	stop(help("Missing argument for vcf file!\n"), call. = FALSE)
}


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# keep one SNP per locus, preferring the higher sample coverage or, if a tie, the first one encountered
# the vcf@fix slot has columns with info per locus; from ipyrad, column 8 has NS=xx;DP=xx
# NS is the number of samples for that SNP
ns <- lapply(vcf@fix[, 8], function(x) {
	strsplit(strsplit(x, split = ";")[[1]][1], split = "=")[[1]][2]
	}
	)
ns <- unlist(ns)
loc <- vcf@fix[, 1]	# locus ids (chrom) are in column 1
id <- vcf@fix[, 3]	# SNP ids (locus_pos) are in column 3
# choose the best SNP per unique locus
# (based on idea from https://github.com/ksil91/Ostrea_PopStructure/blob/master/Scripts/subsetSNPs.py)
id_choose <- rep(NA, length(unique(loc)))
index <- 1
this_locus <- loc[1]
best_id <- id[1]
best_ns <- as.numeric(ns[1])
for (snp in seq_len(length(id))) {
	locus <- loc[snp]
	if (locus != this_locus) {
		id_choose[index] <- best_id
		index <- index + 1
		this_locus <- locus
		best_ns <- 0
		best_id <- id[snp]
	}
	if (as.numeric(ns[snp]) > best_ns) {
		best_ns <- as.numeric(ns[snp])
		best_id <- id[snp]
	}
}
id_choose[index] <- best_id
# filter the vcf for only those SNPs
vcf <- vcf[vcf@fix[, 3] %in% id_choose]
cat("Retained", index, "SNPs (one per locus) from the VCF file\n")


# NOTE: dartR importing of VCF via genind does NOT preserve reference SNP
# Don't use: geni <- vcfR2genind(vcf), followed by: genl <- gi2gl(geni)
# Optionally could use from the dartR package: genl <- gl.read.vcf(vcf)
# but this doesn't allow the above filtering for one SNP per locus

# convert to genlight (will lose loci that are not biallelic)
genl <- vcfR2genlight(vcf)
# assign ploidy for dartR to recognize as SNP data
ploidy(genl) <- 2
# print the object to screen
#genl


# assign correct sample names and population labels to the genlight object if needed
# ALSO: filter the genlight to only retain samples present in the sample table
# format1 = id	sample_name	pop
# format2 = id	pop
# (format2 is for when the id = desired sample name)
len_original <- length(genl@ind.names)
original_names <- genl@ind.names
if (ncol(sample_table) != 2) {
	genl@ind.names <- sample_table$V2[match(genl@ind.names, sample_table$V1)]
	populations <- sample_table$V3[match(genl@ind.names, sample_table$V2)]
} else {
	genl@ind.names <- sample_table$V1[match(genl@ind.names, sample_table$V1)]
	populations <- sample_table$V2[match(genl@ind.names, sample_table$V1)]
}
if (length(genl@ind.names[is.na(genl@ind.names)]) > 0) {		# if there are any NAs (non-matches)
	original_names <- original_names[!is.na(genl@ind.names)]	# to make sure vcf data matches
	populations <- populations[!is.na(genl@ind.names)]
	genl <- genl[!is.na(genl@ind.names), ]
}
len_new <- length(genl@ind.names)
cat("Retained", length(genl@loc.names), "SNPs from the VCF file\n")
cat("Removed", len_original - len_new, "samples not present in the input sample table\n")
cat("Retained", length(genl@ind.names), "samples\n")

# assign population labels for dartR
pop(genl) <- populations


# run a compliance check to get proper headings and stats
genl <- gl.compliance.check(genl, verbose = 0)


# read in the reference and alt alleles from the vcf file and assign to the @loc.all slot of the genlight
alleles <- data.frame(vcf@fix)[, c("ID", "REF", "ALT")]
alleles <- alleles[grep(",", alleles$ALT, invert = TRUE), ]		# exclude multi-allelic loci
alleles$all <- with(alleles, paste0(REF, "/", ALT))

# now fill the loc.all slot with values
genl@loc.all <- alleles$all[match(original_names, alleles$ID)]


# plot an initial visualization of the data
pdf(paste0(output, "_initial.pdf"), paper = "A4")
cat("Plotting an initial visualization of the data\n")
gl.plot(genl, verbose = 0)
invisible(dev.off())



####### Outgroups
# determine names of outgroup samples and remove them if requested
if (remove_outgroups) {
	outgroups <- genl@ind.names[grep("Z", genl@ind.names)]
	if (length(outgroups) > 0) {
		cat("Dropping", length(outgroups), "outgroup samples\n")
		original_names <- original_names[grep("Z", genl@ind.names, invert = TRUE)]
		genl <- gl.drop.ind(genl, ind.list = outgroups, verbose = 0)
		cat("Retained", length(genl@ind.names), "samples\n")
	}
}



####### Replicates
# determine names of replicate pairs
reps <- grep("_R", genl@ind.names, value = TRUE)
if (length(reps) > 0) {
	fewer_names <- grep("_R", genl@ind.names, value = TRUE, invert = TRUE)
	samps <- fewer_names[pmatch(sub("_R.*", "", reps), fewer_names)]

	# run reproducibility, originally based on ideas from script/paper by Kate Farquaharson
	cat("Running reproducibility assessment\n")
	pairs <- cbind(reps, samps)
	pairs <- pairs[rowSums(is.na(pairs)) < 1, ]			# remove rows with missing pairs
	npairs <- nrow(pairs)
	sub_genl <- genl[match(c(pairs), genl@ind.names)]		# create a genlight object with only rep pairs

	# for each locus/SNP, for each pair determine whether they differ
	# then enter an error value for that locus for that pair
	# then compare all the errors vs non-errors for pairs in that locus
	# add a value to a new dataframe for that locus as reproducibility
	# NOTE: I have slightly altered the original (I think) to count when one read is NA and the other isn't as an error
	# NOTE: this approach will create a ratio that is HIGHLY dependent on how many reps are present
	# so the filter applied should be carefully adjusted; I will add a report on what one error = in reproducibility
	df_repro <- data.frame(locus = sub_genl@loc.names, RepAvg = double(length(sub_genl@loc.names)),
				stringsAsFactors = FALSE)
	calls <- as.matrix(sub_genl)
	for (locus_ind in seq_len(length(sub_genl@loc.names))) {	# for each locus
		lcalls <- calls[, locus_ind]
		errors <- vector("list", npairs)
		for (row in 1: npairs) {		# for each replicate pair
			call_a <- lcalls[pairs[[row, 1]]]	# call for indA
			call_b <- lcalls[pairs[[row, 2]]]	# call for indB
			if (is.na(call_a)) {		# if there is no read for indA, assign a number for comparison
				call_a <- 5
			}
			if (is.na(call_b)) {		# if there is no read for indB, assign a number for comparison
				call_b <- 5
			}
			if (call_a != call_b) {
				errors[[row]] <- 1
			} else {
				errors[[row]] <- 0
			}
		}
		# now, sum the errors, then divide by number of reps for reproducibility metric for that SNP
		# NOTE: this assumes that two missing reads are a correct call (wrong); same as in original
		repro <- 1 - (sum(unlist(errors)) / npairs)
		# add to the dataframe
		df_repro[locus_ind, "RepAvg"] <- as.double(repro)
	}

	# add the values to the genlight object
	genl@other$loc.metrics$RepAvg <- df_repro$RepAvg[match(genl@loc.names, df_repro$locus)]

	# report results of reproducibility assessment
	cat("Reproducibility summary:\n")
	summary(df_repro$RepAvg)
	repro_lower <- length(df_repro$RepAvg[df_repro$RepAvg < 1])
	repro_one <- 1 - (1 / npairs)
	cat("There are", repro_lower, "of", nrow(df_repro), "SNPs with less than 100% reproducibility\n")
	cat("A single error reduces reproducibility to", repro_one, "\n")
} else {
	cat("No replicates for reproducibility assessment\n")
}


# remove the replicates from the dataset, if requested
if (length(reps) > 0 && remove_replicates) {
	cat("Dropping", length(reps), "replicate samples\n")
	original_names <- original_names[grep("_R", genl@ind.names, invert = TRUE)]
	genl <- gl.drop.ind(genl, ind.list = reps, verbose = 0)
	cat("Retained", length(genl@ind.names), "samples\n")
}


# calculate average read depth per locus from metadata and attach it to the genlight object
# NOTE: in ipyrad output (VCF 4.0), there is no "AD"=allele depth, so would have to parse out the data used by Kym
# I don't think it matters, since dartR is only using an average depth for all alleles anyway
depths <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
depths <- depths[, !is.na(match(original_names, colnames(depths)))]	# only keep for samples in the genl
depths[depths == 0] <- NA		# turn all zero depth values to NA for ignoring
avgdepths <- data.frame(round(rowMeans(depths, na.rm = TRUE), 1))
avgdepths <- subset(avgdepths, !(is.na(avgdepths)))

# ignore rows that are not in the genlight object
avgdepths <- subset(avgdepths, rownames(avgdepths) %in% genl@loc.names)

# now fill the rdepth column with values
genl@other$loc.metrics$rdepth <- avgdepths[match(genl@loc.names, rownames(avgdepths)), ]

# now remove any loci that have an NA in read depth
# NOTE: I'm not sure why this sometimes happens, possibly because the genlight has loci that had zero read depth?
rm_list <- genl@loc.names[is.na(genl@other$loc.metrics$rdepth)]
if (length(rm_list) > 0) {
	cat("Dropping", length(rm_list), "loci without read depth information\n")
	genl <- gl.drop.loc(genl, loc.list = rm_list, verbose = 0)
	cat("Retained", length(genl@loc.names), "loci\n")
}


# drop monomorphic loci
loc_list <- filt_mono(genl)
if (class(loc_list) == "array") {
	cat("Dropping", length(loc_list), "monomorphic loci\n")
	genl <- gl.drop.loc(genl, loc.list = loc_list, verbose = 0)
	cat("Retained", length(genl@loc.names), "loci\n")
}


# run gl.filter.monomorphs to set flag, then re-calculate stats
genl <- gl.filter.monomorphs(genl, verbose = 0)
cat("Recalculating metrics\n")
genl <- gl.recalc.metrics(genl, verbose = 0)


# check for dartR compliance
cat("Checking dartR compliance\n")
genl <- gl.compliance.check(genl, verbose = 0)
#genl



#######
# Filter the dataset based on loci metrics and run a few quick analyses
#######


# Visualize prior to further filtering on loci metrics
pdf(paste0(output, "_hist_pre.pdf"), paper = "A4")
hist.callrate <- hist(genl@other$loc.metrics$CallRate) # proportion of samples called
hist.het <- hist(genl@other$loc.metrics$FreqHets) # frequency of heterozygotes
hist.homref <- hist(genl@other$loc.metrics$FreqHomRef) # homozygous for ref allele
hist.homsnp <- hist(genl@other$loc.metrics$FreqHomSnp) # homozygous for alt allele
hist.maf <- hist(genl@other$loc.metrics$maf) # minor allele frequency
outliers_pre <- gl.report.rdepth(genl, verbose = 0)
invisible(dev.off())


# Apply filters on loci metrics
filter <- genl
if (length(reps) > 0) {
	filter <- gl.filter.reproducibility(filter, threshold = repro_thresh, verbose = 3)
}
filter <- gl.filter.callrate(filter, method = "loc", threshold = miss_loc_thresh,
			plot = FALSE, verbose = 3)
filter <- gl.filter.callrate(filter, method = "ind", threshold = miss_ind_thresh,
			plot = FALSE, verbose = 3)
filter <- gl.filter.maf(filter, threshold = minor_allele_thresh, verbose = 3)
# NOTE: rdepth was not working because there were NA values in the read depth column (not sure why)
filter <- gl.filter.rdepth(filter, lower = lower_rdepth, upper = upper_rdepth, verbose = 3)


# remove monomorphic loci again, if present
loc_list <- filt_mono(filter)
if (class(loc_list) == "array") {
	cat("Dropping", length(loc_list), "monomorphic loci\n")
	filter <- gl.drop.loc(filter, loc.list = loc_list, verbose = 0)
	cat("Retained", length(genl@loc.names), "loci\n")
}


# Visualize post filtering
data <- filter
pdf(paste0(output, "_hist_post.pdf"), paper = "A4")
hist.callrate <- hist(data@other$loc.metrics$CallRate) # proportion of samples called
hist.het <- hist(data@other$loc.metrics$FreqHets) # frequency of heterozygotes
hist.homref <- hist(data@other$loc.metrics$FreqHomRef) # homozygous for ref allele
hist.homsnp <- hist(data@other$loc.metrics$FreqHomSnp) # homozygous for alt allele
hist.maf <- hist(data@other$loc.metrics$maf) # minor allele frequency
outliers_post <- gl.report.rdepth(data, verbose = 0)
invisible(dev.off())


## Output some basic diagrams
# create a pcoa
pdf(paste0(output, "_pcoa.pdf"), paper = "A4")
pcoa <- gl.pcoa(data, nfactors = 3)
pcoa_coords <- gl.pcoa.plot(pcoa, data, labels = "pop", plot.out = FALSE, xaxis = 1, yaxis = 2)
invisible(dev.off())

# write pcoa scores
# Optionally could use write.table(pcoa$scores, paste0(output, "_pcoascores.txt"), quote = FALSE, sep = "\t")

# create a NJ tree
pdf(paste0(output, "_nj.pdf"), paper = "A4")
njtree <- njs(dist(as.matrix(data)))		# njs allows missing data
# NOTE: may want to change the distance matrix to something with DNA substitution model
plot(njtree, typ = "unrooted", main = "nj tree", cex.main = 0.6, show.tip.label = FALSE, edge.width = 0)
tiplabels(data$pop, frame = "none", cex = 0.5)
invisible(dev.off())



#######
# Output files for further downstream analyses
#######


# output a list of loci names so that they can be extracted from fasta output for phylogenetic analysis
# NOTE: this needs to be tested to see if fasta output uses the same loci names
# could use: write.table(data@loc.names, paste0(output, "_loci.txt"), quote = FALSE)


# export a structure file (NOTE: will not be in pop order)
if (output_structure) {
	cat("Outputting Structure file\n")
	gl2structure(data, outfile = paste0(output, "_structure.str"), addcolumns = data$pop,
		exportMarkerNames = FALSE, outpath = getwd(), verbose = 0)
}


# export a faststructure file (NOTE: doesn't add ind/pop data)
if (output_faststructure) {
	cat("Outputting Faststructure file\n")
	gl2faststructure(data, outfile = paste0(output, "_faststructure.str"),
		outpath = getwd(), verbose = 0)
}


# export a binary nexus file of variant sites for SNAPP
if (output_snapp) {
	cat("Outputting SNAPP file\n")
	gl2snapp(data, outfile = paste0(output, "_snapp.nex"), outpath = getwd(), verbose = 0)
}


# export a nucleotide nexus file for SVD-Quartets (NOTE: this needs @loc.all slot filled)
if (output_svdquartets) {
	cat("Outputting SVD-Quartets file\n")
	gl2svdquartets(data, outfile = paste0(output, "_svd.nex"), method = 2,
		outpath = getwd(), verbose = 0)
}


# export a treemix file (NOTE: output is already gzipped)
if (output_treemix) {
	cat("Outputting TREEMIX file\n")
	gl2treemix(data, outfile = paste0(output, "_treemix.txt.gz"), outpath = getwd(), verbose = 0)
}
