
##########
# Authors: Rachel Binks and Ben Anderson
# Date: June - July 2021
# Description: filter a VCF fiile of SNPs from Stacks or ipyrad using the dartR package
# NOTE: this only works for loci with two alleles; also, ideally, the data will already be one SNP per locus
# but this script has a built in filter to remove SNPs from the same chromosome slot
##########


# set parameters needed for output (may adjust after a run)

remove_outgroupZ <- TRUE	# set to remove outgroup samples which have a "Z" in their names
remove_replicates <- TRUE	# set to remove replicate samples which have an "_R" at the end of the common text
rep_thresh <- 0.8		# set the reproducibility threshold (NOTE: sensitive to number of reps included)
missing_loc_thresh <- 0.5	# set the call rate by locus threshold (min proportion of samples with a locus)
missing_ind_thresh <- 0.2	# set the call rate by individual threshold (min proportion of loci per individual)
minor_allele_thresh <- 0.02	# set the minor allele threshold (min frequency of allele)
lower_rdepth <- 5		# set lower read depth threshold (min read depth for a SNP)
upper_rdepth <- 200		# set upper read depth threshold (max read depth for a SNP)
output_structure <- FALSE	# set whether to output a Structure format file
output_faststructure <- FALSE	# set whether to output a Faststructure format file
output_snapp <- FALSE		# set whether to output a binary nexus for SNAPP
output_svdquartets <- FALSE	# set whether to output a nucleotide nexus for SVD-Quartets
output_treemix <- FALSE		# set whether to output a TREEMIX file



#######
# Define functions
#######

# a helper function for when the script is called incorrectly or without arguments
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


# a function to check for monomorphic loci (no variation) and return a list
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
suppressMessages(library(vcfR))
suppressMessages(library(dartR))
suppressMessages(library(ape))


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)

vcf <- read.vcfR(vcf_file, verbose = FALSE)


# NOTE: dartR importing of VCF via genind does NOT preserve reference SNP
#geni <- vcfR2genind(vcf)	# convert to genind
#genl <- gi2gl(geni)		# convert to genlight object
# OR
#genl <- gl.read.vcf(vcf)	# convert to genlight object after reading VCF
# OR
genl <- vcfR2genlight(vcf)	# convert to genlight
ploidy(genl) <- 2		# needed for dartR to recognize as SNP data
genl


# assign correct sample names and population labels to the genlight object if needed
# format1 = id	sample_name	pop
# format2 = id	pop
# (format2 is for when the id = desired sample name)
if (ncol(sample_table) != 2) {
	indNames(genl) <- sample_table$V2[match(indNames(genl), sample_table$V1)]
	pop(genl) <- sample_table$V3[match(indNames(genl), sample_table$V2)]
} else {
	pop(genl) <- sample_table$V2[match(indNames(genl), sample_table$V1)]
}


# remove SNPs if there are multiple per "chromosome" = locus
# from https://stackoverflow.com/questions/8041720/randomly-select-on-data-frame-for-unique-rows
if (length(duplicated(genl@chromosome)[duplicated(genl@chromosome)]) > 0) {
	df <- genl@chromosome
	loci <- sample(1: length(df))
	df <- df[loci]
	dups <- duplicated(df)
	df <- df[!dups]
	loci <- loci[!dups]
	loci <- loci[sort(loci, index.return = TRUE)$ix]
	loci_include <- genl@loc.names[loci]
	loci_exclude <- genl@loc.names[-loci]
	cat("Removing", length(loci_exclude), "SNPs from the same loci (random)\n")
	# NOTE: slicing it manually leads to problems with @other slot dimensions, even though it is faster
	genl <- genl[, match(loci_include, genl@loc.names)]
	# NOTE: as long as the compliance check comes after this step, should be no problem
	#genl <- gl.drop.loc(genl, loc.list = loci_exclude, verbose = 0)
	genl
}


# run a compliance check to get proper headings and stats
genl <- gl.compliance.check(genl, verbose = 0)


# read in the reference and alt alleles from the vcf file and assign to the @loc.all slot of the genlight
alleles <- data.frame(vcf@fix)[, c("ID", "REF", "ALT")]
alleles <- alleles[grep(",", alleles$ALT, invert = TRUE), ]		# exclude multi-allelic loci
alleles$all <- with(alleles, paste0(REF, "/", ALT))

# now fill the loc.all slot with values
genl@loc.all <- alleles$all[match(locNames(genl), alleles$ID)]


# plot an initial visualization of the data
pdf(paste0(output, "_initial.pdf"), paper = "A4")
cat("Plotting an initial visualization of the data\n")
gl.plot(genl, verbose = 0)
dev.off()


# calculate average read depth per locus from metadata and attach it to the genlight object
# NOTE: in ipyrad output (VCF 4.0), there is no "AD"=allele depth, so would have to parse out the data used by Kym 
# I don't think it matters, since dartR is only using an average depth for all alleles anyway
depths <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
depths[depths == 0] <- NA		# turn all zero depth values to NA for ignoring
avgdepths <- data.frame(round(rowMeans(depths, na.rm = TRUE), 1))
avgdepths <- subset(avgdepths, !(is.na(avgdepths)))

# ignore rows that are not in the genlight object
avgdepths <- subset(avgdepths, rownames(avgdepths) %in% locNames(genl))

# now fill the rdepth column with values
genl@other$loc.metrics$rdepth <- avgdepths[match(locNames(genl), rownames(avgdepths)), ]

# now remove any loci that have an NA in read depth
# NOTE: I'm not sure why this happens, but it may be because the genlight has loci that had zero read depth(?)
rm_list <- locNames(genl)[is.na(genl@other$loc.metrics$rdepth)]
if (length(rm_list) > 0) {
	cat("Dropping", length(rm_list), "loci without read depth information\n")
	genl <- gl.drop.loc(genl, loc.list = rm_list, verbose = 0)
}



####### Outgroups
# determine names of outgroup samples and remove them if requested,
if (remove_outgroupZ) {
	outgroups <- indNames(genl)[grep("Z", indNames(genl))]
	cat("Dropping", length(outgroups), "outgroup samples\n")
	genl <- gl.drop.ind(genl, ind.list = outgroups, verbose = 0)
}



####### Replicates
# determine names of replicate pairs
reps <- grep("_R", indNames(genl), value = TRUE)
fewer_names <- grep("_R", indNames(genl), value = TRUE, invert = TRUE)
samps <- fewer_names[pmatch(sub("_R.*", "", reps), fewer_names)]


# run reproducibility, originally based on ideas from script/paper by Kate Farquaharson
cat("Running reproducibility assessment\n")
pairs <- cbind(reps, samps)
pairs <- pairs[rowSums(is.na(pairs)) < 1, ]			# remove rows with missing pairs
npairs <- nrow(pairs)
sub_genl <- genl[match(c(pairs), indNames(genl))]		# create a genlight object with only rep pairs

# for each locus/SNP, for each pair determine whethere they differ
# then enter an error for that locus for that pair
# then compare all the errors vs non-errors for pairs in that locus
# add a value to a new dataframe for that locus as reproducibility
# NOTE: I have slightly altered the original (I think) to count when one read is NA and the other isn't as an error
# NOTE: this approach will create a ratio that is HIGHLY dependent on how many reps are present
# so the filter applied should be carefully adjusted; I will add a report on what one error = in reproducibility
df_repro <- data.frame(locus = locNames(sub_genl), RepAvg = double(length(locNames(sub_genl))),
			stringsAsFactors = FALSE)
calls <- as.matrix(sub_genl)
for (locus_ind in 1: length(locNames(sub_genl))) {	# for each locus
	lcalls <- calls[, locus_ind]
	errors <- vector("list", npairs)
	for (row in 1: npairs) {		# for each replicate pair
		callA <- lcalls[pairs[[row, 1]]]	# call for indA
		callB <- lcalls[pairs[[row, 2]]]	# call for indB
		if (is.na(callA)) {		# if there is no read for indA, assign a number for comparison
			callA <- 5
		}
		if (is.na(callB)) {		# if there is no read for indB, assign a number for comparison
			callB <- 5
		}
		if (callA != callB) {
			errors[[row]] <- 1
		} else {
			errors[[row]] <- 0
		}
	}
	# now, sum the errors, then divide by number of reps for reproducibility metric for that SNP
	# NOTE: this assumes that two missing reads are a correct call (wrong); same as in original
	repro <- 1 - (sum(unlist(errors)) / npairs) 
	# add to the dataframe
	#df_repro[locus_ind, "locus"] <- locNames(sub_genl)[locus_ind]
	df_repro[locus_ind, "RepAvg"] <- as.double(repro)	
}

# add the values to the genlight object
genl@other$loc.metrics$RepAvg <- df_repro$RepAvg[match(locNames(genl), df_repro$locus)]

# report results of reproducibility assessment
cat("Reproducibility summary:\n")
summary(df_repro$RepAvg)
repro_lower <- length(df_repro$RepAvg[df_repro$RepAvg < 1])
repro_one <- 1 - (1 / npairs)
cat("There are", repro_lower, "of", nrow(df_repro), "SNPs with less than 100% reproducibility\n")
cat("A single error reduces reproducibility to", repro_one, "\n")


# remove the replicates from the dataset, if requested
if (remove_replicates) {
	cat("Dropping", length(reps), "replicate samples\n")
	genl <- gl.drop.ind(genl, ind.list = reps, verbose = 0)
}


# drop monomorphic loci
loc_list <- filt_mono(genl)
if (class(loc_list) == "array") {
	cat("Dropping", length(loc_list), "monomorphic loci\n")
	genl <- gl.drop.loc(genl, loc.list = loc_list, verbose = 0)
}


# run gl.filter.monomorphs to set flag, then re-calculate stats
genl <- gl.filter.monomorphs(genl, verbose = 0)
cat("Recalculating metrics\n")
genl <- gl.recalc.metrics(genl, verbose = 0)


# check for dartR compliance
#genl$other$loc.metrics <- as.data.frame(genl$other$loc.metrics)
cat("Checking dartR compliance\n")
genl <- gl.compliance.check(genl, verbose = 0)
#names(genl@other$loc.metrics)
genl



#######
# Filter the dataset based on loci metrics and run a few quick analyses
#######


# Visualize prior to further filtering on loci metrics
pdf(paste0(output, "_histograms_prefilter.pdf"), paper = "A4")
hist.callrate = hist(genl@other$loc.metrics$CallRate)	# proportion of samples called
hist.het = hist(genl@other$loc.metrics$FreqHets)	# frequency of heterozygotes
hist.homref = hist(genl@other$loc.metrics$FreqHomRef)	# homozygous for ref allele
hist.homsnp = hist(genl@other$loc.metrics$FreqHomSnp)	# homozygous for alt allele
hist.maf = hist(genl@other$loc.metrics$maf)		# minor allele frequency
outliers_pre <- gl.report.rdepth(genl, verbose = 0)
dev.off()


# Apply filters on loci metrics
filter <- genl
filter <- gl.filter.reproducibility(filter, threshold = rep_thresh, verbose = 3)
filter <- gl.filter.callrate(filter, method = "loc", threshold = missing_loc_thresh,
			plot = FALSE, verbose = 3)
filter <- gl.filter.callrate(filter, method = "ind", threshold = missing_ind_thresh, 
			plot = FALSE, verbose = 3)
filter <- gl.filter.maf(filter, threshold = minor_allele_thresh, verbose = 3)
# NOTE: rdepth was not working because there were NA values in the read depth column (not sure why)
filter <- gl.filter.rdepth(filter, lower = lower_rdepth, upper = upper_rdepth, verbose = 3)


# remove monomorphic loci again, if present
loc_list <- filt_mono(filter)
if (class(loc_list) == "array") {
	cat("Dropping", length(loc_list), "monomorphic loci\n")
	filter <- gl.drop.loc(filter, loc.list = loc_list, verbose = 0)
}


# Visualize post filtering
data <- filter
pdf(paste0(output, "_histograms_postfilter.pdf"), paper = "A4")
hist.callrate = hist(data@other$loc.metrics$CallRate)	# proportion of samples called
hist.het = hist(data@other$loc.metrics$FreqHets)	# frequency of heterozygotes
hist.homref = hist(data@other$loc.metrics$FreqHomRef)	# homozygous for ref allele
hist.homsnp = hist(data@other$loc.metrics$FreqHomSnp)	# homozygous for alt allele
hist.maf = hist(data@other$loc.metrics$maf)		# minor allele frequency
outliers_post <- gl.report.rdepth(data, verbose = 0)
dev.off()


# Output some basic diagrams
pdf(paste0(output, "_pcoa_nj.pdf"), paper = "A4")

# create a pcoa
pcoa <- gl.pcoa(data, nfactors = 3)
gl.pcoa.scree(pcoa)
gl.pcoa.plot(pcoa, data, labels = "pop", xaxis = 1, yaxis = 2)

# create a NJ tree
njtree <- njs(dist(as.matrix(data)))		# njs allows missing data
# NOTE: may want to change the distance matrix to something with DNA substitution model
plot(njtree, typ = "unrooted", main = "nj tree", cex.main = 0.6, show.tip.label = FALSE, edge.width = 0)
tiplabels(data$pop, frame = "none", cex = 0.5)

dev.off()


# write pcoa scores
write.table(pcoa$scores, paste0(output, "_pcoascores.txt"), quote = FALSE, sep = "\t")



#######
# Output files for further downstream analyses
#######


# output a list of loci names so that they can be extracted from fasta output for phylogenetic analysis
# NOTE: this needs to be tested to see if fasta output uses the same loci names
write.table(locNames(data), paste0(output, "_loci_names.txt"), quote = FALSE)


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
