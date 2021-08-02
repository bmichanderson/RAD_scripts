
##########
# Authors: Ben Anderson, based partly on McCartney-Melstad et al. 2019 (https://github.com/atcg/clustOpt)
# Date: July-Aug 2021
# Description: summarize various metrics for optimizing parameters for ipyrad assembly of RAD data
##########


# a helper function for when the script is called without arguments
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to summarize outputs from tests/queries of metrics used",
			"for evaluating clustering thresholds in RAD datasets\n")
		cat("Usage: Rscript MM_ipyrad_summary.R -o output_prefix -v vcf_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-e\tThe seq error stats file\n")
		cat("\t-h\tThe heterozygosity stats file\n")
		cat("\t-l\tThe locus and snps recovered stats file\n")
		cat("\t-p\tThe paralogs stats file\n")
		cat("\t-t\tThe bootstraps stats file\n")
		cat("\t--merror\tThe Mastretta-Yanes error rates file\n")
		cat("\t--mdist\tThe Mastretta-Yanes population distances file\n")
		cat("\t--mpca\tThe Mastretta-Yanes PCA var explained table\n")
		cat("\t-s\tThe samples file with the format: sample <tab> lat <tab> lon\n")
		cat("\t-v\tA list of vcf files to be analysed, one per line with format:\n")
		cat("\t\tclust_num <tab> relative_path_to_vcf\n")
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

	outpre <- "output"

	estats_present <- FALSE
	hstats_present <- FALSE
	lstats_present <- FALSE
	pstats_present <- FALSE
	tstats_present <- FALSE
	merror_present <- FALSE
	mdist_present <- FALSE
	mpca_present <- FALSE
	samples_present <- FALSE
	vcf_present <- FALSE

	for (index in 1:length(args)) {
		if (args[index] == "-o") {
			outpre <- args[index + 1]
		} else if (args[index] == "-e") {
			estats_present <- TRUE
			estats_file <- args[index + 1]
		} else if (args[index] == "-h") {
			hstats_present <- TRUE
			hstats_file <- args[index + 1]
		} else if (args[index] == "-l") {
			lstats_present <- TRUE
			lstats_file <- args[index + 1]
		} else if (args[index] == "-p") {
			pstats_present <- TRUE
			pstats_file <- args[index + 1]
		} else if (args[index] == "-t") {
			tstats_present <- TRUE
			tstats_file <- args[index + 1]
		} else if (args[index] == "--merror") {
			merror_present <- TRUE
			merror_file <- args[index + 1]
		} else if (args[index] == "--mdist") {
			mdist_present <- TRUE
			mdist_file <- args[index + 1]
		} else if (args[index] == "--mpca") {
			mpca_present <- TRUE
			mpca_file <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samp_file <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present = TRUE
			vcf_list_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}


########
# Collate and plot various stats files
########

# read in the files when present and plot
if (estats_present) {
	# plot the variation in sequencing error rate for clustering threshold
	pdf(file = paste0(outpre, "_seqerror.pdf"))
	error_table <- read.table(estats_file, sep = "\t", header = FALSE)
	boxplot(error_table$V2 ~ error_table$V1,
		main = "Estimated sequencing error rates across samples",
		ylab = "Error rate", xlab = "Clustering threshold (%)")
	invisible(dev.off())
}

if (hstats_present) {
	# plot the heterozygosity stats per sample per clustering threshold
	pdf(file = paste0(outpre, "_heterozygosity.pdf"))
	hets_table <- read.table(hstats_file, sep = "\t", header = FALSE)
	boxplot(hets_table$V2 * 100 ~ hets_table$V1, 
		main = "Heterozygosity across samples\n(heterozygous sites compared to total sites)", 
		ylab = "Heterozygous sites (%)", xlab = "Clustering threshold (%)")
	invisible(dev.off())
}

if (lstats_present) {
	# plot the number of total loci and SNPs recovered
	pdf(file = paste0(outpre, "_loci_snps.pdf"))
	loci_snps_table <- read.table(lstats_file, sep = "\t", header = FALSE)
	plot(loci_snps_table$V1, loci_snps_table$V2, main = "Total number of loci recovered", 
		ylab = "Number of loci", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(loci_snps_table$V1), max(loci_snps_table$V1), by = 1))
	plot(loci_snps_table$V1, loci_snps_table$V3, main = "Total number of SNPs recovered", 
		ylab = "Number of SNPs", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(loci_snps_table$V1), max(loci_snps_table$V1), by = 1))
	invisible(dev.off())
}

if (pstats_present) {
	# calculate and plot the proportion of paralogs filtered
	pdf(file = paste0(outpre, "_paralogs.pdf"))
	paralogs_table <- read.table(pstats_file, sep = "\t", header = FALSE)
	paralogs_table$V4 <- paralogs_table$V2 / paralogs_table$V3
	plot(paralogs_table$V1, paralogs_table$V4, 
		main = "Proportion paralogous loci filtered\n(using max_shared_Hs_locus)", 
		ylab = "Proportion paralogous loci", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(paralogs_table$V1), max(paralogs_table$V1), by = 1))
	invisible(dev.off())
}

if (tstats_present) {
	# assess/plot the variation in bootstrap support for ML trees
	pdf(file = paste0(outpre, "_bootstraps.pdf"))
	supp_table <- read.table(tstats_file, sep = "\t", header = FALSE)
	boxplot(supp_table$V2 ~ supp_table$V1, 
		main = "Average ufbootstrap support",
		ylab = "Average ufbootstrap support (10 IQ-Tree reps)", xlab = "Clustering threshold (%)")
	invisible(dev.off())
}

if (merror_present) {
	# plot Mastretta-Yanes errors
	pdf(file = paste0(outpre, "_merror.pdf"))
	myerror_table <- read.table(merror_file, sep = "\t", header = FALSE)
	boxplot(myerror_table$V2 ~ myerror_table$V1,
		main = "Locus error rates\n(loci present in one rep but not the other, relative to total in either/both)",
		ylab = "Locus error rate (%)", xlab = "Clustering threshold (%)")
	boxplot(myerror_table$V3 ~ myerror_table$V1,
		main = "Allele error rates\n(loci present in both that differ, relative to total in both)",
		ylab = "Allele error rate (%)", xlab = "Clustering threshold (%)")
	boxplot(myerror_table$V4 ~ myerror_table$V1,
		main = "SNP error rates\n(SNPs called in both that differ, relative to total called in both)",
		ylab = "SNP error rate (%)", xlab = "Clustering threshold (%)")
	invisible(dev.off)
}

if (mdist_present) {
	# plot Mastretta-Yanes population Euclidean distances
	pdf(file = paste0(outpre, "_mdist.pdf"))
	mydist_table <- read.table(mdist_file, sep = "\t", header = FALSE)
	boxplot(100 * mydist_table$V3 ~ mydist_table$V1,
		main = "Average intrapopulation Euclidean distances",
		ylab = "Average intrapopulation Euclidean distance (normalized %)", xlab = "Clustering threshold (%)")
	invisible(dev.off)
}

if (mpca_present) {
	# plot the percent explained by the first four PCs
	pdf(file = paste0(outpre, "_mpca.pdf"))
	mypca_table <- read.table(mpca_file, sep = "\t", header = FALSE)
	mypca_table[, 6] <- rowSums(mypca_table[, 2:5])
	plot(mypca_table$V1, mypca_table$V6, 
		main = "Variation included in the first four axes of the PCoA",
		ylab = "Variation explained by first four axes of the PCoA (%)", xlab = "Clustering threshold (%)",
		xaxt = "n", pch = 16)
	axis(1, at = seq(min(mypca_table$V1), max(mypca_table$V1), by = 1))
	#!!!!!!!!!!!!!!!!!!!
	# May want to add plotting of PCoA
	# Need to read in scores, perhaps
	#!!!!!!!!!!!!!!!!!!!
	invisible(dev.off)
}




########
# Run calculations and plot stats for VCF files
########

if (vcf_present) {
	if (! samples_present) {
		stop(help("\n*****\nSpecify a samples file with location information\n*****\n"), call. = FALSE)
	}

	# read in the list table and samples table
	vcf_list <- read.table(vcf_list_file)
	latlon_table <- read.table(samp_file, sep = "\t", header = FALSE)

	# load the libraries
	suppressMessages(library(SNPRelate))
	suppressMessages(library(geosphere))

	# set a couple parameters
	axes_to_use <- 8		# the number of PCoA axes to consider for sum of variation explained
	missing_rate <- 0.9		# the maximum amount of missing data per SNP (e.g. 0.9 == present in >10%)


	# get a first estimate of genetic distance from the midpoint of the vcf files
	snpgdsVCF2GDS(vcf_list[round(nrow(vcf_list)/2), 2], "temp.gds", verbose = FALSE, method = "biallelic.only")
	gds <- snpgdsOpen("temp.gds")
	# create an identity-by-state estimation
	MIDPOINTibs_est <- snpgdsIBS(gds, autosome.only = FALSE, num.thread = 2, 
				missing.rate = missing_rate, verbose = FALSE)
	# clost the file and remove it
	snpgdsClose(gds)
	file.remove("temp.gds")


	# cycle through the VCF files
	for (rownumber in 1:nrow(vcf_list)) {
		# read in the input file and convert
		snpgdsVCF2GDS(vcf_list[rownumber, 2], "temp.gds", verbose = FALSE, method = "biallelic.only")
		gds <- snpgdsOpen("temp.gds")

		####
		# 1 Percent variation explained by first few PCS
		####
		# perform PCoA and calculate the proportion of variance explained by the first few axes
		pca <- snpgdsPCA(gds, num.thread = 2, verbose = FALSE, missing.rate = missing_rate, 
				autosome.only = FALSE)
		sumvar <- sum(pca$varprop[1: axes_to_use] * 100)
		# add to the table and report
		vcf_list[rownumber, 3] <- sumvar
		cat("Used", length(pca$snp.id), "SNPs to generate a PCoA; the first", axes_to_use, "axes explained", 
			sumvar, "% of the variation\n")

		####
		# 2 Relationship between missing data and variation
		####
		# calculate missing SNPs between samples
		# capture the genotypes in a dataframe - 012 (other numbers = missing)
		mymat <- read.gdsn(index.gdsn(gds, "genotype"))
		mydf <- as.data.frame(mymat)
		sample_names <- read.gdsn(index.gdsn(gds, "sample.id"))
		rownames(mydf) <- sample_names
		# limit to a SNP set with no more than a certain proportion missing data
		snpset <- snpgdsSelectSNP(gds, autosome.only = FALSE, missing.rate = missing_rate)
		sub_mat <- mymat[, snpset]
		sub_df <- mydf[, snpset]
		#sub_mat <- mymat
		#sub_df <- mydf
		# create a matrix and determine the missing data relationships between samples
		missing_matrix <- matrix(nrow = length(rownames(sub_df)), ncol = length(rownames(sub_df)))
		rownames(missing_matrix) <- rownames(sub_df)
		colnames(missing_matrix) <- rownames(sub_df)
		counter <- 1
		for (sample in 1: length(rownames(sub_df))) {
			if (counter %% 10 == 0) {
				cat("Processed", counter, "samples of", length(rownames(sub_df)), "\n")
			}
			counter <- counter + 1
			for (other_sample in 1: length(rownames(sub_df))) {
				if (sample == other_sample) {
					prop_miss <- NA
				} else {
					call1miss2 <- sub_mat[sample, ] == c(0,1,2) & sub_mat[other_sample, ] != c(0,1,2)
					miss1call2 <- sub_mat[sample, ] != c(0,1,2) & sub_mat[other_sample, ] == c(0,1,2)
					total_miss <- sum(c(call1miss2, miss1call2), na.rm = TRUE)
					call1call2 <- sub_mat[sample, ] == c(0,1,2) & sub_mat[other_sample, ] == c(0,1,2)
					total_poss <- sum(call1call2, na.rm = TRUE) + total_miss
					prop_miss <- total_miss / total_poss
				}
				# enter the value in the correct location in the missing_matrix
				missing_matrix[sample, other_sample] <- prop_miss
			}		
		}
		# create identity-by-state estimation and hierarchical clustering dendrogram (using same missing rate)
		#ibs_mat <- snpgdsIBS(gds, autosome.only = FALSE, missing.rate = missing_rate,
		#			num.thread = 2, verbose = FALSE)
		#ibs_est <- snpgdsHCluster(snpgdsIBS(gds, autosome.only = FALSE, 
		#			missing.rate = missing_rate, num.thread = 2, verbose = FALSE))
		# calculate correlation between missing data and genetic divergence
		ibs_est <- MIDPOINTibs_est
		miss_corr <- cor(c(missing_matrix[lower.tri(missing_matrix, diag = FALSE)]),
				c(1 - ibs_est$ibs[lower.tri(ibs_est$ibs, diag = FALSE)]), method = "pearson")
		# add to the table and report
		vcf_list[rownumber, 4] <- miss_corr
		cat("Used", ncol(sub_df), "SNPs to assess proportion missing data; max:", 
			max(missing_matrix, na.rm = TRUE), "min:", min(missing_matrix, na.rm = TRUE), 
			"and correlation with divergence:", miss_corr, "\n")


		####
		# 3 Relationship between geographic distance and variation
		####
		# limit to only ingroup samples, since this assumes within species
		# NOTE: if we are doing species delimitation, this may be a poor metric across all samples
		sample_names <- read.gdsn(index.gdsn(gds, "sample.id"))
		out_index <- grep("Z", sample_names)
		sampset <- sample_names[-out_index]

		# create identity-by-state estimation for divergence
		ibs_mat <- snpgdsIBS(gds, autosome.only = FALSE, num.thread = 2, missing.rate = missing_rate,
					sample.id = sampset, verbose = FALSE)
		# compute a geographic distance matrix based on samples present in the analysis
		geo_mat <- matrix(nrow = length(ibs_mat$sample.id), ncol = length(ibs_mat$sample.id))
		rownames(geo_mat) <- ibs_mat$sample.id
		colnames(geo_mat) <- ibs_mat$sample.id
		for (sample in 1: length(ibs_mat$sample.id)) {
			for (other_sample in 1: length(ibs_mat$sample.id)) {
				if (sample == other_sample) {
					dist <- NA
				} else {
					lat1 <- latlon_table[, 2][latlon_table[, 1] == ibs_mat$sample.id[sample]]
					lon1 <- latlon_table[, 3][latlon_table[, 1] == ibs_mat$sample.id[sample]]
					lat2 <- latlon_table[, 2][latlon_table[, 1] == ibs_mat$sample.id[other_sample]]
					lon2 <- latlon_table[, 3][latlon_table[, 1] == ibs_mat$sample.id[other_sample]]
					dist <- distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine) / 1000
				}
				# add to the correct matrix cell
				geo_mat[sample, other_sample] <- dist
			}
		}
		# calculate the correlation between divergence and distance
		dist_corr <- cor(c(1 - ibs_mat$ibs[lower.tri(ibs_mat$ibs, diag = FALSE)]),
				c(geo_mat[lower.tri(geo_mat, diag = FALSE)]), method = "pearson")
		# calculate the slope of the relationship between divergence and distance
		model <- lm(1 - ibs_mat$ibs[lower.tri(ibs_mat$ibs, diag = FALSE)] ~ 
				geo_mat[lower.tri(geo_mat, diag = FALSE)])
		slope <- model$coefficients[2]
		# add to the table and report
		vcf_list[rownumber, 5] <- dist_corr
		vcf_list[rownumber, 6] <- 100 * 100 * slope
		cat("Used", length(ibs_mat$snp.id), "SNPs to assess the relationship between divergence and geography\n",
			"The correlation and slope for clust", vcf_list[rownumber, 1], "is:", 
			dist_corr, "and", 100 * 100 * slope, "% SNP divergence per 100 km\n\n")


		# clost the file and remove it
		snpgdsClose(gds)
		file.remove("temp.gds")
	}


	# Plot the results
	# 1
	pdf(file = paste0(outpre, "_varPCoA.pdf"))
	plot(vcf_list$V1, vcf_list$V3, main = paste0("Cumulative variance in the first ", axes_to_use,
		" principal components"), ylab = "Cumulative variance (%)", xlab = "Clustering threshold (%)", 
		xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	invisible(dev.off())

	# 2
	pdf(file = paste0(outpre, "_missing.pdf"))
	plot(vcf_list$V1, vcf_list$V4, 
		main = "Correlation between missing data and genetic distance\n(estimated from middle parameter)",
	        ylab = "Pearson correlation between missing data and dissimilarity", 
		xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	invisible(dev.off())

	# 3
	pdf(file = paste0(outpre, "_geodist.pdf"))
	plot(vcf_list$V1, vcf_list$V5, main = "Correlation between genetic distance and geographic distance",
	        ylab = "Pearson correlation between dissimilarity and distance", xlab = "Clustering threshold (%)", 
		xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	plot(vcf_list$V1, vcf_list$V6, main = "Slope of the relationship between SNP divergence and distance",
	        ylab = "Increase in SNP divergence (%) per 100 km", xlab = "Clustering threshold (%)", 
		xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	invisible(dev.off())
}


