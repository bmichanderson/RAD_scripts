
##########
# Authors: Ben Anderson, based partly on McCartney-Melstad et al. 2019 (https://github.com/atcg/clustOpt)
# Date: July-Aug 2021
# Description: summarize various metrics for optimizing parameters for ipyrad assembly of RAD data
##########


# set a couple parameters for calculations
axes_to_use <- 8		# the number of PCoA axes to consider for sum of variation explained
missing_rate <- 0.96		# the maximum proportion of missing data per SNP (e.g. 0.9 is present in >=10%)
				# a value of 1.0 means all SNPs, a value of 0.96 means SNPs in >= 4% of samples


# a helper function for when the script is called without arguments
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to summarize outputs from tests/queries of metrics used",
			"for evaluating clustering thresholds in RAD datasets\n")
		cat("Usage: Rscript MM_ipyrad_summary.R -o output_prefix -v vcf_file\n")
		cat("Options (none are required except for -v):\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-e\tThe seq error stats file\n")
		cat("\t-h\tThe heterozygosity stats file\n")
		cat("\t-l\tThe locus and snps recovered stats file\n")
		cat("\t--pall\tThe paralogs1 stats file (based on full assembly)\n")
		cat("\t--pind\tThe paralogs2 stats file (per sample)\n")
		cat("\t-t\tThe bootstraps stats file\n")
		cat("\t--merr\tThe Mastretta-Yanes error rates file\n")
		cat("\t--mdis\tThe Mastretta-Yanes population distances file\n")
		cat("\t--mparis\tThe Paris et al count of loci and SNPs found in >= 80% of individuals\n")
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
	p1stats_present <- FALSE
	p2stats_present <- FALSE
	tstats_present <- FALSE
	merror_present <- FALSE
	mdist_present <- FALSE
	mparis_present <- FALSE
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
		} else if (args[index] == "--pall") {
			p1stats_present <- TRUE
			p1stats_file <- args[index + 1]
		} else if (args[index] == "--pind") {
			p2stats_present <- TRUE
			p2stats_file <- args[index + 1]
		} else if (args[index] == "-t") {
			tstats_present <- TRUE
			tstats_file <- args[index + 1]
		} else if (args[index] == "--merr") {
			merror_present <- TRUE
			merror_file <- args[index + 1]
		} else if (args[index] == "--mdis") {
			mdist_present <- TRUE
			mdist_file <- args[index + 1]
		} else if (args[index] == "--mparis") {
			mparis_present <- TRUE
			mparis_file <- args[index + 1]
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

if (any(p1stats_present, p2stats_present)) {
	# calculate and plot the proportion of paralogs filtered
	pdf(file = paste0(outpre, "_paralogs.pdf"))
	if (p1stats_present) {
		paralogs_table <- read.table(p1stats_file, sep = "\t", header = FALSE)
		paralogs_table$V4 <- paralogs_table$V2 / paralogs_table$V3
		plot(paralogs_table$V1, paralogs_table$V4, 
			main = "Proportion paralogous loci filtered\n(using max_shared_Hs_locus)", 
			ylab = "Proportion paralogous loci", 
			xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
		axis(1, at = seq(min(paralogs_table$V1), max(paralogs_table$V1), by = 1))
	}
	if (p2stats_present) {
		paralogs_table <- read.table(p2stats_file, sep = "\t", header = FALSE)
		paralogs_table$V4 <- paralogs_table$V3 / paralogs_table$V2
		boxplot(paralogs_table$V4 ~ paralogs_table$V1, 
			main = paste0("Proportion paralogous loci filtered\n(using maxAlleles)"),
			ylab = "Proportion paralogous loci", xlab = "Clustering threshold (%)")
	}
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
		main = paste0("Locus error rates\n",
			"(loci present in one rep but not the other, relative to total in either/both)"),
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
		ylab = "Average intrapopulation Euclidean distance (normalized %)", 
		xlab = "Clustering threshold (%)")
	invisible(dev.off)
}

if (mparis_present) {
	# plot Paris et al counts of loci and SNPs in >= 80% of individuals
	pdf(file = paste0(outpre, "_mparis.pdf"))
	count_table <- read.table(mparis_file, sep = "\t", header = FALSE)
	plot(count_table$V1, count_table$V2, main = "Total number of loci recovered in >= 80% of individuals", 
		ylab = "Number of loci", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(count_table$V1), max(count_table$V1), by = 1))
	plot(count_table$V1, count_table$V3, main = "Total number of SNPs recovered in >= 80% of individuals", 
		ylab = "Number of SNPs", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(count_table$V1), max(count_table$V1), by = 1))
	invisible(dev.off())
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

	# cycle through the VCF files
	for (rownumber in 1:nrow(vcf_list)) {
		# read in the input file and convert
		snpgdsVCF2GDS(vcf_list[rownumber, 2], "temp.gds", verbose = FALSE, method = "biallelic.only")
		gds <- snpgdsOpen("temp.gds")
		cat("Read in a VCF file with", length(read.gdsn(index.gdsn(gds, "snp.id"))), "SNPs\n")

		# exclude samples not present in the sample file
		sample_names <- read.gdsn(index.gdsn(gds, "sample.id"))
		matches <- latlon_table$V1[match(sample_names, latlon_table$V1)]
		samp_include <- sample_names[!is.na(matches)]
		cat("Excluding", length(sample_names) - length(samp_include),
			"samples not present in the sample file\n")

		# exclude low coverage and monomorphic SNPs based on possibly a reduced sample set
		loci_include <- snpgdsSelectSNP(gds, sample.id = samp_include, 
							missing.rate = missing_rate, remove.monosnp = TRUE, 
							autosome.only = FALSE, verbose = FALSE)

		# capture how many SNPs are left
		num_snps <- length(loci_include)
		cat("Found", num_snps, "SNPs after filters for low coverage",
			"(<", 1 - missing_rate,") and monomorphic SNPs\n")
		vcf_list[rownumber, 3] <- num_snps


		####
		# 1 Percent variation explained by first few PCS
		####
		# perform PCoA and calculate the proportion of variance explained by the first few axes
		pca <- snpgdsPCA(gds, snp.id = loci_include, sample.id = samp_include, 
					num.thread = 2, autosome.only = FALSE, verbose = FALSE)
		sumvar <- sum(pca$varprop[1: axes_to_use] * 100)
		# add to the table and report
		vcf_list[rownumber, 4] <- sumvar
		cat("Used", length(pca$snp.id), "SNPs to generate a PCoA; the first", axes_to_use, "axes explained", 
			sumvar, "% of the variation\n")


		####
		# create identity-by-state estimation and dendrogram
		ibs_est <- snpgdsIBS(gds, sample.id = samp_include, snp.id = loci_include,  
					num.thread = 2, autosome.only = FALSE, verbose = FALSE)
		ibs_dend <- snpgdsHCluster(ibs_est)


		# (this next step may or may not be desired, but reduces computation time for missing)
		# randomly exclude original SNPs if there are multiple per "chromosome" = locus
		# from https://stackoverflow.com/questions/8041720/randomly-select-on-data-frame-for-unique-rows
		chromosomes <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
		loci_names <- read.gdsn(index.gdsn(gds, "snp.id"))
		snps_index <- loci_names %in% loci_include
		chromosomes <- chromosomes[snps_index]
		loci_names <- loci_names[snps_index]
		if (length(duplicated(chromosomes)[duplicated(chromosomes)]) > 0) {
			df <- chromosomes
			loci <- sample(1: length(df))
			df <- df[loci]
			dups <- duplicated(df)
			df <- df[!dups]
			loci <- loci[!dups]
			loci <- loci[sort(loci, index.return = TRUE)$ix]
			loci_include <- loci_names[loci]
			loci_exclude <- loci_names[-loci]
			cat("Excluding", length(loci_exclude), "SNPs and keeping", 
				length(loci_include), "SNPs, one (random) per locus\n")
		}


		####
		# 2 Relationship between missing data and genomic dissimilarity
		####
		# calculate missing SNPs between samples (this is ~ equivalent to loci if previous step was run)
		# capture the genotypes in a dataframe - 012 (other numbers = missing)
		mymat <- read.gdsn(index.gdsn(gds, "genotype"))
		sample_names <- read.gdsn(index.gdsn(gds, "sample.id"))
		rownames(mymat) <- sample_names
		loci_names <- read.gdsn(index.gdsn(gds, "snp.id"))
		colnames(mymat) <- loci_names
		# limit to the selected samples and SNP set
		sub_mat <- mymat[rownames(mymat) %in% samp_include, colnames(mymat) %in% loci_include]
		# create a matrix and determine the missing data relationships between samples
		# missing data is computed as: missed in either (not both) / (missed in either + missed in neither)
		missing_matrix <- matrix(nrow = nrow(sub_mat), ncol = nrow(sub_mat))
		rownames(missing_matrix) <- rownames(sub_mat)
		colnames(missing_matrix) <- rownames(sub_mat)
		counter <- 1
		cat("Processing samples for pairwise missing data across all SNPs\n")
		for (samp_ind in 1: nrow(sub_mat)) {
			if (counter %% 10 == 0) {
				cat("Processed", counter, "samples of", nrow(sub_mat), "\n")
			}
			counter <- counter + 1
			for (samp2_ind in 1: nrow(sub_mat)) {
				if (samp_ind == samp2_ind) {
					prop_miss <- NA
				} else {
					calls1 <- sub_mat[samp_ind, ] %in% c(0,1,2)
					calls2 <- sub_mat[samp2_ind, ] %in% c(0,1,2)
					call1miss2 <- calls1 & !calls2
					miss1call2 <- !calls1 & calls2
					total_miss <- sum(c(call1miss2, miss1call2), na.rm = TRUE)
					call1call2 <- calls1 & calls2 
					total_poss <- total_miss + sum(call1call2, na.rm = TRUE)
					prop_miss <- total_miss / total_poss
				}
				# enter the value in the correct location in the missing_matrix
				missing_matrix[samp_ind, samp2_ind] <- prop_miss
			}		
		}
		# estimate correlation between missing rate and genetic dissimilarity
		miss_corr <- cor(c(missing_matrix[lower.tri(missing_matrix, diag = FALSE)]),
				c(ibs_dend$dist[lower.tri(ibs_dend$dist, diag = FALSE)]), method = "pearson")
		# add to the table and report
		vcf_list[rownumber, 5] <- miss_corr
		cat("Used", ncol(sub_mat), "SNPs to assess proportion missing data:\n",
			min(missing_matrix, na.rm = TRUE), "(min) to", max(missing_matrix, na.rm = TRUE), "(max)\n", 
			"Correlation with divergence:", miss_corr, "\n")


		####
		# 3 Relationship between geographic distance and genomic dissimilarity
		####
		# ideally limit to only ingroup samples, since this assumes within species
		# NOTE: if we are doing species delimitation, this may be a poor metric across all samples
		# compute a geographic distance matrix based on samples present in the analysis
		geo_mat <- matrix(nrow = length(ibs_est$sample.id), ncol = length(ibs_est$sample.id))
		rownames(geo_mat) <- ibs_est$sample.id
		colnames(geo_mat) <- ibs_est$sample.id
		for (samp_ind in 1: length(ibs_est$sample.id)) {
			for (samp2_ind in 1: length(ibs_est$sample.id)) {
				if (samp_ind == samp2_ind) {
					dist <- NA
				} else {
					samp <- ibs_est$sample.id[samp_ind]
					samp2 <- ibs_est$sample.id[samp2_ind]
					lat1 <- latlon_table[, 2][latlon_table[, 1] == samp]
					lon1 <- latlon_table[, 3][latlon_table[, 1] == samp]
					lat2 <- latlon_table[, 2][latlon_table[, 1] == samp2]
					lon2 <- latlon_table[, 3][latlon_table[, 1] == samp2]
					dist <- distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine) / 1000
				}
				# add to the correct matrix cell
				geo_mat[samp_ind, samp2_ind] <- dist
			}
		}
		# calculate the correlation between dissimilarity (1 - similarity) and distance
		dist_corr <- cor(c(1 - ibs_est$ibs[lower.tri(ibs_est$ibs, diag = FALSE)]),
				c(geo_mat[lower.tri(geo_mat, diag = FALSE)]), method = "pearson")
		# calculate the slope of the relationship between dissimilarity and distance (use linear model)
		model <- lm(1 - ibs_est$ibs[lower.tri(ibs_est$ibs, diag = FALSE)] ~ 
				geo_mat[lower.tri(geo_mat, diag = FALSE)])
		slope <- model$coefficients[2]
		# add to the table and report
		vcf_list[rownumber, 6] <- dist_corr
		vcf_list[rownumber, 7] <- 100 * 100 * slope
		cat("Used", length(ibs_est$snp.id), 
			"SNPs to assess the relationship between genomic dissimilarity and geography\n",
			"The correlation and slope for clust", vcf_list[rownumber, 1], "is:", 
			dist_corr, "and", 100 * 100 * slope, "% increased genomic dissimilarity per 100 km\n\n")


		# close the file and remove it
		snpgdsClose(gds)
		file.remove("temp.gds")
	}


	# Plot the results
	# 1
	pdf(file = paste0(outpre, "_varPCoA.pdf"))
	plot(vcf_list$V1, vcf_list$V4, main = paste0("Cumulative variance in the first ", axes_to_use,
		" principal components"), ylab = "Cumulative variance (%)", xlab = "Clustering threshold (%)", 
		xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	invisible(dev.off())

	# 2
	pdf(file = paste0(outpre, "_missing.pdf"))
	plot(vcf_list$V1, vcf_list$V5, 
		main = "Correlation between missing data and genomic dissimilarity",
		ylab = "Pearson correlation between missing data and dissimilarity", 
		xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	# also plot total number of SNPs observed in the VCFs (will be different from overall summary stats)
	plot(vcf_list$V1, vcf_list$V3, main = "Total number of SNPs recovered", 
		ylab = "Number of SNPs", xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	invisible(dev.off())

	# 3
	pdf(file = paste0(outpre, "_geodist.pdf"))
	plot(vcf_list$V1, vcf_list$V6, main = "Correlation between genomic dissimilarity and geographic distance",
		ylab = "Pearson correlation between dissimilarity and geographic distance", 
		xlab = "Clustering threshold (%)", xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	plot(vcf_list$V1, vcf_list$V7, main = "Slope of the relationship between genomic dissimilarity and distance",
		ylab = "Increase in genomic dissimilarity (%) per 100 km", xlab = "Clustering threshold (%)", 
		xaxt = "n", pch = 16)
	axis(1, at = seq(min(vcf_list$V1), max(vcf_list$V1), by = 1))
	invisible(dev.off())
}

