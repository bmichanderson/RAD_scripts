
#############
# Author: B. Anderson based partly on O'Leary et al. 2018 DOI: 10.1111/mec.14792
# Date: Sep 2021
# Description: plot summaries of VCF filtering effects
#############


## define a helper function for calls
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to plot summaries (histograms, etc.) after filtering steps with vcftools\n")
		cat("Usage: Rscript filter_plot.R -p prefix\n")
		cat("Option:\n")
		cat("\t-p\tThe file name prefix used when running vcftools\n")
        cat("\t\tNOTE: this script expects the following files to be present in the working directory:\n")
        cat("\t\t*.idepth *.imiss *.ldepth.mean *.lmiss *.het\n")
	} else {
		cat(help_message)
	}
}


## define a function to calculate the mode and one to estimate from continuous distribution
# from https://stackoverflow.com/questions/2547402/how-to-find-the-statistical-mode
calc_mode <- function(vector_x) {
    uniq_vals <- unique(vector_x)
    val_tab <- tabulate(match(vector_x, uniq_vals))
    uniq_vals[val_tab == max(val_tab)]
}
est_mode <- function(vector_x) {
    dens <- density(vector_x)
    dens$x[which.max(dens$y)]
}


## parse command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
    prefix_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-p") {
			prefix <- args[index + 1]
            prefix_present <- TRUE
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
    }
}
if (!prefix_present) {
    stop(help(), call. = FALSE)
}


## load files and put into single dataframes
# individuals
idepth <- read.table(paste0(prefix, ".idepth"), header = TRUE, stringsAsFactors = FALSE)
imiss <- read.table(paste0(prefix, ".imiss"), header = TRUE, stringsAsFactors = FALSE)
het <- read.table(paste0(prefix, ".het"), header = TRUE, stringsAsFactors = FALSE)
indiv <- data.frame(idepth)
indiv$MISS <- imiss$F_MISS[match(indiv$INDV, imiss$INDV)]
indiv$Fis <- het$F[match(indiv$INDV, het$INDV)]
# sites
ldepth <- read.table(paste0(prefix, ".ldepth.mean"), header = TRUE, stringsAsFactors = FALSE)
lmiss <- read.table(paste0(prefix, ".lmiss"), header = TRUE, stringsAsFactors = FALSE)
site <- data.frame(ldepth)
site$MISS <- lmiss$F_MISS[(site$CHR == lmiss$CHR) && (site$POS == lmiss$POS)]
counts <- table(site$CHR)
count_df <- as.data.frame(counts)
count_df$MEAN_DEPTH <- ldepth$MEAN_DEPTH[match(count_df$Var1, ldepth$CHR)]
# note this is the first SNP depth, not the entire locus


## start a pdf for plotting
pdf(file = paste0(prefix, "_plots.pdf"))


## histograms
# depth
mean_val <- mean(indiv$MEAN_DEPTH, na.rm = TRUE)
mode_val <- est_mode(indiv$MEAN_DEPTH)
hist(indiv$MEAN_DEPTH, main = paste("Individual Read Depth\nMean:", round(mean_val, 2), "Mode:", round(mode_val, 2)),
    xlab = "Mean read depth per individual", ylab = "Number of individuals")
abline(v = mean_val, col = "forestgreen", lty = 2)

mean_val <- mean(site$MEAN_DEPTH, na.rm = TRUE)
mode_val <- est_mode(site$MEAN_DEPTH)
if (max(site$MEAN_DEPTH) > 1000) {
    set_max <- 1000
} else {
    set_max <- max(site$MEAN_DEPTH)
}
hist(site$MEAN_DEPTH, xlim = c(0, set_max),
    main = paste("Site Read Depth\nMean:", round(mean_val, 2), "Mode:", round(mode_val, 2),
            "Max:", max(site$MEAN_DEPTH)),
    xlab = "Mean read depth per site", ylab = "Number of sites")
abline(v = mean_val, col = "forestgreen", lty = 2)

# missingness
mean_val <- mean(indiv$MISS, na.rm = TRUE)
mode_val <- est_mode(indiv$MISS)
hist(indiv$MISS, xlim = c(0, 1),
    main = paste("Individual Missing Data\nMean:", round(mean_val, 2), "Mode:", round(mode_val, 2)),
    xlab = "Proportion missing data per individual", ylab = "Number of individuals")
abline(v = mean_val, col = "forestgreen", lty = 2)

mean_val <- mean(site$MISS, na.rm = TRUE)
mode_val <- est_mode(site$MISS)
hist(site$MISS, xlim = c(0, 1),
    main = paste("Site Missing Data\nMean:", round(mean_val, 2), "Mode:", round(mode_val, 2)),
    xlab = "Proportion missing data per site", ylab = "Number of sites")
abline(v = mean_val, col = "forestgreen", lty = 2)

# heterozygosity
mean_val <- mean(indiv$Fis, na.rm = TRUE)
mode_val <- est_mode(indiv$Fis)
hist(indiv$Fis, xlim = c(0, 1),
    main = paste("Individual Heterozygosity (Fis)\nMean:", round(mean_val, 2), "Mode:", round(mode_val, 2)),
    xlab = "Heterozygosity per individual (Fis)", ylab = "Number of individuals")
abline(v = mean(indiv$Fis, na.rm = TRUE), col = "forestgreen", lty = 2)

# sites per locus
mean_val <- mean(c(counts), na.rm = TRUE)
mode_val <- calc_mode(c(counts))
hist(c(counts), breaks = max(c(counts)),
    main = paste("Sites per locus\nMean:", round(mean_val, 2), "Mode:", mode_val),
    xlab = "Number of sites per locus", ylab = "Number of loci")
abline(v = mean(c(counts), na.rm = TRUE), col = "forestgreen", lty = 2)


## plots
plot(indiv$MEAN_DEPTH, indiv$MISS,
    xlab = "Mean read depth per individual", ylab = "Proportion missing data per individual")
abline(v = mean(indiv$MEAN_DEPTH, na.rm = TRUE), col = "forestgreen", lty = 2)
abline(h = mean(indiv$MISS, na.rm = TRUE), col = "blue", lty = 2)

plot(indiv$Fis, indiv$MISS,
    xlab = "Heterozygosity per individual", ylab = "Proportion missing data per individual")
abline(v = mean(indiv$Fis, na.rm = TRUE), col = "forestgreen", lty = 2)
abline(h = mean(indiv$MISS, na.rm = TRUE), col = "blue", lty = 2)

plot(indiv$Fis, indiv$MEAN_DEPTH,
    xlab = "Heterozygosity per individual", ylab = "Mean read depth per individual")
abline(v = mean(indiv$Fis, na.rm = TRUE), col = "forestgreen", lty = 2)
abline(h = mean(indiv$MEAN_DEPTH, na.rm = TRUE), col = "blue", lty = 2)

plot(site$MEAN_DEPTH, site$MISS,
    xlab = "Mean read depth per site", ylab = "Proportion missing data per site")
abline(v = mean(site$MEAN_DEPTH, na.rm = TRUE), col = "forestgreen", lty = 2)
abline(h = mean(site$MISS, na.rm = TRUE), col = "blue", lty = 2)

plot(c(counts), count_df$MEAN_DEPTH,
    xlab = "Number of sites per locus", ylab = "Mean read depth for first site")
abline(v = mean(c(counts), na.rm = TRUE), col = "forestgreen", lty = 2)
abline(h = mean(site$MEAN_DEPTH, na.rm = TRUE), col = "blue", lty = 2)


## stop plotting
invisible(dev.off())

## report details of the dataset
cat("\nThe", prefix, "dataset contains", nrow(indiv), "individuals,",
    length(unique(site$CHR)), "loci and", nrow(site), "sites\n\n")
