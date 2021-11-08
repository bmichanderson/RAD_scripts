
##########
# Author: Ben Anderson
# Date: Oct 2021
# Description: run a PCA on a VCF file and plot; now can alternatively plot from a covariance matrix from PCAngsd
##########


# set whether running interactively (won't work with command arguments...)
interactive <- FALSE


# load required libraries
suppressMessages(library(adegenet))
suppressMessages(library(vcfR))


# set colours (needs to be at least as long as number of pops or remainder will be black)
draw_colours <- c(
    "forestgreen",      # 22 8B 22 / 34 139 34
    "darkslateblue",    # 48 3D 8B / 72 61 139
    "lightskyblue",     # 87 CE FA / 135 206 250
    "darkorange",       # FF 8C 00 / 255 140 0
    "mediumorchid",     # BA 55 D3 / 186 85 211
    "lawngreen",        # 7C FC 00 / 124 252 0
	"steelblue",		# 46 82 B4 / 70	130	180
	"aquamarine",		# 7F FF D4 / 127 255 212
	"peru",				# CD 85 3F / 80 52 25
	"sandybrown",		# F4 A4 60 / 244 164 96
	"khaki",			# F0 E6 8C / 94 90 55
	"lavender"			# E6 E6 FA / 90 90 98
)


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to run a PCA and plot it using an input VCF file or a covariance matrix file\n",
			"It is assumed that the SNPs in the VCF are already one per locus and biallelic\n")
		cat("Usage: Rscript PCA_plot.R -o output -s sample_file [-v vcf_file, -c covmat_file]\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-s\tThe file with tab-delimited sample IDs and pops, one per line, for labelling/colouring the plot\n")
		cat("\t-v\tThe VCF file to be analysed (takes priority)\n")
		cat("\t-c\tThe covariance matrix file (not used if VCF present)\n")
	} else {
		cat(help_message)
	}
}

# a function to plot a PCA
plot_pca <- function(pca_scores, pca_var, x = 1, y = 2, pca_colour, legend, legend_colour, method) {
	par(mar = c(5.1, 4.1, 4.1, 6.1))
	plot(pca_scores[, x], pca_scores[, y],
		xlab = paste0("PC", x, " (", pca_var[x], "%)"),
		ylab = paste0("PC", y, " (", pca_var[y], "%)"),
		main = paste0("Principal Component Analysis\n(PC", y, " vs PC", x, ") using ", method),
		col = adjustcolor(pca_colour, alpha.f = 0.8),
		pch = 19,
		panel.first = {
			grid()
			abline(h = 0, v = 0)
		})
	legend("topright", legend,
			xpd = TRUE, inset = c(-0.2, 0),
			col = legend_colour, pch = 19)
	par(mar = c(5.1, 4.1, 4.1, 2.1))
}


# Read in and format the data

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
	covmat_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			output <- args[index + 1]
		} else if (args[index] == "-s") {
			samples_present <- TRUE
			samples_file <- args[index + 1]
		} else if (args[index] == "-v") {
			vcf_present <- TRUE
			vcf_file <- args[index + 1]
		} else if (args[index] == "-c") {
			covmat_present <- TRUE
			covmat_file <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}
if (! samples_present) {
	stop(help("Missing argument for sample/pops file!\n"), call. = FALSE)
}
if (all(! vcf_present, ! covmat_present)) {
	stop(help("Missing argument for vcf file or covariance matrix file!\n"), call. = FALSE)
}


# read in the input files
sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
if (vcf_present) {
	vcf <- read.vcfR(vcf_file, verbose = FALSE)
	cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
		length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")
} else if (covmat_present) {
	covmat <- as.matrix(read.table(covmat_file))
}


# Processing
if (vcf_present) {
	# convert to genlight
	genl <- vcfR2genlight(vcf)

	# check that the samples match
	if (length(indNames(genl)) != nrow(sample_table)) {
		stop(help("Number of samples in the VCF does not match number in the pops file\n"), call. = FALSE)
	}

	for (sample in sample_table$V1) {
		if (! sample %in% indNames(genl)) {
			stop(help("Sample", sample, "not found in VCF\n"), call. = FALSE)
		}
	}

	# assign population labels to the genlight object
	populations <- sample_table$V2[match(indNames(genl), sample_table$V1)]
	pop(genl) <- populations

	# create PCA
	pca <- glPca(genl, nf = 5, parallel = TRUE)

	# calculate percent variation explained for each PC
	pca_var <- round(100 * pca$eig / sum(pca$eig), digits = 1)

	# set the method used text for plotting
	method <- paste0("glPca on ", nLoc(genl), " SNPs")

	# set what the pca_scores comprise
	pca_scores <- pca$scores

} else if (covmat_present) {
	# capture pop labels
	populations <- sample_table$V2

	# convert to eigen
	eig <- eigen(covmat)

	# calculate percent variation explained for the eigenvalues
	pca_var <- round(100 * eig$values / sum(eig$values), digits = 1)

	# set the method used text for plotting
	method <- paste0("a covariance matrix based on allele frequencies")

	# set what the pca_scores comprise
	pca_scores <- eig$vectors
}


# assign colours to the individuals/populations
indiv_colours <- rep("black", length(populations))
pop_colours <- rep("black", length(unique(populations)))
ind <- 1
for (pop in unique(populations)) {
	if (ind > length(draw_colours)) {
		break
	}
	indiv_colours[grep(pop, populations)] <- draw_colours[ind]
	pop_colours[grep(pop, unique(populations))] <- draw_colours[ind]
	ind <- ind + 1
}


# Plot PCAs
pdf(paste0(output, "_pca.pdf"), paper = "A4")
plot_pca(pca_scores, pca_var, x = 1, y = 2, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pca(pca_scores, pca_var, x = 1, y = 3, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pca(pca_scores, pca_var, x = 1, y = 4, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pca(pca_scores, pca_var, x = 2, y = 3, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pca(pca_scores, pca_var, x = 3, y = 4, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
invisible(dev.off())


# To interact with points and labels for exploring use plotly
if (interactive) {
	suppressMessages(library("plotly"))
	data <- as.data.frame(pca_scores)
	if (vcf_present) {
		x_vals <- data$PC1
		y_vals <- data$PC2
		z_vals <- data$PC3
		text_field <- rownames(data)
	} else {
		x_vals <- data$V1
		y_vals <- data$V2
		z_vals <- data$V3
		text_field <- sample_table$V1
	}
	plot_ly(data = data, type = "scatter", mode = "markers",
			x = ~x_vals, y = ~y_vals, text = text_field, size = 15,
			color = populations, colors = pop_colours)
	plot_ly(data = data, type = "scatter3d", mode = "markers",
			x = ~x_vals, y = ~y_vals, z = ~z_vals,
			text = text_field, size = 15,
			color = populations, colors = pop_colours)
}
