
##########
# Author: Ben Anderson
# Date: Oct 2021
# Description: run a PCA on a VCF file and plot
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
		cat("A script to run a PCA and plot it using an input VCF file\n",
			"It is assumed that the SNPs in the VCF are already one per locus and biallelic\n")
		cat("Usage: Rscript PCA_plot.R -o output -s sample_file -v vcf_file\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-s\tThe file with tab-delimited sample IDs and pops, one per line, for labelling/colouring the plot\n")
		cat("\t-v\tThe VCF file to be analysed\n")
	} else {
		cat(help_message)
	}
}

# a function to plot a PCA
plot_pca <- function(pca_scores, pca_var, x = 1, y = 2, pca_colour, legend, legend_colour, num_snps) {
	par(mar = c(5.1, 4.1, 4.1, 6.1))
	plot(pca_scores[, x], pca_scores[, y],
		xlab = paste0("PC", x, " (", pca_var[x], "%)"),
		ylab = paste0("PC", y, " (", pca_var[y], "%)"),
		main = paste0("Principal Component Analysis\n(PC", y, " vs PC", x, ") using ", num_snps, " SNPs"),
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
cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
	length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")


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

# create PCA
pca <- glPca(genl, nf = 5, parallel = TRUE)


# calculate percent variation explained for each PC
pca_var <- round(100 * pca$eig / sum(pca$eig), digits = 1)


# Plot PCAs
pdf(paste0(output, "_pca.pdf"), paper = "A4")
plot_pca(pca$scores, pca_var, x = 1, y = 2, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, num_snps = nLoc(genl))
plot_pca(pca$scores, pca_var, x = 1, y = 3, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, num_snps = nLoc(genl))
plot_pca(pca$scores, pca_var, x = 1, y = 4, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, num_snps = nLoc(genl))
plot_pca(pca$scores, pca_var, x = 2, y = 3, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, num_snps = nLoc(genl))
plot_pca(pca$scores, pca_var, x = 3, y = 4, pca_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, num_snps = nLoc(genl))
invisible(dev.off())


# To interact with points and labels for exploring use plotly
if (interactive) {
	suppressMessages(library("plotly"))
	data <- as.data.frame(pca$scores)
	plot_ly(data = data, x = ~PC1, y = ~PC2,
			text = rownames(data), size = 15,
			color = genl$pop, colors = pop_colours)
	plot_ly(data = data, x = ~PC1, y = ~PC3,
			text = rownames(data), size = 15,
			color = genl$pop, colors = pop_colours)
	plot_ly(data = data, x = ~PC1, y = ~PC4,
			text = rownames(data), size = 15,
			color = genl$pop, colors = pop_colours)
	plot_ly(data = data, x = ~PC2, y = ~PC3,
			text = rownames(data), size = 15,
			color = genl$pop, colors = pop_colours)
	plot_ly(data = data, x = ~PC3, y = ~PC4,
			text = rownames(data), size = 15,
			color = genl$pop, colors = pop_colours)
}
