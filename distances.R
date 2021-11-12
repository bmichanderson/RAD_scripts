
##########
# Author: Ben Anderson
# Date: Nov 2021
# Description: calculate distances between samples from a VCF file,
#				output distance matrix as Nexus and plot NJ tree and PCoA
##########


# set whether running interactively (won't work with command arguments...)
interactive <- FALSE


# load required libraries
suppressMessages(library(adegenet))
suppressMessages(library(ape))
suppressMessages(library(pofadinr))
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
	"tomato",			# FF 63 47 / 100 39 28
	"khaki",			# F0 E6 8C / 94 90 55
	"deeppink",			# FF 14 93 / 100 8 58
	"slategrey"			# 70 80 90 / 44 50 56
)


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate distances from an input VCF file, then output the matrix and plot NJ and PCoA\n")
		cat("Usage: Rscript distances.R -d dist_method -o output -v vcf_file -s samples_file\n")
		cat("Options:\n")
		cat("\t-d\tThe distance method 1: GENPOFAD [default] or 2: Euclidean\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-s\tSamples and populations as tab-delimited sample IDs and pops, one per line\n")
	} else {
		cat(help_message)
	}
}

# a function to plot a PCoA
plot_pcoa <- function(pcoa_points, pcoa_var, x = 1, y = 2, pcoa_colour, legend, legend_colour, method) {
	par(mar = c(5.1, 4.1, 4.1, 6.1))
	plot(pcoa_points[, x], pcoa_points[, y],
		xlab = paste0("PC", x, " (", pcoa_var[x], "%)"),
		ylab = paste0("PC", y, " (", pcoa_var[y], "%)"),
		main = paste0("Principal Coordinate Analysis\n(PC", y, " vs PC", x, ") using ", method),
		col = adjustcolor(pcoa_colour, alpha.f = 0.8),
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
	dist_method <- 1
	output <- "output"
	vcf_present <- FALSE
	samples_present <- FALSE
	for (index in seq_len(length(args))) {
		if (args[index] == "-d") {
			dist_method <- as.integer(args[index + 1])
		} else if (args[index] == "-o") {
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
if (any(c(! vcf_present, ! samples_present))) {
	stop(help("Missing argument for vcf file and/or samples file!\n"), call. = FALSE)
}


# read in the input files
if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
}
vcf <- read.vcfR(vcf_file, verbose = FALSE)
cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
	length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")


# calculate distances between samples
if (dist_method == 1) {
	dnabin <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE)
	individuals <- rownames(dnabin)
	## GENPOFAD (allows better comparison between hets and homo; ambiguity codes)
	distance <- dist.snp(dnabin, model = "GENPOFAD")
	dist_suffix <- "_distGENPOFAD.nex"
	method <- paste0("GENPOFAD distances from ", dim(dnabin)[2], " SNPs")
} else if (dist_method == 2) {
	genl <- vcfR2genlight(vcf)
	individuals <- indNames(genl)
	## Euclidean
	distance <- dist(as.matrix(genl))
	dist_suffix <- "_distEuclidean.nex"
	method <- paste0("Euclidean distances from ", nLoc(genl), " SNPs")
} else {
	stop(help("Distance method needs to be specified as 1 or 2!\n"), call. = FALSE)
}


# Output the distance matrix
taxa <- individuals
taxa_block <- paste0("BEGIN TAXA;\n\tDIMENSIONS NTAX=", length(taxa), ";\n\t",
					"TAXLABELS ", paste(taxa, collapse = " "), ";\nEND;\n")
dist_block <- paste0("BEGIN DISTANCES;\n\tFORMAT\n\t\tTRIANGLE=BOTH\n\t\tDIAGONAL\n\t\t",
					"LABELS=LEFT\n\t;\n\tMATRIX\n")
outfile <- file(paste0(output, dist_suffix), open = "w")
writeLines("#NEXUS", con = outfile)
writeLines(taxa_block, con = outfile)
writeLines(dist_block, con = outfile)
write.table(as.matrix(distance), file = outfile, col.names = FALSE,
			append = TRUE, quote = FALSE)
writeLines("\t;\nEND;\n", con = outfile)
close(outfile)


# Assign colours for plotting
populations <- sample_table$V2[match(individuals, sample_table$V1)]
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


# Plot NJ trees
pdf(paste0(output, "_nj.pdf"), paper = "A4")
bionjtree <- bionjs(distance)
plot(bionjtree, type = "unrooted", main = paste0("BioNJ Tree using ", method),
	show.tip.label = TRUE, font = 1, cex = 0.25, edge.width = 0.75,
	underscore = TRUE, lab4ut = "axial", tip.color = indiv_colours)
invisible(dev.off())


# Calculate and plot PCoA
pdf(paste0(output, "_pcoa.pdf"), paper = "A4")
pcoa <- cmdscale(distance, k = 4, eig = TRUE)
pcoa_var <- round(100 * pcoa$eig / sum(pcoa$eig), digits = 1)
pcoa_points <- pcoa$points
plot_pcoa(pcoa_points, pcoa_var, x = 1, y = 2, pcoa_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pcoa(pcoa_points, pcoa_var, x = 1, y = 3, pcoa_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pcoa(pcoa_points, pcoa_var, x = 1, y = 4, pcoa_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pcoa(pcoa_points, pcoa_var, x = 2, y = 3, pcoa_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
plot_pcoa(pcoa_points, pcoa_var, x = 3, y = 4, pcoa_colour = indiv_colours,
		legend = unique(populations), legend_colour = pop_colours, method = method)
invisible(dev.off())


# To interact with points and labels for exploring use plotly
if (interactive) {
	suppressMessages(library("plotly"))
	data <- as.data.frame(pcoa_points)
	x_vals <- data$V1
	y_vals <- data$V2
	z_vals <- data$V3
	text_field <- sample_table$V1
	plot_ly(data = data, type = "scatter", mode = "markers",
			x = ~x_vals, y = ~y_vals, text = text_field, size = 15,
			color = populations, colors = pop_colours)
	plot_ly(data = data, type = "scatter3d", mode = "markers",
			x = ~x_vals, y = ~y_vals, z = ~z_vals,
			text = text_field, size = 15,
			color = populations, colors = pop_colours)
}
