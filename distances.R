
##########
# Author: Ben Anderson
# Date: Nov 2021
# Description: calculate distances between samples from a VCF file, output distance matrices as Nexus and plot NJ trees
##########


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
	"sandybrown",		# F4 A4 60 / 244 164 96
	"khaki",			# F0 E6 8C / 94 90 55
	"lavender"			# E6 E6 FA / 90 90 98
)


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to calculate distances from an input VCF file\n")
		cat("Usage: Rscript distances.R -o output -v vcf_file [-s samples_file]\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-v\tThe VCF file to be analysed\n")
		cat("\t-s\tColour the NJ tree by these tab-delimited sample IDs and pops, one per line [optional]\n")
	} else {
		cat(help_message)
	}
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
	vcf_present <- FALSE
	samples_present <- FALSE
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
if (! vcf_present) {
	stop(help("Missing argument for vcf file!\n"), call. = FALSE)
}


# read in the input files
if (samples_present) {
	sample_table <- read.table(samples_file, sep = "\t", header = FALSE)
}
vcf <- read.vcfR(vcf_file, verbose = FALSE)
cat("Read in a VCF with", ncol(vcf@gt) - 1, "samples,",
	length(unique(vcf@fix[, 1])), "loci and", nrow(vcf@fix), "SNPs\n")


# convert to genlight and DNAbin objects
genl <- vcfR2genlight(vcf)
dnabin <- vcfR2DNAbin(vcf, consensus = TRUE, extract.haps = FALSE)


# calculate distances between samples
## Euclidean
dist_euclidean <- dist(as.matrix(genl))
## GENPOFAD (allows better comparison between hets and homo; ambiguity codes)
dist_genpofad <- dist.snp(dnabin, model = "GENPOFAD")


# Output the distance matrices
taxa <- rownames(as.matrix(genl))
taxa_block <- paste0("BEGIN TAXA;\n\tDIMENSIONS NTAX=", length(taxa), ";\n\t",
					"TAXLABELS ", paste(taxa, collapse = " "), ";\nEND;\n")
dist_block <- paste0("BEGIN DISTANCES;\n\tFORMAT\n\t\tTRIANGLE=BOTH\n\t\tDIAGONAL\n\t\t",
					"LABELS=LEFT\n\t;\n\tMATRIX\n")
## Euclidean
outfile <- file(paste0(output, "_distEuc.nex"), open = "w")
writeLines("#NEXUS", con = outfile)
writeLines(taxa_block, con = outfile)
writeLines(dist_block, con = outfile)
write.table(as.matrix(dist_euclidean), file = outfile, col.names = FALSE,
			append = TRUE, quote = FALSE)
writeLines("\t;\nEND;\n", con = outfile)
close(outfile)
## GENPOFAD
outfile <- file(paste0(output, "_distGENP.nex"), open = "w")
writeLines("#NEXUS", con = outfile)
writeLines(taxa_block, con = outfile)
writeLines(dist_block, con = outfile)
write.table(as.matrix(dist_genpofad), file = outfile, col.names = FALSE,
			append = TRUE, quote = FALSE)
writeLines("\t;\nEND;\n", con = outfile)
close(outfile)


# If sample info is present, assign colours
if (samples_present) {
	populations <- sample_table$V2[match(indNames(genl), sample_table$V1)]
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
	tip_cols <- indiv_colours
} else {
	tip_cols <- rep("black", length(taxa))
}


# Plot NJ trees
pdf(paste0(output, "_nj.pdf"), paper = "A4")
## Euclidean
bionjtree <- bionjs(dist_euclidean)
plot(bionjtree, type = "unrooted", main = "BioNJ Tree using Euclidean distances",
	show.tip.label = TRUE, font = 1, cex = 0.25, edge.width = 0.75,
	underscore = TRUE, lab4ut = "axial", tip.color = tip_cols)
## GENPOFAD
bionjtree <- bionjs(dist_genpofad)
plot(bionjtree, type = "unrooted", main = "BioNJ Tree using GENPOFAD distances",
	show.tip.label = TRUE, font = 1, cex = 0.25, edge.width = 0.75,
	underscore = TRUE, lab4ut = "axial", tip.color = tip_cols)
invisible(dev.off())
