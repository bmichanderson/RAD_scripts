
##########
# Author: B.M. Anderson
# Date: Nov 2021
# Modified: Feb 2026 (adjusted figure size and added output to text file)
# Description: assess/plot the amount of missing data per sample in a VCF (arg1)
##########


# parse the command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	stop("Missing argument for VCF file", call. = FALSE)
} else {
	vcf_file <- args[1]
}


# load the library and analyse the VCF
suppressMessages(library(vcfR))
myvcf <- read.vcfR(vcf_file, verbose = FALSE)
gt <- extract.gt(myvcf, element = "GT")
missing <- apply(gt, MARGIN = 2, function(x) { sum(is.na(x)) })
missing <- 100 * missing / nrow(myvcf)

# set width relative to samples
mywidth <- max(40, length(missing) / 5)

# create PDF
pdf("missing.pdf", width = mywidth, height = 20)
par(mar = c(20, 4, 4, 2))
par(fig = c(0, 0.9, 0, 1))
barplot(missing, las = 2, cex.axis = 2, ylim = c(0, max(missing) + 0.1 * max(missing)))
par(fig = c(0.85, 1, 0, 1), new = TRUE)
boxplot(missing, axes = FALSE, ylim = c(0, max(missing) + 0.1 * max(missing)))
mtext("Percent missing data per sample", side = 3,
	outer = TRUE, line = -4, cex = 4)
invisible(dev.off())

# output the data to text file
connect <- file("missing.txt", open = "w")
write.table(missing, file = connect, col.names = FALSE, quote = FALSE)
close(connect)

# summarise
summary(missing)
