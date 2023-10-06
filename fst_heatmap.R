
##########
# Author: Ben Anderson
# Date: Mar 2022
# Modified: Oct 2023 (changed orientation and legend)
# Description: create a population heatmap based on an input Fst matrix
##########


# Define functions

# a helper function for errors or no args
help <- function(help_message) {
	if (missing(help_message)) {
		cat("A script to create a heatmap from an input Fst matrix\n")
		cat("Usage: Rscript fst_heatmap.R -o out_pref -f fst_matrix -p pops_file -c color_pal\n")
		cat("Options:\n")
		cat("\t-o\tThe output file name prefix [default output]\n")
		cat("\t-f\tThe Fst matrix\n")
		cat("\t-p\tPopulations, one per line in the order desired for output [optional]\n")
		cat("\t-c\tColour palette to use [default \"greens\"]\n")
	} else {
		cat(help_message)
	}
}

# a function to plot a heatmap
# for making legend, see https://stackoverflow.com/a/13355440 and https://stackoverflow.com/a/70522655
heatmapper <- function(fstmat, palette = "greens", ...) {
	# heatmap
	image(fstmat,
		col = hcl.colors(n = 100, palette = palette, rev = TRUE),
		axes = FALSE,
		main = expression("Pairwise F"[ST]),
		useRaster = TRUE)
	axis(1, at = seq(1, 0, length.out = ncol(fstmat)),
		labels = colnames(fstmat), las = 2, lwd.ticks = 0)
	axis(2, at = seq(0, 1, length.out = ncol(fstmat)),
		labels = colnames(fstmat), las = 2, lwd.ticks = 0)
	# legend
	subx <- grconvertX(c(0.7, 0.9), from = "user", to = "ndc")
	suby <- grconvertY(c(0.6, 1), from = "user", to = "ndc")
	op <- par(fig = c(subx, suby),
		mar = c(1, 1, 1, 0),
		new = TRUE)
	legend_colours <- as.raster(hcl.colors(n = 100, palette = palette))
	range_min <- min(fstmat, na.rm = TRUE)
	range_max <- max(fstmat, na.rm = TRUE)
	if (range_max - range_min < 0.2) {
		mult <- 100
		step <- 0.02
		rnum <- 2
	} else {
		mult <- 10
		step <- 0.1
		rnum <- 1
	}
	legend_seq <- seq(ceiling(range_min * mult) / mult,
		floor(range_max * mult) / mult, by = step)
	legend_labels <- format(round(legend_seq, rnum), nsmall = rnum)
	plot(x = c(0, 2), y = c(0, 1), type = "n",
		axes = FALSE, xlab = "", ylab = "", main = "")
	axis(side = 4, at = legend_seq / range_max, pos = 1, labels = FALSE,
		col = 0, col.ticks = 1)
	mtext(legend_labels, side = 4, line = -0.5, at = legend_seq / range_max, las = 2)
	rasterImage(legend_colours, xleft = 0, ybottom = 0, xright = 1, ytop = 1)
	par(op)
}


# Read in and format the data

# parse the command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
	stop(help(), call. = FALSE)
} else {
	catch_args <- vector("list")
	i <- 1
	out_pref <- "output"
	fst_present <- FALSE
	pops_present <- FALSE
	col_pal <- "greens"
	for (index in seq_len(length(args))) {
		if (args[index] == "-o") {
			out_pref <- args[index + 1]
		} else if (args[index] == "-p") {
			pops_present <- TRUE
			pops_file <- args[index + 1]
		} else if (args[index] == "-f") {
			fst_present <- TRUE
			fst_file <- args[index + 1]
		} else if (args[index] == "-c") {
			col_pal <- args[index + 1]
		} else {
			catch_args[i] <- args[index]
			i <- i + 1
		}
	}
}

if (! fst_present) {
	stop(help("Missing argument for Fst matrix!\n"), call. = FALSE)
} else if (! pops_present) {
	cat("No pops file specified, so plotting all in the given order\n")
}


# read in the input files; filter if pops_present
fst_mat <- as.matrix(read.csv(fst_file, sep = " ", header = TRUE))
if (pops_present) {
	pops_table <- read.table(pops_file, header = FALSE)
	fst_mat <- fst_mat[rownames(fst_mat) %in% pops_table[, 1], colnames(fst_mat) %in% pops_table[, 1]]
}


# convert the Fst matrix to square (account for multiple input orientations)
# this will allow us to re-order it if pops_present
if (is.na(fst_mat[1, 2])) {	# lower triangle matrix
	cat("Assuming lower trianglular Fst matrix input\n")
	fst_mat[upper.tri(fst_mat)] <- t(fst_mat)[upper.tri(fst_mat)]
} else if (is.na(fst_mat[2, 1])) {		# upper triangle matrix
	cat("Assuming upper trianglular Fst matrix input\n")
	fst_mat[lower.tri(fst_mat)] <- t(fst_mat)[lower.tri(fst_mat)]
} else {
	cat("Assuming square Fst matrix input\n")
}
diag(fst_mat) <- NA


# order the Fst matrix by pop if pops_present
if (pops_present) {
	myfsts <- fst_mat
	ord <- pops_table[, 1]
	myfsts <- myfsts[ord, ord]
} else {
	myfsts <- fst_mat
}


# set the bottom to NA and flip horizontally for plotting
myfsts[lower.tri(myfsts)] <- NA
myfsts <- myfsts[, seq(ncol(myfsts), 1)]


# start creating a pdf
pdf(paste0(out_pref, "_Fst.pdf"), width = 7, height = 7)

# plot the heatmap
par(mar = c(5, 4, 4, 5) + 0.1)
heatmapper(myfsts, palette = col_pal)

# stop creating the pdf
invisible(dev.off())


# similarly, a png
png(paste0(out_pref, "_Fst.png"), width = 20, height = 20, units = "cm", res = 600)
par(mar = c(5, 4, 4, 5) + 0.1)
heatmapper(myfsts, palette = col_pal)
invisible(dev.off())


# output the Fst matrix if it was sorted
if (pops_present) {
	write.table(myfsts, file = paste0(out_pref, "_sorted_Fst.txt"),
		quote = FALSE, row.names = TRUE)
}
