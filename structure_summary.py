#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Sep 2021
# Description: summarize results from an ipyrad Structure run using the ipyrad.sif container
# NOTE: Run this script using the ipyrad.sif container (e.g. singularity exec ipyrad.sif python structure_summary.py ...)
##########################


import sys
import argparse
import itertools
import ipyrad.analysis as ipa
import toyplot
from toyplot import pdf


# set colours to be used (needs to be as long as the range of K values that need to be plotted or it will re-use)
draw_colours = [
    'forestgreen',      # 22 8B 22 / 34 139 34
    'darkslateblue',    # 48 3D 8B / 72 61 139
    'lightskyblue',     # 87 CE FA / 135 206 250
    'darkorange',       # FF 8C 00 / 255 140 0
    'mediumorchid',     # BA 55 D3 / 186 85 211
    'lawngreen',        # 7C FC 00 / 124 252 0
	'steelblue',		# 46 82 B4 / 70	130	180
	'aquamarine',		# 7F FF D4 / 127 255 212
	'olive',			# 80 80 00 / 128 128 0
	'sandybrown',		# F4 A4 60 / 244 164 96
]


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to summarize results from Structure runs via the ipyrad.sif container')


# add arguments to parse
parser.add_argument('-s', type = str, dest = 'snps_file', help = 'The *.snps.hdf5 file containing the data used for the runs')
parser.add_argument('-n', type = str, dest = 'name', help = 'The name that was used for the run (prefix of the existing and to be created output files)')
parser.add_argument('-p', type = str, dest = 'pops_file', help = 'The populations file in the tab-separated form "sampleID    pop", one per line')
parser.add_argument('-k', type = str, dest = 'k_vals', help = 'The min and max K values run, space-delimited in a single set of quotes " "')
parser.add_argument('-w', type = str, dest = 'workdir', help = 'The working directory where the results are')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
snps_file = args.snps_file
name = args.name
pops_file = args.pops_file
k_vals = [int(i) for i in args.k_vals.split()]
workdir = args.workdir

if any([not name, not pops_file, not k_vals, not workdir, not snps_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)


# load the pops file and create the pop dictionary
samples = []
pops = []
pop_dict = {}
with open(pops_file, 'r') as sample_input:
	for line in sample_input:
		sampleID = line.rstrip().split()[0]
		pop = line.rstrip().split()[1]
		samples.append(sampleID)
		if pop not in pops:
			pops.append(pop)
			pop_dict[pop] = [sampleID]
		else:
			pop_dict[pop].append(sampleID)


# load the structure object
struct = ipa.structure(
	name = name,
	workdir = workdir,
	load_only = True,
	data = snps_file,
	imap = pop_dict,
)


# generate the Evanno table and plot
etable = struct.get_evanno_table(range(k_vals[0], k_vals[-1] + 1))
canvas = toyplot.Canvas(width = 400, height = 300)
# log probability
axes = canvas.cartesian(ylabel = 'Estimated Ln(Prob) Mean')
axes.plot(etable.estLnProbMean, color = draw_colours[0], marker = 'o')
axes.scatterplot(etable.estLnProbMean + etable.estLnProbStdev, color = 'black', marker = 'x')
axes.scatterplot(etable.estLnProbMean - etable.estLnProbStdev, color = 'black', marker = 'x')
axes.y.spine.style = {'stroke': draw_colours[0]}
# delta K
axes = axes.share('x', ylabel = 'deltaK', ymax = etable.deltaK.max() + etable.deltaK.max() * 0.1)
axes.plot(etable.deltaK, color = draw_colours[1], marker = 'o')
axes.y.spine.style = {'stroke': draw_colours[1]}
axes.x.ticks.locator = toyplot.locator.Explicit(range(len(etable.index)), etable.index)
axes.x.label.text = 'K (N ancestral populations)'
toyplot.pdf.render(canvas, fobj = name + '_evanno.pdf')


# make Clumpp figures
colour_cycler = itertools.cycle(draw_colours)
colours = [next(colour_cycler)]

for k_val in range(k_vals[0], k_vals[-1] + 1):
	table = struct.get_clumpp_table(k_val, max_var_multiple = 10)
	onames = list(itertools.chain(*pop_dict.values()))
	table = table.loc[onames]
	canvas = toyplot.Canvas(width = 4000, height = 2000)
	axes = canvas.cartesian(bounds = ('10%', '90%', '10%', '50%'))
	colours.append(next(colour_cycler))
	axes.bars(table, color = colours)
	ticklabels = [i for i in table.index.tolist()]
	axes.x.ticks.locator = toyplot.locator.Explicit(labels = ticklabels)
	axes.x.ticks.labels.angle = -60
	axes.x.ticks.show = True
	axes.x.ticks.labels.offset = 20
	axes.x.ticks.labels.style = {'font-size': '14px'}
	toyplot.pdf.render(canvas, fobj = name + '_K' + str(k_val) + '.pdf')

