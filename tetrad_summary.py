#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Sep 2021
# Modified: Oct 2021
# Description: summarize results from an ipyrad tetrad run using the tetrad.sif container
# NOTE: Run this script using the tetrad.sif container (e.g. singularity exec tetrad.sif python tetrad_summary.py ...)
##########################


import sys
import argparse
import ipyrad.analysis as ipa
import toytree
import toyplot
from toyplot import pdf
import itertools


# set colours to be used (needs to be as long as the number of groups or it will re-use)
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
parser = argparse.ArgumentParser(description = 'A script to summarize results from tetrad runs via the tetrad.sif container')


# add arguments to parse
parser.add_argument('-s', type = str, dest = 'snps_file', help = 'The *.snps.hdf5 file containing the data used for the runs')
parser.add_argument('-n', type = str, dest = 'name', help = 'The name that was used for the run (prefix of the existing and to be created output files)')
parser.add_argument('-p', type = str, dest = 'pops_file', help = 'The populations file in the tab-separated form "sampleID    pop", one per line')
parser.add_argument('-w', type = str, dest = 'workdir', help = 'The working directory where the results are')
parser.add_argument('-o', type = str, dest = 'outgroup', help = 'The string present in outgroup samples IDs')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
snps_file = args.snps_file
name = args.name
pops_file = args.pops_file
workdir = args.workdir
outgroup = args.outgroup

if any([not name, not pops_file, not workdir, not snps_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)
if not outgroup:
	outgroup = 'Z'


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


# create a colours dictionary
colour_cycler = itertools.cycle(draw_colours)
pop_colours = {}
for pop in pops:
	pop_colours[pop] = next(colour_cycler)
sample_colours = {}
for sample in samples:
	for pop, list in pop_dict.items():
		if sample in list:
			sample_colours[sample] = pop_colours[pop]
			break


# load the tetrad object
tet = ipa.tetrad(
	name = name,
	workdir = workdir,
	load_only = True,
	data = snps_file,
	imap = pop_dict,
)


## plot trees

# main tree
tre = toytree.tree(tet.trees.tree).root([i for i in samples if outgroup in i])
node_sizes = []
support_vals = tre.get_node_values('support')
for index in range(tre.nnodes):
	if support_vals[index]:
		node_sizes.append(25)
	else:
		node_sizes.append(0)
tip_colours = []
tip_labels = tre.get_tip_labels()
for index in range(tre.ntips):
	tip_colours.append(sample_colours[tip_labels[index]])
tre_plot = tre.draw(
	height = 3000,
	width = 1500,
	tip_labels_style = {'font-size': '14px'},
	tip_labels_colors = tip_colours,
	node_labels = 'support',
	node_labels_style = {'font-size': '12px'},
	node_sizes = node_sizes,
	node_style = {'fill': 'lightgray'},
	use_edge_lengths = False,
)
toyplot.pdf.render(tre_plot[0], fobj = name + '_tree.pdf')

# bootstrap majority rule consensus tree
tre = toytree.tree(tet.trees.cons).root([i for i in samples if outgroup in i])
node_sizes = []
support_vals = tre.get_node_values('support')
for index in range(tre.nnodes):
	if support_vals[index]:
		node_sizes.append(25)
	else:
		node_sizes.append(0)
tip_colours = []
tip_labels = tre.get_tip_labels()
for index in range(tre.ntips):
	tip_colours.append(sample_colours[tip_labels[index]])
tre_plot = tre.draw(
	height = 3000,
	width = 1500,
	tip_labels_style = {'font-size': '14px'},
	tip_labels_colors = tip_colours,
	node_labels = 'support',
	node_labels_style = {'font-size': '12px'},
	node_sizes = node_sizes,
	node_style = {'fill': 'lightgray'},
	use_edge_lengths = False,
)
toyplot.pdf.render(tre_plot[0], fobj = name + '_constree.pdf')

# cloud tree
mtre = toytree.mtree(tet.trees.boots)
mtre.treelist = [i.root([i for i in samples if outgroup in i]) for i in mtre.treelist]
mtre_plot = mtre.draw_cloud_tree(
	height = 3000,
	width = 1500,
	tip_labels_style = {'font-size': '14px'},
	use_edge_lengths = False,
)
toyplot.pdf.render(mtre_plot[0], fobj = name + '_cloudtree.pdf')
