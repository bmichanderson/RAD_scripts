#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Oct 2021
# Description: summarize/visualize results from New Hybrids runs
##########################


import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



# set colours to be used
draw_colours = [
    'forestgreen',      # 22 8B 22 / 34 139 34
    'darkslateblue',    # 48 3D 8B / 72 61 139
    'darkorange',       # FF 8C 00 / 255 140 0
    'mediumorchid',     # BA 55 D3 / 186 85 211
    'lawngreen',        # 7C FC 00 / 124 252 0
	'lightskyblue',     # 87 CE FA / 135 206 250
]


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to summarize/visualize results from New Hybrids runs')


# add arguments to parse
parser.add_argument('-p', type = str, dest = 'prob_file', help = 'The probabilities file for class membership for a single run, typically \"aa-PofZ.txt\"')
parser.add_argument('-o', type = str, dest = 'output', help = 'Name of the output pdf without extension [default \"output\"]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
prob_file = args.prob_file
out_pre = args.output

if not prob_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# load the prob file as a dataframe and rename columns to the classes
prob_df = pd.read_csv(prob_file, sep = '\t')
prob_df.rename(columns = {	'1.000/0.000/0.000/0.000': 'Pure0',
							'0.000/0.000/0.000/1.000': 'Pure1',
							'0.000/0.500/0.500/0.000': 'F1',
							'0.250/0.250/0.250/0.250': 'F2',
							'0.500/0.250/0.250/0.000': 'Bx0',
							'0.000/0.250/0.250/0.500': 'Bx1'},
							inplace = True)


# create a bar chart of assignment probabilities to the six classes
fig, ax = plt.subplots()
ax = prob_df.plot.bar('IndivName', prob_df.columns[range(2, 8)], stacked = True,
					color = draw_colours, xlabel = '', ylabel = 'Probabilities',
					width = 1, edgecolor = 'black', figsize = (11, 8.5))
for spine in ax.spines:
    ax.spines[spine].set_visible(False)
ax.legend(loc = 'upper right', bbox_to_anchor = (1.15, 1), frameon = False)
ax.set_yticks(np.arange(0, 1.1, 0.1))
plt.savefig(out_pre + '.pdf', bbox_inches = 'tight')
