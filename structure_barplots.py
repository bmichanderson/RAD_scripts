#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Nov 2021
# Modified: Mar 2022
# Description: create barplots from Q matrices from Structure-like runs and output from CLUMPAK
##########################


import sys
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# set colours to be used (these can be adjusted)
draw_colours = [
    '#228b22', 	# forestgreen
    '#483d8b', 	# darkslateblue
    '#87cefa',	# lightskyblue
    '#ba55d3',	# mediumorchid
	'#7fffd4',	# aquamarine
	'#ff8c00',	# darkorange
	'#cd853f',	# peru
	'#ff6347',	# tomato
	'#f0e68c',	# khaki
	'#ff1493',	# deeppink
	'#708090',	# slategrey
    '#7cfc00',	# lawngreen
	'#4682b4',	# steelblue
]


# set colours to greyscale for neutral plotting (could comment out)
draw_colours = [
	'1',
	'0',
	'0.5',
	'0.75',
	'0.25',
	'0.88',
	'0.62',
	'0.38',
	'0.12',
]


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to create barplots for Structure-like Q matrices')


# add arguments to parse
parser.add_argument('-o', type = str, dest = 'out_pre', help = 'The prefix for the output pdf [default \"output\"')
parser.add_argument('-p', type = str, dest = 'pops_file', help = 'The populations file in the tab-delimited form ' +
					'"sampleID    pop", one per line; samples must be in the same order as in the Q matrix')
parser.add_argument('-q', type = str, dest = 'Q_file', help = 'The Q matrix file containing the whitespace-delimited assignment proportions ' +
					'without headers or any line information; i.e. a table with K columns and as many rows as samples')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
out_pre = args.out_pre
pops_file = args.pops_file
Q_file = args.Q_file

if any([not pops_file, not Q_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)
if not out_pre:
	out_pre = 'output'


# load the pops file, keeping the order
sample_df = pd.read_csv(pops_file, sep = '\t', header = None)
sample_df.rename(columns = {0: 'Sample', 1: 'Pop'},	inplace = True)


# load the Q matrix
Q_df = pd.read_csv(Q_file, delim_whitespace = True, header = None)


# add the Q values to the sample dataframe
num_apops = len(Q_df.columns)
for K in range(num_apops):
	sample_df.loc[:, 'Ancestral_pop' + str(K + 1)] = Q_df[K]


# sort by pop
plot_df = sample_df.sort_values(['Pop', 'Sample'])


# NOTE: The next K sorting bit is unnecessary if using output from CLUMPAK
# Otherwise, you may have to manually order columns in the q files for consistency

# sort the K columns to consistently colour
# while the columns aren't sorted yet
# for pop in pops
#	which column has the highest value for individual 1
#		select that column as next in the sort order if not already there
# NOTE: it may still not be sorted (pop not majority in any first ind)
# if so, compare column max values, and continue until all assigned
#K = 1
#sort_order = []
#while K <= num_apops:
#	for pop in set(plot_df['Pop']):
#		df1 = plot_df[plot_df['Pop'] == pop]
#		apop = df1[df1.columns[2:]].idxmax(axis = 1).iloc[0]
#		if apop not in sort_order:
#			sort_order.append(apop)
#			K = K + 1
#	maxs = plot_df.iloc[:, 2:].max().sort_values(ascending = False)
#	for apop in list(maxs.index):
#		if apop not in sort_order:
#			sort_order.append(apop)
#			K = K + 1
#plot_df = plot_df[['Sample', 'Pop'] + sort_order]


# create the barplot
fig, ax = plt.subplots()
ax = plot_df.plot.bar('Sample', plot_df.columns[range(2, 2 + num_apops)], stacked = True,
					color = draw_colours[: num_apops], xlabel = '', ylabel = '',
					width = 1, edgecolor = 'black', figsize = (30, 5),
					legend = None)
for spine in ax.spines:
    ax.spines[spine].set_visible(False)
plt.tick_params(axis = 'y', left = False, labelleft = False)
plt.tick_params(axis = 'x', bottom = False)
plt.savefig(out_pre + '.svg', bbox_inches = 'tight', format = 'svg')
#plt.savefig(out_pre + '.pdf', bbox_inches = 'tight')
