#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: Dec 2023
# Modified: Feb 2026 (adjusted threshold options and screen output format)
# Description: report a clone threshold for a set of sample comparisons (from vcf_similarity.py),
#	based on replicate pairs in a text file; also report what sample comparisons (not reps) are below it
##########################


import sys
import argparse
import numpy
import pandas


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to report clonality from a comparison table and replicate pairs')


# add arguments to parse
parser.add_argument(type = str, dest = 'comp_file', help = 'The comparisons file with three columns (no header): ' + \
	'individual 1 (string), individual 2 (string), distance (float)')
parser.add_argument('-r', type = str, dest = 'reps', help = 'A text file of replicate pairs, tab separated, one pair per line')
parser.add_argument('-t', type = str, dest = 'thresh', help = 'Threshold type to use: "d" double the average (default), ' + \
	'"s" two standard deviations above the mean, or "q" 1.5 * IQR above Q3')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
comp_file = args.comp_file
rep_file = args.reps
thresh = args.thresh


if any([not comp_file, not rep_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)

if not thresh:
	thresh = 'd'


# read in the comparisons file
complist = []
with open(comp_file, 'r') as infile:
	for line in infile:
		parts = line.strip().split()
		complist.append(parts)


# read in the replicate pair combinations
replist = []
replicates = []
with open(rep_file, 'r') as infile:
	for line in infile:
		items = line.strip().split()
		replist.append(items)
		replicates.append(items[0])
		replicates.append(items[1])


# assess distances for replicate pairs
repdists = []
for reppair in replist:
	for comp in complist:
		if all([reppair[0] in comp, reppair[1] in comp]):
			repdists.append(comp[2])
			break

dist_vals = [float(item) for item in repdists]
mydf = pandas.Series(dist_vals)

first_quartile = mydf.quantile(0.25)
third_quartile = mydf.quantile(0.75)
interquartile = third_quartile - first_quartile
outlier_threshold = third_quartile + (1.5 * interquartile)

avg = numpy.mean(dist_vals)
maxrep = max(dist_vals)
minrep = min(dist_vals)
stdev = numpy.std(dist_vals)
avg_threshold = 2 * avg
stdev_threshold = avg + (2 * stdev)


print('Average distance between rep pairs: ' + '{:.4f}'.format(avg))
print('Minimum distance between rep pairs: ' + '{:.4f}'.format(minrep))
print('Maximum distance between rep pairs: ' + '{:.4f}'.format(maxrep))
print('')
print('Threshold based on double the average: ' + '{:.4f}'.format(avg_threshold))
print('Threshold based on two standard deviations above the mean: ' + '{:.4f}'.format(stdev_threshold))
print('Threshold based on interquartile range outliers: ' + '{:.4f}'.format(outlier_threshold))


if thresh == 'd':
	threshold = avg_threshold
elif thresh == 's':
	threshold = stdev_threshold
elif thresh == 'q':
	threshold = outlier_threshold
else:
	print('Unrecognised option for threshold argument -t\n')
	parser.print_help(sys.stderr)
	sys.exit(1)


# report which comparisons are below the threshold (putative clones)
highcomps = []
lowcomps = []
rephighcomps = []
replowcomps = []
for comp in complist:
	if any([[comp[0], comp[1]] in replist, [comp[1], comp[0]] in replist]): 	# a rep pair
		if float(comp[2]) > threshold:	# a high comparison
			rephighcomps.append(comp)
	elif all([comp[0] in replicates, comp[1] in replicates]):	# involves only replicates (not paired in list)
		if float(comp[2]) < threshold:	# a low comparison
			replowcomps.append(comp)
	else:	# any combo
		if float(comp[2]) < threshold:	# a low comparison
			lowcomps.append(comp)
		else:	# a high comparison
			highcomps.append(comp)

if len(rephighcomps) > 0:
	print('')
	print('Replicate combos above threshold:')
	for item in sorted(rephighcomps, key = lambda x: x[2]):
		print('\t'.join([item[0], item[1], '{:.8f}'.format(float(item[2]))]))

if len(replowcomps) > 0:
	print('')
	print('Unlisted replicate combos below threshold:')
	for item in sorted(replowcomps, key = lambda x: x[2]):
		print('\t'.join([item[0], item[1], '{:.8f}'.format(float(item[2]))]))

if len(lowcomps) > 0:
	print('')
	print('Combos below threshold (possible clones):')
	for item in sorted(lowcomps, key = lambda x: x[2]):
		print('\t'.join([item[0], item[1], '{:.8f}'.format(float(item[2]))]))

maxhigh = max([float(item[2]) for item in highcomps])
minhigh = min([float(item[2]) for item in highcomps])
lowfive = sorted(highcomps, key = lambda x: x[2])[0:5]

print('')
print('Maximum distance between samples including non-reps: ' + '{:.4f}'.format(maxhigh))
print('Minimum distance above threshold including non-reps: ' + '{:.4f}'.format(minhigh))
print('')
print('Lowest five comparisons above threshold including non-reps:')
for item in lowfive:
	print('\t'.join([item[0], item[1], '{:.8f}'.format(float(item[2]))]))
