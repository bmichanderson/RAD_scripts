#!/usr/bin/env python3

##########################
# Author: B. Anderson
# Date: Dec 2023
# Description: report a clone threshold for a set of sample comparisons (from vcf_similarity.py),
#	based on replicate pairs in a text file; also report what sample comparisons (not reps) are below it
##########################


import sys
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to report clonality from a comparison table and replicate pairs')


# add arguments to parse
parser.add_argument(type = str, dest = 'comp_file', help = 'The comparisons file with three columns (no header): ' + \
	'individual 1 (string), individual 2 (string), distance (float)')
parser.add_argument('-r', type = str, dest = 'reps', help = 'A text file of replicate pairs, tab separated, one pair per line')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
comp_file = args.comp_file
rep_file = args.reps

if any([not comp_file, not rep_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)


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

avg = sum([float(item) for item in repdists]) / len(repdists)
maxrep = max([float(item) for item in repdists])

print('Maximum distance between rep pairs: ' + '{:.4f}'.format(maxrep))
print('Average distance between rep pairs: ' + '{:.4f}'.format(avg))
print('Threshold: ' + '{:.4f}'.format(2 * avg))


# report which comparisons are below 2*avg (putative clones)
highcomps = []
lowcomps = []
rephighcomps = []
replowcomps = []
for comp in complist:
	if any([[comp[0], comp[1]] in replist, [comp[1], comp[0]] in replist]): 	# a rep pair
		if float(comp[2]) > (2 * avg):	# a high comparison
			rephighcomps.append(comp)
	elif all([comp[0] in replicates, comp[1] in replicates]):	# involves only replicates (not paired in list)
		if float(comp[2]) < (2 * avg):	# a low comparison
			replowcomps.append(comp)
	else:	# any combo
		if float(comp[2]) < (2 * avg):	# a low comparison
			lowcomps.append(comp)
		else:	# a high comparison
			highcomps.append(comp)

if len(rephighcomps) > 0:
	print('')
	print('Replicate combos above threshold:')
	for item in sorted(rephighcomps, key = lambda x: x[2]):
		print(item)

if len(replowcomps) > 0:
	print('')
	print('Unlisted replicate combos below threshold:')
	for item in sorted(replowcomps, key = lambda x: x[2]):
		print(item)

if len(lowcomps) > 0:
	print('')
	print('Combos below threshold (possible clones):')
	for item in sorted(lowcomps, key = lambda x: x[2]):
		print(item)

maxhigh = max([float(item[2]) for item in highcomps])
minhigh = min([float(item[2]) for item in highcomps])
lowfive = sorted(highcomps, key = lambda x: x[2])[0:5]

print('')
print('Maximum distance between samples including non-reps: ' + '{:.4f}'.format(maxhigh))
print('Minimum distance above threshold including non-reps: ' + '{:.4f}'.format(minhigh))
print('Lowest five comparisons above threshold including non-reps:')
for item in lowfive:
	print(item)
