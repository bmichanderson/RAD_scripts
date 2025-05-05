#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: Oct 2021
# Modified: April 2022; May 2025 (cleaned up; added optional samples file input for converting
# 	tip labels and population table for removing autapomorphic SNPs for SNAPP/ER)
# Description: convert a VCF file (VCF 4.0) to a Nexus file for SplitsTree or SNAPP/ER
##########################


import sys
import argparse
import os
import random


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file to a Nexus input format required by SplitsTree or SNAPP/ER;' +
	' it is up to the user to ensure the VCF has already been filtered to the SNPs of interest (e.g. biallelic, one per locus)')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to convert')
parser.add_argument('-o', type = str, dest = 'out_format', help = 'The type of output desired: \"snapp\" or \"splits\" [default]')
parser.add_argument('-r', type = int, dest = 'random_down', help = 'Randomly downsample to this number of SNP loci')
parser.add_argument('-s', type = str, dest = 'sample_table', help = 'File with sample labels and desired display labels, tab separated and one per line')
parser.add_argument('-a', type = str, dest = 'autapo', help = 'File with sample labels and population labels, tab separated and one per line;' +
	' if provided, this will turn on filtering of autapomorphic SNPs (only variable in single populations)')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
out_format = args.out_format
random_down = args.random_down
sample_table = args.sample_table
autapo = args.autapo

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if out_format:
	if all([out_format.lower() != 'splits', out_format.lower() != 'snapp']):
		parser.print_help(sys.stderr)
		sys.exit(1)
else:
	out_format = 'splits'


# process the samples table if present
if sample_table:
	sample_dict = {}
	with open(sample_table, 'r') as sample_file:
		for line in sample_file:
			pieces = line.strip().split('\t')
			sample_dict[pieces[0]] = pieces[1]


# process the pops table if present
if autapo:
	pops = []
	pop_dict = {}
	with open(autapo, 'r') as pop_file:
		for line in pop_file:
			pieces = line.strip().split('\t')
			sample = pieces[0]
			pop = pieces[1]
			if pop not in pops:
				pops.append(pop)
				pop_dict[pop] = [sample]
			else:
				pop_dict[pop].append(sample)


# process the VCF to grab sample labels and genotypes for each SNP locus
with open(vcf_file, 'r') as vcf:
	snps = []
	count_snps = 0
	for line in vcf:
		if line.startswith('#'):		# a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
				if autapo:
					pop_labels = []
					for sample_label in sample_labels:
						for pop in pops:
							if sample_label in pop_dict[pop]:
								pop_labels.append(pop)
		else:
			calls = line.rstrip().split()[9:]
			genotypes = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) 0/0 then split that into the two alleles 0 and 0 as a list
					Note: some genotypes may be phased ('|'), so account for that
				'''
				genotypes.append(call.split(':')[0].replace('|', '/').split('/'))
			snps.append(genotypes)
			count_snps = count_snps + 1
	print('Read in a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNP loci')


# if there was a samples table, translate the labels
if sample_table:
	for index, sample in enumerate(sample_labels):
		sample_labels[index] = sample_dict.get(sample, sample)		# don't replace if the sample isn't in the dictionary


# remove monomorphic SNP loci
keep_indices = []
count_mono = 0
index = 0
for genotypes in snps:
	gts = []
	mono = True
	for gt in genotypes:
		if gt == ['.', '.']:		# ignore missing
			continue
		else:
			if gt not in gts:
				gts.append(gt)
		if len(gts) > 1:
			mono = False
			break
	if not mono:
		keep_indices.append(index)
	else:
		count_mono = count_mono + 1
	index = index + 1
if count_mono > 0:
	snps = [snps[i] for i in keep_indices]
	print('Dropped ' + str(count_mono) + ' monomorphic SNPs and retained ' + str(len(snps)))


# if there was a populations table provided, drop autapomorphic SNPs
if autapo:
	count_autapo = 0
	keep_indices = []
	index = 0
	for genotypes in snps:
		comparisons = []
		autapomorphic = True
		for pop in set(pop_labels):
			pop_genos = [''.join(sorted(genotypes[ind])) for ind, item in enumerate(pop_labels) if item == pop]
			pop_genos = [item for item in pop_genos if item != '..']
			for item in set(pop_genos):
				comparisons.append(item)
			genos = [item for item in set(comparisons)]
			if len(set(comparisons)) > 2:
				autapomorphic = False
				break
			elif len(set(comparisons)) == 2:
				if all([len([item for item in comparisons if item == genos[0]]) > 1,
					len([item for item in comparisons if item == genos[1]]) > 1]):
					autapomorphic = False
					break
				else:
					continue
			else:
				continue
		if not autapomorphic:
			keep_indices.append(index)
		else:
			count_autapo = count_autapo + 1
		index = index + 1
	if count_autapo > 0:
		snps = [snps[i] for i in keep_indices]
		print('Dropped ' + str(count_autapo) + ' autapomorphic SNPs and retained ' + str(len(snps)))


# randomly downsample the SNP loci to the desired number
if random_down:
	if len(snps) > random_down:
		snps = random.sample(snps, random_down)
		print('Randomly downsampled to ' + str(random_down) + ' SNPs')


# create the output character lines to write to the Nexus file
lines_out = []
for index, sample_label in enumerate(sample_labels):
	gts = []
	for genotypes in snps:
		gt = genotypes[index]
		if out_format == 'snapp':			# 0 homo, 1 het, 2 homo alt
			if '.' in gt:	# missing data
				gts.append('?')
			else:
				gts.append(str(int(gt[0]) + int(gt[1])))
		elif out_format == 'splits':		# each allele as a binary character
			if '.' in gt:	# missing data
				gts.append('??')
			else:
				gts.append(str(gt[0]))
				gts.append(str(gt[1]))
	line_out = ''.join(gts)
	lines_out.append(line_out)


# set output parameters
if out_format == 'snapp':
	datatype = 'SNP'
	missing = '?'
	symbols = '012'
	num_chars = len(snps)
	data_block = ('BEGIN DATA;\n\tDIMENSIONS NTAX=' + str(len(sample_labels)) +
		' NCHAR=' + str(num_chars) + ';\n\t' +
		'FORMAT\n\t\tDATATYPE=' + datatype + '\n\t\tMISSING=' + missing + '\n\t\t' +
		'SYMBOLS=\"' + symbols + '\";\n\tMATRIX\n')
elif out_format == 'splits':
	datatype = 'STANDARD'
	missing = '?'
	symbols = '01'
	num_chars = len(snps) * 2
	data_block = ('BEGIN TAXA;\n\tDIMENSIONS NTAX=' + str(len(sample_labels)) + ';\n\t' +
		'TAXLABELS ' + ' '.join(sample_labels) + ';\nEND;\n' +
		'BEGIN CHARACTERS;\n\tDIMENSIONS NCHAR=' + str(num_chars) + ';\n\t' +
		'FORMAT\n\t\tDATATYPE=' + datatype + '\n\t\tMISSING=' + missing + '\n\t\t' +
		'LABELS=NO\n\t\tSYMBOLS=\"' + symbols + '\";\n\tMATRIX\n')

# create the Nexus file
with open(os.path.splitext(vcf_file)[0] + '_' + out_format + '.nex', 'w') as outfile:
	outfile.write('#NEXUS\n')
	outfile.write(data_block)
	if out_format == 'snapp':
		for index, line_out in enumerate(lines_out):
			outfile.write('\t\t' + sample_labels[index] + '\t' + line_out + '\n')
	elif out_format == 'splits':
		for line_out in lines_out:
			outfile.write('\t\t' + line_out + '\n')
	outfile.write('\t;\nEND;\n')


# report completion
print('Converted a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNPs to Nexus format with ' +
	str(len(lines_out)) + ' samples and ' + str(num_chars) + ' characters')
