#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Oct 2021
# Description: convert a VCF file (VCF 4.0) to a Nexus file, for SplitsTree or SNAPP/ER
##########################


import sys
import argparse
import random

# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file to a Nexus input format required by SplitsTree or SNAPP;' +
						' it is up to the user to ensure the VCF has the SNPs of interest (e.g. biallelic, one per locus)')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to convert')
parser.add_argument('-o', type = str, dest = 'out_format', help = 'The type of output desired: \"snapp\" or \"splits\" [default]')
parser.add_argument('-r', type = int, dest = 'random_down', help = 'Randomly downsample to this number of SNP loci')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
out_format = args.out_format
random_down = args.random_down

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if out_format:
	if all([out_format.lower() != 'splits', out_format.lower() != 'snapp']):
		parser.print_help(sys.stderr)
		sys.exit(1)
else:
	out_format = 'splits'


# process the VCF to grab sample labels and genotypes for each SNP locus
with open(vcf_file, 'r') as vcf:
	snps = []
	count_snps = 0
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
		else:
			calls = line.rstrip().split()[9:]
			genotypes = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) 0/0 then split that into the two alleles 0 and 0 as a list
				'''
				genotypes.append(call.split(':')[0].split('/'))
			snps.append(genotypes)
			count_snps = count_snps + 1
	print('Read in a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNP loci')


# only retain non-monomorphic SNP loci
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
	new_snps = [snps[i] for i in keep_indices]
	snps = new_snps
	print('Dropped ' + str(count_mono) + ' monomorphic SNPs and retained ' + str(len(snps)) + ' loci after filtering for monomorphic')


# randomly downsample the SNP loci to the desired number
if random_down:
	if len(snps) > random_down:
		new_snps = random.sample(snps, random_down)
		snps = new_snps
		print('Randomly downsampled to ' + str(random_down) + ' SNP loci')


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
	datatype = 'INTEGERDATA'
	missing = '?'
	symbols = '012'
	num_chars = len(snps)
elif out_format == 'splits':
	datatype = 'STANDARD'
	missing = '?'
	symbols = '01'
	num_chars = len(snps) * 2


# create the Nexus file
with open(vcf_file + '.nex', 'w') as outfile:
	outfile.write('#NEXUS\n')
	taxa_block = ('BEGIN TAXA;\n\tDIMENSIONS NTAX=' + str(len(sample_labels)) + ';\n\t' +
				'TAXLABELS ' + ' '.join(sample_labels) + ';\nEND;\n')
	char_block = ('BEGIN CHARACTERS;\n\tDIMENSIONS NCHAR=' + str(num_chars) + ';\n\t' +
				'FORMAT\n\t\tDATATYPE=' + datatype + '\n\t\tMISSING=' + missing + '\n\t\t' +
				'LABELS=NO\n\t\tSYMBOLS=\"' + symbols + '\";\n\tMATRIX\n')
	outfile.write(taxa_block)
	outfile.write(char_block)
	for line_out in lines_out:
		outfile.write('\t\t\t' + line_out + '\n')
	outfile.write('\t;\nEND;\n')


# report completion
print('Converted a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNP loci to Nexus format with ' +
	str(len(lines_out)) + ' samples and ' + str(num_chars) + ' characters')
