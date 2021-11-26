#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Nov 2021
# Description: create input files for HyDe scripts for hybrid detection
##########################


import sys
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to create input files for HyDe')


# add arguments to parse
parser.add_argument('-v', type = str, dest = 'vcf_file',  help = 'The VCF file to obtain sequences from;'
				+ ' SNPs are assumed to be unlinked and biallelic')
parser.add_argument('-s', type = str, dest = 'sample_file', help = 'The tab-delimited file of sample name and pop/taxon, one per line')
parser.add_argument('-o', type = str, dest = 'out_prefix', help = 'The output prefix [default \"out\"]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
sample_file = args.sample_file
out_prefix = args.out_prefix

if not all([vcf_file, sample_file]):
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_prefix:
	out_prefix = 'out'


# create an ambiguity dictionary
amb_dict = {
	'CT': 'Y',
	'TC': 'Y',
	'AG': 'R',
	'GA': 'R',
	'AT': 'W',
	'TA': 'W',
	'CG': 'S',
	'GC': 'S',
	'GT': 'K',
	'TG': 'K',
	'AC': 'M',
	'CA': 'M'
}

# process the sample file into a tuple list (sample, pop/taxon)
sample_list = []
with open(sample_file, 'r') as infile:
	for line in infile:
		parts = line.strip().split()
		sample_list.append((parts[0], parts[1]))
# sort it by pop/taxon
sample_list.sort(key = lambda x: x[1])


# process the VCF to extract a single base representation for each SNP
with open(vcf_file, 'r') as vcf:
	snps = []
	count_snps = 0
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
		else:
			parts = line.rstrip().split()
			base1 = parts[3]
			base2 = parts[4]
			calls = parts[9:]
			bases = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) 0/0 then split that into the two alleles 0 and 0
				'''
				genotype = ''.join(call.split(':')[0].split('/'))
				if genotype == '00':
					bases.append(base1)
				elif genotype == '11':
					bases.append(base2)
				elif any([genotype == '10', genotype == '01']):
					bases.append(amb_dict[base1 + base2])
				elif genotype == '..':		# missing data
					bases.append('?')
				else:
					print('Encountered a genotype that does not fit expected pattern!')
			snps.append(bases)
			count_snps = count_snps + 1
	print('Read in a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNPs')


# create the output lines to write to the file in the order of the sorted samples
lines_out = []
for sample_tuple in sample_list:
	sample = sample_tuple[0]
	this_index = -9
	for index, sample_label in enumerate(sample_labels):
		if sample_label == sample:
			this_index = index
			break
	if this_index < 0:
		print('Sample not present in VCF! Exiting...')
		sys.exit(1)
	bases = []
	for snp in snps:
		bases.append(snp[this_index])
	line_out = sample + '\t' + ''.join(bases)
	lines_out.append(line_out)


# create the output files
with open(out_prefix + '_data.txt', 'w') as outfile:
	for line_out in lines_out:
		outfile.write(line_out + '\n')

with open(out_prefix + '_map.txt', 'w') as outfile:
	for sample_tuple in sample_list:
		outfile.write(sample_tuple[0] + '\t' + sample_tuple[1] + '\n')


# report completion
print('Wrote output for ' + str(len(sample_list)) + ' samples')
