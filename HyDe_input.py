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
				+ ' SNPs are assumed to be unlinked')
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


# process the VCF to extract a single base representation for each SNP per individual
with open(vcf_file, 'r') as vcf:
	snps = []
	count_snps = 0
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
		else:
			parts = line.rstrip().split()
			all_bases = []
			all_bases.append(parts[3])	# the ref allele
			if ',' in parts[4]:		# a multi-allele
				for allele in parts[4].split(','):
					all_bases.append(allele)
			else:
				all_bases.append(parts[4])
			calls = parts[9:]
			bases = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					or for a multi-allele A T,C with genotype 1/2:50:30,0,20,0 for a TC with total 50 depth
					We want to grab the genotype (GT) 0/0 then split that into the two alleles 0 and 0
					or for a multi-allele 0/3 then split into 0 and 3, but have to index
				'''
				genotype = ''.join(call.split(':')[0].split('/'))
				allele1 = genotype[0]
				allele2 = genotype[1]
				if genotype == '..':		# missing data
					bases.append('?')
				elif allele1 == allele2:
					bases.append(all_bases[int(allele1)])
				else:
					bases.append(amb_dict[all_bases[int(allele1)] + all_bases[int(allele2)]])
			snps.append(bases)
			count_snps = count_snps + 1
	print('Read in a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNPs')


# create the output lines to write to the file in the order of the sorted samples
lines_out = []
samples_out = []
sample_excl = 0
for sample_tuple in sample_list:
	sample = sample_tuple[0]
	this_index = -9
	for index, sample_label in enumerate(sample_labels):
		if sample_label == sample:
			this_index = index
			break
	if this_index < 0:
		print('Sample ' + sample + ' is not present in the VCF')
		sample_excl = sample_excl + 1
	else:
		bases = []
		samples_out.append(sample_tuple)
		for snp in snps:
			bases.append(snp[this_index])
		line_out = sample + '\t' + ''.join(bases)
		lines_out.append(line_out)


# determine how many unique taxa:
num_taxa = len(set([item[1] for item in samples_out]))


# create the output files
with open(out_prefix + '_data.txt', 'w') as outfile:
	for line_out in lines_out:
		outfile.write(line_out + '\n')

with open(out_prefix + '_map.txt', 'w') as outfile:
	for sample_tuple in samples_out:
		outfile.write(sample_tuple[0] + '\t' + sample_tuple[1] + '\n')


# report completion
print('Wrote output: ' + str(len(samples_out)) + ' samples, ' +
	str(num_taxa) + ' taxa and ' + str(count_snps) + ' SNPs')
if sample_excl > 0:
	print('Did not write output for ' + str(sample_excl) + ' sample(s) in the samples file')
