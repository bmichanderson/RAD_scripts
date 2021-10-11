#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Oct 2021
# Description: convert a VCF file (VCF 4.0) to the New Hybrids input format
##########################


import sys			# allows access to command line arguments
import argparse
import random


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file to the input format required by New Hybrids;' +
												' it is up to the user to ensure the VCF has the SNPs of interest (e.g. already one per locus)')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to convert')
parser.add_argument('-s', type = str, dest = 'samples_file', help = 'A samples file of sample identifiers (same as in the VCF) that you want to keep, ' +
                                                                'one per line; if not provided, all samples will be kept')
parser.add_argument('-r', type = int, dest = 'random_down', help = 'Randomly downsample to this number of SNP loci')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
samples_file = args.samples_file
random_down = args.random_down


# process the samples file, if present
if samples_file:
	samples = []
	with open(samples_file, 'r') as input_file:
		for line in input_file:
			samples.append(line.rstrip())


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


# potentially filter samples
if samples_file:
	keep_indices = []
	new_snps = []
	new_samples = []
	for index, sample_label in enumerate(sample_labels):
		if sample_label in samples:
			keep_indices.append(index)
			new_samples.append(sample_label)
	for genotypes in snps:
		new_snps.append([genotypes[i] for i in keep_indices])
	snps = new_snps
	print('Retained ' + str(len(snps[0])) + ' samples after dropping those not in the samples file')


# only retain non-monomorphic SNP loci
keep_indices = []
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
	index = index + 1
new_snps = [snps[i] for i in keep_indices]
snps = new_snps
print('Retained ' + str(len(snps)) + ' loci after filtering for monomorphic')


# randomly downsample the SNP loci to the desired number
if random_down:
	if len(snps) > random_down:
		new_snps = random.sample(snps, random_down)
		snps = new_snps
		print('Randomly downsampled to ' + str(random_down) + ' SNP loci')


# create the output lines as ordered in the samples file if present
lines_out = []
write_index = 1
if samples_file:
	for sample in samples:
		if sample in new_samples:
			this_index = new_samples.index(sample)
			gts = []
			for genotypes in snps:
				gt = genotypes[this_index]
				if '.' in gt:	# missing data
					gts.append('0 ')
				else:
					gts.append(''.join(str(int(i) + 1) for i in gt))
			line_out = str(write_index) + ' n  ' + sample + '\t' + ' '.join(gts)
			lines_out.append(line_out)
			write_index = write_index + 1
		else:
			print('Sample ' + sample + ' was not present in the VCF but was present in your samples file')
else:
	for index, sample_label in enumerate(sample_labels):
		gts = []
		for genotypes in snps:
			gt = genotypes[index]
			if '.' in gt:	# missing data
				gts.append('0 ')
			else:
				gts.append(''.join(str(int(i) + 1) for i in gt))
		line_out = str(write_index) + ' n  ' + sample_label + '\t' + ' '.join(gts)
		lines_out.append(line_out)
		write_index = write_index + 1


# create the New Hybrids file
with open(vcf_file + '.newhybs', 'w') as outfile:
	outfile.write('NumIndivs ' + str(len(lines_out)) + '\n')
	outfile.write('NumLoci ' + str(len(snps)) + '\n')
	outfile.write('Digits 1\nFormat Lumped\n\n')
	for line_out in lines_out:
		outfile.write(line_out + '\n')


# report completion
print('Converted a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNP loci to New Hybrids format with ' +
	str(len(lines_out)) + ' samples and ' + str(len(snps)) + ' SNP loci')
