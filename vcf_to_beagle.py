#!/usr/bin/env python

##########################
# Author: B. Anderson, with inspiration for some formulas from https://github.com/darencard/RADpipe/blob/master/genotype_from_VCF.py
# Date: Oct 2021
# Description: convert a VCF file with PL fields into Beagle format
##########################


import sys
import os
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file with PL fields to a Beagle format')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to convert')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)


# process the VCF to grab sample labels and alleles and PL values for each SNP locus
with open(vcf_file, 'r') as vcf:
	phredls = []
	snp_names = []
	alleles = []
	count_snps = 0
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
		else:
			pieces = line.rstrip().split()
			snp_name = pieces[2]
			allele_ref = pieces[3]
			allele_alt = pieces[4]
			pls = []
			calls = pieces[9:]
			for call in calls:
				'''	The format expected is GT:DP:CATG:PL, e.g. 0/0:85:0,85,0,0:0,25.51,13.2
					We want to grab the Phred-scaled genotype likelihood field (PL) 0,25.51,13.2
					then split that into the three values as a list
				'''
				pls.append(call.split(':')[-1].split(','))
			phredls.append(pls)
			alleles.append((allele_ref, allele_alt))
			snp_names.append(snp_name)
			count_snps = count_snps + 1
	print('Read in a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNP loci')


# create the output lines to write to the Beagle file
header_line = ('Marker\tRef\tAlt\t' + ''.join(3 * (sample + '\t') for sample in sample_labels))
lines_out = []
for snp_index, phredl in enumerate(phredls):
	line_out = (str(snp_names[snp_index]) + '\t' + '\t'.join(alleles[snp_index]) + '\t')
	val_list = []
	for sample_index, sample_label in enumerate(sample_labels):
		this_pl = phredl[sample_index]
		if this_pl == ['.']:		# missing data
			# Anders Albrechtsen recommended using three numbers, not zeroes
			# I think we should keep it in line with the "normalized" below
			# This is a bit tricky, but we are saying any genotype is equally likely
			norm_00 = 1 / 3
			norm_01 = 1 / 3
			norm_11 = 1 / 3
		else:
			like_00 = float(10 ** (float(this_pl[0]) / -10))
			like_01 = float(10 ** (float(this_pl[1]) / -10))
			like_11 = float(10 ** (float(this_pl[2]) / -10))
			sum_likes = float(like_00 + like_01 + like_11)
			norm_00 = float(like_00 / sum_likes)
			norm_01 = float(like_01 / sum_likes)
			norm_11 = float(like_11 / sum_likes)
		norms = [norm_00, norm_01, norm_11]
		val_list.append('\t'.join('{0:.4f}'.format(item) for item in norms))
	line_out = line_out + '\t'.join(val_list)
	lines_out.append(line_out)


# create and write to the output file
with open(os.path.basename(vcf_file) + '.beagle', 'w') as outfile:
	outfile.write(header_line + '\n')
	for line_out in lines_out:
		outfile.write(line_out + '\n')
