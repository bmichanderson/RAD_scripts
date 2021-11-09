#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Oct 2021
# Description: convert a VCF file without PL fields into one with (for Tiger to update with values)
##########################


import sys
import os
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file without PL fields to one with')


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


# process the VCF
with open(vcf_file, 'r') as vcf, open('modPL_' + os.path.basename(vcf_file), 'w') as outfile:
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				outfile.write('##FORMAT=<ID=PL,Number=G,Type=Float,Description=\"Phred-scaled genotype likelihoods\">\n')
				outfile.write(line)
			else:
				outfile.write(line)
		else:
			pieces = line.rstrip().split()
			line_out = pieces[:8]
			line_out.append(pieces[8] + ':PL')
			for sample_piece in pieces[9:]:
				line_out.append(sample_piece + ':0.0,0.0,0.0')
			outfile.write('\t'.join(line_out) + '\n')

# report completion
print('Done')
