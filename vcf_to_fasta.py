#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Dec 2021
# Description: convert a VCF file (VCF 4.0) to a fasta alignment
##########################


import sys
import argparse
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file to a fasta alignment')


# add arguments to parse
parser.add_argument('-v', type = str, dest = 'vcf_file', help = 'The VCF file to convert')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)


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


# process the VCF to grab genotypes for each SNP locus
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
			genotypes = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) 0/0 then split that into the two alleles 0 and 0 as a list
					Then convert that into bases
				'''
				genotype = call.split(':')[0].split('/')
				allele1 = genotype[0]
				allele2 = genotype[1]
				if allele1 == '.':		# missing data
					genotypes.append(['?', '?'])
				else:
					genotypes.append([all_bases[int(allele1)], all_bases[int(allele2)]])
			snps.append(genotypes)
			count_snps = count_snps + 1
	print('Processed a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) + ' SNPs')


# grab the sequences and convert to Seq Records
new_records = []
for index, sample_label in enumerate(sample_labels):
	bases1 = []
	bases2 = []
	sequence = []
	for genotypes in snps:
		genotype = genotypes[index]
		if genotype[0] == genotype[1]:
			sequence.append(genotype[0])
		else:
			sequence.append(amb_dict[genotype[0] + genotype[1]])
	new_sequence = Seq(''.join(sequence))
	new_record = SeqRecord(new_sequence, id = sample_label, description = sample_label)
	new_records.append(new_record)


# create the alignment and write the fasta
myalign = MultipleSeqAlignment(new_records)
with open(os.path.basename(vcf_file).replace('.vcf', '.fasta'), 'w') as outfile:
	AlignIO.write(myalign, outfile, 'fasta')


# report completion
print('Converted a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) +
	' SNPs to fasta of length ' + str(myalign.get_alignment_length()) + ' bp')
