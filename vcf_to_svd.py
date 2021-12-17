#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Dec 2021
# Description: convert a VCF file (VCF 4.0) to a Nexus file for SVDQuartets
##########################


import sys
import argparse
import os


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert a VCF file to a Nexus input format required by SVDQuartets;' +
						' it is up to the user to ensure the VCF has the SNPs of interest (e.g. one per locus)')


# add arguments to parse
parser.add_argument('-v', type = str, dest = 'vcf_file', help = 'The VCF file to convert')
parser.add_argument('-s', type = str, dest = 'sample_file',
					help = 'The tab-delimited file of sample name and pop/taxon, one per line [required if specifying species]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
sample_file = args.sample_file

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


# process the sample file (if present) into a tuple list (sample, pop/taxon)
if sample_file:
	sample_list = []
	with open(sample_file, 'r') as infile:
		for line in infile:
			parts = line.strip().split()
			sample_list.append((parts[0], parts[1]))


# process the VCF to grab sample labels and genotypes for each SNP locus
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


# if there wasn't a sample list provided, create one with pop = sample (2 alleles per pop/ind)
if not sample_file:
	sample_list = []
	for sample in sample_labels:
		sample_list.append((sample, sample))


# create the output data lines to write to the Nexus files
lines_out = []
align_lines_out = []
samples_out = []
align_samples_out = []
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
		bases1 = []
		bases2 = []
		sequence = []
		# sample names can't have dashes in PAUP, so substitute with periods
		sample = sample.replace('-', '.')
		samples_out.append([sample + '_A', sample_tuple[1]])
		samples_out.append([sample + '_B', sample_tuple[1]])
		align_samples_out.append(sample)
		for genotypes in snps:
			genotype = genotypes[this_index]
			bases1.append(genotype[0])
			bases2.append(genotype[1])
			if genotype[0] == genotype[1]:
				sequence.append(genotype[0])
			else:
				sequence.append(amb_dict[genotype[0] + genotype[1]])
		line_out1 = sample + '_A' + '\t' + ''.join(bases1)
		line_out2 = sample + '_B' + '\t' + ''.join(bases2)
		lines_out.append(line_out1)
		lines_out.append(line_out2)
		align_lines_out.append(sample + '\t' + ''.join(sequence))


# determine if any samples in the VCF were not in the sample list
samples_missed = 0
for sample in sample_labels:
	if sample not in [item[0] for item in sample_list]:
		print('Sample ' + sample + ' is present in the VCF but not the sample list (so not included)')
		samples_missed = samples_missed + 1


# set output parameters
datatype = 'NUCLEOTIDE'
missing = '?'
num_chars = len(snps)
num_taxa = len(samples_out)
align_num_taxa = len(align_samples_out)
taxa_labels = [item[0] for item in samples_out]
species = []
pops_out = set([item[1] for item in samples_out])
for pop in pops_out:
	members = [item[0] for item in samples_out if item[1] == pop]
	pop = pop.replace('-', '.')		# quick substitution in case pops have dashes
	species.append(pop + ' : ' + ' '.join(members))


# create the Nexus files (one for SVDQuartets, one as an alignment for e.g. IQTREE)
with open(os.path.basename(vcf_file) + '.nex', 'w') as outfile:
	outfile.write('#NEXUS\n\n')
	taxa_block = ('BEGIN TAXA;\n\tDIMENSIONS NTAX=' + str(num_taxa) + ';\n\t' +
				'TAXLABELS ' + ' '.join(taxa_labels) + ';\nEND;\n')
	data_block = ('BEGIN DATA;\n\tDIMENSIONS NTAX=' + str(num_taxa) + ' NCHAR=' + str(num_chars) + ';\n\t' +
				'FORMAT DATATYPE=' + datatype + ' MISSING=' + missing + ';\n\tMATRIX\n')
	sets_block = ('BEGIN SETS;\nTAXPARTITION svdspecies =\n\t' + ',\n\t'.join(species) + ';\nEND;\n')
	outfile.write(taxa_block + '\n')
	outfile.write(data_block)
	for line_out in lines_out:
		outfile.write(line_out + '\n')
	outfile.write('\t;\nEND;\n\n')
	outfile.write(sets_block)

with open(os.path.basename(vcf_file) + '_align.nex', 'w') as outfile:
	outfile.write('#NEXUS\n\n')
	taxa_block = ('BEGIN TAXA;\n\tDIMENSIONS NTAX=' + str(align_num_taxa) + ';\n\t' +
				'TAXLABELS ' + ' '.join(align_samples_out) + ';\nEND;\n')
	data_block = ('BEGIN DATA;\n\tDIMENSIONS NTAX=' + str(align_num_taxa) + ' NCHAR=' + str(num_chars) + ';\n\t' +
				'FORMAT DATATYPE=' + datatype + ' MISSING=' + missing + ';\n\tMATRIX\n')
	outfile.write(taxa_block + '\n')
	outfile.write(data_block)
	for line_out in align_lines_out:
		outfile.write(line_out + '\n')
	outfile.write('\t;\nEND;\n\n')


# report completion
if sample_excl > 0:
	print('Did not include ' + str(sample_excl) + ' samples missing from the VCF')
if samples_missed > 0:
	print('Did not include ' + str(samples_missed) + ' samples missing from the sample list')
print('Converted a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_snps) +
	' SNPs to Nexus format with ' +	str(len(pops_out)) + ' species, ' + str(len(lines_out)) +
	' alleles and ' + str(num_chars) + ' characters')
print('Also created the corresponding alignment for ' + str(len(align_samples_out)) +
	' individuals and ' + str(num_chars) + ' variable sites')
