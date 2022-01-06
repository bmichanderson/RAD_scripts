#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Dec 2021
# Description: filter a VCF file by samples and/or snps
##########################


import sys
import argparse
import statistics


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to filter a VCF file')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to filter')
parser.add_argument('-o', type = str, dest = 'outpre', help = 'Output file prefix [default \"output\"]')
parser.add_argument('-s', type = str, dest = 'samps', help = 'A file with sample names to be dropped, one per line')
parser.add_argument('--bial', type = str, dest = 'bial', help = 'Whether a locus must be biallelic to be kept (\"yes\" or \"no\" [default])')
parser.add_argument('--mincov', type = float, dest = 'mincov', help = 'The minimum cover as a proportion of samples a locus must be found in to be kept')
parser.add_argument('--mac', type = int, dest = 'mac', help = 'The minimum minor allele count (how many times the minor allele must be called to keep the locus)')
parser.add_argument('--maxd', type = int, dest = 'maxd', help = 'The maximum call depth found anywhere in a locus to keep it')
parser.add_argument('--maxmd', type = int, dest = 'maxmd', help = 'The maximum mean call depth for a locus to keep it')
parser.add_argument('--mind', type = int, dest = 'mind', help = 'The minimum call depth found anywhere in a locus to keep it')
parser.add_argument('--minmd', type = int, dest = 'minmd', help = 'The minimum mean call depth for a locus to keep it')



# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
out_pre = args.outpre
sample_file = args.samps
bial = args.bial
mincov = args.mincov
minmac = args.mac
maxd = args.maxd
maxmd = args.maxmd
mind = args.mind
minmd = args.minmd

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# parse the samples file, if present
if sample_file:
	samples_remove = []
	with open(sample_file, 'r') as samples:
		for sample in samples:
			samples_remove.append(sample.rstrip())


# parse and filter the VCF
with open(vcf_file, 'r') as vcf, open(out_pre + '.vcf', 'w') as outfile:
	in_snps = 0
	in_samples = 0
	out_snps = 0
	for line in vcf:
		if line.startswith('#'):		# a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				pieces = line.rstrip().split()
				sample_labels = pieces[9:]
				in_samples = len(sample_labels)
				if sample_file:			# if we need to remove any
					indices = []
					for index, sample in enumerate(sample_labels):
						if sample not in samples_remove:
							indices.append(index)
					newline = '\t'.join(pieces[:9] + [pieces[9:][index] for index in indices])
					sample_labels = [sample_labels[index] for index in indices]
					outfile.write(newline + '\n')
				else:
					outfile.write(line)
			else:
				outfile.write(line)
		else:
			in_snps = in_snps + 1
			pieces = line.rstrip().split()
			if bial:
				if bial.lower() == 'yes':
					if len(pieces[4].split(',')) > 1:
						continue
			if sample_file:		# if samples were removed
				calls = [pieces[9:][index] for index in indices]
			else:
				calls = pieces[9:]
			depths = []
			genotypes = []
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the genotype (GT) and depth (DP)
				'''
				gt = call.split(':')[0].split('/')		# a list, e.g. ['0', '0']
				genotypes.append(gt)		
				if gt != ['.', '.']:		# if not an N (shouldn't be needed, but currently is)
					depths.append(int(call.split(':')[1]))
			# calculate stats depending on what filters need to be applied
			keep_locus = True
			while keep_locus:
				if mincov:
					count_present = len([geno for geno in genotypes if geno != ['.', '.']])
					if count_present / len(sample_labels) < mincov:
						keep_locus = False
				if minmac:
					count_0 = ''.join(''.join(geno) for geno in genotypes).count('0')
					count_1 = ''.join(''.join(geno) for geno in genotypes).count('1')
					if any([count_0 < minmac, count_1 < minmac]):
						keep_locus = False
				if maxd:
					if max(depths) > maxd:
						keep_locus = False
				if maxmd:
					if sum(depths) > 0:
						if statistics.mean([depth for depth in depths if depth > 0]) > maxmd:
							keep_locus = False
					else:
						keep_locus = False
				if mind:
					if min(depths) < mind:
						keep_locus = False
				if minmd:
					if sum(depths) > 0:
						if statistics.mean([depth for depth in depths if depth > 0]) < minmd:
							keep_locus = False
					else:
						keep_locus = False
				# write the locus if it should be kept
				if keep_locus:
					newline = '\t'.join(pieces[:9] + calls)
					outfile.write(newline + '\n')
					out_snps = out_snps + 1
					break
	print('Read a VCF file with ' + str(in_samples) + ' samples and ' + str(in_snps) + ' SNP loci' +
		' and filtered it to ' + str(len(sample_labels)) + ' samples and ' + str(out_snps) + ' SNP loci')