#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Sep 2021
# Modified: Nov 2021 (added a sample exclude option); Dec 2023 (increased summary and added consensus)
# Description: extract loci from an ipyrad .loci file based on an input list
##########################


import sys
import argparse
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# create an ambiguity dictionary
amb_dict = {
	'AC': 'M',
	'AG': 'R',
	'AT': 'W',
	'CG': 'S',
	'CT': 'Y',
	'GT': 'K',
	'ACG': 'V',
	'ACT': 'H',
	'AGT': 'D',
	'CGT': 'B',
	'ACGT': 'N'
}


# create a function for computing the consensus sequence of a list of aligned strings
def make_consensus(seq_list, locus_num):
	con_list = []
	for position in range(len(seq_list[0])):
		align_slice = [sequence[position] for sequence in seq_list]
		nbases = [nbase for nbase in align_slice if nbase in ['-', 'N']]
		bases = [base for base in align_slice if base not in ['-', 'N']]
		base_list = []
		if len(bases) != 0:		# a residue is there
			if bases.count('A') > 0:
				base_list.append('A')
			if bases.count('C') > 0:
				base_list.append('C')
			if bases.count('G') > 0:
				base_list.append('G')
			if bases.count('T') > 0:
				base_list.append('T')

			if  any([len(base_list) == 0, len(nbases) > len(bases)]):		# no bases above threshold or more missing data than data
				base = 'N'
			elif len(base_list) == 1:
				base = base_list[0]
			else:
				base = amb_dict[''.join(base_list)]

			con_list.append(base)
		else:
			con_list.append('N')

	return(SeqRecord(Seq(''.join(con_list)), id = str(locus_num), name = str(locus_num), description = str(locus_num)))


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract loci from an ipyrad .loci file using an input list of locus numbers. ' +
	'The loci will be saved as individual fasta files in the current directory with prefix "locus_"')


# add arguments to parse
parser.add_argument('loci_file', type = str, help = 'The .loci file to extract loci from')
parser.add_argument('-c', action = 'store_true', help = 'Whether to output a file with consensus sequences of the loci (optional)')
parser.add_argument('-l', type = str, dest = 'loci_list', help = 'A text file with a list of desired loci numbers, one per line (optional). ' +
	'If not provided, will extract all loci')
parser.add_argument('-s', type = str, dest = 'samp_list', help = 'A text file with a list of sample names to exclude, one per line (optional). ' +
	'If not provided, will include all samples')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
consens = args.c
loci_file = args.loci_file
loci_list = args.loci_list
samp_list = args.samp_list


# check for consensus and create a list to hold them
if consens:
	consensus_list = []


# read in the list of loci (if present) and assign to a list
if loci_list:
	locus_nums = []
	with open(loci_list, 'r') as loci:
		for locus in loci:
			locus_nums.append(locus.rstrip())


# read in the list of samples (if present) and assign to a list
if samp_list:
	samp_names = []
	with open(samp_list, 'r') as samples:
		for sample in samples:
			samp_names.append(sample.rstrip())
	print('Will exclude ' + str(len(samp_names)) + ' samples from extracted loci')


# extract the loci from the .loci file
with open(loci_file, 'r') as locfile:
	this_locus = []
	found_loci = 0
	processed_loci = 0
	empty_loci = 0
	locus_lengths = []
	for line in locfile:
		if line.startswith('//'):        # this indicates the final line of a locus
			processed_loci = processed_loci + 1
			locus_num = line.rstrip().split('|')[1]
			if loci_list:
				if locus_num in locus_nums:
					found_loci = found_loci + 1
					if len(this_locus) > 0:
						with open('locus_' + str(locus_num) + '.fasta', 'w') as outfile:
							for entry in this_locus:
								outfile.write('>' + entry.rstrip().split()[0] + '\n' + entry.rstrip().split()[1] + '\n')
						locus_lengths.append(len(this_locus[0].rstrip().split()[1]))
						if consens:
							text_seqs = [entry.rstrip().split()[1] for entry in this_locus]
							conseq = make_consensus(text_seqs, locus_num)
							consensus_list.append(conseq)
					else:
						empty_loci = empty_loci + 1
			else:
				found_loci = found_loci + 1
				if len(this_locus) > 0:
					with open('locus_' + str(locus_num) + '.fasta', 'w') as outfile:
						for entry in this_locus:
							outfile.write('>' + entry.rstrip().split()[0] + '\n' + entry.rstrip().split()[1] + '\n')
					locus_lengths.append(len(this_locus[0].rstrip().split()[1]))
					if consens:
						text_seqs = [entry.rstrip().split()[1] for entry in this_locus]
						conseq = make_consensus(text_seqs, locus_num)
						consensus_list.append(conseq)
				else:
					empty_loci = empty_loci + 1
			this_locus = []
		else:
			if samp_list:
				if line.split()[0] not in samp_names:
					this_locus.append(line)
			else:
				this_locus.append(line)


# output the consensus if requested
with open('loci_consensus.fasta', 'w') as outfile:
	for conseq in consensus_list:
		SeqIO.write(conseq, outfile, 'fasta')


# summarize the results
if loci_list:
	print('Processed ' + str(processed_loci) + ' loci and extracted ' + str(found_loci) + ' loci of the targeted ' + str(len(locus_nums)))
else:
	print('Processed ' + str(processed_loci) + ' loci and extracted ' + str(found_loci) + ' loci')
if empty_loci > 0:
	print('Did not export ' + str(empty_loci) + ' loci that only included excluded taxa')

mean = sum(locus_lengths)/len(locus_lengths)
stdev = statistics.pstdev(locus_lengths)

print('For ' + str(len(locus_lengths)) + ' retained loci, the mean length was ' + str('{:.2f}'.format(mean)) +
	' bp (min '+ str(min(locus_lengths)) + ', max ' + str(max(locus_lengths)) + ', stdev ' + str('{:.2f}'.format(stdev)) + ')')
