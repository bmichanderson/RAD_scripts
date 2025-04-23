#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: Sep 2021
# Modified: Nov 2021 (added a sample exclude option); Dec 2023 (increased summary and added consensus);
#	Dec 2023 (extended to Stacks allelic output); Mar 2025 (decreased counter reporting; adjusted consensus calculations; improved Stacks efficiency)
# Description: extract loci from an ipyrad .loci file (or Stacks populations.samples.fa) based on an input list
# Note: this script will convert the allelic data from Stacks (if run that way) to a single consensus per sample
##########################


import sys
import argparse
import statistics
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# create a function for computing the consensus sequence of a list of aligned sequences
# because the degenerate_consensus is behaving oddly, go back to the old approach
# note: it was inserting a 'V' whenever the input position had only 'N' or non-ACGT
# the current approach (below) will at least avoid this, though any lone non-ambiguity will be consensus
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

def make_consensus(seq_list, locus_num):
#	mymotif = motifs.create(seq_list)
#	consensus = mymotif.degenerate_consensus
	con_list = []
	for position in range(len(seq_list[0])):
		align_slice = [entry[position] for entry in seq_list]
		bases = [base for base in align_slice if base not in ['-', 'N']]
		base_list = []
		if len(bases) != 0:
			if bases.count('A') > 0:
				base_list.append('A')
			if bases.count('C') > 0:
				base_list.append('C')
			if bases.count('G') > 0:
				base_list.append('G')
			if bases.count('T') > 0:
				base_list.append('T')

			if len(base_list) == 1:
				con_list.append(base_list[0])
			elif len(base_list) > 1:
				con_list.append(amb_dict[''.join(base_list)])
			else:		# position only contains ambiguities
				# for simplicity, just return the first one (this isn't a true consensus)
				con_list.append(bases[0])

		elif 'N' in align_slice:
			con_list.append('N')
		else:
			con_list.append('-')

	return(SeqRecord(Seq(''.join(con_list)), id = str(locus_num), name = str(locus_num), description = str(locus_num)))


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract loci from an ipyrad .loci file using an input list of locus numbers. ' +
	'The loci will be saved as individual fasta files in the current directory with prefix "locus_". ' +
	'Now also works for Stacks input (use arg --stacks)')


# add arguments to parse
parser.add_argument('loci_file', type = str, help = 'The .loci file to extract loci from')
parser.add_argument('-c', action = 'store_true', help = 'Flag for whether to output a file with consensus sequences of the loci (optional)')
parser.add_argument('-l', type = str, dest = 'loci_list', help = 'A text file with a list of desired loci numbers, one per line (optional). ' +
	'If not provided, will extract all loci')
parser.add_argument('-s', type = str, dest = 'samp_list', help = 'A text file with a list of sample names to exclude, one per line (optional). ' +
	'If not provided, will include all samples')
parser.add_argument('--stacks', action = 'store_true', help = 'Flag for whether the input file is from Stacks (optional) [default: ipyrad]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)
args = parser.parse_args()
consens = args.c
loci_file = args.loci_file
loci_list = args.loci_list
samp_list = args.samp_list
stacks = args.stacks


# read in the list of loci (if present) and assign to a list
if loci_list:
	locus_nums = []
	with open(loci_list, 'r') as loci:
		for locus in loci:
			locus_nums.append(locus.rstrip())
	print('Will attempt to extract ' + str(len(locus_nums)) + ' loci')


# read in the list of samples (if present) and assign to a list
if samp_list:
	samp_names = []
	with open(samp_list, 'r') as samples:
		for sample in samples:
			samp_names.append(sample.rstrip())
	print('Will exclude ' + str(len(samp_names)) + ' samples from extracted loci')


# check for consensus and create a list to hold them
if consens:
	consensus_list = []
	print('Will generate a consensus sequence for each locus')


# extract the loci from the .loci file
with open(loci_file, 'r') as locfile:
	counter = 1
	found_loci = 0
	processed_loci = 0
	empty_loci = 0
	locus_lengths = []

	###### ipyrad
	if not stacks:		# the ipyrad .loci file has special formatting
		this_locus = []
		for line in locfile:
			if line.startswith('//'):        # this indicates the final line of a locus
				processed_loci = processed_loci + 1
				if processed_loci > 1000 * counter:
					print('Processed ' + str(1000 * counter) + ' loci')
					counter = counter + 1
				locus_num = line.strip().split('|')[1]
				if loci_list:
					if locus_num in locus_nums:
						found_loci = found_loci + 1
						if len(this_locus) > 0:
							with open('locus_' + str(locus_num) + '.fasta', 'w') as outfile:
								for entry in this_locus:
									SeqIO.write(entry, outfile, 'fasta')
							locus_lengths.append(len(this_locus[0].seq))
							if consens:
								seqs = [entry.seq for entry in this_locus]
								conseq = make_consensus(seqs, locus_num)
								consensus_list.append(conseq)
						else:
							empty_loci = empty_loci + 1
				else:
					found_loci = found_loci + 1
					if len(this_locus) > 0:
						with open('locus_' + str(locus_num) + '.fasta', 'w') as outfile:
							for entry in this_locus:
								SeqIO.write(entry, outfile, 'fasta')
						locus_lengths.append(len(this_locus[0].seq))
						if consens:
							seqs = [entry.seq for entry in this_locus]
							conseq = make_consensus(seqs, locus_num)
							consensus_list.append(conseq)
					else:
						empty_loci = empty_loci + 1
				this_locus = []
			else:
				parts = line.strip().split()
				if samp_list:
					if parts[0] not in samp_names:
						this_locus.append(SeqRecord(Seq(parts[1]), id = str(parts[0]), name = str(parts[0]), description = str(parts[0])))
				else:
					this_locus.append(SeqRecord(Seq(parts[1]), id = str(parts[0]), name = str(parts[0]), description = str(parts[0])))

	###### Stacks
	else:		# the Stacks file is already fasta (note: two seqs per sample = alleles); assuming ordered by locus
		this_locus_num = ''
		store = False
		seq_list = []
		for line in locfile:
			if line.startswith('#'):
				continue
			elif line.startswith('>'):
				parts = line.strip().split()		# format: ">CLocus_1_Sample_100_Locus_1_Allele_0 [sampleID]"
				locus_num = parts[0].split('_')[1]
				sample = parts[1].strip('[').strip(']')
				if loci_list:
					if locus_num in locus_nums:		# a target locus
						store = True
					else:
						store = False
				else:
					store = True
			else:
				if locus_num != this_locus_num:
					processed_loci = processed_loci + 1
					if processed_loci > 1000 * counter:
						print('Processed ' + str(1000 * counter) + ' loci')
						counter = counter + 1

					# process the stored seq_list to ensure each sample is represented by only one sequence
					# generate loci output and store consensus if requested
					if len(seq_list) > 0:
						found_loci = found_loci + 1
						seq_out_list = []
						sample_list = sorted(list(set([item[0] for item in seq_list])))
						for samp in sample_list:
							seq_recs = [item[1] for item in seq_list if item[0] == samp]
							if len(seq_recs) == 1:	# only represented by one sequence
								print('Sample ' + str(samp) + ' has one allele for locus ' + str(this_locus_num))
								seq_out_list.append(seq_recs[0])		
							elif len(seq_recs) > 2:	# shouldn't happen
								print('Sample ' + str(samp) + ' has more than two alleles for locus ' + str(this_locus_num))
							else:		# normal two alleles
								sample_seq = make_consensus([entry.seq for entry in seq_recs], samp)
								seq_out_list.append(sample_seq)
						with open('locus_' + str(this_locus_num) + '.fasta', 'w') as outfile:
							for entry in seq_out_list:
								SeqIO.write(entry, outfile, 'fasta')
						locus_lengths.append(len(seq_out_list[0].seq))
						if consens:
							conseq = make_consensus([entry.seq for entry in seq_out_list], this_locus_num)
							consensus_list.append(conseq)

					# reset the locus_num and seq_list
					this_locus_num = locus_num
					seq_list = []

				# store the line (sequence)
				if store:
					if samp_list:
						if sample not in samp_names:		# a desirable sample
							seq_list.append([sample,
								SeqRecord(Seq(line.strip()), id = str(sample), name = str(sample), description = str(sample))])
					else:
						seq_list.append([sample,
							SeqRecord(Seq(line.strip()), id = str(sample), name = str(sample), description = str(sample))])

		# run once more for last locus
		if len(seq_list) > 0:
			found_loci = found_loci + 1
			seq_out_list = []
			sample_list = sorted(list(set([item[0] for item in seq_list])))
			for samp in sample_list:
				seq_recs = [item[1] for item in seq_list if item[0] == samp]
				if len(seq_recs) == 1:	# only represented by one sequence
					print('Sample ' + str(samp) + ' has one allele for locus ' + str(this_locus_num))
					seq_out_list.append(seq_recs[0])				
				elif len(seq_recs) > 2:	# shouldn't happen
					print('Sample ' + str(samp) + ' has more than two alleles for locus ' + str(this_locus_num))
				else:		# normal two alleles
					sample_seq = make_consensus([entry.seq for entry in seq_recs], samp)
					seq_out_list.append(sample_seq)
			with open('locus_' + str(this_locus_num) + '.fasta', 'w') as outfile:
				for entry in seq_out_list:
					SeqIO.write(entry, outfile, 'fasta')
			locus_lengths.append(len(seq_out_list[0].seq))
			if consens:
				conseq = make_consensus([entry.seq for entry in seq_out_list], this_locus_num)
				consensus_list.append(conseq)


# output the consensus if requested
if consens:
	with open('loci_consensus.fasta', 'w') as outfile:
		for conseq in consensus_list:
			SeqIO.write(conseq, outfile, 'fasta')


# summarize the results
if loci_list:
	print('Processed ' + str(processed_loci) + ' loci and extracted ' + str(found_loci) + ' loci of the targeted ' + str(len(locus_nums)))
else:
	print('Processed ' + str(processed_loci) + ' loci and extracted ' + str(found_loci) + ' loci')
if empty_loci > 0:
	print('Did not export ' + str(empty_loci) + ' loci that only included excluded samples')

mean = sum(locus_lengths)/len(locus_lengths)
stdev = statistics.pstdev(locus_lengths)

print('For ' + str(len(locus_lengths)) + ' retained loci, the mean length was ' + str('{:.2f}'.format(mean)) +
	' bp (min '+ str(min(locus_lengths)) + ', max ' + str(max(locus_lengths)) + ', stdev ' + str('{:.2f}'.format(stdev)) + ')')
