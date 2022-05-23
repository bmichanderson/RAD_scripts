#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: May 2022
# Description: create a sample set of population sequences for testing
##########################


import sys
import argparse
import random
import copy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to create a sample set of population sequences for testing')


# add arguments to parse
parser.add_argument(type = str, dest = 'outpre', help = 'The prefix for output files')
parser.add_argument('-l', type = int, dest = 'seqlen', help = 'The length of the sequences (default: 500)')
parser.add_argument('-p', type = int, dest = 'num_pops', help = 'The number of pops to create (default: 4)')
parser.add_argument('-i', type = int, dest = 'inds', help = 'The number of individuals per pop (default: 4)')
parser.add_argument('-v', type = float, dest = 'var', help = 'The proportion of bp that can have SNPs (default: 0.1)')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()
seqlen = args.seqlen
num_pops = args.num_pops
inds = args.inds
var = args.var
outpre = args.outpre

if not seqlen:
	seqlen = 500

if not num_pops:
	num_pops = 4

if not inds:
	inds = 4

if not var:
	var = 0.1


# determine the lengths of sequence bits
var_bits = round(var * seqlen)
qvar_bits = round(var_bits / 4)	# 1/4 of the SNPs
hvar_bits = round(var_bits / 2)	# 1/2 of the SNPs
rvar_bits = var_bits - qvar_bits - hvar_bits	# the rest of the SNPs
con_bits = seqlen - var_bits


# create the bp options
bps = ['A', 'T', 'C', 'G']


# create the ambiguity dictionaries
amb_dict = {
	'GT': 'K',
	'TG': 'K',
	'AC': 'M',
	'CA': 'M',
	'AG': 'R',
	'GA': 'R',
	'CG': 'S',
	'GC': 'S',
	'AT': 'W',
	'TA': 'W',
	'CT': 'Y',
	'TC': 'Y'
}

amb_dict2 = {
	'K': ['G', 'T'],
	'M': ['A', 'C'],
	'R': ['A', 'G'],
	'S': ['C', 'G'],
	'W': ['A', 'T'],
	'Y': ['C', 'T']
}

# set up the list of sequences
seqlists = list()


# iterate through and create the sequences
for pop_index in range(num_pops):
	for ind in range(inds):
		# we'll make 1/4 of the SNPs fixed (between pops if <=4 pops)
		myseq = bps[pop_index % 4] * qvar_bits
		seqlists.append(list(myseq))


# iterate again and add bases
seqindex = 0
for pop_index in range(num_pops):
	# we'll create a basic sequence to change
	base_seq = random.choices(bps, k = hvar_bits)
	# create an alternative allele at 1/2 of locations
	alt_seq = copy.deepcopy(base_seq)
	indices = random.sample(range(len(alt_seq)), k = round(len(alt_seq) / 2))
	for index in indices:
		alt_seq[index] = random.choice(list(set(bps).difference(alt_seq[index])))
	for ind in range(inds):
		# choose which of the two seqs to start with
		this_seq = copy.deepcopy(random.choice([base_seq, alt_seq]))
		# randomly change up to 1/2 of alts to hets (so up to 1/4 of original)
		choices = round(len(this_seq) / random.choice([4, 5, 6, 7, 8]))
		new_indices = random.sample(indices, k = choices)
		for index in new_indices:
			base1 = base_seq[index]
			base2 = alt_seq[index]
			this_seq[index] = amb_dict[base1 + base2]
		# add random bases for the remaining SNPs (including missing)
		this_seq.extend(random.choices(list(''.join(bps) + '?'), k = rvar_bits))
		# extend the sequence in the list
		seqlists[seqindex].extend(this_seq)
		seqindex = seqindex + 1


# add constant DNA characters to all seqs
constant_seq = random.choices(bps, k = con_bits)
for seq in seqlists:
	seq.extend(constant_seq)


# turn the sequences into SeqRecords for writing to file
seqrecs = list()
for index, seq in enumerate(seqlists):
	seqrec = SeqRecord(Seq(''.join(seq)), id = 'ind' + str(index + 1), description = 'ind' + str(index + 1))
	seqrecs.append(seqrec)


# write the sequences to fasta
with open(outpre + '.fasta', 'w') as out1:
	SeqIO.write(seqrecs, out1, 'fasta')


# write an Arlequin type file and a Nexus file and a single line fasta
# we need to write two haplotypes per ind
with open(outpre + '.arp', 'w') as out2, open(outpre + '.nex', 'w') as out3, open(outpre + '.fsa', 'w') as out4:
	profile_block = '[Profile]\n' + \
		'\tTitle="Fst Test"\n' + \
		'\tNbSamples=' + str(num_pops) + '\n' + \
		'\tGenotypicData=1\n' + \
		'\tMissingData="?"\n' + \
		'\tDataType=DNA\n' + \
		'\tLocusSeparator=NONE\n' + \
		'\tGameticPhase=0\n'
	out2.write(profile_block)
	out2.write('[Data]\n  [[Samples]]\n')

	out3.write('#NEXUS\n')
	data_block = ('BEGIN DATA;\n\tDIMENSIONS NTAX=' + str(num_pops * inds) +
				' NCHAR=' + str(seqlen) + ';\n\t' +
				'FORMAT\n\t\tDATATYPE=DNA\n\t\tMISSING=?\n\t\t' +
				';\n\tMATRIX\n')
	out3.write(data_block)

	ind_index = 0
	for pop in range(num_pops):
		out2.write('\tSampleName="pop' + str(pop + 1) + '"\n')
		out2.write('\tSampleSize=' + str(inds * 2) + '\n')
		out2.write('\tSampleData={\n')
		for ind in range(inds):
			this_seq = seqlists[ind_index]
			seq1 = list()
			seq2 = list()
			for bp in this_seq:
				if any([bp in bps, bp == '?']):
					seq1.append(bp)
					seq2.append(bp)
				elif bp in ['K', 'M', 'R', 'S', 'W', 'Y']:
					seq1.append(amb_dict2[bp][0])
					seq2.append(amb_dict2[bp][1])
				else:
					print('PROBLEM!')

			out2.write('\t\tind' + str(ind_index + 1) + ' 2\t' + ''.join(seq1) + '\n')
			out2.write('\t\t\t' + ''.join(seq2) + '\n')

			out3.write('\t\tind' + str(ind_index + 1) + 'a\t' + ''.join(seq1) + '\n')
			out3.write('\t\tind' + str(ind_index + 1) + 'b\t' + ''.join(seq2) + '\n')

			out4.write('>ind' + str(ind_index + 1) + 'a\n')
			out4.write(''.join(seq1) + '\n')
			out4.write('>ind' + str(ind_index + 1) + 'b\n')
			out4.write(''.join(seq2) + '\n')

			ind_index = ind_index + 1
		out2.write('}\n')
	
	out3.write('\t;\nEND;\n')


# write a populations file
with open(outpre + '_pops.tab', 'w') as out5:
	ind_index = 1
	for pop in range(num_pops):
		for ind in range(inds):
			out5.write('ind' + str(ind_index) + '\tpop' + str(pop + 1) + '\n')
			ind_index = ind_index + 1


# and one for the two haplotypes
with open(outpre + '_pops2.tab', 'w') as out6:
	ind_index = 1
	for pop in range(num_pops):
		for ind in range(inds):
			out6.write('ind' + str(ind_index) + 'a\tpop' + str(pop + 1) + '\n')
			out6.write('ind' + str(ind_index) + 'b\tpop' + str(pop + 1) + '\n')
			ind_index = ind_index + 1
