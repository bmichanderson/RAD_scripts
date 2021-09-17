#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Sep 2021
# Description: extract loci from an ipyrad .loci file based on an input list
##########################


import sys			# allows access to command line arguments
import argparse


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to extract loci from an ipyrad .loci file using an input list of locus numbers. ' +
                                                'The loci will be saved as individual fasta files in the current directory with prefix "locus_"')


# add arguments to parse
parser.add_argument('loci_file', type = str, help = 'The .loci file to extract loci from')
parser.add_argument('-l', type = str, dest = 'loci_list', help = 'A text file with a list of loci numbers, one per line (required)')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

loci_file = args.loci_file
loci_list = args.loci_list

if not loci_list:
	parser.print_help(sys.stderr)
	sys.exit(1)


# read in the list of loci and assign to a list
locus_nums = []
with open(loci_list, 'r') as loci:
    for locus in loci:
        locus_nums.append(locus.rstrip())


# extract the loci from the .loci file
with open(loci_file, 'r') as locfile:
    this_locus = []
    found_loci = 0
    processed_loci = 0
    for line in locfile:
        if line.startswith('//'):        # this indicates the final line of a locus
            processed_loci = processed_loci + 1
            locus_num = line.rstrip().split('|')[1]
            if locus_num in locus_nums:
                found_loci = found_loci + 1
                with open('locus_' + str(locus_num) + '.fasta', 'w') as outfile:
                    for entry in this_locus:
                        outfile.write('>' + entry.rstrip().split()[0] + '\n' + entry.rstrip().split()[1] + '\n')
            this_locus = []
        else:
            this_locus.append(line)


# summarize the results
print('Processed ' + str(processed_loci) + ' loci and extracted ' + str(found_loci) + ' loci of the targeted ' + str(len(locus_nums)))
