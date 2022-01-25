#!/usr/bin/env python

##########################
# Author: B. Anderson
# Date: Jan 2022
# Description: convert an IQTree Nexus concordance factor tree to Newick
##########################


import sys
import argparse
from Bio import Phylo


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to convert an IQTree Nexus concordance factor tree to Newick')


# add arguments to parse
parser.add_argument('-t', type = str, dest = 'tree_file', help = 'The IQTree Nexus file to convert')
parser.add_argument('-o', type = str, dest = 'out_pre', help = 'The output prefix [default \"output\"]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

tree_file = args.tree_file
out_pre = args.out_pre

if not tree_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'


# process the treefile to relabel internal nodes to the site concordance factors
# a typical comment looks like:
# [&sCF="55.74",sCF/sDF1/sDF2="55.74/19.87/24.39",sCF_N="12.29",sCF_N/sDF1_N/sDF2_N="12.29/4.2/5.15",sDF1="19.87",sDF1_N="4.2",sDF2="24.39",sDF2_N="5.15",sN="21.639"]
# we want the second one (sCF/sDF1/sDF2)
tree = Phylo.read(tree_file, 'nexus')
for node in tree.get_nonterminals():
	if node.comment:
		node.name = str(node.comment).split(',')[1].split('=')[1].strip('"')

# write the output tree
Phylo.write(tree, open(out_pre + '.tre', 'w'), 'newick')
