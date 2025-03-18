#!/usr/bin/env python3

##########################
# Author: B.M. Anderson
# Date: Nov 2021
# Modified: Mar 2022 (report stats; corrected for ipyrad behaviour); Mar 2025 (added standard deviation and per sample)
# Description: capture and plot read depth information from a VCF file (VCF 4.0)
##########################


import sys
import argparse
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statistics


# instantiate the parser
parser = argparse.ArgumentParser(description = 'A script to capture and plot read depth by site from a VCF file')


# add arguments to parse
parser.add_argument('vcf_file', type = str, help = 'The VCF file to convert')
parser.add_argument('-o', type = str, dest = 'output', help = 'Name of the output pdf without extension [default \"output\"]')
parser.add_argument('-m', type = str, dest = 'maxd', help = 'Maximum potential read depth to graph [default 500]')


# parse the command line
if len(sys.argv[1:]) == 0:		# if there are no arguments
	parser.print_help(sys.stderr)
	sys.exit(1)

args = parser.parse_args()

vcf_file = args.vcf_file
out_pre = args.output
maxd = args.maxd

if not vcf_file:
	parser.print_help(sys.stderr)
	sys.exit(1)

if not out_pre:
	out_pre = 'output'
if not maxd:
	maxd = 500
else:
	maxd = int(maxd)


# define a function to set a max for graphing
def use_max(my_list):
	return(min(maxd, max(my_list)))


# process the VCF to grab sample labels and depths for each SNP locus
with open(vcf_file, 'r') as vcf:
	sample_list = []
	snps = []
	count_site = 0
	bad_site = 0
	for line in vcf:
		if line.startswith('#'):        # a header INFO line
			if line.startswith('#CHROM'):		# the line with sample names
				sample_labels = line.rstrip().split()[9:]
				for sample in sample_labels:
					sample_list.append([])
		else:
			calls = line.rstrip().split()[9:]
			depths = []
			index = 0
			for call in calls:
				'''	The format is GT:DP:CATG, so 0/0:85:0,85,0,0 for a homozygous AA with 85 depth
					We want to grab the depth (DP)
				'''
				call_pieces = call.split(':')
				gt = call_pieces[0].split('/')
				if gt != ['.', '.']:		# if not an N (shouldn't be needed, but currently is)
					depth = int(call_pieces[1])
					depths.append(depth)
					if depth > 0:
						sample_list[index].append(depth)
				index = index + 1

			# check if there is at least one depth > 0 at this site
			# if so, append the depths to the snp list
			if len([item for item in depths if item > 0]) > 0:
				snps.append(depths)
				count_site = count_site + 1
			else:
				bad_site = bad_site + 1
	print('Read in a VCF file with ' + str(len(sample_labels)) + ' samples and ' + str(count_site) + ' SNP sites')
	if bad_site > 0:
		print('There were also ' + str(bad_site) + ' SNP sites with no depths > 0 !')


# calculate the mean and max depth per site
# the mean is based only on present data (no missing or 0 depth)
# also collate the individual SNP depths for reporting overall standard deviation
means = []
maxs = []
for site in snps:
	means.append(round(statistics.mean([item for item in site if item > 0])))
	maxs.append(round(max(site)))


# calculate the mean depth across all sites per sample
# also calculate the overall depths for standard deviation
sample_means = []
all_depths = []
for sample_entry in sample_list:
	sample_means.append(round(statistics.mean(sample_entry)))
	all_depths.extend(sample_entry)


# plot histograms
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, figsize = (10, 10))
plt.subplots_adjust(hspace = 0.2)
ax1.hist(all_depths, density = False, range = [0, use_max(all_depths)],
	bins = round(use_max(all_depths) / 5))
ax1.set_title('Depth per individual SNP')
ax1.set_ylabel('SNPs')
ax2.hist(means, density = False, range = [0, use_max(means)],
	bins = round(use_max(means) / 5))
ax2.set_title('Mean depth per site')
ax2.set_ylabel('Site')
ax3.hist(maxs, density = False, range = [0, use_max(maxs)],
	bins = round(use_max(maxs) / 5))
ax3.set_title('Maximum depth per site')
ax3.set_ylabel('Site')
ax3.set_xlabel('Depth')
ax3.set_xlim(0, use_max(maxs) + 0.01 * use_max(maxs))
plt.savefig(out_pre + '_depth.pdf')


# plot depth per sample
array_list = []
for sample_entry in sample_list:
	array_list.append(numpy.array(sample_entry))

target_width = max(len(sample_labels) / 5, 18)
plt.figure(2, figsize = (target_width, target_width / 2))
plt.xticks(rotation = 90)
plt.boxplot(array_list, labels = sample_labels, showfliers = False)
plt.savefig(out_pre + '_depth_boxplot.png', dpi = 300)


# report mean depth per sample
print('Mean depth per sample:')
for index, mean in enumerate(sample_means):
	print(str(sample_labels[index]) + '\t' + str(mean))


# report statistics
print('Minimum mean depth per site: ' + str(min(means)))
print('Maximum mean depth per site: ' + str(max(means)))
print('Overall mean depth across sites: ' + str(round(statistics.mean(all_depths), 2)))
print('Overall stdev depth across sites: ' + str(round(statistics.stdev(all_depths), 2)))
