#!/bin/bash

##################################
# Author: B. Anderson
# Date: Jul-Aug 2021
# Description: run a series of scripts to report metrics from replicate ipyrad runs to optimize parameters
##################################


# Script paths for the analysis/plotting script and the Mastretta-Yanes error script
MM_script="/home/benjamin/RAD_scripts/MM_ipyrad_summary.R"
MY_script="/home/benjamin/RAD_scripts/MY_replicate_error.R"


######### Data setup
# Output from the runs should be in the working directory, in folders named clust_{num}
# An example for getting the data from Pawsey and putting in the current directory might be:
#
# for num in {88..99}; do if [ ! -d clust_"$num" ]; then mkdir clust_"$num"; fi && \
# rsync -vrt --progress anderson@hpc-data.pawsey.org.au:/scratch/pawsey0220/anderson/MM_test/MMtest/rep_clust_"$num"_outfiles/*.vcf clust_"$num"/; done
#
# for num in {88..99}; do rsync -vrt --progress anderson@hpc-data.pawsey.org.au:/scratch/pawsey0220/anderson/MM_test/MMtest/rep_clust_"$num"_*/*.txt clust_"$num"/; done
#
# for num in {88..99}; do rsync -vrt --progress anderson@hpc-data.pawsey.org.au:/scratch/pawsey0220/anderson/MM_test/MMtest/rep_clust_"$num"_tree_outfiles/*.phy clust_"$num"/; done
#
# The text files from the ipyrad runs are important, as well as the vcf files
# The .phy files are useful, but would require also generating runs of trees before being able to plot
#
# The naming convention for the ouput might be "rep_clust" then _(paramvalue); it should be consistent
#########


# function for when the script is called incorrectly or without arguments
usage()
{
	echo -e "This script is for running a series of metrics on replicate ipyrad runs\n"
	echo -e "Usage:\topt_ipyrad.sh -n name_prefix -o out_prefix -p param_range -r reps_present -t tree_folder\n"
	echo -e "Options:\n" \
		"\t-n\tNaming convention for the ipyrad output; i.e. the text prefixing _paramvalue* in output\n" \
		"\t-o\tOutput prefix to name files produced during the run (default \"output\")\n" \
		"\t-p\tThe range of the clustering parameter, e.g. \"88 99\"\n" \
		"\t-r\tWhether replicates are present (yes [default] or no)\n" \
		"\t-s\tSamples file in the format: sample(same as vcf label) <tab> pop <tab> lat <tab> lon\n" \
		"\t-t\tIf you have phylogenetic trees, indicate the name of the folder they are in\n" \
		"\t\t(They should be in the format: value_rep.treefile, e.g. 88_4.treefile\n"
	exit 1
}

if [ $# -eq 0 ]; then		# if no arguments
	usage
fi


# Parse arguments

have_trees="False"
reps_present="yes"

FILES=()

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-n)
	name_prefix="$2"
	shift
	shift
	;;
	-o)
	out_prefix="$2"
	shift
	shift
	;;
	-p)
	param_range="$2"
	shift
	shift
	;;
	-r)
	reps_present="$2"
	shift
	shift
	;;
	-s)
	samp_file="$(readlink -f $2)"		# the full path to the samples file
	shift
	shift
	;;
	-t)
	tree_folder="$2"
	have_trees="True"
	shift
	shift
	;;
	*)
	FILES+=("$1")		# capture arguments without options in an array
	shift
	;;
esac
done


# Check options

if [ -z "$name_prefix" ]; then
	usage
fi

if [ -z "$out_prefix" ]; then
	out_prefix="output"
fi

if [ -z "$param_range" ]; then
	echo -e "\n\n*** Specify a clustering parameter range in the form: \"start end\" **\n\n"
	usage
fi

if [ -z "$samp_file" ]; then
	usage
fi


echo -e "\nStarting at $(date)\n"



########
# Create text files across parameter values to plot with the R script
########

# 1 Proportion inferred paralogs -- could potentially be measured in different ways (?)
# A) requires the final *_stats.txt file from ipyrad step 7
# We grab fields 3 and 4 of the first table, which correspond to number flagged for hets out of total after filter
# B) uses the s5 stats file for individual samples
# We need to grab the filter by hets and max alleles (fields 4 and 5) as well as the total (field 2)
# minus the amount filtered by depth (field 3)
for num in $(seq $param_range)
do
	grep "max_shared_het" clust_"$num"/"$name_prefix"_"$num"_stats.txt | \
	tr -s " " "\t" | cut -f 3,4 | sed "s/^/$num\t/" >> "$out_prefix"_paralogs1.tab && \
	tail -n +2 clust_"$num"/s5_consens*.txt | \
	tr -s " " "\t" | cut -f 2,3,4,5 | sed "s/^/$num\t/" >> "$out_prefix"_paralogs2.tab
done


# 2 Heterozygosity
# (requires the s5 stats file)
for num in $(seq $param_range)
do
	tail -n +2 clust_"$num"/s5_consens*.txt | \
	tr -s " " "\t" | cut -f 10 | sed "s/^/$num\t/" >> "$out_prefix"_heterozygosity.tab
done


# 3 Number of loci and SNPs recovered
# (requires the final *_stats.txt file from ipyrad step 7
for num in $(seq $param_range)
do
	echo "$num" >> temp1 && \
	grep "total_filtered_loci" clust_"$num"/"$name_prefix"_"$num"_stats.txt | \
	tr -s " " "\t" | cut -f 4 >> temp2 && \
	grep "snps matrix size" clust_"$num"/"$name_prefix"_"$num"_stats.txt | \
	tr -s " " "\t" | cut -f 5 | tr -d ")," >> temp3
done

paste temp1 temp2 temp3 > "$out_prefix"_loci_snps.tab && rm temp1 temp2 temp3


# 4 Estimated sequencing error rate
# (requires the s4 stats file)
for num in $(seq $param_range)
do
	tail -n +2 clust_"$num"/s4_joint*.txt | \
	tr -s " " "\t" | cut -f 3 | sed "s/^/$num\t/" >> "$out_prefix"_seqerror.tab
done


# Plot them with the R script

Rscript "$MM_script" -p1 "$out_prefix"_paralogs1.tab -p2 "$out_prefix"_paralogs2.tab \
	-h "$out_prefix"_heterozygosity.tab -l "$out_prefix"_loci_snps.tab \
	-e "$out_prefix"_seqerror.tab -o "$out_prefix"




########
# Create a list of the vcf files and submit to the R script for processing
########

if [ -f "$out_prefix"_vcf_list.txt ]; then
	rm "$out_prefix"_vcf_list.txt
fi

for num in $(seq $param_range)
do
	echo -e "${num}\tclust_${num}/${name_prefix}_${num}.vcf" >> "$out_prefix"_vcf_list.txt
done


# create the input samples file from the submitted one
cut -f 1,3,4 "$samp_file" > temp_samp_file.tab

# If wanting to exclude outgroups for this step, can submit a different sample file:
# e.g. grep -v "Z" temp_samp_file.tab > new_sampe_file.tab




# Run the analyses and plot results
echo -e "\nRunning analyses of vcf files to determine variation in PCoAs and relationships between" \
	"missing data and genetic distance, and geographic distance and genetic distance\n"
Rscript "$MM_script" -v "$out_prefix"_vcf_list.txt -o "$out_prefix" -s temp_samp_file.tab




########
# Phylogenetic trees
########

# For this step, we need a set of trees to calculate bootstrap supports
# To get those, first (independently from this script) filter the phylip files for missing data
#
#for num in $(seq $param_range)
#do
#	cd clust_"$num" && \
#	~/scripts/clean_alignment.py -f phylip -p 50 "$name_prefix"_"$num"_tree.phy && \
#	cd ..
#done
#
# Then, run trees on Pawsey with the filtered phylips
#
# Upload data:
#for clust in $(seq $param_range)
#do
#	rsync -vrt clust_"$clust"/*clean.phy \
#	anderson@hpc-data.pawsey.org.au:/scratch/pawsey0220/anderson/MM_test/MMtest/"$name_prefix"_"$clust"_tree_outfiles/
#done
#
# Run 10 trees per parameter value (on supercomputer):
#for clust in $(seq $param_range)
#do
#	for rep in {0..9}
#	do
#		sbatch ~/scripts/iqtree.sbatch -o no \
#		-s /scratch/pawsey0220/anderson/MM_test/MMtest/"$name_prefix"_"$clust"_tree_outfiles/"$name_prefix"_"$clust"_tree_clean.phy \
#		-c "--prefix ${clust}_${rep} --ufboot 1000 -m GTR+I+G"
#	done
#done
#
# Move all the trees and associated output files into a single folder in the ipyrad run directory
#
# Download the data:
# rsync -vrt --progress anderson@hpc-data.pawsey.org.au:/scratch/pawsey0220/anderson/MM_test/iqtrees .
#
# Trees should be named {num}_{rep}.treefile

# If there is a folder of trees, collect the average bootstrap values for plotting, then submit to the R script
if [ "$have_trees" == "True" ]; then
	echo -e "\nCompiling bootstrap supports and plotting\n"
	cd "$tree_folder"
	for num in $(seq $param_range)
	do
		for rep in {0..9}
		do
			avg_support=$(grep -o -e ")[0-9]*:" "$num"_"$rep".treefile | \
			tr -d ")",":" | \
			awk '{n += $1; l += 1}; END{print n / l}')
			echo -e "$num\t$avg_support" >> temp_"$num"_bootstraps.tab
		done
	done
	cat temp*_bootstraps.tab > ../"$out_prefix"_bootstraps.tab && rm temp*_bootstraps.tab
	cd ..
	Rscript "$MM_script" -t "$out_prefix"_bootstraps.tab -o "$out_prefix"
fi



########
# Run Mastretta-Yanes error estimation and Euclidean distances
########


# create the input samples file from the submitted one
cut -f 1,2 "$samp_file" > temp_samp_file.tab

# If wanting to exclude outgroups for this step, can submit a different sample file:
# e.g. grep -v "Z" temp_samp_file.tab > new_sampe_file.tab

echo -e "\nRunning Mastretta-Yanes et al. error estimation (if reps present) and Euclidean distances\n"

# Submit each vcf file to the script, along with the sample file, then move the output to unique named files
for num in $(seq $param_range)
do
	# submit vcf and sample file and run script
	Rscript "$MY_script" -s temp_samp_file.tab \
	-v "clust_${num}/${name_prefix}_${num}.vcf" -o "temp" -r "$reps_present"

	# make output unique
	if [ "$reps_present" == "yes" ]; then
		tail -n +2 temp_error_table.tab | sed "s/^/$num\t/" > "error_table_${num}.tab"
		rm temp_error_table.tab
	fi
	tail -n +2 temp_dist_table.tab | sed "s/^/$num\t/" > "dist_table_${num}.tab"
	rm temp_dist_table.tab

	# also grab the PCoA information
#	for file in temp*pca*.tab
#	do
#		mv "$file" "${file/temp/$num}"
#	done
#	tail -n +2 "${num}_pca_table.tab" | sed "s/^/$num\t/" > "pca_table_${num}.tab"		# omit header
#	rm "${num}_pca_table.tab"
done

# concatenate the files together and submit to the R script for plotting
if [ "$reps_present" == "yes" ]; then
	cat error_table_*.tab > "$out_prefix"_error_tables.tab
	rm error_table_*.tab
fi
cat dist_table_*.tab > "$out_prefix"_dist_tables.tab
#cat pca_table_*.tab > "$out_prefix"_pca_tables.tab
rm dist_table_*.tab #pca_table_*.tab

if [ "$reps_present" == "yes" ]; then
	Rscript "$MM_script" --merror "$out_prefix"_error_tables.tab --mdist "$out_prefix"_dist_tables.tab \
		-o "$out_prefix" -s temp_samp_file.tab #--mpca "$out_prefix"_pca_tables.tab
else
	Rscript "$MM_script" --mdist "$out_prefix"_dist_tables.tab \
		-o "$out_prefix" -s temp_samp_file.tab #--mpca "$out_prefix"_pca_tables.tab
fi

rm temp_samp_file.tab

echo -e "\nFinished at $(date)!\n"
