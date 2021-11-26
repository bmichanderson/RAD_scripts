#! /bin/bash

###############
# Author: B. Anderson
# Date: Sep 2021
# Modified: Oct 2021 (simplified filtering and removed F2); Nov 2021 streamlined and made sequential optional
# Description: filter a VCF
# NOTE: some ideas from O'Leary et al. 2018 DOI: 10.1111/mec.14792
# NOTE: requires vcftools in the PATH
# WARNING: it removes some files it creates, so best run when the directory does not have anything named with out_prefix
###############


# Script path (for previous versions where stats plotting was of interest)
#plot_script="/home/benjamin/RAD_scripts/filter_plot.R"


# Set default parameters
mincov="0.9"				# minimum proportion coverage of samples for a SNP
maxm="0.4"					# maximum proportion missing data for a sample
mind="6"					# minimum depth for a genotype
minmd="10"					# minimum mean depth for a locus
maxd="500"					# maximum depth for a genotype
mac="3"						# alleles must occur >= this many times to be kept

# Set the optional iterative filters
iter_snp=(0.2 0.4 0.6 0.8)	# the increasing mincov filters (needs to be the same length as iter_indv)
iter_indv=(0.9 0.7 0.5 0.3)	# the decreasing maxm filters


# Define functions

# function for when the script is called incorrectly or without arguments
usage()
{
	echo -e "This script is for filtering a VCF file\n"
	echo -e "Usage:\tfilter.sh -o out_prefix -v vcf_file -r run_iter" \
			"[--mincov --maxm --mind --minmd --maxd --mac]\n"
	echo -e "Options:\n" \
		"\t-o\tOutput prefix to name files produced during the run [default \"out\"]\n" \
		"\t-v\tVCF file\n" \
		"\t-d\tRun depth filter, \"y\" [default] or \"n\"\n" \
		"\t-r\tRun iterative final filter, \"y\" or \"n\" [default]\n" \
		"\t--mincov\tminimum proportion sample coverage for a SNP [default 0.9]\n" \
		"\t--maxm\tmaximum proportion missing for a sample [default 0.4]\n" \
		"\t--mind\tminimum depth for a genotype call [default 6]\n" \
		"\t--minmd\tminimum mean depth for a locus [default 10]\n" \
		"\t--maxd\tmaximum depth for a genotype call [default 500]\n" \
		"\t--mac\tminimum minor allele count [default 3]\n"
	exit 1
}

if [ $# -eq 0 ]; then		# if no arguments
	usage
fi


# function to calculate stats with vcftools (arguments: vcf, prefix)
calc_stats()
{
    vcftools --vcf $1 --out $2 --depth && \
    vcftools --vcf $1 --out $2 --site-mean-depth && \
    vcftools --vcf $1 --out $2 --missing-indv && \
    vcftools --vcf $1 --out $2 --missing-site && \
    vcftools --vcf $1 --out $2 --het
}


# function to filter by depth (arguments: vcf, prefix, mind, minmd, maxd, mac)
filter_depth()
{
    vcftools --vcf $1 --out "$2"dp --minDP $3 --min-meanDP $4 --maxDP $5 --recode --recode-INFO-all && \
    #calc_stats "$2"dp.recode.vcf "$2"dp && \
    #Rscript $plot_script -p "$2"dp |& tee -a "$out_prefix"_log.txt && \
    vcftools --vcf "$2"dp.recode.vcf --out "$2"dpmac --mac $6 --recode --recode-INFO-all && \
    rm "$2"dp.*
    #calc_stats "$2"dpmac.recode.vcf "$2"dpmac && \
    #Rscript $plot_script -p "$2"dpmac |& tee -a "$out_prefix"_log.txt
}


# function for filtering SNPs by sample cover (arguments: vcf, prefix, mincov)
filter_snp()
{
    vcftools --vcf $1 --out $2 --max-missing $3 --recode --recode-INFO-all
}


# function for filtering samples by missingness (arguments: vcf, prefix, maxm)
filter_indv()
{
    vcftools --vcf $1 --out $2 --missing-indv && \
    awk -F"\t" -v max="$3" 'NR > 1 && $5 >= max' "$2".imiss | cut -f 1 > "$2".to_remove && \
    vcftools --vcf $1 --out $2 --remove "$2".to_remove --recode --recode-INFO-all && \
	rm "$2".imiss "$2".log "$2".to_remove
}


# Parse arguments
FILES=()
while [[ $# -gt 0 ]]
do
key="$1"
case $key in
	-o)
	out_prefix="$2"
	shift
	shift
	;;
	-v)
	vcf_file="$2"
	shift
	shift
	;;
	-d)
	run_depth="$2"
	shift
	shift
	;;
	-r)
	run_iter="$2"
	shift
	shift
	;;
	--mincov)
	mincov="$2"
	shift
	shift
	;;
	--maxm)
	maxm="$2"
	shift
	shift
	;;
	--mind)
	mind="$2"
	shift
	shift
	;;
	--minmd)
	minmd="$2"
	shift
	shift
	;;
	--maxd)
	maxd="$2"
	shift
	shift
	;;
	--mac)
	mac="$2"
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
if [ -z "$vcf_file" ]; then
	usage
fi
if [ -z "$out_prefix" ]; then
	out_prefix="out"
fi
if [ -z "$run_depth" ]; then
	run_depth="y"
fi
if [ -z "$run_iter" ]; then
	run_iter="n"
fi


# Run the filtering

# could insert a stats calculation at some point
#if [ "$run_stats" == "y" ]; then
	#calc_stats $vcf_file $out_prefix && \
	#Rscript $plot_script -p $out_prefix |& tee -a "$out_prefix"_log.txt && \
	#rm "$out_prefix".*
#fi

if [ "$run_depth" == "n" ]; then
	# F1 straight to strict filter, SNP first (no depth)
	filter_snp $vcf_file "$out_prefix"_snp $mincov && \
	filter_indv "$out_prefix"_snp.recode.vcf "$out_prefix"_snp_indv $maxm && \
	#calc_stats "$out_prefix"_snp_indv.recode.vcf "$out_prefix"_F1 && \
	#Rscript $plot_script -p "$out_prefix"_F1 |& tee -a "$out_prefix"_log.txt && \
	mv "$out_prefix"_snp_indv.recode.vcf "$out_prefix"_F1_out.vcf && \
	rm "$out_prefix"_snp* #"$out_prefix"_F1.*
else
	# F3 depth first, then to strict filter, SNP first
	filter_depth $vcf_file "$out_prefix" $mind $minmd $maxd $mac && \
	filter_snp "$out_prefix"dpmac.recode.vcf "$out_prefix"_depth_snp $mincov && \
	filter_indv "$out_prefix"_depth_snp.recode.vcf "$out_prefix"_depth_snp_indv $maxm && \
	#calc_stats "$out_prefix"_depth_snp_indv.recode.vcf "$out_prefix"_F3 && \
	#Rscript $plot_script -p "$out_prefix"_F3 |& tee -a "$out_prefix"_log.txt && \
	mv "$out_prefix"_depth_snp_indv.recode.vcf "$out_prefix"_F3_out.vcf && \
	rm "$out_prefix"_depth_snp* #"$out_prefix"_F3.*
fi

# run an iterative filter if requested
if [ "$run_iter" == "y" ]; then
	if [ "$run_depth" == "y" ]; then
		mv "$out_prefix"dpmac.recode.vcf temp_in.vcf && \
		rm "$out_prefix"dpmac.*
	else
		cp $vcf_file temp_in.vcf
	fi
	# F4 depth first (or not), then iterative filters, then to strict filter, indv first
	for iter in $(seq 1 ${#iter_snp[@]})
	do
		filter_snp temp_in.vcf "$out_prefix"_"$iter"_snp ${iter_snp[iter - 1]} && \
		filter_indv "$out_prefix"_"$iter"_snp.recode.vcf "$out_prefix"_"$iter"_snp_indv ${iter_indv[iter - 1]} && \
		mv "$out_prefix"_"$iter"_snp_indv.recode.vcf temp_in.vcf
	done
	rm "$out_prefix"_*_snp*
	filter_indv temp_in.vcf "$out_prefix"_iter_indv $maxm && \
	filter_snp "$out_prefix"_iter_indv.recode.vcf "$out_prefix"_iter_indv_snp $mincov && \
	#calc_stats "$out_prefix"_iter_indv_snp.recode.vcf "$out_prefix"_F4 && \
	#Rscript $plot_script -p "$out_prefix"_F4 |& tee -a "$out_prefix"_log.txt && \
	mv "$out_prefix"_iter_indv_snp.recode.vcf "$out_prefix"_F4_out.vcf && \
	rm "$out_prefix"_iter_indv* temp_in.vcf #"$out_prefix"_F4.* 
elif [ "$run_depth" == "y" ]; then
	rm "$out_prefix"dpmac.*	
fi
