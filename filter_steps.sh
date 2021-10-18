#! /bin/bash

###############
# Author: B. Anderson
# Date: Sep 2021
# Modified: Oct 2021 (simplified filtering and removed F2)
# Description: run a successive series of filtering steps on VCF files
# NOTE: based on O'Leary et al. 2018 DOI: 10.1111/mec.14792
# NOTE: it requires R scripts, as well as vcftools in the PATH
# WARNING: it removes some files it creates, so best run when the directory does not have anything named with out_prefix
###############


# Script paths
plot_script="/home/benjamin/RAD_scripts/filter_plot.R"
indv_miss_script="/home/benjamin/RAD_scripts/indv_miss.R"


# Set parameters
min_cover="0.9"            # minimum proportion coverage of samples for a SNP
max_indv_miss="0.2"        # maximum proportion missing data for a sample
minDP="6"                  # minimum depth for a genotype
min_meanDP="10"            # minimum mean depth for a locus (across *all* individuals)
minor_count="3"            # alleles must occur >= this many times to be kept
iter_snp=(0.2 0.4 0.6 0.8)     # the progressive filters (needs to be the same length as iter_indv)
iter_indv=(0.9 0.7 0.5 0.3)


# function for when the script is called incorrectly or without arguments
usage()
{
	echo -e "This script is for running a series of filtering steps on a VCF file\n"
	echo -e "Usage:\tfilter_steps.sh -o out_prefix -r run_setting -v vcf_file\n"
	echo -e "Options:\n" \
		"\t-o\tOutput prefix to name files produced during the run (default \"raw\")\n" \
		"\t-r\tThe run setting:\t1 (default) = generate stats for just the raw,\n" \
        "\t\t\t\t\t2 = successively filter and generate stats\n" \
		"\t-v\tThe raw VCF file\n"
	exit 1
}

if [ $# -eq 0 ]; then		# if no arguments
	usage
fi


# function to calculate stats with vcftools (arguments: vcf prefix)
calc_stats()
{
    vcftools --vcf $1 --out $2 --depth && \
    vcftools --vcf $1 --out $2 --site-mean-depth && \
    vcftools --vcf $1 --out $2 --missing-indv && \
    vcftools --vcf $1 --out $2 --missing-site && \
    vcftools --vcf $1 --out $2 --het
}


# function for initial filtering for depth (arguments: vcf prefix minDP min_meanDP minor_count)
filter_depth()
{
    vcftools --vcf $1 --out "$2"mDPmmDP --minDP $3 --min-meanDP $4 --recode --recode-INFO-all && \
    calc_stats "$2"mDPmmDP.recode.vcf "$2"mDPmmDP && \
    Rscript $plot_script -p "$2"mDPmmDP |& tee -a "$out_prefix"_log.txt && \
    vcftools --vcf "$2"mDPmmDP.recode.vcf --out "$2"mDPmmDPmac --mac $5 --recode --recode-INFO-all && \
    rm "$2"mDPmmDP.* && \
    calc_stats "$2"mDPmmDPmac.recode.vcf "$2"mDPmmDPmac && \
    Rscript $plot_script -p "$2"mDPmmDPmac |& tee -a "$out_prefix"_log.txt
}


# function for filtering by SNP cover (arguments vcf prefix mincov)
filter_snp()
{
    vcftools --vcf $1 --out $2 --max-missing $3 --recode --recode-INFO-all
}


# function for filtering by individual missingness (arguments vcf prefix maxindvmiss)
filter_indv()
{
    vcftools --vcf $1 --out $2 --missing-indv && \
    Rscript $indv_miss_script -l $3 -p $2 && \
    vcftools --vcf $1 --out $2 --remove "$2".to_remove --recode --recode-INFO-all
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
	-r)
	run_setting="$2"
	shift
	shift
	;;
	-v)
	vcf_file="$2"
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
	out_prefix="raw"
fi
if [ -z "$run_setting" ]; then
    run_setting="1"
    echo -e "Run setting not specified, so just summarizing stats for raw VCF\n"
elif [ "$run_setting" != "2" ] && [ "$run_setting" != "1" ]; then
    usage
fi


## Run the filtering
# calculate stats for raw
if [ "$run_setting" == "1" ]; then
    calc_stats $vcf_file $out_prefix && \
    Rscript $plot_script -p $out_prefix |& tee -a "$out_prefix"_log.txt && \
    rm "$out_prefix".*
fi
# run further filtering if requested
if [ "$run_setting" == "2" ]; then
    # F1 straight to strict filter, SNP first
    filter_snp $vcf_file "$out_prefix"_snp $min_cover && \
    filter_indv "$out_prefix"_snp.recode.vcf "$out_prefix"_snp_indv $max_indv_miss && \
    calc_stats "$out_prefix"_snp_indv.recode.vcf "$out_prefix"_F1 && \
    Rscript $plot_script -p "$out_prefix"_F1 |& tee -a "$out_prefix"_log.txt && \
    mv "$out_prefix"_snp_indv.recode.vcf "$out_prefix"_F1_out.vcf && \
    rm "$out_prefix"_snp* "$out_prefix"_F1.*

    # F3 depth first, then to strict filter, SNP first
    filter_depth $vcf_file "$out_prefix" $minDP $min_meanDP $minor_count && \
    filter_snp "$out_prefix"mDPmmDPmac.recode.vcf "$out_prefix"_depth_snp $min_cover && \
    filter_indv "$out_prefix"_depth_snp.recode.vcf "$out_prefix"_depth_snp_indv $max_indv_miss && \
    calc_stats "$out_prefix"_depth_snp_indv.recode.vcf "$out_prefix"_F3 && \
    Rscript $plot_script -p "$out_prefix"_F3 |& tee -a "$out_prefix"_log.txt && \
    mv "$out_prefix"_depth_snp_indv.recode.vcf "$out_prefix"_F3_out.vcf && \
    rm "$out_prefix"_depth_snp* "$out_prefix"_F3.*

    # F4 depth first (already done), then to iterative filters, then to strict filter, indv first
    mv "$out_prefix"mDPmmDPmac.recode.vcf temp_in.vcf && \
    rm "$out_prefix"mDPmmDPmac.*
    for iter in $(seq 1 ${#iter_snp[@]})
    do
        filter_snp temp_in.vcf "$out_prefix"_"$iter"_snp ${iter_snp[iter - 1]} && \
        filter_indv "$out_prefix"_"$iter"_snp.recode.vcf "$out_prefix"_"$iter"_snp_indv ${iter_indv[iter - 1]} && \
        mv "$out_prefix"_"$iter"_snp_indv.recode.vcf temp_in.vcf
    done
    rm "$out_prefix"_*_snp*
    filter_indv temp_in.vcf "$out_prefix"_iter_indv $max_indv_miss && \
    filter_snp "$out_prefix"_iter_indv.recode.vcf "$out_prefix"_iter_indv_snp $min_cover && \
    calc_stats "$out_prefix"_iter_indv_snp.recode.vcf "$out_prefix"_F4 && \
    Rscript $plot_script -p "$out_prefix"_F4 |& tee -a "$out_prefix"_log.txt && \
    mv "$out_prefix"_iter_indv_snp.recode.vcf "$out_prefix"_F4_out.vcf && \
    rm "$out_prefix"_iter_indv* "$out_prefix"_F4.* temp_in.vcf
fi
