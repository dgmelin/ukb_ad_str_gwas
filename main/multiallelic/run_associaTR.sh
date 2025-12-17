#!/bin/bash

default_wrkbase=multiallelics/istrs/custom_filter/ # or whatever the working directory is

usage() {
    echo "Usage: $0 [-c chr] [-w working dir] [-h]"
    exit 1
}

wrk_base=""
chr=""

while getopts ":c:w:h" opt; do
    case $opt in
        c) chr="$OPTARG" ;;
        w) wrk_base="$OPTARG" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done
shift $((OPTIND -1))

tmpdir=$SCRATCH/chr$chr
mkdir -p $tmpdir

if [ -z "$wrk_base" ]; then
    wrk_base=$default_wrkbase
    echo "Info: wrkbase is not set, using default for chr${chr}: $wrk_base"
fi
echo "running in $wrk_base"

infile=$wrk_base/prepare/chr${chr}_annotated.vcf.gz
outdir=${wrk_base}/associatrs/chr${chr}
combine_covars=base/scripts/combine_covars.py # This can be found in helpers

mkdir -p $outdir

echo "extract sample ids for covar-sorting"
bcftools query -l $infile > $tmpdir/sorted.ids

echo "combine covars"
covars=$wrk_base/sample_infos/white_british_all.covars
phenos=$wrk_base/sample_infos/white_british_all.phe

python $combine_covars -v -s -f $tmpdir/sorted.ids -c $covars -p $phenos -o $outdir/${chr}_covars.npy

echo "run association"
associaTR $outdir/$chr.assoc.tsv $infile ad_risk_by_proxy $outdir/${chr}_covars.npy --same-samples --vcftype hipstr

rm $tmpdir/sorted.ids
rmdir $tmpdir
