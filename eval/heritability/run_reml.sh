#!/bin/bash
gcta=tools/gcta-1.94.4-linux-kernel-3-x86_64/gcta64


while getopts ":i:g:o:p:h" opt; do
    case $opt in
        i) infile="$OPTARG" ;;
        c) covars="$OPTARG" ;;
        p) pheno="$OPTARG" ;;
        o) outfile="$OPTARG" ;;
        h) echo "Usage: $(basename "$0") [-i infile] [-g groupfile] [-p pheno] [-o outfile]" >&2; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Missing argument for -$OPTARG" >&2; exit 1 ;;
    esac
done

$gcta --reml \
    --mgrm $infile \
    --pheno $pheno --qcovar $covars.quantitative --covar $covars.discrete \
    --out $outfile --threads 10 --prevalence 0.05 --reml-maxit 10000
