#!/bin/bash

while getopts ":i:o:h" opt; do
    case $opt in
        i) infile="$OPTARG" ;;
        o) outfile="$OPTARG" ;;
        h) echo "Usage: $(basename "$0") [-i infile] [-o outfile]" >&2; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Missing argument for -$OPTARG" >&2; exit 1 ;;
    esac
done

#point this to your gcta binary
gcta=tools/gcta-1.94.4-linux-kernel-3-x86_64/gcta64

outdir=$(dirname "$outfile")
mkdir -p $outdir

$gcta --bfile $infile --ld-score-region 200 --out $outfile --threads 10


