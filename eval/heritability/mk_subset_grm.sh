#!/bin/bash

while getopts ":i:g:o:h" opt; do
    case $opt in
        i) infile="$OPTARG" ;;
        g) groupfile="$OPTARG" ;;
        o) outfile="$OPTARG" ;;
        h) echo "Usage: $(basename "$0") [-i infile] [-g groupfile] [-o outfile]" >&2; exit 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Missing argument for -$OPTARG" >&2; exit 1 ;;
    esac
done

group_id=$1

gcta=tools/gcta-1.94.4-linux-kernel-3-x86_64/gcta64

echo "running group $groupfile for $infile."

$gcta --bfile $infile --extract $groupfile  --make-grm --out $outfile --threads 10
