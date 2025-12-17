#!/bin/bash

key=$1

if [[ -z "$key" ]]; then
    echo "Usage: $0 {strs|snps|combined}" >&2
    exit 1
fi

case "$key" in
    strs|snps|combined) ;;
    *)
        echo "Error: key must be one of strs, snps, combined" >&2
        exit 1
        ;;
esac

plink_ressources="--memory 32000 --threads 2"
basedir=data/heritability/prep/filtered/$key
outfile=$basedir/all_merged

for chr in {1..22}; do
    infile=$basedir/${chr}_${key}_filtered
    echo $infile >> $basedir/merge_list.txt
done

infile=$basedir/merge_list.txt

plink --merge-list $infile --make-bed --out $outfile $plink_ressources
