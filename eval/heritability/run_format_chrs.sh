#!/bin/bash

chr=$1

plink_ressources="--memory 32000 --threads 2"

workdir=data/heritability
ids_file=$workdir/final_ids.formatted
infile=$workdir/merged_beds/white_british_${chr}_snpfinal

outdir=$workdir/prep/filtered/

# run SNPs
outfile=$outdir/snps/${chr}_snps_filtered
mkdir -p $outdir/snps

plink --bfile $infile \
    --keep $ids_file \
    --allow-extra-chr \
    --chr ${chr} \
    --make-bed \
    $plink_ressources \
    --out $outfile

# run STRs
infile=$workdir/merged_beds/white_british_${chr}_strfinal
outfile=$outdir/strs/${chr}_strs_filtered

mkdir -p $outdir/strs
plink --bfile $infile \
    --keep $ids_file \
    --allow-extra-chr \
    --chr ${chr} \
    --make-bed \
    $plink_ressources \
    --out $outfile

# run combined
infile=$workdir/merged_beds/${chr}_merged
outfile=$outdir/combined/${chr}_combined_filtered
mkdir -p $outdir/combined

plink --bfile $infile \
    --keep $ids_file \
    --allow-extra-chr \
    --chr ${chr} \
    --make-bed \
    $plink_ressources \
    --out $outfile


