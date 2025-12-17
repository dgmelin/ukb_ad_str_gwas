#!/bin/bash

chr=$1

tmpdir=tmp/chr${chr}
mkdir -p $tmpdir

data_in=/path/to/imputed/data
# eg data_in=data/imputed/full/$chr\_merged.vcf.gz
wrk_base=multiallelics/istrs/custom_filter/ # or whatever the working directory is
outdir=${wrk_base}/detailled/prepare
mkdir -p $outdir

# Script that annotates the regions correctly
annotate_script=$wrk_base/annotate_imputed.py
# path to regions to be extracted
# in our case this was 
regions=$wrk_base/regions/chr${chr}_regions.bed

# first extract the regions from vcf
if [ ! -f "${data_in}.tbi" ]; then
    echo "Indexing VCF for chromosome $chr"
    tabix -p vcf "$data_in"
fi

echo "Extracting regions for chromosome $chr"
bcftools view -R $regions $data_in -Oz -o $outdir/chr${chr}_extracted.vcf.gz

echo "Extracting multiallelic and biallelic variants for chromosome $chr"
bcftools view -m3 $data_in -Oz -o $outdir/chr${chr}_multiallelic.vcf.gz
bcftools view -m2 -M2 $data_in -Oz -o $outdir/chr${chr}_biallelic.vcf.gz


# prefer multiallelic over biallelic records
echo "Filtering biallelic variants for chromosome $chr"
bcftools query -f'%CHROM\t%POS\n' $outdir/chr${chr}_multiallelic.vcf.gz | sort -u > $tmpdir/multi_sites.txt

bcftools view -T ^$tmpdir/multi_sites.txt $outdir/chr${chr}_biallelic.vcf.gz -Oz -o $outdir/chr${chr}_biallelic_filtered.vcf.gz
# 
# rejoin specifically
echo "Rejoin biallelic variants for chromosome $chr"
bcftools norm -m+any -Oz -o $outdir/chr${chr}_biallelic_rejoined.vcf.gz $outdir/chr${chr}_biallelic_filtered.vcf.gz
echo "index rejoined"
bcftools index $outdir/chr${chr}_multiallelic.vcf.gz &
bcftools index $outdir/chr${chr}_biallelic_rejoined.vcf.gz

wait

echo "merging filtered biallelic variants and mutliallelics for chromosome $chr"
bcftools concat -a -Oz -o $outdir/chr${chr}_combined.vcf.gz $outdir/chr${chr}_biallelic_rejoined.vcf.gz $outdir/chr${chr}_multiallelic.vcf.gz

# prepare using my script
echo "Annotating merged VCF for chromosome $chr"
python $annotate_script $outdir/chr${chr}_combined.vcf.gz $outdir/chr${chr}_annotated.vcf.gz

echo "Preparation for chromosome $chr completed"