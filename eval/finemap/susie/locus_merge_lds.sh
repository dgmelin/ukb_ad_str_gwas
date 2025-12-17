#!/bin/bash
id=$1

plink_resources="--memory 32000 --threads 1 --seed 42"

chr=$(echo $id | cut -f1 -d_ | sed 's/chr//')

echo "running for locus $id on chr $chr"

#preparation
base=data/finemapping/susie/ld_matrices/
out=$base/out/$id
mkdir -p $out

#take variant data that was used to generate association results
snp_in=data/snps/$chr/white_british_all_varqced_2
str_in=data/strs/$chr/white_british_all_varqced_2

varfile=$base/varlists/$id\_surroundings
var_id=$(echo $id | sed 's/_/:/g')

tmp=tmp/$id
mkdir $tmp

# extract vars around chr
plink2 --pfile $snp_in --allow-extra-chr --extract $varfile.snplist --make-bed --out $tmp/snps_extracted
plink2 --pfile $str_in --allow-extra-chr --extract $varfile.strlist --make-bed --out $tmp/strs_extracted

#merge
plink --bfile $tmp/strs_extracted \
  --bmerge $tmp/snps_extracted \
  --allow-extra-chr \
  --make-bed \
  --out $out/merged \
    $plink_resources

# calulate lds
plink --bfile $out/merged \
  --allow-extra-chr \
  --r square \
  --write-snplist \
  --out $out/$id\_merged_lds \
    $plink_resources

plink --bfile $out/merged \
  --allow-extra-chr \
  --r2 square \
  --out $out/$id\_merged_lds2 \
    $plink_resources

plink --bfile $out/merged \
  --allow-extra-chr \
  --ld-snp $var_id \
  --ld-window 500e6 \
  --ld-window-kb 500 \
  --ld-window-r2 1e-30 \
  --r --out $out/$id\_single_r \
  $plink_resources

cat $out/$id\_single_r.ld | sed -E 's/^[[:space:]]+//g' | sed -E 's/[[:space:]]+/\t/g' > $out/$id\_single_r.ld.tsv

plink --bfile $tmp/snps_extracted \
  --allow-extra-chr \
  --r square \
  --out $out/$id\_snps_lds \
    $plink_resources

plink --bfile $tmp/strs_extracted \
  --allow-extra-chr \
  --r square \
  --out $out/$id\_strs_lds \
    $plink_resources

# cleanup
rm -rf $tmp