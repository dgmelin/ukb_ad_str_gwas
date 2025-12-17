# Running COJO Analysis
## tools 
- plink V1 and V2
- gcta cojo: https://yanglab.westlake.edu.cn/software/gcta/#COJO


## Steps
calculate frequencies for the respective datasets

```bash
snp_data=/path/to/snp/data
str_data=/path/to/str/data
outpath=/path/to/output/

mkdir -p $outpath

for chr in {1..22}; 
do
    plink2 --pfile $snp_data/${chr}/white_british_all_varqced_2 \
    --allow-extra-chr \
    --freq cols=ref,alt,freq,nobs \
    --out $outpath/$chr\_snp_freqs

    plink2 --pfile $str_data/${chr}/white_british_all_varqced_2 \
    --allow-extra-chr \
    --freq cols=ref,alt,freq,nobs \
    --out $outpath/$chr\_str_freqs
done
```

Also, convert plinks pfiles to bfiles to be able to merge concat SNPs and STRs into a combined dataset: 
```bash
#process SNPs
plink2 --pfile $snp_data/${chr}/white_british_all_varqced_2 \
--allow-extra-chr --make-bed --out $outpath/$chr\_snp_convert

#process STRs
plink2 --pfile $str_data/${chr}/white_british_all_varqced_2 \
--allow-extra-chr --make-bed --out $outpath/$chr\_str_convert

#combine using plink1
plink --bfile $outpath/$chr\_snp_convert --allow-extra-chr --bmerge $outpath/$chr\_snp_convert --make-bed --out $outpath/$chr\_merged 
```

Next, combine summary stats and allele frequencies into COJO-format using `to_cojo.awk`

```bash
snp_freqs=$outpath/${chr}_snp_freqs.afreq
str_freqs=$outpath/${chr}_str_freqs.afreq
$snp_assocs=$snp_data/$chr/$cohort\_assoc.ad_risk_by_proxy.glm.linear
$str_assocs=$str_data/$chr/$cohort\_assoc.ad_risk_by_proxy.glm.linear

printf "snp\t"
awk -f to_cojo.awk $snp_freqs $snp_assocs > $outpath/${chr}_snp.ma
printf "str\t"
awk -f to_cojo.awk $str_freqs $str_assocs > $outpath/${chr}_str.ma
printf "fin\n"  

cat $outpath/${chr}_str.ma > $outpath/${chr}_combined.ma
tail -n +2 $outpath/${chr}_snp.ma >> $outpath/${chr}_combined.ma
```

...and finally, run the actual cojo: 
```bash
mkdir -p $outpath/results
gcta64 --bfile $outpath/${chr}_merged \
    --chr $chr --maf 0.01 --cojo-slct --cojo-p 5e-5 \
    --cojo-file $outpath/${chr}_combined.ma --out $outpath/results/${chr}_cojo
```