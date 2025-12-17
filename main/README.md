# Main GWAS Analysis
This firectory contains the file(s) to run our main association analysis. 
Starting point are plink files after imputation (-QC). 
Variant QC steps are not necessary as they will be performed along the way. 
**Requiered Files**: 
- plink-files (*.pvar, *.pgen, *.psam) containing genotype data
- info/phenotype file containing infos about sex and phenotype for each sample
- covars file containing the covariates (currently must be placed at out)

Runs can be parallelized, I ran them on a chromosome level. 
To run for all samples on the white_british - cohort across all chrs, run eg: 


```bash
infos=/path/to/sample/infos
source=/path/to/data
out_base=/path/to/output/dir
cohort=White_british
for chr in {1..22}; do
    out_dir=${out}/$chr
    mkdir -p $out_dir; 
    cp $infos/$cohort.covars $out_dir
    bash run_group_analysis.sh \
        -o $out_dir/$cohort \
        -i $infos/$cohort.infos \
        -p $source/$cohort.$chr \
        -m 6000 -t 1; 
done
```

After all runs finished, results were combined simply with: 

```bash
head -n 1 $out_base/1/$cohort\_assoc.ad_risk_by_proxy.glm.linear > white_british_combined_chr_all.glm.linear; 
for chr in {1..22}; do
    tail -n +2 $out_base/$chr/$cohort\_assoc.ad_risk_by_proxy.glm.linear | grep ADD > white_british_combined_chr_all.glm.linear; 
done
```


