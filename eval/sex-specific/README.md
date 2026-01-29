# Sex Stratification
We ran separate sex-stratified association analyses to investigate potential sex-specific risks. 
We performed additional sex-stratified analyses and a genome-wide interaction analysis using Genotype x Sex as the interaction term.  
In order to run both analyses, we only need slight changes in the command for PLINK2. 
Based on the files after variant qc that we created in the second to last step for the main analyses (see `main/run_group_analysis.sh`), we now run: 
```bash
$in_path=/path/to/varqced/data
$out_path=/path/to/sex_strat/analysis
chr=22
# Sex specific: 
plink2 \
    --pfile ${in_path}/white_british_all_varqced_2 \
    --keep-females \
    --covar ${in_path}/white_british_all.covars \
    --glm --out ${out_path}/chr$chr\_females

plink2 \
    --pfile ${in_path}/white_british_all_varqced_2 \
    --keep-males \
    --covar ${in_path}/white_british_all.covars \
    --glm --out ${out_path}/chr$chr\_males

# sex-interaction
$out_path=/path/to/sex_strat/interaction
plink2 \
    --pfile ${in_path}/white_british_all_varqced_2 \
    --covar ${in_path}/white_british_all.covars \
    --glm interaction --parameters 1-19 \
    --out ${out_path}/chr$chr\_interaction

# Note, take parameters since we only want ADDxSex. 
# Looking at plinks documentation, we can use params 1=ADD, 2-18=sex+seq+pc1-15
# Next would be ADDxcovar1. This is sex, since it is included in the main plink infos
```