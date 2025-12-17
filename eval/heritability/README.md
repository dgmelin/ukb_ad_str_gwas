# Heritability analysis
To assess the contribution of STRs to heritability of our data, we ran [GCTA-LDMS-GREML](https://yanglab.westlake.edu.cn/software/gcta/#GREMLinWGSorimputeddata) analysis.
Here, we mainly follow the steps provided by the authors. 
For computational reasons, we ran the analysis in multiple batches since calculation of the pairwise GRM relationship matrix for all samples together is estimated to exceed reasonalble memory limits. 
We also did not include samples with proxy-risk for AD in the analysis but diagnosed cases and controls (i.e. risk score $< 1$, no affected parents) only.
All 2,947 AD cases were included, along with 15,000 randomly selected controls.
The controls were split into 16 non-overlapping batches. The remaining controls were put into a 17th batch, which was filled up with already used controls to ensure all batches are the same size.

We run heritability analysis for SNPs and STRs alone and a combined SNP-STR dataset.  

## Preparation
To be able to merge SNP and STR data using plink, they need to be converted into plink-1 fomat (bed, bim, fam). 

### Steps for GCTA LDMS-GREML
Exemplary for the `combined` data on any batch (`$batch_id`)
#### 1. Filter chromosomes for extra-chrs and extract samples using `run_format_chrs.sh`
This can be done via `bash run_format_chrs.sh $chr`. 
In the script, file paths need to be adapted accordingly. 

#### 2. Merge Chromosomes using `run_merge_formatted.sh`
```bash
bash run_merge_formatted.sh combined
```

#### 3. Calculate segment based LD-Score using `run_format_ldms.sh`
```bash
bash run_format_ldms.sh \ 
    -i data/batch$batch_id/combined/all_merged \
    -o data/batch$batch_id/combined/ldscored.score.ld
```

#### 4. Stratify vars using `strat_snps.R`
This script is adapted from the GCTA-Website ([link](https://yanglab.westlake.edu.cn/software/gcta/#GREMLinWGSorimputeddata), Step 2, Option 1).
```bash
Rscript strat_snps.R data/batch$batch_id/combined/ldscored.score.ld
```

#### 5. Crate GRMs for stratified vars `mk_subset_grm.sh`
Run this via 
```bash 
bash run_format_ldms.sh \
    -i data/batch$batch_id/combined/all_merged
    -g data/batch$batch_id/combined/ld_group1\
    -o data/batch$batch_id/combined/ldgrm_group1
```

#### 6. run final analysis with `run_reml.sh`
```bash
infile=data/batch$batch_id/combined/multi_GRMs.txt
# prepare all grm files
for i in {1..4}; 
do 
    readlink -f data/batch$batch_id/combined/ldgrm_group$i
done > $infile

bash run_reml.sh \
    -i $infile \
    -o heritability/results/combined_$batch_id
    -p data/infos/white_british.pheno \
    -c data/infos/white_british.covars 

```