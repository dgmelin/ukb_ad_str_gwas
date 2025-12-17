# Imputation Quality Evaluation
Here, we want to assess how well the genotypes from imputation (iSTRs) match the popSTR genotypes coming from genotyping (gSTRs).
The main challenge is that we do not have a 1:1 mapping between the variants, as imputation can yield different alleles (e.g., different repeat units) or slightly shifted positions. 
Even if variants match by position, they might not match in terms of alleles.
We chose to evaluate imputation quality in two ways:
1. On the large scale, variants coming from WGS and imputation were matched by genomic position and length of reference and alternative alleles. 
2. On a smaller scale, we manually matched variants passing at least suggestive significance in either gSTRs or iSTRs. These variants were matched with respect to their true allelic correspondence, independent of position or allele representation.

All matching variants that we defined to be matching in either way were then evaluated by comparison of  
a. Beta Values  
b. Allele Frequencies
c. Allele Lengths across samples

## 1. Large-scale evaluation by position and allele matching
Matching of variants is described in the notebook `iSTR_vs_gSTR.ipynb`, section **'match variants'**.
In order to calculate actual allele lengths for the popSTR data, the actual repeat motifs need to be extracted from the raw data. 
The logic for this can be adapted from the `correlate_genotypes.py`.get_geno_lens().

## 2. Manual matching of significant variants
Re manual matching, see Manuscript for now. 

## Comparison Metrics
Beta values and allele frequencies were compared on a variant level. 
Allele lengths were compared per variant across all samples.  

### Correlation of Allele Lengths
The evaluation of genotype concordance is implemented in `run_geno_correlation.sh` and the associated python script `correlate_genotypes.py`.
The main idea is: 
1. Extract genotypes from final data using plinks --export A and either `--snp varname` or `--extract varlist.txt` (one variant per line) for both, the popSTR and the imputed data. 
This will most likely result in two tables that we would 
2. need to merge (using python and pandas) and 
3. to correlate. 
4. Resulting score can be saved in a text file or something. 

Pseudocode:
```bash
while read -r col1 col2 col3 rest; do
    varA="$col2"
    varB="$col3"
    # Do something with $varA and $varB
    plink --snp "$varA" --bfile popSTR_data --recode A --out popSTR_geno
    plink --snp "$varB" --bfile imputed_data --recode A --out imputed_geno
    # Merge and correlate using a Python script
    python correlate_genotypes.py popSTR_geno.raw imputed_geno.raw > correlation_result
done < input.txt
```

