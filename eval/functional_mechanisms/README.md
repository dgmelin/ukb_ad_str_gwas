# Functional mechanisms
We ran eQTL, meQTL and eQTM analyses to indetify potential functional mechanisms and to check for potenital associations for our top STR hits with gene expression. 

## meQTL and eQTL analysis
Code for meQTL and eQTL is mostly identical since both is using the [Matrix eQTL package](https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html). Please also refer to the examples there. 
The analyses mostly differ only in the input files for expression data. 
To run qtl analyses in bash, do

```bash
Rscript runQTL.R \
    str_genotypes.txt \ % str variants genotype
    str_loc.txt \ % variants location information
    expr.txt \ % expression (meQTL or RNA)
    gene_loc.txt \ % gene location information
    covars.txt \ % covairates file
    out_cis.txt \ % output file for cis
    out_tra.txt % output for trans
```

## QTM
Please refer to the `eQTM_mapping_OPTIMA.R` script for further infos. 


## Credits
Special thanks to 
- Olena Ohlei and M. Muaaz Aslam for running the analyses and providing the code for the QTL analyses.  
- Marit P. Junge for running the QTM analyses and for providing the respective code. 