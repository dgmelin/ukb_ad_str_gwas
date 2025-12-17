# Susie - Analysis
Here, we run SusieR analysis on combined SNP and STR summary statistics. 
## Prep
First, we extracted variants +- 500kb around the leading str variant positon at each locus for both, SNPs and STRs using `prepare_loci.py`.

Next, we run `locus_merge_lds.sh` to calculate ld-matrices for each locus on snp data, str data and combined data.

## Main
To run susie, we use the `run_susie.R` script. This iterates across all loci in the leading loci list created by the python script. Based on the ld-matrices created before it runs [susie_rss](https://stephenslab.github.io/susieR/reference/susie_rss.html) and creates the respective pip values for all variants. 
It also creates some initial plots, similar to the ones improved using python/seaborn in the next step. 

## Eval
We run the script `eval_susie.py` to create final plots and tables for all loci. As I'm at most a R-beginner, I prefer to run this in python. 