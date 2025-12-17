# Conditioning Analysis
We wanted to assess the effects of SNP results on STRs and vice versa. 
Here, we repeat association studies for STR and SNP regions, while additionaly conditioning for the respective leading variant int he user study. 
In detail: 
1. We adjusted the STR signal for the effects of the respective lead SNP (±500kb) per STR locus derived from the SNP-based AD GWAS by Bellenguez et al. 
2. We performed reciprocal analyses where the SNP association statistics (our results) in any of the top STR loci were conditioned on leading STR genotype. 

# Steps
1. Extract SNP-STR variant pairs from our STR-GWAS results and the Bellenguez Summary Stats. 
> Check the preparation section of `condition.ipynb`
2. Run actual conditioning analysis
> Run the `run_condition.sh` script. 
> This can be done using, for example
```bash
for ln in $(tail snp_str.ids); do
    # Run Conditioning STRs on SNPs
    snp_id=$(cut -f 3)
    # extract chr from SNP id
    chr=$(echo $id | cut -f 1 -d'_' )
    out=condition/strs/$id
    #select source from where to extract the variant
    source=/path/to/snp/data
    bash $run_script -i $id -c $chr -a STR -s $source -o $out -m 64000 -t 2
    
    # Same for SNPs on STRs
    str_id = $(cut -f 1)
    # extract chr from SNP id
    chr=$(echo $id | cut -f 1 -d'_' )
    out=condition/strs/$id
    #select source from where to extract the variant
    source=/path/to/STR/data
    bash $run_script -i $id -c $chr -a SNP -s $source -o $out -m 64000 -t 2
```
3. Evaluate and Plot results
> Check the eval section of `condition.ipynb`