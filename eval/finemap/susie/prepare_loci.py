import pandas as pd
import os

def check_mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


locus_range = 5e5  # 500kb
str_df = pd.read_csv('data/strs_biallelic_results.csv', sep="\t")
snp_df = pd.read_csv('data/snp_results.csv', sep="\t")
out_path = 'data/finemapping/susie/varlists'
check_mkdir(out_path)
loci = []
for i, row in str_df[str_df.islocmax == 1].iterrows():
    chr = row['#CHROM']
    min, max = row.POS-locus_range, row.POS+locus_range 
    str_range_vars = str_df[(str_df['#CHROM'] == chr) & (str_df.POS >= min) & (str_df.POS <= max)].ID
    outfile = f'{out_path}/{row.ID.replace(":", "_")}_surroundings.strlist'
    str_range_vars.to_csv(outfile, index=False, header=False, sep='\t')
    snp_range_vars = snp_df[(snp_df['#CHROM'] == chr) & (snp_df.POS >= min) & (snp_df.POS <= max)].ID
    outfile = f'{out_path}/{row.ID.replace(":", "_")}_surroundings.snplist'
    snp_range_vars.to_csv(outfile, index=False, header=False, sep='\t')
    locus_df = pd.concat([snp_df[snp_df.ID.isin(snp_range_vars)], str_df[str_df.ID.isin(str_range_vars)]], ignore_index=True)
    locus_df.to_csv(f'{out_path}/{row.ID.replace(":", "_")}_surroundings.tsv', sep='\t', index=False)
    loci.append(row.ID)
    
with open(f'{out_path}/loci.list', 'w') as loci_file:
    for locus in loci:
        loci_file.write(f"{locus}\n")