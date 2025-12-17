import pandas as pd
import numpy as np

def get_maxes_idxs(df, winsize=250000):
    keep=[]
    for i, r in df.iterrows():
        chr = r['#CHROM']
        pos = r.POS
        cands = df[(df['#CHROM'] == chr) 
                & (abs(df.POS-pos)<winsize)]
        min_idx = cands.P.idxmin()
        if min_idx == i: 
            keep.append(i)

    return keep

biallelic_df = pd.read_csv('/path/to/imputed_str_biallelic.csv', sep="\t").astype({'POS': int, '#CHROM': int})
max_idxs = get_maxes_idxs(biallelic_df[biallelic_df.P<1e-5], winsize=250000)
biallelic_df['islocmax'] = np.where(biallelic_df.index.isin(max_idxs), 1, 0)
mask = biallelic_df.islocmax == 1
hits = biallelic_df[mask][['#CHROM', 'POS', 'P']].sort_values(by=['#CHROM', 'POS']).reset_index(drop=True)
margin = 500000
hits['#CHROM'] = 'chr' + hits['#CHROM'].astype(str)
hits['upper'] = hits['POS'] + margin
hits['lower'] = hits['POS'] - margin
hits[['#CHROM', 'lower', 'upper']].to_csv('path/to/data/multiallelic/extract_regions.bed', sep="\t", header=False, index=False)