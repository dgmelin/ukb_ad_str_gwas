import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


data_path = 'data/finemapping/susie'
out_path = f'{data_path}/eval'

gw_t = 0.05/334605 

def check_mkdir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)

color_palette = sns.color_palette("tab10", 10)
color_dict = {i+1: color_palette[i] for i in range(10)}
style_dict = {'STR' : '*', 'SNP' : 'o'}
ld_palette = sns.color_palette("viridis", as_cmap=True)

def run():
    varpath = f'{data_path}/varlists'

    with open(f'{varpath}/loci.list', 'r') as loci_file:
        loci = [line.strip() for line in loci_file.readlines()]

    check_mkdir(f'{out_path}')

    for locus in loci:

        l = locus.replace(":", "_")
        locus_file = f'{data_path}/finemapping_results/{l}_susie_outputs.tsv'
        if os.path.exists(locus_file):
            locus_df = pd.read_csv(
                locus_file, sep="\t" ).drop(columns='islocmax').replace(np.nan, 0)
        else:
            print(f'skipping {l} as it does not exist')
            continue

        print(f'Processing {locus}')
        locus_df = pd.read_csv(
            f'{data_path}/finemapping_results/{l}_susie_outputs.tsv', sep="\t" ).drop(columns='islocmax').replace(np.nan, 0)

        locus_df['R2 to lead variant'] = locus_df['ld_to_lead_matrix'] ** 2
        locus_df = locus_df.drop(columns=['ld_to_lead_matrix'])

        fig ,axes = plt.subplots(nrows=2,sharex=True,figsize=(15,7),height_ratios=(4,1))
        col_to_plot = "-log10"
        sns.scatterplot(data=locus_df, x="POS", y=col_to_plot, ax=axes[0], hue='R2 to lead variant', s=50, edgecolors=None, style='source', palette=ld_palette, markers=style_dict)

        sns.scatterplot(data=locus_df.loc[locus_df["cs"]>0], x="POS", y=col_to_plot, ax=axes[0],
                    facecolors="none", edgecolors="red", marker="o", s=80, linewidth=1, label="Variants in credible sets")

        axes[0].set_xlabel("position")
        axes[0].set_ylabel(col_to_plot)
        axes[0].legend()
        axes[0].set_title(f'Susie results for {locus}')
        sns.scatterplot(locus_df, x='POS', y='pip', ax=axes[1], hue='cs', palette=color_dict, style='source', markers=style_dict, legend=True, edgecolors='black')

        axes[1].set_title('Posterior Incusion Probabilities (PIP)')
        lgd = axes[0].legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.01, fontsize=13)
        lgd2 = axes[1].legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.01, fontsize=13)

        out_suf = 'suggestive' if locus_df.P.min() > gw_t else 'genomewide'
        check_mkdir(f'{out_path}/{out_suf}')

        fig.savefig(f'{out_path}/{out_suf}/{l}_susie_finemap.png', dpi=300, bbox_inches='tight', bbox_extra_artists=(lgd,lgd2))
        plt.close()

        locus_df['locus'] = l*locus_df.shape[0]

        cs_hits = locus_df[locus_df.cs > 0]
        if cs_hits.shape[0] > 0:
            cs_hits.sort_values(by="cs").to_csv(
                f'{out_path}/{out_suf}/{l}_susie_cred_sets.tsv', sep="\t", index=False)
        else: 
            best_hits = cs_hits.sort_values(by="pip", ascending=False).head(5)
            best_hits.to_csv(
                f'{out_path}/{out_suf}/{l}_susie_best_hits.tsv', sep="\t", index=False)
            

if __name__ == '__main__':
    run()