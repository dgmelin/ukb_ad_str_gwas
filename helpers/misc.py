# NOTE: 
# This is a collection of functions that are used in multiple places throughout the repository
# Refer to their appropriate sections for examples and usage

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# %% Settings and constants
custom_palette = sns.color_palette()[:4] 
manhatti_palette = {}
for c in range(1, 23):
    manhatti_palette[c] = custom_palette[c % len(custom_palette)]

font_kwargs={
    'fontsize':18
}

suggestive_t=1e-5
gw_t = 0.05/334605 

scatter_palette={
    'SNP':'tab:blue', 
    'STR':'tab:orange', 
    'adjusted':'tomato', 
    'unadjusted':'teal'
}

marker_palette={
    'lead STR after conditioning\non Bellengeuz\' SNP' : "s",
    'lead STR after conditioning on lead SNP' : "s",
    'lead STR' : "D",
    'lead SNP after conditioning\non lead STR' : "s",
    'lead SNP in Bellenguez/Jansen' : "D",
    'lead SNP' : "D"
}

# %% Data processing functions
def get_associaTR_results(file_path, is_meta=False, cut=True, cut_threshold=10, p_name='p_ad_risk_by_proxy', rename_chr=True):
    df = pd.read_csv(file_path, sep="\t").rename(columns={
        p_name:'P'
    })
    initial = df.shape[0]
    df = df.dropna(axis=0)
    df = df[df['locus_filtered'] == 'False']
    after = df.shape[0]
    print(f'drop {initial-after} rows because of NA values')
    df = df.astype({'P':float})
    df.loc[df.P<sys.float_info.min, 'P']=sys.float_info.min
    min_ct = df[df.P == sys.float_info.min].shape[0]
    if min_ct > 0:
        print(f'set {min_ct} P values to maximum possible min of {sys.float_info.min}')
    df['-log10'] = -np.log10(df.P)

    if cut:
        df=df[df['-log10'] < cut_threshold]
    if rename_chr:
        df['chrom'] = df['chrom'].str.replace('chr', '').astype(int)
    return df.sort_values(by=['chrom', 'pos']).reset_index(drop=True)
    
def get_results(file_path, is_meta=False, cut=True, cut_threshold=10):
    df = pd.read_csv(file_path, sep="\t")
    initial = df.shape[0]
    df = df.dropna(axis=0)
    after = df.shape[0]
    print(f'drop {initial-after} rows because of NA values')
    df = df.astype({'P':float})
    df.loc[df.P<sys.float_info.min, 'P']=sys.float_info.min

    min_ct = df[df.P == sys.float_info.min].shape[0]
    if min_ct > 0:
        print(f'set {min_ct} P values to maximum possible min of {sys.float_info.min}')
    df['-log10'] = -np.log10(df.P)

    if not is_meta:
        df = df[df['TEST']=='ADD']

    if cut:
        df=df[df['-log10'] < cut_threshold]
    
    return df

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

# %% Plotting
def manhatti(df, title, sort_keys=['#CHROM', "POS"], hue_key='#CHROM', style_key=None, y_steps=10, legend=None, 
             sort_key='-log10', sort_asc=True, plot_y_thresh=True, do_full=True, suggest_t=-np.log10(suggestive_t), p_sig_t=-np.log10(gw_t), palette=manhatti_palette):
    plot_df = df.copy(deep=True)
    plt.figure(figsize=(30,10))
    if do_full:
        chr_maxxs = plot_df.groupby('#CHROM')['POS'].max()
        chroms_sorted = sorted(chr_maxxs.keys())
        chr_starts = {chroms_sorted[0]: 0}
        for idx, c in enumerate(chroms_sorted[1:], start=1):
            prev_c = chroms_sorted[idx - 1]
            chr_starts[c] = chr_starts[prev_c] + chr_maxxs[prev_c] + 1
        plot_df['i'] = [p + chr_starts[c] for p,c in zip(plot_df['POS'], plot_df['#CHROM'])]
    else:
        plot_df = df.sort_values(['numeric_chr', sort_keys[1]])
        plot_df = plot_df.reset_index(drop=True); 
        plot_df['i'] = plot_df.index

    plot_df = plot_df.sort_values(by=sort_key, ascending=sort_asc)
    plot = sns.scatterplot(plot_df, x='i', y='-log10', hue=hue_key, legend=legend, linewidth=0.3, style=style_key, palette=palette)
    labels_df=plot_df.groupby(sort_keys[0])['i'].median()
    plot.set_xlabel('chr')
    plot.set_xticks(labels_df)
    plot.set_xticklabels(labels_df.index)
    max_x = plot_df.i.max()
    max_x = max_x + 0.01*max_x
    min_x = plot_df.i.min() - 0.01*max_x
    plt.xlim([min_x, max_x])

    max_y=int(max(plot_df['-log10'].max(), y_steps))+1
    step = int(max_y/y_steps)
    
    plot.set_yticks(range(1, max_y, step))
    if plot_y_thresh:
        plt.hlines(y=[p_sig_t, suggest_t], xmin=min_x, xmax=max_x, colors=['tab:red', 'k'], linestyles='dashed')

    plot.grid(axis='y')
    plot.set_title(title)
    plt.grid(axis='y')
    plt.title(title)
    plt.tight_layout()

    return plot.figure

def broken_manhatti(df, title, sort_keys=['#CHROM', "POS"], hue_key='#CHROM', style_key=None, legend=None, 
             sort_key='-log10', sort_asc=True, plot_y_thresh=True, do_full=True, suggest_t=5, p_sig_t=-np.log10(gw_t), plot_break=(30,100), kwargs={}):
    plot_df = df.copy(deep=True)
    f, (outlier_ax, plot_ax) = plt.subplots(2, 1, sharex=True, figsize=(30,10), height_ratios=[1,8])
    if do_full:
        chr_maxxs = plot_df.groupby('#CHROM')['POS'].max()
        chroms_sorted = sorted(chr_maxxs.keys())
        chr_starts = {chroms_sorted[0]: 0}
        for idx, c in enumerate(chroms_sorted[1:], start=1):
            prev_c = chroms_sorted[idx - 1]
            chr_starts[c] = chr_starts[prev_c] + chr_maxxs[prev_c] + 1
        plot_df['i'] = [p + chr_starts[c] for p,c in zip(plot_df['POS'], plot_df['#CHROM'])]
    else:
        plot_df = df.sort_values(sort_keys)
        plot_df = plot_df.reset_index(drop=True); 
        plot_df['i'] = plot_df.index

    plot_df = plot_df.sort_values(by=sort_key, ascending=sort_asc)
    plot = sns.scatterplot(plot_df, x='i', y='-log10', hue=hue_key, palette=manhatti_palette, legend=None, linewidth=0.3, style=style_key, ax=plot_ax, **kwargs)
    
    labels_df=plot_df.groupby(sort_keys[0])['i'].agg(['median', 'min', 'max'])
    # labels_df['mid'] = (labels_df['min'] + labels_df['max']) / 2 
    # labels_df['mid'] = labels_df.apply(lambda row: row['mid']/2 if row['mid'] == row['max'] else row['mid'], axis=1)
    
    labels_df['mid'] = (labels_df['max'].shift() + labels_df['max']) / 2
    labels_df.loc[labels_df['mid'].isna(), 'mid'] = labels_df['max'].iloc[0] / 2

    plot_ax.set_xticks(labels_df['max'])
    plot_ax.set_xticklabels([], minor=False)

    plot_ax.set_xticks(labels_df['mid'], minor=True)
    plot_ax.set_xticklabels(labels_df.index, minor=True)
    plot_ax.xaxis.grid(True)
    
    max_x = plot_df.i.max()
    max_x = max_x + 0.01*max_x
    min_x = plot_df.i.min() - 0.01*max_x
    plt.xlim([min_x, max_x])
    plot_ax.set_yticks(range(0, 30, 5))

    plot.set_xlabel('Chromosome', fontdict=font_kwargs)
    plot.set_ylabel('-log10(p)', fontdict=font_kwargs)

    if plot_y_thresh:
        plt.hlines(y=[p_sig_t, suggest_t], xmin=min_x, xmax=max_x, colors=['tab:red', 'k'], linestyles='dashed')
    plt.grid(axis='y')
        # plt.tight_layout()

    sns.scatterplot(plot_df, x='i', y='-log10', hue=hue_key, palette=manhatti_palette, legend=legend, linewidth=0.3, style=style_key, ax=outlier_ax)
    y_max_outliers = outlier_ax.get_ylim()[1]
    y_shift = y_max_outliers + y_max_outliers*0.1
    outlier_ax.set_ylim(bottom=plot_break[1], top=y_shift)  # outliers only
    plot_ax.set_ylim(0, plot_break[0])  # most of the data

    outlier_ax.spines["bottom"].set_visible(False)
    plot_ax.spines["top"].set_visible(False)

    # outlier_ax.xaxis.tick_top()
    outlier_ax.tick_params(top=False, bottom=False)

    d = .25  # proportion of vertical to horizontal extent of the slanted line
    ax_kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    outlier_ax.plot([0, 1], [0, 0], transform=outlier_ax.transAxes, **ax_kwargs)
    outlier_ax.set_ylabel('')
    plot_ax.plot([0, 1], [1, 1], transform=plot_ax.transAxes, **ax_kwargs)



    plt.subplots_adjust(hspace=0.05)
    outlier_ax.set_title(title, **font_kwargs)
    plt.tight_layout()

    return plot.figure


def plt_qq(df, title, do_save=False, save_path='data/qq_filtered_apoe/qq.png'):
    plt.figure(figsize=(10,10))
    obs=df['-log10'].sort_values()
    expect = np.sort(-np.log10(np.random.uniform(size=obs.shape)))
    qq_df = pd.DataFrame({'expected':expect, 'observed':obs})
    sns.regplot(qq_df, x='expected', y='observed', ci=None, line_kws=dict(color="C1", linewidth=0.99, linestyle='dashed'))
    z = stats.norm.ppf(df['P']/2)
    infl_coeff = np.median(z**2)/stats.chi2.ppf(0.5, 1)
    figtext = r'$\lambda$' + '={:.4f}'.format(infl_coeff)
    plt.text(0.2, 0.8, figtext, transform=plt.gca().transAxes)
    plt.title(title)
    plt.tight_layout()
    plt.xlabel('Expected -log10(p)')
    plt.ylabel('Observed -log10(p)')
    if do_save:
        plt.savefig(save_path)


def locus_zoom_review(df, chr, pos, plt_range=5e5, 
               title='', plot_y_thresh=True, highligh_lead=True, suggest_t=5, p_sig_t=-np.log10(1.6*10**(-7)),
               palette=None, style_palette=marker_palette, hue_key='type', kwargs={}, style_key=None, highlight_x=None
    ): 

    if plt_range == 0:
        plt_df = df[(df['#CHROM']==chr)]
    else:
        plt_df = df[(abs(df.POS-pos)<plt_range) & (df['#CHROM']==chr)]
    plt_df = plt_df.sort_values(by='-log10' )
    ax = sns.scatterplot(plt_df, x='POS', y='-log10', hue=hue_key, palette=palette, style=style_key, markers=style_palette, **kwargs)
    
    if plot_y_thresh:
        max_x = plt_df.POS.max()
        min_x = plt_df.POS.min()
        plt.hlines(y=[p_sig_t, suggest_t], xmin=min_x, xmax=max_x, colors=['tab:red', 'k'], linestyles='dashed')
    
    if highligh_lead:
        max_y = plt_df['-log10'].max()
        max_y = ax.get_ylim()[1]*0.97
        if highlight_x is None:
            if plt_df is None or plt_df.empty:
                highlight_x = max_x - plt_range/2
            else: 
                highlight_x = plt_df.iloc[0].POS
                if highlight_x > max_x or highlight_x < min_x:
                    highlight_x = max_x - plt_range/2

        plt.vlines(x=highlight_x, ymin=0, ymax=max_y, alpha=0.1, color='orange', linewidth=30)
    ax.set_ylabel('-log10(p)', **font_kwargs)
    ax.set_xlabel('Genomic coordinate (hg38)', **font_kwargs)

    xticks_df = pd.DataFrame({'POS':[pos - plt_range, pos, pos+plt_range]}).astype('int64')
    ax.set_xticks(xticks_df.POS.to_list())
    ax.set_xticklabels(xticks_df.POS)
    
    plt.title(title)
    plt.grid(axis='y')
    return ax
    
def strip_compare(df, palette=scatter_palette, show_t=True, show_legend=False, title="", plot_break=(20,300), axis_steps = 4):
    f, (plot_ax, outlier_ax) = plt.subplots(1, 2, sharey=True, width_ratios=[8,1])
    g = sns.stripplot(df, y='ID', x='-log10', hue='type', legend=show_legend, palette=palette, ax=plot_ax)
    sns.stripplot(df, y='ID', x='-log10', hue='type', legend=False, palette=palette, ax=outlier_ax)
    if show_t:
        y_min, y_max = plt.ylim()
        if y_min < y_max:
            tmp = y_min
            y_min = y_max
            y_max = tmp
        plot_ax.vlines(x=[-np.log10(gw_t), 5], ymax=y_max, ymin=y_min, colors=['tab:red', 'k'], linestyles='dashed')

    outlier_ax.set_xlim(left=plot_break[1])  # outliers only
    plot_ax.set_xlim(0, plot_break[0])  # most of the data
    outlier_ax.spines["left"].set_visible(False)
    plot_ax.spines["right"].set_visible(False)

    # outlier_ax.xaxis.tick_top()
    outlier_ax.tick_params(left=False, right=False)


    step_size = int(plot_break[0] / axis_steps)
    plot_ax.set_xticks(range(0, plot_break[0] + 1, step_size))


    d = .5  # proportion of vertical to horizontal extent of the slanted line
    ax_kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    outlier_ax.plot([0, 0], [0, 1], transform=outlier_ax.transAxes, **ax_kwargs)
    outlier_ax.set_ylabel('')
    plot_ax.plot([1, 1], [1, 0], transform=plot_ax.transAxes, **ax_kwargs)

    g.set_xlabel('-log10(p)', **font_kwargs)
    outlier_ax.set_xlabel("")
    g.set_title(title)
    l = g.get_legend()
    
    sns.move_legend(
        g, 
        loc="lower left",
        bbox_to_anchor=(-.5, -.10),
        # ncol=2, 
        title=None, frameon=True,
    )
    return df[['#CHROM', 'type', 'POS', 'ID', 'P', 'BETA', '-log10']]