from matplotlib import pyplot as plt
import pandas as pd
import argparse
import seaborn as sns

var_header = [
    "CHROM", "POS", "ID", "REF", "ALT", "FILTER", "INFO"
]

def get_options():    
    parser = argparse.ArgumentParser(description="Join tables script")
    parser.add_argument("-i", "--imputed", type=str, required=True, help="Inputed data prefix (expects prefix.raw and prefix.var)")
    parser.add_argument("-g", "--genotyped", type=str, required=True, help="Inputed data prefix (expects prefix.raw and prefix.var)")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output prefix")
    parser.add_argument("-s", "--silence", action='store_true', help="Suppress output messages")
    parser.add_argument("-p", "--no_plot", action='store_true', help="Suppress plot generation")
    return parser.parse_args()

def log(message, verbose=True):
    if verbose:
        print(message)

def run(args):
    v = not args.silence

    #load data
    imputed_df, imputed_var_id = load_data(f'{args.imputed}.raw', sep="\t")
    imputed_df['IID'] = [x.split('_')[0] for x in imputed_df['IID'].astype(str)]
    log(f"loaded imputed data with shape {imputed_df.shape}", v)
    genotyped_df, genotyped_var_id = load_data(f'{args.genotyped}.raw', sep="\t")
    log(f"loaded genotyped data with shape {genotyped_df.shape}", v)

    #process genotypes
    geno_specs = pd.read_csv(f'{args.genotyped}.var', sep="\t", names=var_header)[["ID", "REF", "ALT", "INFO"]]
    impu_specs = pd.read_csv(f'{args.imputed}.var', sep="\t", names=var_header)[["ID", "REF", "ALT"]]

    geno_lens = get_len_dict(geno_specs, source='genotyped')
    impu_lens = get_len_dict(impu_specs, source='imputed')

    imputed_df['sum_len'] = imputed_df[imputed_var_id].map(impu_lens)
    genotyped_df['sum_len'] = genotyped_df[genotyped_var_id].map(geno_lens)

    merged_df = pd.merge(
        imputed_df.astype({'IID': str}), genotyped_df.astype({'IID': str}), 
        on="IID", how='right', suffixes=('_imputed', '_genotyped')
    )

    log(f"merged data to shape {merged_df.shape}", v)

    merged_df.to_csv(f'{args.output}.tsv', index=False, sep="\t")

    log(f"Calculate pearson correlation", v)
    corr = merged_df[['sum_len_imputed', 'sum_len_genotyped']].corr(method='pearson').iloc[0,1]
    n=merged_df[['sum_len_imputed', 'sum_len_genotyped']].dropna().shape[0]
    with open(f'{args.output}_correlation.txt', 'w') as f:
        f.write('imputed_variant_id\tgenotyped_variant_id\tcorr\tn\n')
        f.write(f'{imputed_var_id}\t{genotyped_var_id}\t{corr}\t{n}\n')
        f.close()

    log(f'{imputed_var_id} vs {genotyped_var_id} pearson correlation:\n{corr}\n')
    log(f"Correlation result saved to {args.output}_correlation.txt", v)

    if args.no_plot:
        return
    
    # Plot
    counts = merged_df.groupby(["sum_len_imputed","sum_len_genotyped"]).size().reset_index(name="count")
    sns.scatterplot(
        data=counts, x="sum_len_imputed", y="sum_len_genotyped", size="count", sizes=(10, 200), alpha=0.5
    ).set_title(f'Correlation: {corr:.4f}')

    # Save plot
    plt.savefig(f'{args.output}_correlation.png')
    plt.close()


def load_data(file_path, sep="\t"):
    df = pd.read_csv(file_path, sep=sep)
    var_id = df.keys()[-1]
    df = df[["IID", var_id]] 
    return df, var_id

def get_geno_lens(df):
    ref_len = len(df.iloc[0]['REF'])
    info = df.iloc[0]['INFO'] 
    info_dict = {k: v for k, v in (x.split('=') for x in info.split(';') if '=' in x)}
    if 'Motif' in info_dict:
        mot_len = len(info_dict['Motif'])
    else:
        mot_len = ref_len
    alt = df.iloc[0]['ALT']
    af = float(alt.replace('<', '>').replace('>', ''))
    al = af * mot_len

    return ref_len, al

def get_len_dict(df, source):
    if source == 'genotyped':
        rl, al = get_geno_lens(df)
    elif source == 'imputed':
        rl = len(df.iloc[0]['REF'])
        al = len(df.iloc[0]['ALT'])
    else:
        raise ValueError("source must be 'genotyped' or 'imputed'")
    
    return {
        0 : 2 * rl,
        1 : rl + al,
            2 : 2 * al
        }

if __name__ == "__main__":
    run(get_options())

