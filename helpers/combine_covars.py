import pandas as pd
import argparse

#defaults 
str_keys=["#IID"]
snp_keys=["#FID", "IID"]

def get_options():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--initialCovars', type=str, default=None)
    parser.add_argument('-s','--snpCovars', type=str, default=None)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--create_dummies', action='store_true')
    parser.add_argument('-o', '--outfile', type=str, default='covars')

    return parser.parse_args()

def main():
    options=get_options()
    initial = pd.read_csv(options.initialCovars, sep='\t')
    snp_covars = pd.read_csv(options.snpCovars, sep='\t')

    v = options.verbose
    if v:
        msg = f'loaded two files with dimensions:\ninitial:\t{initial.shape}\nsnpCovars:\t{snp_covars.shape}'
        print(msg)

    if options.create_dummies:
        dummy_cols = snp_covars.keys()[1:]
        dummy_df = pd.get_dummies(snp_covars, columns=dummy_cols)

        if v:
            print(f'created dummy variables, snpCovars have size {dummy_df.shape} now')
        
        snp_covars = dummy_df
        
    snp_covars = snp_covars.rename(columns={'IID':'#IID'})
    combined = initial.merge(snp_covars, on='#IID')
    if v:
        print(combined.head())

    combined.to_csv(options.outfile, sep='\t', index=False, na_rep='NA')

if __name__=='__main__':
    main()