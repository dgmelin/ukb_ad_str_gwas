#!/bin/bash

#read input file from command line
# Format should be iSTR-ID gSTR-ID 
input=$1

istr_base=data/iSTR/white_british/all/
gstr_base=data/gSTR/white_british/white_british_varqced_2

base=data/imputation_quality
out_base=${base}/out
#Make this point to the correlate_genotypes.py script
py_base=${base}/correlate_genotypes.py

tail -n +2 "$input" | while read -r col1 col2 rest; do
    geno_var="$col1"
    impu_var="$col2"

    idi=$(echo "$impu_var" | cut -d: -f-2 | sed 's/:/_/')
    idg=$(echo "$geno_var" | sed 's/:/_/')

    out_dir="$out_base/${idi}_vs_${idg}"
    mkdir -p $out_dir

    #find relevant files
    chr=$(echo "$impu_var" | cut -d: -f1 | sed 's/chr//')
    istr=$istr_base/$chr/white_british_all_varqced_2
    gstr=$gstr_base/$chr/white_british_varqced_

    #extract genotypes
    plink2 --snp "$geno_var" --pfile $gstr --recode A --out "$out_dir/gstr"
    plink2 --snp "$impu_var" --pfile $istr --recode A --out "$out_dir/istr"

    # strip additional infos and only extract variant line
    grep -F "$geno_var" $gstr.pvar > $out_dir/gstr.var
    grep -F "$impu_var" $istr.pvar > $out_dir/istr.var

    #run correlation
    python $py_base -i $out_dir/istr -g $out_dir/gstr -o $out_dir/${idi}_vs_${idg}
done
