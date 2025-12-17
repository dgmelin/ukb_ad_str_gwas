#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -p plinkPrefix -i info_path -o outPath"
   echo -e "\t-i variant ID, mandatory"
   echo -e "\t-c variant chr, mandatory"
   echo -e "\t-s source, to extrat the variant from, mandatory"
   echo -e "\t-a association analysis, either SNP or STR, default STR"
   echo -e "\t-v covariates file, defaults to data/white_british_all.covars"
   echo -e "\t-o output directory, default pwd"
   echo -e "\t-m plink memory usage, in MB (see plink doc). Optional, default 4000"
   echo -e "\t-t plink threads, optional. default 1"
   exit 1 # Exit script after printing help
}

while getopts "i:o:a:c:s:h:m:v:t" opt
do
   case "$opt" in
      i ) id="$OPTARG" ;;
      c ) chr="$OPTARG" ;;
      s ) source="$OPTARG" ;;
      a ) assoc="$OPTARG" ;;
      o ) out="$OPTARG" ;;
      m ) pmemory="$OPTARG" ;;
      t ) pthreads="$OPTARG" ;;
      v ) covars_file="$OPTARG" ;;
      h ) helpFunction ;;
      ? ) helpFunction ;; # also Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$id" ] || [ -z "$chr" ] || [ -z $source ]
then
   echo "required input params missing.";
   helpFunction
fi

if [ -z "$assoc" ] 
then
    assoc='STR'
fi

if [ -z "$out" ] 
then
    out=$(pwd)
fi

if [ -z "$pthreads" ] 
then
    pthreads=1
fi

if [ -z "$pmemory" ] 
then
    pmemory=4000
fi

if [ -z "$covars_file" ] 
then
    covars_file=data/white_british_all.covars
fi

p_assoc="dummy_to_replace"
if [ $assoc == "STR" ]
then 
    # Enter own path to STR association data here
    p_assoc="data/assoc/str_data/white_british/all/${chr}/white_british_all_varqced_2"
else
    p_assoc="data/assoc/snp_data/white_british/all/${chr}/white_british_all_varqced_2"
fi

py_file=scripts/combine_covars.py

out_prefix=${out}/$id

printf "extracting %s from %s to %s\n" $id $source $out_prefix

# Export relevant variant to condition for
plink2 --pfile $source  --snp $id --export A --out $out_prefix  --allow-extra-chr --memory $pmemory --threads $pthreads

#Format it to integrate into covariates
cut -f 2,7 $out_prefix.raw > $out_prefix.edit
python $py_file -v -s $out_prefix.edit -i $covars_file -o $out_prefix.covars

# Run assoc with additional covariate
plink2 --pfile $p_assoc \
    --glm --covar $out_prefix.covars --out $out_prefix\_assoc --allow-extra-chr --memory $pmemory --threads $pthreads 


