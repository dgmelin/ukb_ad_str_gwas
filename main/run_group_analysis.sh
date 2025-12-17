#!/bin/bash
# real example call from starting dir: 
   #/data/liga_starad/ukb
# bash run_var_qc.sh -o experiments/arraySnps/all -p imputation/data/initial_snps/merged -i sample_infos/all_samples/ -m 16000 -t 2

helpFunction()
{
   echo ""
   echo "Usage: $0 -p plinkPrefix -i info_path -o outPath"
   echo -e "\t-p prefix for initial plink pfile inputs, required"
   echo -e "\t-i path to info/phenotype file , required"
   echo -e "\t-o output prefix, required"
   echo -e "\t-m plink memory usage, in MB (see plink doc). Optional, default 4000"
   echo -e "\t-t plink threads, optional. default 1"
   exit 1 # Exit script after printing help
}

while getopts "p:i:o:h:m:t" opt
do
   case "$opt" in
      p ) pfiles="$OPTARG" ;;
      i ) infos="$OPTARG" ;;
      o ) out="$OPTARG" ;;
      m ) pmemory="$OPTARG" ;;
      t ) pthreads="$OPTARG" ;;
      h ) helpFunction ;;
      ? ) helpFunction ;; # also Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$pfiles" ] || [ -z "$infos" ] || [ -z "$out" ]
then
   echo "required input params missing.";
   helpFunction
fi

if [ -z "$pthreads" ] 
then
    pthreads=1
fi

if [ -z "$pmemory" ] 
then
    pmemory=4000
fi

# Begin script in case all parameters are correct
printf "running on \nplinkdir %s,\n infos %s \noutput %s\n" $pfiles $infos $out

# prepare cohort_file
cut -f -2 $infos > $out\_sample_infos
cut -f 1,3 $infos > $out.phenos
grep "control" $infos | cut -f 1 > $out\_ctrl.ids
#cp $infos $out/

# (re-) do variant qc
plink2 --pfile $pfiles --make-pgen --out $out\_update_1 \
   --update-sex $out\_sample_infos --pheno $out.phenos \
   --maf 0.01 --geno 0.02 --sort-vars\
   --keep $out\_sample_infos \
   --memory $pmemory --threads $pthreads; \

plink2 --pfile $out\_update_1 --keep $out\_ctrl.ids --hwe 5e-06 --write-snplist 'allow-dups' -out $out\_hwe_ctrls --memory $pmemory --threads $pthreads; \
plink2 --pfile $out\_update_1 --extract $out\_hwe_ctrls.snplist --make-pgen --out $out\_varqced_2 --memory $pmemory --threads $pthreads; \

# combine with covars? 
# -> combine before and place them at out directory?
 

# run association analysis
plink2 --pfile $out\_varqced_2 --covar $out.covars --glm --out $out\_assoc --memory $pmemory --threads $pthreads