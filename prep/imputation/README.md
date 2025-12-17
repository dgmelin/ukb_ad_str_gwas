# Imputation
This directory contains the necessary steps to run the imputation of STRs based on SNP data from the UK Biobank and the SNP-STR haplotype panel from here
https://github.com/gymrek-lab/EnsembleTR.  
Please note that this is not intended to be run as-is. Instead, it should work as a reference to outline the necessary steps to perform the imputation rather than providing the fully functional code. 

We also need the following tools installed:
- [plink2](https://www.cog-genomics.org/plink/2.0/)
- [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard)
- [bcftools](http://www.htslib.org/)
- [vcftools](https://vcftools.github.io/index.html)
- [Java](https://www.java.com/en/download/) 
- [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
- [Conform-gt](https://faculty.washington.edu/browning/conform-gt/conform-gt.html)



## Steps (brief)
- Lift SNP data to hg38
- conform ids from SNP data and reference panel
- impute using beagle
- run post-imputation steps to clean data

## Steps (more detailled)
### Preparation
We start with BULK data plink files on nexus (Field ID 22418) and perform initial QC steps to increase reliability of imputation. 
That is, we filter for genotyping efficiency and minor allele frequency: 
```bash
prep=data/imputation/prep
plink2 -bfile /path/to/bulk/data/chrAB --geno 0.02 --maf 0.01 --make-pgen $prep/chrAB_prefiltered
```
Next, in order to use the data for imputation, plink files need to be converted into vcf format. 
Here, we want to make sure that the reference-allele is appearing first in the vcf output using the `ref-first` flag and `--fa` reference genome:
```bash
gen_ref=data/genome_references/human_g1k_v37.fasta
plink2 --pfile $prep/chrAB_prefiltered --export vcf-4.2 bgz ref-first --fa ${gen_ref} --out $prep/chrAB_converted 
```
Then, we lift the vcf file from hg37 to hg38 using LiftoverVcf from gatk. 
In our case, we needed to decompress our vcf data and at a `chr` prefix to all variants. 
```bash
#references for liftover
chain_file=data/genome_references/hg19ToHg38.over.chain.gz
ref_file=data/genome_references/hg38.fa

#unzip and add prefix
gzip -cd $prep/chrAB_converted | sed 's/^[^#]/chr&/g' > $tmpfile

output_prefix=$prep/chrAB_lifted
gatk LiftoverVcf -I $tmpfile \
    -O $output_prefix.vcf \
    --REJECT $output_prefix.rejected.vcf \
    -CHAIN $chain_file \
    -R $ref_file \
    --DISABLE_SORT true \
    --RECOVER_SWAPPED_REF_ALT true &> $output_prefix.log
```

### Conforming IDs and Imputation
Next, we need to conform the SNP ids between our SNP data and the reference panel. 
This is the last step before we can run the actual imputation. 
Both steps are run subsequently. 
```bash
conform_jar=data/imputation/conform-gt.24May16.cec.jar
beagle_jar=data/imputation/beagle.27Jan18.7e1.jar

#This is a placeholder
chr=AB
output_prefix=$prep/chr${chr}_conformed  

#This is the path to the SNP-STR-panel for the respective chromosome
ref_panel=data/imputation/panel/chr${chr}_SNP_TRs.vcf.gz

input_file=$prep/chr${chr}_lifted.vcf

echo "started chr%s conforming:\t" $chr
java -Xmx60g -jar $conform_jar \
    gt=${input_file} \
    ref=$ref_panel \
    chrom=chr$chr \
    match=POS \
    out=${output_prefix}_conformed &> ${output_prefix}_conformed.log
printf "fin\n"

output_prefix=data/imputation/imputed/chr${chr}
printf "started chr%s imputation:\t" $chr
java -Xmx60g -jar $beagle_jar \
        gt=${output_prefix}_conformed.vcf.gz \
        ref=$ref_panel \
        out=${output_prefix}_imputed &> ${output_prefix}_imputed.log
printf "fin\n"
```

### Post-Imputation
#### Initial filtering and QC
After imputation, we perform some initial filtering steps. First, we remove all SNPs that have been introduced through the imputation using vcftools.
We also add a new info field to the vcf header that describes the 'END' position of the STRs. 
This line says `##INFO=<ID=END,Number=1,Type=Integer,Description="ending position">` and is contained in the separate file `headinfo.txt`. 
```bash 
wrk_dir=data/imputation/data
add_info=$wrk_dir/headinfo

tmpdir=$wrk_dir/tmpdata
mkdir -p $tmpdir
tmp=$tmpdir/chr${chr}

out=data/imputation/imputed/chr${chr}

# remove SNP data
/usr/bin/time vcftools --gzvcf $f --keep-only-indels --recode --recode-INFO-all --stdout | bgzip -c > $tmp.STRs.vcf.gz
tabix $tmp.STRs.vcf.gz

head_info=$wrk_dir/headinfo.txt
 # reheader
cat <(bcftools view -h $tmp.STRs.vcf.gz | head -n 5) <(cat $head_info) <(bcftools view -h $tmp.STRs.vcf.gz | tail -n +6) > $tmpdir/header
bcftools reheader -h $tmpdir/header -o $tmp.rehead.vcf.gz $tmp.STRs.vcf.gz
```

Next, we perform multiallelic splitting and filter for STRs that have an imputation info score of at least 0.7. 
```bash
# split multiallelic
bcftools norm -m-any $tmp.rehead.vcf.gz -Oz -o $tmp.split.vcf.gz; 

# filter
echo "fin split, filter"
bcftools filter -e 'INFO/DR2<0.7' -Oz -o $tmp.split.filter.vcf.gz $tmp.split.vcf.gz 
```
Finally, we give the STRs more descriptive IDs that contain chr:pos:ref:alt information and index the final vcf file. 
```bash
# give more descriptive ids
cat <(bcftools view -h ${tmp}.split.filter.vcf.gz) \
    <(bcftools view -H ${tmp}.split.filter.vcf.gz | \
    awk 'BEGIN{OFS="\t";FS="\t"}{$3=$1":"$2":"$4":"$5; print}') | \
bgzip -c > ${out}.prepared.vcf.gz;

echo "fin rename, do index and fin"
bcftools index ${out}.prepared.vcf.gz

# clean up tmpdir
rm -r $tmpdir
```

Lastly, we convert the vcf into plink format to be able to run the association studies:
```bash
plink2 --vcf ${out}.prepared.vcf.gz --make-pgen --out ${out}_prepared
```

#### Note: 
To run this on the large scale for our UKB cohort, we neded to split the data/samples into smaller chunks for imputation.
This is not included here, but can be done by including a subsetting step as the first step before conforming/imputation. 
We also need a file that contains a list of sample IDs to include in the subset, i.e `sample_ids_042.txt`:
```bash
# subset samples
sample_file=$prep/sample_ids_042.txt
bcftools view -S $sample_file -Oz -o $tmp.subset.vcf.gz $input_file
```
After imputation, the smaller chunks can be merged back together using bcftools merge:
```bash
bcftools merge $out*.prepared.vcf.gz -Oz -o ${out}_mergedSTRs.vcf.gz
```
