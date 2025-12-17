# Run Meta
## Prep
First, we want to prepare the outputs of the association studies accordingly. 
The expected columns differ for plink2s association outputs and the columns expected by plinks --meta:
```
Source is 
#CHROM(1) POS(2) ID(3) REF(4) ALT(5) A1(6) TEST(7) OBS_CT(8) BETA(9) SE(10) T_STAT(11) P(12) ERRCODE(13)
Target is
CHR BP SNP A2 A1 OBS_CT BETA SE (Z_STAT) P -> (or similar) 
```
```bash
awk -f format_to_meta.awk /path/to/white_british/association/results > white_british.meta
awk -f format_to_meta.awk /path/to/other_white/association/results > other_white.meta
``` 
## Run
Next, we run the meta study using plink1

```bash
plink --meta-analysis white_british.meta other_white.meta + logscale
```
Since the output contains many spaces insead of tabs, we fix it to make it readable. 
```bash
cp plink.meta all_white.meta
for i in {1..10}; 
    do echo $i; 
    sed -i 's/  / /g' all_white.meta; 
done; 
sed -i 's/^ //g' all_white.meta; 
sed -i 's/ /\t/g' all_white.meta
```