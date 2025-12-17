# Multiallelic Approach
We added a multiallelic analysis branch to ensure that splitting of multiallelic variants does not introduce false positive findings. 
Here, we took variants +- 500kb around leading STR variants at at least suggestively significant loci. 

**Preparation**<br>
To do so, we prepared a bed-file that contain the regions of interest (see `extract_regions.py`)  
Next,  absed on these regions the script `run_prepare_extract.sh` is preparing the data coming from imputation such that it can be analyzed using the associaTR tool.
It perfomrs the steps mentioned in the methods section and can be run via `bash run_prepare_extract.sh $chr` for the respective chromosome.
First, we extract the regions of interest from the initial files. Next, biallelic records and multiallelic records are extracted. For positions where a multiallelic variants exists, the multiallelic one is taken. Otherwise, biallelic variants at the same position are rejoined and finally, remaining variants are merged together.  

Lastly, we use the `annotate_imputed.py` script to format the entries in the way that the files are accepted by associaTR. 

**Analysis**
To run the actual association analysis, run the `run_associaTR.sh` script via `bash run_associaTR.sh -c $chr -w /path/to/working/directory`. This working directory should contain the files prepared in the preparation step, as well as sample info files. 



