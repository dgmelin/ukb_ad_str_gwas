# STR GWAS on AD in UKB
This repository contains the code accompanying the paper '**GWAS on short tandem repeats identifies novel genetic mechanisms in Alzheimer’s disease**' ([doi: https://doi.org/10.1101/2025.03.13.25323833](https://doi.org/10.1101/2025.03.13.25323833)). 
This code is intended to illustrate the main steps and methods used to generate the results presented in the paper, reflecting the structure and logic of the analyses.
However, please note that it is provided as a reference implementation and is not ready for direct use "out of the box". 
Before applying the code to your own data, you will need to adapt paths and parameters to fit your specific environment.  

This repo is structured into three main sections: 
- [prep](./prep/): contains the steps we performed prior to the main analyses, 
- [main](./main/): contains the steps we performed to run the actual GWAS, 
- [eval](./eval/): contains the evaluation step we we performed post GWAS, such as evaluation of main results, finemapping.

If you encounter any issues or if any questions are raised, I'm happy to help. 

### Associated Data
* **Imputation Panel:** Available on Zenodo [10.5281/zenodo.8365671](https://doi.org/10.5281/zenodo.8365671)
* **GWAS Summary Statistics:** Available on Zenodo [10.5281/zenodo.17908177](https://doi.org/10.5281/zenodo.17908177)

## License and Data Access

The code is provided as is, without any guarantee of completeness or functionality outside the context of the original study. 
The code is licensed under the GNU General Public License v3.0. 
Researchers are welcome to reuse and modify it for their own work.
We kindly ask that you cite the accompanying manuscript if you make use of this code in your own work.

**Strict Data Separation**  
Please note that the open-source license applies exclusively to the code in this repository. 
It does not apply to the underlying biological or clinical data.
Due to data protection regulations, the raw and individual-level genetic data from the UK Biobank and the Oxford Project to Investigate Memory and Ageing (OPTIMA) cannot be shared publicly.

- No individual-level participant data, restricted phenotypes, or raw genotypes are hosted in this repository.
- Researchers wishing to run these scripts on the original study cohorts must independently apply for data access through the UK Biobank Access Management System and the respective OPTIMA data management committees.
- Any use of the original data is strictly governed by the legal terms of the respective Material Transfer Agreements (MTAs) and Data Use Agreements (DUAs), which explicitly prohibit unauthorized distribution or attempts to re-identify participants.

For access to the publicly available GWAS summary statistics and the SNP-STR imputation panel generated in this study, please see the Associated Data section above for the respective Zenodo repository links and the original publication. 
