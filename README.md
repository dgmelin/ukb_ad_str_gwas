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

## Disclaimer

The code is provided as is, without any guarantee of completeness or functionality outside the context of the original study.
Researchers are welcome to reuse and modify it for their own work.
We kindly ask that you cite the accompanying manuscript if you make use of this code in your own work.