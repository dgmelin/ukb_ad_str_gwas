library(MatrixEQTL)
## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

args <- commandArgs(trailingOnly=TRUE)

n_args <- 7

if (length(args) < n_args) {
  stop(paste0("Error: Too few arguments supplied (", length(args), ").\n",
              "Expected ", n_args, " arguments:\n",
              "1. SNP_file_name\n2. snps_location_file_name\n3. expression_file_name\n",
              "4. gene_location_file_name\n5. covariates_file_name\n",
              "6. output_file_name_cis\n7. output_file_name_tra"))
} else if (length(args) > n_args) {
  stop(paste0("Error: Too many arguments supplied (", length(args), ").\n",
              "Expected exactly ", n_args, " arguments."))
}

# Genotype file name
SNP_file_name            <- args[1] # Path to genotype file
snps_location_file_name  <- args[2] # Path to file with variant locations

# Gene expression file name
expression_file_name     <- args[3] # Path to gene expression file
gene_location_file_name  <- args[4] # Path to file with gene locations
# Covariates file name
# Set to character() for no covariates
covariates_file_name     <- args[5] # Path to covariates file

# Output file name
output_file_name_cis     <- args[6]
output_file_name_tra     <- args[7]

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1; # save all cis output
pvOutputThreshold_tra = 5e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6; # 1Mb

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = FALSE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = TRUE);