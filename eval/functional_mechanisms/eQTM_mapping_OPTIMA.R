## This script is used to perform an eQTM mapping with known AD-CpGs. First, 
## the CpGs and transcripts that have a selected distance to each other are 
## determined. These are then used to calculate linear regressions. Correlation 
## plots are also created.

#.libPaths(.libPaths()[grepl("rtrack_env", .libPaths())])


library(data.table)
library(ggplot2)
library(rlang)
library(performance)
library(rtracklayer)
library(dplyr)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)
cis <- args[1] ## Maximum distance between CpGs and transcripts
cis <- as.numeric(cis)
out_path <- args[2] ## TODO change path to save results

print(cis)


DNAm <-read.csv("path/to/M.csv",
	header = TRUE, row.names = 1, check.names = FALSE) ## TODO Change Path to DNAm data
cpgs <- read.csv("path/to/CpGs_gmelin.csv") ## TODO Change Path to CpGs to be tested


DNAm <- DNAm[is.element(rownames(DNAm), cpgs[,1]),]

RNA <- read.csv("path/to/G.csv",
	header = TRUE, row.names = 1, check.names = FALSE) ## TODO Change Path to mRNA data
pheno <- read.csv("path/to/C.csv", 
	header = TRUE, row.names = 1, check.names = FALSE) ## TODO Change Path to covariates. We used: age, sex, AD disease status, study, RIN, PMI, DNAm-PCs and mRNA-PCs



annotation_g <- fread("path/to/G.bed6")
annotation_m <- fread("path/to/M.bed6")


## prepare dataframe for pairs of CpGs and transcripts for the eQTM Mapping
eQTMs <- data.frame(mt_id = character(),
		gt_id = character())

## Select CpGs and transcripts with the selected distance to each other
for(i in 1:23){
  anno_g <- subset(annotation_g, chrom == i)
  anno_m <- subset(annotation_m, chrom == i)

  dnam <- anno_m[is.element(anno_m$name, row.names(DNAm)),]
  mrna <- anno_g[is.element(anno_g$name, row.names(RNA)),]


  if(nrow(dnam) != 0){
    for(j in 1:nrow(dnam)){
      for(n in 1:nrow(mrna)){
      
        ## CpG is located in the transcript
        if(mrna$chromStart[n] < dnam$chromStart[j] & dnam$chromStart[j] < mrna$chromEnd[n]){
          new_row <- data.frame(mt_id = dnam$name[j], gt_id = mrna$name[n])
          eQTMs <- rbind(eQTMs, new_row)
        }
      
        ## CpG is located in front of the transcript
        if(dnam$chromStart[j] < mrna$chromStart[n]){
          if(mrna$chromStart[n] - dnam$chromStart[j] < cis){
            new_row <- data.frame(mt_id = dnam$name[j], gt_id = mrna$name[n])
            eQTMs <- rbind(eQTMs, new_row)
          }
        }
      
        ## CpG is located behind the transcript
        if(dnam$chromStart[j] > mrna$chromEnd[n]){
          if(dnam$chromStart[j] - mrna$chromEnd[n] < cis){
            new_row <- data.frame(mt_id = dnam$name[j], gt_id = mrna$name[n])
            eQTMs <- rbind(eQTMs, new_row)
          }
        }
      }
    }
  }
}



print(c("Number of tests:", nrow(eQTMs)))


## add gene id and name
anno_g <- fread("path/to/GTF_GRCh38_all_252989.txt") ## TODO change path to file with gene ID and name
eQTMs <- merge(eQTMs, anno_g[,c("transcript_id", "gene_id", "gene_name")], by.x = "gt_id", by.y = "transcript_id") ## TODO change column names

## Prepare dataframe for results
results <- data.frame(dnam = eQTMs$mt_id,
		rna = eQTMs$gt_id,
		gene_id = eQTMs$gene_id,
		gene_name = eQTMs$gene_name,
		outlier = numeric(length = nrow(eQTMs)),
		est = numeric(length = nrow(eQTMs)),
		std_error = numeric(length = nrow(eQTMs)),
		p = numeric(length = nrow(eQTMs)),
		q = numeric(length = nrow(eQTMs)),
		r = numeric(length = nrow(eQTMs)),
		rse = numeric(length = nrow(eQTMs)),
		vif = numeric(length = nrow(eQTMs)))


## function to remove outliers iteratively 
remove_outlier <- function(df, column, threshold = 4) {
  
  repeat {
    mean_val <- mean(df[[column]], na.rm = TRUE)
    sd_val <- sd(df[[column]], na.rm = TRUE)
  
    new_df <- df[abs(df[[column]] - mean_val) <= threshold * sd_val, ]
    
    if (nrow(new_df) == nrow(df)) {
      break
    }
    
    df <- new_df
  }
  
  return(df)
}


## For loop to remove the outliers, calculate the linear regressions and create the correlation plots		 
print("Start linear regressions")
for (i in 1:nrow(results)) {
  
  ## extract cpg and transcript for regression
  dnam <- as.data.frame(t(DNAm[results$dnam[i],]))
	dnam$id <- row.names(dnam)
  mrna <- as.data.frame(t(RNA[results$rna[i],]))
  mrna$id <- row.names(dnam)
	pheno$id <- row.names(pheno)
  
	
	## merge DNAm, mRNA and covariate data data for regression and plotting
  plot_data <- merge(dnam, mrna, by= "id")
  plot_data <- merge(plot_data, pheno, by = "id")
  
	## remove outliers
	n_with_outlier <- nrow(plot_data)
  plot_data <- remove_outlier(plot_data, 2) ## remove outliers based on DNAm
  plot_data <- remove_outlier(plot_data, 3) ## remove outliers based on mRNA
	n_without_outlier <- nrow(plot_data)

	
	## z-standatization for better comparability of the different data sets 
	if(!all(plot_data[,3] == 0)){ ## Z-standadization is not possible if all values are 0
	
	  plot_data[,2] <- scale(plot_data[,2])
	  plot_data[,3] <- scale(plot_data[,3])

	  ## linear regression with covariates 
	  lm_model <- lm(as.formula(paste(names(plot_data)[3], "~", names(plot_data)[2],
	                                "+age+sex+status+pmi_h+rin+PC1.x+PC2+PC3+PC1.y")), ## TODO change covariates
	                 data = plot_data)
  
	
	  results$est[i] <- summary(lm_model)$coefficients[2,1]
    results$std_error[i] <- summary(lm_model)$coefficients[2,2]
    results$p[i] <- summary(lm_model)$coefficients[2,4]
	  results$rse[i] <- summary(lm_model)$sigma ## residual standard error
	  results$vif[i] <- check_collinearity(lm_model)$VIF[1] ## variance inflation factor
	  results$outlier[i] <- n_with_outlier - n_without_outlier ## number of outliers
	results$r[i] <- round(cor(plot_data[2], plot_data[3], method = "pearson"),4) 	
	
	  ## create correlation plots
	  if(!is.na(results$p[i])){
		  plot_data$status <- ifelse(plot_data$status == 1,"AD","HC")
		  R <- round(cor(plot_data[2], plot_data[3], method = "pearson"),4) ## calculate Pearson correlation coefficient

		  cor_plot <- ggplot(plot_data, aes_string(x = names(plot_data)[2], y = names(plot_data)[3])) +
		    geom_point(aes(color = status)) +	## Color of the dots indicates whether AD or HC
		    geom_smooth(method = "lm", se = FALSE, color = "black")+ ## regression line
		    annotate("text", x = min(plot_data[,2])+0.75, y = max(plot_data[,3])+0.25, label = paste("R=", R))+ ## pearson correlation coefficient
		    theme_minimal() +
		    labs(color = "Disease status") +
		    scale_color_manual(
		      values = c("AD" = "#813513", "HC" = "#32584B"))
		  
		  ggsave(
		    filename = paste0("plots/OPTIMA_corplot_", as.character(cis), "_", ## TODO: change file name
		                      names(plot_data)[2], "_", names(plot_data)[3], ".png"), 
		    plot = cor_plot,
		    width = 6, height = 4)
		  }
	  }
  }

## remove invalid p-values so that an FDR correction can be performed
results <- results[!is.na(results$p),] 
results <- results[results$p != 0,]

## perform FDR correction
results$q <- p.adjust(results$p, method = "fdr")

## Add mRNA or lncRNA
results <- results %>% 
  left_join(
    anno_g %>%                             
      select(transcript_id, transcript_biotype) %>% 
      distinct(transcript_id, .keep_all = TRUE),  
    by = c("rna" = "transcript_id")
  )


## Save results
write.csv(results, out_path, quote = FALSE, row.names = FALSE)


















