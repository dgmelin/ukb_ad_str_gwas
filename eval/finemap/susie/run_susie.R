library(data.table)
library(dplyr)
library(susieR)
library(ggplot2)
library(viridis)
library(ggnewscale)

data_path <- "data/finemapping/susie"

check_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

is_symmetric_tolerance <- function(M, tol = 1e-8) {
  max(abs(M - t(M))) < tol
}

is_symmetric_exact <- function(M) {
  identical(M, t(M))
}

force_symmetry <- function(M) {
  (M + t(M)) / 2
}

varpath <- file.path(data_path, "varlists")
loci <- readLines(file.path(varpath, "loci.list"))
  
result_path <- file.path(data_path, "finemapping_results")
check_mkdir(result_path)
  

for (locus in loci) {
  cat("Processing", locus, "\n")
  l <- gsub(":", "_", locus)
    
  # Load LD matrices
  R_df <- as.matrix(fread(file.path(data_path, "ld_scores", paste0(l, "_merged_lds.ld")), header = FALSE))

  if (!is_symmetric_exact(R_df)){
    if (!is_symmetric_tolerance(R_df)){
      warning("LD matrix is really not symmetric")
    } else {
      R_df <- force_symmetry(R_df)
    }
  }
    
  # Load main variant table
  locus_df <- fread(file.path(varpath, paste0(l, "_surroundings.tsv")))
  locus_df <- locus_df %>% distinct(ID, .keep_all = TRUE)
  
  # Order main results table by varlist
    # This should make sure that this matrix is in the same order as the ld-matrix
  varlist <- fread(file.path(data_path, "ld_scores", paste0(l, "_merged_lds.snplist")),
                     header = FALSE, col.names = "ID")
  locus_df <- locus_df %>%
  mutate(id_cat = factor(ID, levels = varlist$ID, ordered = TRUE)) %>%
  arrange(id_cat) %>%
  mutate(cs = 0, pip=0)
    
  # Load LD-to-lead scores
    # This is mainly to double check whther vars in the matrix and results_df are in the correct order
  ld_table <- fread(file.path(data_path, "ld_scores", paste0(l, "_single_r.ld.tsv"))) %>%
  select(ID = SNP_B, ld_to_lead = R)
    
  locus_df <- left_join(locus_df, ld_table, by = "ID")
    
  # Add column for ld_to_lead_matrix (i.e. LD from matrix vs vector)
  lead_idx <- which(locus_df$ID == locus)
  if (length(lead_idx) == 0) {
    warning(paste("Lead SNP", locus, "not found in ID column"))
    next
  }
    
  lead_ld_vector <- R_df[lead_idx, ]
  locus_df$ld_to_lead_matrix <- lead_ld_vector
  locus_df$ld_diff <- locus_df$ld_to_lead - locus_df$ld_to_lead_matrix
    
  if (sum(locus_df$ld_diff > 1e-30, na.rm = TRUE) > 0) {
    warning_n <- sum(locus_df$ld_diff > 1e-30, na.rm = TRUE)
    print(sprintf("Warning: %d variants with ld_diff > 1e-30", warning_n))
    fwrite(locus_df, file.path(result_path, paste0(l, "_susie_finemap_warning.tsv")), sep = "\t")
  }
    
    # Run SuSiE
  fit <- susie_rss(
    bhat = as.matrix(locus_df$BETA),
    shat = as.matrix(locus_df$SE),
    R = R_df,
    L = 10,
    n = 295551
  )

  cs_list <- susie_get_cs(fit, coverage = 0.5, min_abs_corr = 0.5, Xcorr = R_df)$cs

  pip <- susie_get_pip(fit)
    

  # Set pip values for all variants
  for (i in seq_along(pip)) {
      locus_df$pip[i] <- pip[i]
  }

  if (length(cs_list) == 0) {
    print(paste("Nothing found for", locus))
  } else {
    # Set cs indeces for variants in cred sets
    for (i in seq_along(cs_list)) {
      cs_idx <- cs_list[[i]]
      # print(cs_idx)
      if (!is.null(cs_idx)) {
        locus_df$cs[cs_idx] <- i
        # locus_df$pip[cs_idx] <- pip[cs_idx]
      }
    }
  }
  
  # Output Results
  locus_df %>%
    fwrite(file.path(result_path, paste0(l, "_susie_outputs.tsv")), sep = "\t")
    
  # Plot
  p <- ggplot(locus_df, aes(x = POS, y = `-log10`, color = ld_to_lead_matrix)) +
    geom_point(aes(shape = source), size = 2.5, stroke = 0.3) +
    scale_color_viridis(option = "D", name = "LD to lead") +
    new_scale_color() +
    geom_point(data = filter(locus_df, cs > 0),aes(color = factor(cs), shape = source), size = 2, stroke = 0.5) +
    labs(x = "Position", y = "-log10", title = paste("Credible sets for", locus)) +
    theme_bw(base_size = 13)
    
  ggsave(file.path(result_path, paste0(l, "_susie_finemap_r.png")), plot = p, width = 12, height = 5, dpi = 300)
    
  cat("Number of credible sets:", length(cs_list), "\n")
  }

