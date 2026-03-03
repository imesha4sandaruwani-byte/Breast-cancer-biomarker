library(survival)
library(data.table)

# -----------------------------
# THR-6E 6 genes
# -----------------------------
thr6e_ids <- cancer.signatures.new$`THR-6E`$sig    
thr6e_sym <- cancer.signatures.new$`THR-6E`$symb  

# Map IDs → symbols
id2sym <- setNames(thr6e_sym, thr6e_ids)

# -----------------------------
# Generate all 5-gene combinations
# -----------------------------
gene_combinations <- combn(thr6e_ids, 5, simplify = FALSE)

# -----------------------------
# Loop over combinations and compute p-values
# -----------------------------
results <- list()

for (i in seq_along(gene_combinations)) {
  
  sig_ids <- gene_combinations[[i]]
  sig_sym <- id2sym[as.character(sig_ids)]  
  
  idx <- which(as.character(new.study.med$genes) %in% sig_ids)
  if(length(idx) < 2) next
  
  b <- list(
    genes = new.study.med$genes[idx],
    data  = new.study.med$data[idx, ],
    survival = new.study.med$survival
  )
  
  # Diagnostic checks
  if(any(!is.finite(as.matrix(b$data)))) next
  if(any(apply(b$data, 1, var) < 1e-12)) next
  if(length(unique(binaryPC1Score(b))) < 2) next
  
  # Compute score
  score <- binaryPC1Score(b)
  
  # Cox / log-rank test
  res <- tryCatch(logRankP(score, b$survival), error = function(e) NULL)
  if(is.null(res)) next
  
  # Save results with gene symbols
  results[[i]] <- list(
    Genes = paste(sig_sym, collapse = ","),
    Pvalue = res$survival.index
  )
}

# -----------------------------
# Convert to data.frame
# -----------------------------
df_results <- rbindlist(lapply(results, function(x) data.table(
  Genes = x$Genes,
  Pvalue = x$Pvalue
)))

# Sort by p-value
df_results <- df_results[order(Pvalue)]

# View
print(df_results)

# Save
save(df_results, file="THR6E_all5gene_combinations_symbols.RData")

#####################################################################3
#single gene
# Generate single-gene combinations (each gene separately)
gene_combinations <- combn(thr6e_ids, 1, simplify = FALSE)

# -----------------------------
# Loop through each gene and compute p-values
# -----------------------------
results <- list()

for (i in seq_along(gene_combinations)) {
  
  gene_id <- gene_combinations[[i]]
  gene_sym <- id2sym[as.character(gene_id)]
  
  # Find expression data index
  idx <- which(as.character(new.study.med$genes) %in% gene_id)
  if(length(idx) < 1) next  # skip if gene not found
  
  b <- list(
    genes = new.study.med$genes[idx],
    data  = new.study.med$data[idx, , drop=FALSE],  # keep as matrix
    survival = new.study.med$survival
  )
  
  # Diagnostic checks
  if(any(!is.finite(as.matrix(b$data)))) next
  if(apply(b$data, 1, var) < 1e-12) next  # skip near-zero variance
  if(length(unique(binaryPC1Score(b))) < 2) next  # skip constant score
  
  # Compute score and Cox log-rank p-value
  score <- binaryPC1Score(b)
  res <- tryCatch(logRankP(score, b$survival), error = function(e) NULL)
  if(is.null(res)) next
  
  # Save result
  results[[i]] <- list(
    Gene = gene_sym,
    Pvalue = res$survival.index
  )
}

# -----------------------------
# Remove any NULLs
# -----------------------------
results <- Filter(Negate(is.null), results)

# -----------------------------
# Convert to data.table
# -----------------------------
df_results <- rbindlist(lapply(results, function(x) data.table(
  Gene = x$Gene,
  Pvalue = x$Pvalue
)))

# Sort by p-value
df_results <- df_results[order(Pvalue)]

# View results
print(df_results)

# Save results
save(df_results, file="THR6E_single_gene_results.RData")