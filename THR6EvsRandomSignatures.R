#############Random Gene Set Analysis 1: GO-Filtered Random Gene Sets########

#Map signature genes to their GO terms
bm_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = bm_signature,
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID

gene_to_go_all <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = bm_entrez,
  keytype = "ENTREZID",
  columns = c("GO", "ONTOLOGY", "EVIDENCE")
)
dim(gene_to_go_all)

# Remove genes without GO annotations 
gene_to_go_all <- gene_to_go_all[!is.na(gene_to_go_all$GO), ]
dim(gene_to_go_all)

#Map GO terms back to genes
go_terms <- unique(gene_to_go_all$GO)
length(go_terms)

go_to_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_terms,
  keytype = "GO",
  columns = c("ENTREZID", "SYMBOL")
)$ENTREZID

length(go_to_genes)

#Convert go related genes ENTREZID → SYMBOL 
go_related_gene_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_to_genes,
  keytype = "ENTREZID",
  columns = "SYMBOL"
)$SYMBOL

length(go_related_gene_symbols)
go_related_gene_symbols<- unique(na.omit(go_related_gene_symbols))

# ============================================================
# Remove signature genes from all_genes
# ============================================================
non_sig_genes <- setdiff(all_genes, bm_signature)
length(non_sig_genes)

# ============================================================
# Remove all genes linked to GO terms
# ============================================================
random_genes_pool1 <- setdiff(non_sig_genes, go_related_gene_symbols)
length(random_genes_pool1)

# ---------------------------------------------
# Define random sampling procedure-1 function
# ---------------------------------------------

gene_removal_procedure1<- function(random_genes_pool) {
  
  # Generate maximum number of non-overlapping random sets
  results_list <- list()
  iter <- 1
  sig_size <- 6
  
  tmp_pool <- random_genes_pool
  
  while (length(tmp_pool) >= sig_size) {
    
    sampled <- sample(tmp_pool, sig_size)
    
    results_list[[iter]] <- data.table(
      Iteration = iter,
      Genes = list(sampled)   
    )
    
    tmp_pool <- setdiff(tmp_pool, sampled)
    iter <- iter + 1
  }
  
  final_dt <- rbindlist(results_list)
  return(final_dt)
}



# Number of repetitions
n_reps <- 1000

# Store results for all iterations
all_random_pvalues <- vector("list", n_reps)

for (rep in 1:n_reps) {
  
  # ---- 1) Generate random samples ----
  random_samples <- gene_removal_procedure1(random_genes_pool1)
  
  # Initialize p-value vector for this repetition
  random_pvalues <- numeric(nrow(random_samples))
  
  for(i in 1:nrow(random_samples)) {
    gene_set <- random_samples$Genes[[i]] 
    gene_idx <- which(colnames(rfs_expr) %in% gene_set)
    
    if (length(gene_idx) != length(bm_signature)) {
      random_pvalues[i] <- NA
      next
    }
    
    expr <- rfs_expr[, gene_idx]
    
    # PCA
    pca <- prcomp(expr, center = TRUE, scale. = TRUE)
    pc1 <- pca$x[, 1]
    
    # Split by median
    median_pc1 <- median(pc1, na.rm = TRUE)
    risk_group <- ifelse(pc1 <= median_pc1, "Low", "High")
    
    # Survival
    surv_obj <- Surv(rfs_expr$RFS_MONTHS, rfs_expr$RFS_STATUS)
    
    # MH test
    fit <- survdiff(surv_obj ~ risk_group)
    random_pvalues[i] <- 1 - pchisq(fit$chisq, df = 1)
  }
  
  # Store this repetition's p-values
  all_random_pvalues[[rep]] <- random_pvalues
}

# Combine all into a data frame
pval_df_all <- bind_rows(
  lapply(1:n_reps, function(i) {
    data.frame(
      Rep = i,
      Sample = 1:length(all_random_pvalues[[i]]),
      Pvalue = all_random_pvalues[[i]]
    )
  })
)


#ggplot(pval_df_all, aes(x = Sample, y = -log10(Pvalue))) +
# geom_point(color = "darkgreen") +
#geom_hline(yintercept = -log10(bm_pvalue), color = "red", linetype = "dashed", size = 1) +
#labs(
# x = "Random Sample Index",
#y = "-log10(p-value)",
#title = "Random Gene Sets vs BM Signature"
#  ) +
# theme_minimal() +
#annotate("text", x = max(pval_df_all$Sample)*0.7, y = -log10(bm_pvalue)+0.5,
#        label = "BM signature", color = "red")

ggplot(pval_df_all, aes(x = "", y = -log10(Pvalue))) +
  geom_violin(
    fill = "lightgreen",
    alpha = 0.6,
    trim = FALSE
  ) +
  geom_jitter(
    width = 0.08,
    alpha = 0.08,
    size = 0.4,
    color = "darkgreen"
  ) +
  geom_hline(
    yintercept = -log10(bm_pvalue),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = "Prognostic Performance of GO-Filtered Random Gene Signatures versus THR-6E",
    y = expression(-log[10]*"(MH test p-value)"),
    x = "GO-filtered random 6-gene signatures"
  ) +
  annotate(
    "text",
    x = 1.15,
    y = -log10(bm_pvalue) + 0.4,
    label = "THR-6E signature",
    color = "red",
    size = 4
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.line.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggplot(pval_df_all, aes(x = "", y = -log10(Pvalue))) +
  geom_violin(
    fill = "lightgreen",
    alpha = 0.6,
    trim = FALSE
  ) +
  geom_jitter(
    width = 0.08,
    alpha = 0.05,
    size = 0.3,
    color = "darkgreen"
  ) +
  geom_hline(
    yintercept = -log10(bm_pvalue),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = "",                 # ← keep discrete
    y = -log10(bm_pvalue) + 0.4,
    label = "THR-6E signature",
    color = "red",
    size = 4,
    hjust = -0.15
  ) +
  labs(
    title = "GO-Filtered Random Gene Signatures versus THR-6E",
    y = expression(-log[10]*"(MH test p-value)"),
    x = "GO-filtered random 6-gene signatures"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.line.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


##Bonferroni correction##
m_tests <- nrow(pval_df_all)

pval_df_all_Bon <- pval_df_all %>%
  mutate(
    Pvalue_Bonferroni = p.adjust(Pvalue, method = "bonferroni"),
    NegLog10_Bonferroni = -log10(Pvalue_Bonferroni)
  )

# Plot with raw BM p-value
ggplot(pval_df_all_Bon, aes(x = "", y = NegLog10_Bonferroni)) +
  geom_violin(fill = "lightblue", alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.08, alpha = 0.08, size = 0.4, color = "steelblue") +
  geom_hline(
    yintercept = -log10(bm_pvalue),  # <-- raw BM p-value
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = "Prognostic Performance of GO-filtered Random Gene signatures",
    y = expression(-log[10]*"(Bonferroni-corrected p-value)"),
    x = "GO-filtered random 6-gene signatures"
  ) +
  annotate(
    "text",
    x = 1.15,
    y = -log10(bm_pvalue) + 0.4,
    label = "THR-6E signature",
    color = "red"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.line.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )



#############Random Gene Set Analysis 2: Single-Gene Replacement Strategy########
# ---------------------------------------------
# Create random gene pool
# ---------------------------------------------
random_genes_pool2 <- setdiff(all_genes, bm_signature)
length(random_genes_pool2)

# ---------------------------------------------
# Define random sampling procedure-2 function
# ---------------------------------------------

gene_removal_procedure2 <- function(random_genes_pool2,bm_signature) {
  
  results <- list()
  iter <- 1
  
  for (bm_gene in bm_signature) {
    
    # BM genes except the one being removed
    remaining_bm <- setdiff(bm_signature, bm_gene)
    
    for (pool_gene in random_genes_pool2) {
      
      new_signature <- c(remaining_bm, pool_gene)
      
      results[[iter]] <- data.table(
        Removed_BM_Gene = bm_gene,
        Added_Pool_Gene = pool_gene,
        Genes = list(new_signature)
      )
      
      iter <- iter + 1
    }
  }
  
  rbindlist(results)
}

# ---------------------------------------------
# Generate random samples 
# ---------------------------------------------
##hare pool 1 can be 2
random_samples <- gene_removal_procedure2(random_genes_pool2, bm_signature)
dim(random_samples)

# Initialize p-value vector for this repetition
random_pvalues <- numeric(nrow(random_samples))

for(i in 1:nrow(random_samples)) {
  gene_set <- random_samples$Genes[[i]] 
  gene_idx <- which(colnames(rfs_expr) %in% gene_set)
  
  # Skip if some genes are missing
  if (length(gene_idx) != length(bm_signature)) {
    random_pvalues[i] <- NA
    next
  }
  
  expr <- rfs_expr[, gene_idx]
  
  # PCA
  pca <- prcomp(expr, center = TRUE, scale. = TRUE)
  pc1 <- pca$x[, 1]
  
  # Split by median
  median_pc1 <- median(pc1, na.rm = TRUE)
  risk_group <- ifelse(pc1 <= median_pc1, "Low", "High")
  
  # Survival analysis
  surv_obj <- Surv(rfs_expr$RFS_MONTHS, rfs_expr$RFS_STATUS)
  
  # Log-rank test
  fit <- survdiff(surv_obj ~ risk_group)
  pval <- 1 - pchisq(fit$chisq, df = 1)
  
  random_pvalues[i] <- pval
}


# Combine all into a data frame
pval_df_all2 <- data.frame(
  Sample = 1:length(random_pvalues),
  Pvalue = random_pvalues
)

save(pval_df_all2, file = "pval_df_all_analysis2_withgo.RData")
saveRDS(pval_df_all2, file = "pval_df_all_analysis2_withgo.rds")

pval_df_all2_sorted <- pval_df_all2[order(pval_df_all2$Pvalue), ]
pval_df_all2_sorted$Sample <- 1:nrow(pval_df_all2_sorted)

ggplot(pval_df_all2_sorted, aes(x = Sample, y = -log10(Pvalue))) +
  geom_point(color = "darkgreen", alpha = 0.7, size = 1.8) +
  
  # THR-6E line
  geom_hline(
    yintercept = -log10(bm_pvalue),
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  
  # 1% significance level line
  geom_hline(
    yintercept = -log10(sig_level),
    color = "blue",
    linetype = "dashed",
    size = 1
  ) +
  
  labs(
    x = "Random Gene Signature Index (sorted by significance)",
    y = "-log10(p-value)",
    title = "Prognostic Performance: Randomly Single-Gene Replaced Signatures vs THR-6E"
  ) +
  
  theme_minimal(base_size = 13) +
  
  # BOTH annotations in same plot
  annotate(
    "text",
    x = max(pval_df_all2_sorted$Sample) * 0.7,
    y = -log10(bm_pvalue) + 0.5,
    label = "THR-6E Signature",
    color = "red",
    size = 4
  ) +
  annotate(
    "text",
    x = max(pval_df_all2_sorted$Sample) * 0.7,
    y = -log10(sig_level) + 0.5,
    label = "1% significance level",
    color = "blue",
    size = 4
  )





#####identify the random samples with p-value >= BM signature ########

# Combine p-values with random_samples
random_samples$Pvalue <- random_pvalues

# Identify sets with p-value <= BM signature
above_bm <- random_samples[random_samples$Pvalue <= bm_pvalue, ]

# View the first few rows
head(above_bm)
dim(above_bm)

# Find the index of the highest p-value
idx_max <- which.max(random_samples$Pvalue)

# Extract the top random sample
top_sample <- random_samples[idx_max, ]

# View the information
top_sample

##Bonferroni correction##
m_tests2 <- nrow(pval_df_all2_sorted)
pval_df_all2_sorted <- pval_df_all2_sorted %>%
  mutate(
    Pvalue_Bonferroni = p.adjust(Pvalue, method = "bonferroni")
  )

# BM signature Bonferroni-corrected

ggplot(pval_df_all2_sorted, aes(x = Sample, y = -log10(Pvalue_Bonferroni))) +
  geom_point(color = "darkgreen", alpha = 0.6) +
  
  geom_hline(
    yintercept = -log10(bm_pvalue),   # THR-6E raw p-value
    color = "red",
    linetype = "dashed",
    size = 1
  ) +
  
  geom_hline(
    yintercept = -log10(0.01),
    color = "blue",
    linetype = "dashed",
    size = 1
  ) +
  
  labs(
    x = "Random Gene Signature Index (sorted by significance)",
    y = expression(-log[10]*"(Bonferroni-corrected p-value)"),
    title = "Prognostic Performance: Randomly Single-Gene Replaced Signatures vs THR-6E"
  ) +
  
  annotate(
    "text",
    x = max(pval_df_all2_sorted$Sample)*0.75,
    y = -log10(bm_pvalue) + 0.2,
    label = "THR-6E Signature",
    color = "red"
  ) +
  
  annotate(
    "text",
    x = max(pval_df_all2_sorted$Sample)*0.75,
    y = -log10(0.01) + 0.2,
    label = "1% significance level",
    color = "blue"
  ) +
  
  theme_minimal(base_size = 13)
