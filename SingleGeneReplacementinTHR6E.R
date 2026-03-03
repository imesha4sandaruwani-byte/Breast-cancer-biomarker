#############Single-Gene Replacement Stratergy########
library(data.table)
library(ggplot2)

# ---------------------------------------------
# Define random sampling procedure-2 function
# ---------------------------------------------

gene_removal_procedure2 <- function(random_genes_pool, bm_signature) {
  
  results <- list()
  iter <- 1
  
  for (bm_gene in bm_signature) {
    
    remaining_bm <- setdiff(bm_signature, bm_gene)
    
    for (pool_gene in random_genes_pool) {
      
      new_signature <- c(remaining_bm, pool_gene)
      
      results[[iter]] <- list(
        Removed_BM_Gene = bm_gene,
        Added_Pool_Gene = pool_gene,
        Genes = new_signature
      )
      
      iter <- iter + 1
    }
  }
  
  return(results)
}

# ---------------------------------------------
# Create random gene pool
# ---------------------------------------------
random_genes_pool<- bp_related_data$`THR-6E`$Random_Gene_Pool_id

bm_signature      <- cancer.signatures.new$`THR-6E`$sig

# ---------------------------------------------
# Generate samples 
# ---------------------------------------------
all_samples <- gene_removal_procedure2(random_genes_pool, cancer.signatures.new$`THR-6E`$sig)
length(all_samples)


correction <- TRUE   # or FALSE
sig_level  <- 0.01

pvals <- numeric(length(all_samples))
gene_lists <- vector("list", length(all_samples))

for (i in seq_along(all_samples)) {
  
  sig_genes <- as.character(all_samples[[i]]$Genes)
  gene_idx <- which(as.character(new.study.med$genes) %in% sig_genes)
  
  if (length(gene_idx) < 2) {
    pvals[i] <- NA
    next
  }
  
  e <- list(
    genes    = new.study.med$genes[gene_idx],
    data     = new.study.med$data[gene_idx, ],
    survival = new.study.med$survival
  )
  
  score <- score.method(e)
  res   <- survival.association.method(score, e$survival)
  
  pvals[i] <- res$survival.index
  gene_lists[[i]] <- sig_genes
}


# Apply Bonferroni correction
if (correction) {
  pvals <- pvals * length(pvals)
}


SingleGeneReplaced_ResultsBon <- data.frame(
  Genes  = I(gene_lists),
  Pvalue = pvals
)

save(SingleGeneReplaced_Results_NotCorrected, file = "SingleGeneReplaced_Results.RData")
save(SingleGeneReplaced_ResultsBon, file = "SingleGeneReplaced_Results_BonCorrected.RData")

sum(is.na(SingleGeneReplaced_ResultsBon))

# Sort by p-value (smallest first = most significant)
SingleGeneReplaced_ResultsPvalueBon_sorted <- 
  SingleGeneReplaced_ResultsBon[order(SingleGeneReplaced_ResultsBon$Pvalue), ]

SingleGeneReplaced_ResultsPvalueBon_sorted$Sample <- 
  seq_len(nrow(SingleGeneReplaced_ResultsPvalueBon_sorted))


SingleGeneReplaced_ResultsPvalueBon_sorted <- SingleGeneReplaced_ResultsBon[order(SingleGeneReplaced_ResultsBon$Pvalue), ]
SingleGeneReplaced_ResultsPvalueBon_sorted$Sample <- 1:nrow(SingleGeneReplaced_ResultsPvalueBon_sorted)

bm_pvalue=1.24083e-11 

ggplot(SingleGeneReplaced_ResultsPvalueBon_sorted, aes(x = Sample, y = -log10(Pvalue))) +
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
    x = "Index ",
    y = "-log10(p-value)",
    title = "Prognostic Performance: Randomly Single-Gene Replaced Signatures vs THR-6E (sorted by significance)"
  ) +
  
  theme_minimal(base_size = 13) +
  
  # BOTH annotations in same plot
  annotate(
    "text",
    x = max(SingleGeneReplaced_ResultsPvalueBon_sorted$Sample) * 0.7,
    y = -log10(bm_pvalue) + 0.5,
    label = "THR-6E Signature",
    color = "red",
    size = 4
  ) +
  annotate(
    "text",
    x = max(SingleGeneReplaced_ResultsPvalueBon_sorted$Sample) * 0.7,
    y = -log10(sig_level) + 0.5,
    label = "1% significance level",
    color = "blue",
    size = 4
  )

# Identify random samples with p-value <= BM signature
better_than_bm <- SingleGeneReplaced_ResultsPvalue_sorted[
  SingleGeneReplaced_ResultsPvalue_sorted$Pvalue <= bm_pvalue, 
]

# View first few rows
head(better_than_bm)
dim(better_than_bm)

# Find the random sample with the **lowest p-value**
idx_max <- which.min(SingleGeneReplaced_ResultsPvalue_sorted$Pvalue)

# Extract the top random sample
top_sample <- SingleGeneReplaced_ResultsPvalue_sorted[idx_max, ]

# View details
top_sample