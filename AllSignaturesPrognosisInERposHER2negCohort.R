library(ggplot2)

# Set significance level
sig_level <- 0.01

# Compute p-values for all signatures
real_pvals <- sapply(cancer.signatures.new, function(sig) {
  gene_idx <- which(new.study.med$genes %in% sig$sig)
  if(length(gene_idx) < 2) return(NA)
  e <- list(
    genes = new.study.med$genes[gene_idx],
    data  = new.study.med$data[gene_idx, ],
    survival = new.study.med$survival
  )
  score <- binaryPC1Score(e)
  logRankP(score, e$survival)$survival.index
})


# Clean names: remove ".pvalue" if present
names(real_pvals) <- gsub("\\.pvalue$", "", names(real_pvals))

# Prepare dataframe
plot_df <- data.frame(
  Signature = names(real_pvals),
  pvalue = real_pvals
)
plot_df <- na.omit(plot_df)
plot_df <- plot_df[order(plot_df$pvalue), ]
plot_df$Signature <- factor(plot_df$Signature, levels = plot_df$Signature)

# Add a column for bar color: THR-6E in red
plot_df$color <- ifelse(plot_df$Signature == "THR-6E", "red", "darkgray")

# Plot -log10(p-values)
ggplot(plot_df, aes(x = Signature, y = -log10(pvalue), fill = color)) +
  geom_bar(stat="identity") +
  
  # Significance threshold line
  geom_hline(
    yintercept = -log10(sig_level),
    color = "blue",
    linewidth = 1
  ) +
  
  # Label for threshold
  annotate(
    "text",
    x = length(unique(plot_df$Signature)) + 0.1,
    y = -log10(sig_level) + 0.001,
    label = paste0("log10(", sig_level, ")"),
    color = "blue",
    size = 4,
    fontface = "bold",
    hjust = 1
  ) +
  
  coord_flip() +   # horizontal bars
  scale_fill_identity() +  # use colors from the 'color' column
  labs(
    y = "-log10(p-value)",
    x = "Gene Signature"
  ) +
  theme_minimal(base_size = 14)

# Create BM_Pvalues with signature names and p-values
BM_Pvalues <- plot_df[, c("Signature", "pvalue")]

# View the result
BM_Pvalues
save(BM_Pvalues, file = "BM_Pvalues.RData")

