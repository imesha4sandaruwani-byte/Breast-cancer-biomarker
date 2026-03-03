
#THR- 6E Signature prognosis in ER+/HER2- cohort
sig <- cancer.signatures.new$`THR-6E`
sig_genes <- sig$sig
length(sig_genes) 

# Get indices of these genes in your dataset
gene_idx <- which(new.study.med$genes %in% sig_genes)

# Create the expression + survival list
e <- list(
  genes = new.study.med$genes[gene_idx],
  data  = new.study.med$data[gene_idx, ],
  survival = new.study.med$survival
)
score <- binaryPC1Score(e)   
res <- logRankP(score, e$survival)

# Extract the p-value
res$survival.index   # This is the log-rank p-value
res$hr              # Hazard ratio
res$hr.2.5          # 95% CI lower
res$hr.97.5         # 95% CI upper
library(survival)

# Convert score to a factor
group <- factor(score, levels=c(0,1), labels=c("Low", "High"))

# Create a formula
frml <- as.formula(Surv(e$survival[,1], e$survival[,2]) ~ group)

# Plot KM curve
km(frml)
