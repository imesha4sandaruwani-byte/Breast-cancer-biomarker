#Define the full THR-6E Cox model (baseline)

# PCA on full 6-gene signature
pca_full <- prcomp(
  rfs_expr[, bm_signature, drop = FALSE],
  scale. = TRUE
)

df_full <- data.frame(
  RFS_MONTHS = rfs_expr$RFS_MONTHS,
  RFS_STATUS = rfs_expr$RFS_STATUS,
  THR6E_score = pca_full$x[, 1]
)

# risk score from PCA (full 6 genes)
cox_full <- coxph(
  Surv(RFS_MONTHS, RFS_STATUS) ~ THR6E_score,
  data = df_full
)

summary(cox_full)

#HR
#p-value

loo_results <- list()

for (g in bm_signature) {
  # 1. Remove one gene
  genes_loo <- setdiff(bm_signature, g)
  
  # 2. PCA on remaining 5 genes
  pca_loo <- prcomp(
    rfs_expr[, genes_loo, drop = FALSE],
    scale. = TRUE
  )
  
  df_loo <- data.frame(
    RFS_MONTHS  = rfs_expr$RFS_MONTHS,
    RFS_STATUS = rfs_expr$RFS_STATUS,
    LOO_score = pca_loo$x[, 1]
  )
  
  # Cox model
  cox_loo <- coxph(
    Surv(RFS_MONTHS, RFS_STATUS) ~ LOO_score,
    data = df_loo
  )
  
  s <- summary(cox_loo)
  
  loo_results[[g]] <- data.frame(
    Removed_Gene = g,
    HR = s$coefficients[,"exp(coef)"],
    CI_lower = s$conf.int[,"lower .95"],
    CI_upper = s$conf.int[,"upper .95"],
    p_value = s$coefficients[,"Pr(>|z|)"]
  )
}

loo_results <- do.call(rbind, loo_results)
loo_results
