# ---- All genes (exclude clinical columns) ----
all_genes<- colnames(rfs_expr)[-(1:3)]
head(all_genes)
length(all_genes)

# ---- BM signature ----
bm_signature <-c("CDC20","PIMREG","KIF2C","KIF4A","LMNB2","TPX2")
#PIMREG=FAM64A
bm_idx <- which(colnames(rfs_expr)%in% bm_signature)
bm_idx
bm_expr <- rfs_expr[, bm_idx] 
dim(bm_expr)

## PCA
pca_bm <- prcomp(bm_expr, center = TRUE, scale. = TRUE)
pc1_bm <- pca_bm$x[, 1]

## Risk groups
median_pc1 <- median(pc1_bm, na.rm = TRUE)
risk_group <- ifelse(pc1_bm <= median_pc1, "Low", "High")

# # Survival object
surv_obj <- Surv(rfs_expr$RFS_MONTHS, rfs_expr$RFS_STATUS)

#Mantel–Haenszel test
fit_bm <- survdiff(surv_obj ~ risk_group)
bm_pvalue <- 1 - pchisq(fit_bm$chisq, df = 1)
bm_pvalue

#Kaplan–Meier curve
surv_data <- data.frame(
  RFS_MONTHS = rfs_expr$RFS_MONTHS,
  RFS_STATUS = rfs_expr$RFS_STATUS,
  risk_group = factor(risk_group, levels = c("Low", "High"))
)

fit_km <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ risk_group, data = surv_data)


#Kaplan–Meier Survival Curve: THR-6E Signature
ggsurvplot(
  fit_km,
  data = surv_data,
  pval = TRUE,
  risk.table = TRUE,
  legend.labs = c("Low risk", "High risk"),
  legend.title = "Risk Group",
  xlab = "Time (months)",
  ylab = "Relapse-free survival probability",
  title = "Kaplan–Meier Survival Curve: THR-6E Signature")