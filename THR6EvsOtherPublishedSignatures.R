##Comparison of relapse-free survival prognostic performance across published gene expression signatures using PC1-based median stratification.

library(survival)
library(survminer)
library(dplyr)

compute_signature_km <- function(expr_df, genes, time_col="RFS_MONTHS", status_col="RFS_STATUS") {
  genes <- intersect(genes, colnames(expr_df))
  if(length(genes) < 2) {
    warning("Not enough genes found in dataset.")
    return(NULL)
  }
  
  # PCA
  pca <- prcomp(expr_df[, genes], center=TRUE, scale.=TRUE)
  expr_df$PC1_score <- pca$x[,1]
  
  # Median split
  expr_df$risk_group <- ifelse(expr_df$PC1_score >= median(expr_df$PC1_score), "High", "Low")
  
  # Create Surv object inside the dataframe
  expr_df$SurvObj <- Surv(time = expr_df[[time_col]], event = expr_df[[status_col]])
  
  # Fit KM using the SurvObj in the dataframe
  fit <- survfit(SurvObj ~ risk_group, data = expr_df)
  
  # Log-rank test
  logrank <- survdiff(SurvObj ~ risk_group, data = expr_df)
  pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  
  list(
    fit = fit,
    pval = pval,
    expr_df = expr_df
  )
}

sig_6e <- c("CDC20","PIMREG","KIF2C","KIF4A","LMNB2","TPX2")

res_6e <- compute_signature_km(rfs_expr, sig_6e)
res_6e$pval


ggsurvplot(
  res_6e$fit,
  data = res_6e$expr_df,
  risk.table = TRUE,
  pval = TRUE,
  title = "6e Signature"
)

sig_names <- names(cancer.signatures)

results_all <- lapply(sig_names, function(sig_name) {
  genes <- cancer.signatures[[sig_name]]$symb
  res <- compute_signature_km(rfs_expr, genes)
  
  if(!is.null(res)) {
    data.frame(Signature = sig_name, LogRankP = res$pval)
  } else {
    NULL
  }
})


results_all <- do.call(rbind, results_all)
results_all <- results_all %>% arrange(LogRankP)
results_all


# Add 6e to results_all
results_plot <- rbind(
  data.frame(Signature = "6e", LogRankP = 3.74154e-11),
  results_all
)

# Order by p-value
results_plot <- results_plot %>%
  arrange(LogRankP) %>%
  mutate(
    Signature = factor(Signature, levels = Signature),
    negLog10P = -log10(LogRankP),
    Is6e = Signature == "6e"
  )
library(ggplot2)

ggplot(results_plot, aes(x = Signature, y = negLog10P, fill = Is6e)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(
    values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
    labels = c("Other signatures", "6e signature")
  ) +
  labs(
    x = "",
    y = expression(-log[10]("Log-rank p-value")),
    title = "Comparison of prognostic performance across gene signatures"
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.text.y = element_text(size = 9)
  )


ggplot(results_plot, aes(x = seq_along(Signature), y = negLog10P)) +
  geom_point(aes(color = Is6e), size = 3) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  labs(
    x = "Gene signatures (ranked)",
    y = expression(-log[10]("Log-rank p-value")),
    title = "Prognostic strength of gene signatures"
  ) +
  theme_classic()


# Cox model for 6e
cox_6e <- coxph(
  Surv(RFS_MONTHS, RFS_STATUS) ~ risk_group,
  data = res_6e$expr_df
)

summary(cox_6e)

compute_signature_cox <- function(expr_df, genes,
                                  time_col = "RFS_MONTHS",
                                  status_col = "RFS_STATUS") {
  
  genes <- toupper(trimws(genes))
  colnames(expr_df) <- toupper(trimws(colnames(expr_df)))
  
  genes <- intersect(genes, colnames(expr_df))
  if (length(genes) < 2) return(NULL)
  
  # PCA
  pca <- prcomp(expr_df[, genes], center = TRUE, scale. = TRUE)
  expr_df$PC1 <- pca$x[, 1]
  
  # Median split
  expr_df$risk_group <- factor(
    ifelse(expr_df$PC1 >= median(expr_df$PC1), "High", "Low"),
    levels = c("Low", "High")
  )
  
  # Cox model
  cox <- coxph(
    Surv(expr_df[[time_col]], expr_df[[status_col]]) ~ risk_group,
    data = expr_df
  )
  
  s <- summary(cox)
  
  data.frame(
    HR = s$coefficients[,"exp(coef)"],
    Lower95 = s$conf.int[,"lower .95"],
    Upper95 = s$conf.int[,"upper .95"],
    Pvalue = s$coefficients[,"Pr(>|z|)"]
  )
}

sig_names <- names(cancer.signatures)

cox_results <- lapply(sig_names, function(sig_name) {
  genes <- cancer.signatures[[sig_name]]$symb
  res <- compute_signature_cox(rfs_expr, genes)
  
  if (!is.null(res)) {
    res$Signature <- sig_name
    res
  } else {
    NULL
  }
})

cox_results <- do.call(rbind, cox_results)

# Add 6e
cox_6e_res <- compute_signature_cox(
  rfs_expr,
  c("CDC20","PIMREG","KIF2C","KIF4A","LMNB2","TPX2")
)
cox_6e_res$Signature <- "6e"

cox_results <- rbind(cox_6e_res, cox_results)

cox_results <- cox_results %>%
  arrange(desc(HR)) %>%
  mutate(
    Signature = factor(Signature, levels = Signature),
    Is6e = Signature == "6e"
  )

library(ggplot2)

ggplot(cox_results, aes(x = HR, y = Signature)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = Lower95, xmax = Upper95),
    height = 0.25,
    color = "grey40"
  ) +
  geom_point(
    aes(color = Is6e),
    size = 3
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = c("FALSE" = "black", "TRUE" = "firebrick")
  ) +
  labs(
    x = "Hazard Ratio (log scale)",
    y = "",
    title = "Cox proportional hazards analysis of gene signatures",
    subtitle = "High vs Low risk (PC1 median stratification)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 9)
  )





