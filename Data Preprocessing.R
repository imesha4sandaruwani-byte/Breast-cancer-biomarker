# 1) Load clinical data
clin <- fread("data_clinical_patient.txt",skip=4,header = TRUE, data.table = FALSE)
head(clin)


# 2) Filter ER+/HER2- patients
er_her2_patients <- clin %>%
  filter(ER_IHC == "Positve", HER2_SNP6 == "NEUTRAL") %>%
  select(PATIENT_ID, RFS_MONTHS, RFS_STATUS)

# 3) RFS_STATUS to numeric (0 = Not Recurred, 1 = Recurred)
er_her2_patients$RFS_STATUS <- ifelse(grepl("^1", er_her2_patients$RFS_STATUS), 1,
                                      ifelse(grepl("^0", er_her2_patients$RFS_STATUS), 0, NA))
# 4) Load expression data
expr <- fread("data_mrna_illumina_microarray.txt", data.table = FALSE)

# 5) Keep only columns for ER+/HER2- patients
common_ids <- intersect(er_her2_patients$PATIENT_ID, colnames(expr))
expr_filtered <- expr[, c("Hugo_Symbol", common_ids), drop = FALSE]
expr_filtered <- expr_filtered %>%mutate(across(-Hugo_Symbol, as.numeric))

# 6) Collapse duplicated genes by taking the mean expression
expr_collapsed <- expr_filtered %>%
  group_by(Hugo_Symbol) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 7) Transpose expression to have patients as rows, genes as columns
expr_t <- t(expr_collapsed[,-1])
colnames(expr_t) <- expr_collapsed$Hugo_Symbol
expr_t <- as.data.frame(expr_t)
expr_t$PATIENT_ID <- rownames(expr_t)

# 8) Merge expression with RFS data
er_her2_patients2 <- er_her2_patients %>%filter(PATIENT_ID %in% common_ids)
rfs_expr <- inner_join(er_her2_patients2, expr_t, by = "PATIENT_ID")

# Remove genes with any missing values
expr_clean <- rfs_expr[, c(1:3, which(colSums(is.na(rfs_expr[, -(1:3)])) == 0) + 3)]

#update rfs_expr to the cleaned version
rfs_expr <- expr_clean
rfs_expr$RFS_MONTHS <- as.numeric(rfs_expr$RFS_MONTHS)
rfs_expr$RFS_STATUS <- as.numeric(rfs_expr$RFS_STATUS)

# 9) Save as RDS
saveRDS(rfs_expr, file = "METABRIC_ER_HER2_RFS_expression.rds")

# 10) Load preprocessed data
rfs_expr <- readRDS("METABRIC_ER_HER2_RFS_expression.rds")



dim(rfs_expr)
table(rfs_expr$RFS_STATUS)
head(rfs_expr[, 1:6])