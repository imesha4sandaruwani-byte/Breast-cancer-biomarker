source("./codes/main.R", echo=T)

real_pvals <- sapply(cancer.signatures.new, function(sig){
  
  gene_idx <- which(new.study.med$genes %in% sig$sig)
  
  if(length(gene_idx) < 2) return(NA)
  
  e <- list(
    genes = new.study.med$genes[gene_idx],
    data  = new.study.med$data[gene_idx, ],
    survival = new.study.med$survival
  )
  
  score <- score.method(e)
  
  logRankP(score, e$survival)$survival.index
})

all <- data.frame(p = real_pvals)

rownames(all) <- names(real_pvals)
head(all)
colnames(all)
rownames(all) <- gsub("\\.pvalue$", "", rownames(all))


geneRP1 <- function(GEdata, signatures, study, BP_data, correction = TRUE){
  p.val <- matrix(nrow=nsig, ncol=length(signatures))
  result_all <- list()  
  g.list <- list()
  
  # Removing biomarker genes and GO-terms that are on a level
  set.seed(5)
  for (j in 1:length(signatures)) {
    v <- c()
    v.result <- c()
    sig.size <- sum(GEdata$genes %in% study[[j]])
    
    print(paste("Signature:", signatures[[j]]$name))
    print(paste("Signature size:", sig.size))
    print(paste("Remain pool size:", length(BP_data[[j]]$Random_Gene_Pool_id)))
    
    for(k in 1:nsig){
      a1 <- sample(1:length(BP_data[[j]]$Random_Gene_Pool_id), sig.size)
      b <- selectRandomGene(GEdata, a1 , BP_data[[j]]$Random_Gene_Pool_id)
      
      # Diagnostic check
      if(any(!is.finite(as.matrix(b$data)))) stop("Non-finite values detected")
      if(any(apply(b$data, 1, var) < 1e-12)) stop("Zero or near-zero variance gene detected")
      
      score <- do.call(score.method, list(b))
      x <- do.call(survival.association.method, list(score, b$survival))$survival.index
      v.result <- c(v.result, x)
      
    }
    a <- v.result
    
    # store the pvalues obtained from the 1000 random selection of genes of the same size as the BM. 
    
    
    # Bonferroni correction
    if(correction){
      pvalue <- 0.05/nsig
      p.val[,j] <- a * nsig
    }else{
      pvalue <- 0.05
      p.val[,j] <- a 
    }
    
    v <- c(v, sum(a < pvalue))
    
    g.list[[signatures[[j]]$name]] = BP_data[[j]]$Random_Gene_Pool_id
    result_all[[signatures[[j]]$name]] <- list("Accuracy" = v/nsig)
    
    print(j)
  }
  colnames(p.val) <- names(signatures)
  
  outcome <- list("results" = result_all, "pvalues" = p.val)
  return(outcome)
  
}

# Analysis 
cancer.sig <- study.sig(cancer.signatures.new)

library(GOxploreR)
# without bonferroni correction¨

NotCorrectedResults<- geneRP1(GEdata = new.study.med, signatures = cancer.signatures.new, study =cancer.sig, BP_data =bp_related_data, correction = FALSE)
saveRDS(NotCorrectedResults, "NotCorrectedResults.rds")

# Accuracy
sum.outcomeGRP1(NotCorrectedResults$results, cancer.signatures.new)

# Visualisation of the results
library(beanplot)
tmp <- arrange.sig(cancer.signatures.new, NotCorrectedResults$pvalues)

## save 
pdf("NotCorrectedResults.pdf", w=6, h=12)
myvisplot(tmp)
dev.off()

# with bonferroni correction
CorrectedResults<- geneRP1(GEdata = new.study.med, signatures = cancer.signatures.new, study =cancer.sig, BP_data =bp_related_data, correction = TRUE)
saveRDS(CorrectedResults, "CorrectedResults.rds")

# Accuracy
sum.outcomeGRP1(CorrectedResults$results, cancer.signatures.new)

# Visualisation of the results
library(beanplot)
tmp2 <- arrange.sig(cancer.signatures.new, bon.correction.corrected$pvalues)

## save 
pdf("CorrectedResults.pdf", w=6, h=12)
myvisplot(tmp2)
dev.off()


##quit
#sessionInfo()
#q(save="no")
