
library(amap)
library(survival)

## list of signatures

study.sig <- function(signature){
  sig_cancer <- list()
  # Breast cancer signatures
  for(i in 1:length(signature)){
    sig_cancer[[i]] <- signature[[i]]$sig
  }
  return(sig_cancer)
}

## summarize the outcome association results for gene removal procedure 1

# summarising the results of GRP1 
sum.outcomeGRP1 <- function(result, signatures){
  sum.biomarker <- data.frame(nrow=length(signatures),ncol=2)
  colnames(sum.biomarker) <- c("Signature", "Acc.")
  for (i in 1:length(signatures)) {
    sum.biomarker[i,1] <- signatures[[i]]$name
    sum.biomarker[i,2] <- result[[signatures[[i]]$name]]$Accuracy * 100
  }
  return(sum.biomarker)
}

# summarising the results of GRP2
sum.outcomeGRP2 <- function(result, signatures){
  n <- names(signatures)
  sum.biomarker <- list()
  for (i in 1:length(n)) {
    dat <- matrix(data = NA, nrow = length(result[[n[i]]]$Accuracy), ncol = 7)
    dat[,1] <- result[[n[i]]]$hierarchy
    dat[,2] <- result[[n[i]]]$genesremoved
    cum_sum <- cumsum(result[[n[i]]]$genesremoved)
    dat[,3] <- cum_sum
    dat[,4] <- result[[n[i]]]$genesleft
    dat[,5] <- result[[n[i]]]$gotermsremoved
    dat[,6] <- cumsum(result[[n[i]]]$gotermsremoved)
    dat[,7] <- result[[n[i]]]$Accuracy * 100
    colnames(dat) <- c("Hierachy.level", "Genes.removed","Genes.Cum.sum", "Genes.left","GO-terms.removed","GO-terms.cum.sum", "Acc.")
    dat <- as.data.frame(dat)
    sum.biomarker[[n[i]]] <- dat
  }
  
  ####
  return(sum.biomarker)
}

## Arrange the signatures according to the size

arrange.sig <- function(signatures, p.val){
  all[, "p"] <- log10(all[, "p"])
  gene.outcome <- rbind(all[colnames(p.val), "p"], log10(p.val))
  l <- sapply(signatures, function(z) length(z$sig))
  i <- sort(l, dec=T, index.return=T)$ix
  gene.outcome <- gene.outcome[, i]
  return(gene.outcome)
}

## Select random samples "Our method"

selectRandomGene <- function(e, n, x,...){
  ii <- x[n]
  i <- which(e$genes %in% ii)
  e$data <- e$data[i,]
  e$genes <- e$genes[i]
  return(e)
}

## Visualisation

# green = col=rgb(20,75,0,130,maxColorValue=255)
# dark green = rgb(0,51,0,maxColorValue=255)
# purple = col=rgb(153,10,255, maxColorValue=255)
# orange = rgb(220, 160, 20, maxColorValue=255)
myvisplot <- function(x, ...)
{
  par(mar=c(8, 12, 3, 0.6) + 0.1) # mar=c(bottom, left, top, right)
  par(mgp=c(4, 1, 0))
  
  beanplot(as.data.frame(x)[-1,], range=0, horizontal=TRUE, pch=".", cut = -0.1, ylim = c(min(x) - 1, 3),
           what=c(0,1,1,0),  overallline = "median", 
           xlab=expression("p-value " (log[10])),
           lwd=1.5, las=1, cex.lab=1.5, cex.axis=1.0, col=rgb(153,10,255, maxColorValue=255), border=1, ...)
  
  for (i in 1:ncol(x)) {
     q <- quantile(x[-1, i], 0.03)   #rgb(46,139,87,maxColorValue=255) rgb(0,51,0,maxColorValue=255)
     rect(min(x[-1, i]), i-0.2, q,i+0.3, border=NA, col =rgb(0,51,0,maxColorValue=255))
     print(colnames(x)[i])
     print(q)
   }
  
  for (k in 1:ncol(x)) {
    if((median(x[-1,k]) < log10(0.05)) && (x[1,k] < log10(0.05))){
      points(x[1,k], k, pch=19, col="green", cex=2)
    }else if((median(x[-1,k]) < log10(0.05)) && (x[1,k] > log10(0.05))){
      points(x[1,k], k, pch=19, col="darkgreen", cex=2)
    }else if((median(x[-1,k]) > log10(0.05)) && (x[1,k] > log10(0.05))){
      points(x[1,k], k, pch=19, col="darkred", cex=2)
    }else{
      points(x[1,k], k, pch=19, col="red", cex=2)
    }
    #print(median(x[-1,k]))
  }
  
  abline(v=log10(0.05), lwd=3, col="blue")
  axis(side=3, at=log10(0.05), labels=expression(p=log[10](0.05)),
       tick=T, col="blue", cex.axis=1.5, col.axis="blue")
  return(q)
}

Acc.Lowest.level <- function(x, vbon, sig, ...)
{
  par(mar=c(8, 12, 3, 0.6) + 0.1) # mar=c(bottom, left, top, right)
  par(mgp=c(4, 1, 0))
  
  l <- sapply(sig, function(z) length(z$sig))
  l
  i <- sort(l, dec=T, index.return=T)$ix
  x <- x[, i]
  
  
  beanplot(as.data.frame(x), range=0, horizontal=TRUE, pch=".",
           what=c(0,0,0,0),  overallline = "median",
           xlab=expression("Accuracy(%)"), 
           ylim = c(0,100),
           lwd=1.5, las=1, cex.lab=1.5, cex.axis=1.0, col=rgb(178, 58, 238, maxColorValue=255), 
           border=1, ...)
  
  vbon <- vbon[i]
  for (i in 1:ncol(x)) {
    q <- vbon[i]
    if(q < 5){
      rect(-25, i-0.3, q,i+0.3, border=NA, col=rgb(218,165,32,50,maxColorValue=255))
    }
    else{
      rect(-25, i-0.3, q,i+0.3, border=NA, col=rgb(218,165,32,maxColorValue=255))
      
    }
  }
  
  abline(v= 5, lwd=3, col="red")
  axis(side=3, at=5, labels=expression("5%"),
       tick=T, col="red", cex.axis=1.5, col.axis="red")
}

## The first principal component method
PC1Score <- function(e)
{
  pc1 <- prcomp(t(e$data))$x[,1]
  return(pc1)
}

## score = binary partition based on first principal component
binaryPC1Score <- function(e)
{
  pc1 <- PC1Score(e)
  m <- median(pc1)
  part <- as.numeric(pc1 > m)
  return(part)
}

###########################################################
##measuring association of score with survival
##
## output is a study object with a sample-wise score
## and a study-wise association score

logRankP <- function(score, survival, ...)
{
  c <- summary(coxph(survival ~ score, ...))
  ci <- c$conf.int
  #print(c$coef[1])
  if (c$coef[1]<0) {
    ci <- 1/ci
    x <- ci[3]
    ci[3] <- ci[4]
    ci[4] <- x
    score = 1 - score ##not great formally for non-binary non-normalized scores
    coef <- -c$coef[1]
  }
  
  return(list(survival.index= c$logtest["pvalue"], survival.index.name="p-value",
              hr=ci[1], hr.2.5=ci[3], hr.97.5=ci[4],
              n=c$n, logrank.p=log10(c$logtest["pvalue"]), ebeta=coef, score=score))
}


##run Cox analysis
cox <- function(frml, ...)
{
  c <- summary(coxph(frml, ...))
  ci <- c$conf.int
  if (c$coef[1]<0) {##present HR in term of bad vs. good prognosis
   ci <- 1/ci
    x <- ci[3]
    ci[3] <- ci[4]
    ci[4] <- x
  }
  return(c(rr=ci[1], ci.low=ci[3], ci.high=ci[4], n=c$n, p=c$logtest["pvalue"], ebeta=c$coef[1]))
}

## Function to annotated Kaplan-Meier plot
km <- function(frml, lab="", vl=c(5, 10), display=TRUE,
               cex.rr=1.5, cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
               mar=c(5.1 , 5.3, 4.1, 1.1), sig=NULL, ...)
{
  c <- cox(frml)
  
  # Color: bad prognosis in red
  if (c["ebeta"] > 0) {
    col <- c("blue", "red")
  } else {
    col <- c("red", "blue")
  }
  
  if (display) {
    plot(survfit(frml), 
         col = col, 
         main = lab, 
         xlab = "Time (months)", 
         ylab = "Relapse-free survival probability", 
         mark.time = TRUE, 
         mark = 3,
         cex.axis = cex.axis, 
         cex.lab = cex.lab, 
         cex.main = cex.main, 
         mar = mar, ...)
    
    # Add HR and p-value text (kept as in your original code)
    text(0.25 * cex.rr, 0.07 * cex.rr,
         paste("HR=", signif(c["rr"], 2), " (CI, ", signif(c["ci.low"], 2), "-",
               signif(c["ci.high"], 2), ")\np=", signif(c["p.pvalue"], 2), sep=""), 
         pos = 4, cex = cex.rr)
  }
  
  return(c)
}


