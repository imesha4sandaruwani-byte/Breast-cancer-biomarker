library(GOxploreR)
library(reshape2)
library(ggplot2)

# total number of genes in each study
dat <- data.frame()
for (j in 1:length(cancer.signatures.new)) {
  dat[j,1] <- cancer.signatures.new[[j]]$name
  dat[j,2] <- length(cancer.signatures.new[[j]]$sig)
}
colnames(dat) <- NULL
colnames(dat) <- c("study", "genes")

#----------------------------------

ind <- order(dat$genes)
df <- data.frame(x1 = dat$study, y1 = dat$genes)
df$x1 <- factor(df$x1, levels = df$x1[ind])


theme_dotplot <-
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        legend.position =  "none")

p <- ggplot(df, aes(x=x1, y=y1)) + theme_dotplot +
  coord_flip() +
  ylab("Number of BM genes") +
  xlab("Studies") +
  geom_point(color="blue")

p

## Number of GO-terms in each study
bp_level <- list()
mf_level <- list()
cc_level <- list()
Biomarker.d <- data.frame()
k <- 1
for (i in 1:length(cancer.signatures.new)) {
  bp <- Gene2GOTermAndLevel(gene = cancer.signatures.new[[i]]$sig, organism = "Human", domain = "BP" )
  mf <- Gene2GOTermAndLevel(gene = cancer.signatures.new[[i]]$sig, organism = "Human", domain = "MF" )
  cc <- Gene2GOTermAndLevel(gene = cancer.signatures.new[[i]]$sig, organism = "Human", domain = "CC" )
  
  bp_level[[cancer.signatures.new[[i]]$name]] <- bp$Level
  mf_level[[cancer.signatures.new[[i]]$name]] <- mf$Level
  cc_level[[cancer.signatures.new[[i]]$name]] <- cc$Level
  
  Biomarker.d[k,1] <- cancer.signatures.new[[i]]$name
  Biomarker.d[k,2] <- length(unique(bp$`GO ID`))
  Biomarker.d[k,3] <- length(unique(mf$`GO ID`))
  Biomarker.d[k,4] <- length(unique(cc$`GO ID`))
  k <- k + 1
  print(i)
  
}
colnames(Biomarker.d) <- c("Study", "BP", "MF", "CC")
Biomarker_d <- Biomarker.d
ind <- order(Biomarker.d$BP)
Biomarker.d$Study <- factor(Biomarker.d$Study, levels = Biomarker.d$Study[ind])

p1 <- ggplot(Biomarker.d, aes(x=Study)) + theme_dotplot +
  coord_flip() +
  ylab("Number of GO-terms") +
  xlab("Studies") +
  geom_point(aes(y = BP),color="green") + geom_point(aes(y = MF),color="red") +
  geom_point(aes(y = CC),color="blue")
p1

# Distribution of GO-terms of BP hierarchy levels for each study
final_dat_bp <- data.frame()
name_bp <- names(bp_level)

for(k in 1:length(bp_level)){
  len <- unlist(bp_level[k])
  x <- rep(name_bp[k], length(len))
  dat <- as.data.frame(bp_level[[k]])
  colnames(dat) <- NULL
  #dat <- dat[-length(dat),]
  dat <- cbind(x, dat)
  final_dat_bp <- rbind(final_dat_bp,dat)
  
}


final_dat_bp$x <- factor(final_dat_bp$x)

p2 <- ggplot(final_dat_bp, aes(x = x, y = dat)) + theme_dotplot + 
  ylab("Number of levels") +
  xlab("Studies") + ggtitle("BP") + geom_boxplot(aes(fill = x)) + coord_flip()
p2

# Distribution of GO-terms of MF hierarchy levels for each study
final_dat_mf <- data.frame()
name_mf <- names(mf_level)

for(k in 1:length(mf_level)){
  len <- unlist(mf_level[k])
  x <- rep(name_mf[k], length(len))
  dat <- as.data.frame(mf_level[[k]])
  colnames(dat) <- NULL
  #dat <- dat[-length(dat),]
  dat <- cbind(x, dat)
  final_dat_mf <- rbind(final_dat_mf,dat)
  
}


final_dat_mf$x <- factor(final_dat_mf$x)

p3 <- ggplot(final_dat_mf, aes(x = x, y = dat)) + theme_dotplot + 
  ylab("Number of levels") +
  xlab("Studies") + ggtitle("MF") + geom_boxplot(aes(fill = x)) + coord_flip()
p3


#Distribution of GO-terms of CC hierarchy levels for each study
final_dat_cc <- data.frame()
name_cc <- names(cc_level)

for(k in 1:length(cc_level)){
  len <- unlist(cc_level[k])
  x <- rep(name_cc[k], length(len))
  dat <- as.data.frame(cc_level[[k]])
  colnames(dat) <- NULL
  #dat <- dat[-length(dat),]
  dat <- cbind(x, dat)
  final_dat_cc <- rbind(final_dat_cc,dat)
  
}


final_dat_cc$x <- factor(final_dat_cc$x)

p4 <- ggplot(final_dat_cc, aes(x = x, y = dat)) + theme_dotplot + 
  ylab("Number of levels") +
  xlab("Studies") + ggtitle("CC") + geom_boxplot(aes(fill = x)) + coord_flip()
p4

# Percentage of GO-terms of BP, MF and CC used by each study
library(GO.db)

# Get all GO IDs
all_go_ids <- keys(GO.db, keytype = "GOID")

# Separate by ontology
bp_ids <- all_go_ids[Ontology(all_go_ids) == "BP"]
mf_ids <- all_go_ids[Ontology(all_go_ids) == "MF"]
cc_ids <- all_go_ids[Ontology(all_go_ids) == "CC"]

length(bp_ids)
length(mf_ids)
length(cc_ids)

ind <- order(Biomarker_d$BP)
Biomarker_d$BP <- Biomarker_d$BP/length(bp_ids)
Biomarker_d$MF <- Biomarker_d$MF/length(mf_ids)
Biomarker_d$CC <- Biomarker_d$CC/length(cc_ids)

Biomarker_d$Study <- factor(Biomarker_d$Study, levels = Biomarker_d$Study[ind])

p5 <- ggplot(Biomarker_d, aes(x=Study)) + theme_dotplot +
  coord_flip() +
  ylab("Percentage of GO-terms") +
  xlab("Studies") +
  geom_point(aes(y = BP),color="green") + geom_point(aes(y = MF),color="red") +
  geom_point(aes(y = CC),color="blue")
p5

### Pairwise comparison of the biomarkers

n <- names(cancer.signatures.new)
dat <- matrix(data = NA, ncol = length(n), nrow = length(n))
rownames(dat) <- n
colnames(dat) <- 1:length(n)

sig_cancer <- list()

# Breast cancer signatures
for(i in 1:length(cancer.signatures.new)){
  sig_cancer[[i]] <- cancer.signatures.new[[i]]$sig
}

for (kk in 1:length(sig_cancer)) {
  for (jj in 1:length(sig_cancer)) {
    dat[kk,jj] <- (length(intersect(sig_cancer[[kk]],sig_cancer[[jj]]))/length(sig_cancer[[kk]])) 
  }
}

dat <- as.data.frame(dat)
da <- t(dat)
dat_melt <- melt(da)
Biomarker <- rep(unique(dat_melt$Var2),length(n))
dat_melt$Var1 <- Biomarker

# annotate text 
g.val <- c()
for (e in 1:length(n)) {
  g.val <- c(g.val,sum(dat[e,] != 0))
  
}
g.val <- g.val - 1
p.overlap <- ggplot(dat_melt, aes(Var1,Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "darkgreen") + theme(axis.text.x = element_text(angle = 90))
p.overlap + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  annotate("text", x = 49, y = 1:49, label = as.character(g.val), col=rgb(0,50,0,maxColorValue=255))


## GO-terms of BP overlap between the signatures
sig.go <- list()
dat.go <- matrix(data = NA, ncol = length(n), nrow = length(n))
rownames(dat.go) <- n
colnames(dat.go) <- 1:length(n)

for (i in 1:length(cancer.signatures.new)) {
  bp.go <- Gene2GOTermAndLevel(gene = cancer.signatures.new[[i]]$sig, organism = "Human", domain = "BP" )
  sig.go[[cancer.signatures.new[[i]]$name]] <- unique(bp.go$`GO ID`)
  print(i)
}

for (k in 1:length(sig.go)) {
  for (j in 1:length(sig.go)) {
    dat.go[k,j] <- (length(intersect(sig.go[[k]],sig.go[[j]]))/length(sig.go[[k]])) 
  }
}

dat.go <- as.data.frame(dat.go)
da.go <- t(dat.go)
dat_melt.go <- melt(da.go)
Biomarker.go <- rep(unique(dat_melt.go$Var2),length(n))
dat_melt.go$Var1 <- Biomarker.go

#annotate text
go.val <- c()
for (e in 1:length(n)) {
  go.val <- c(go.val,sum(dat.go[e,] != 0))
  
}
go.val <- go.val - 1
go.overlap <- ggplot(dat_melt.go, aes(Var1,Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "darkred") + theme(axis.text.x = element_text(angle = 90))
go.overlap + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  annotate("text", x = 49, y = 1:49, label = as.character(go.val), col=rgb(0,50,0,maxColorValue=255))



