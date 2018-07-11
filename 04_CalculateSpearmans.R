
setwd("~/Workspace/rnapipeline/tcga-rawcounts-comparison")
load(file = "cancerRawCleanLegacy.RData")
load(file = "healthyRawCleanLegacy.RData")
load(file = "cancerRawClean.RData")
load(file = "healthyRawClean.RData")

### Define Spearmans Function
spearmanCorr <- function(newGenes, legacyGenes){
  dim(newGenes)
  dim(legacyGenes)
  genes.ranked <- data.frame(cbind(rank(newGenes, ties.method = 'average'),
                                   rank(legacyGenes, ties.method = 'average')))
  colnames(genes.ranked) <- c('new', 'legacy')
  rho <- cov(genes.ranked) / (sd(genes.ranked$new) * sd(genes.ranked$legacy))
  n <- length(genes.ranked$new)
  r <- cor(x = genes.ranked$new, y = genes.ranked$legacy, method = 'pearson')
  s <- (n^3 - n) * (1 - r) / 6
  t <- r * sqrt((n - 2) / (1 - r^2))
  p <- 2 * (1-pt(q = t, df = n - 2))
  return(c(rho[[2]], s, p))
}

#### First Normal Data
normal.targets <- normal$targets
normal.targets.legacy <- normal.legacy$targets

head(normal.targets)
head(normal.targets.legacy)

normal.targets <- normal.targets[order(normal.targets$Case), ]
normal.targets.legacy <- normal.targets.legacy[order(normal.targets.legacy$Case), ]
## Do we have them all?
size <- unique(rbind(dim(normal.targets), dim(normal.targets.legacy))[1])
cases <- cbind(as.character(normal.targets$Case), as.character(normal.targets.legacy$Case))
cases <- t(unique(t(cases)))
stopifnot(dim(cases) == c(size[1], 1))
normal.df <- data.frame(normal$Counts)
normal.df.legacy <- data.frame(normal.legacy$Counts)

head(normal.df)
head(normal.df.legacy)

results.df <- data.frame()
for (i in 1:nrow(normal.targets)) {
  results.df <- rbind(results.df, 
                  spearmanCorr(normal.df[, normal.targets[i, "ID"]], 
                  normal.df.legacy[, normal.targets.legacy[i, "ID"]]))
}
colnames(results.df) <- c("rho", "S", "p")
library(ggplot2)
# Basic histogram
df <- data.frame(normal.targets$Case, results.df$rho)
colnames(df) <- c("Case", "Rho")
df <- df[order(df$Rho), ]
df
ggplot(data=df, aes(x=Case, y=Rho, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point() 
mean(results.df$rho)
write.table(df,  file="normal_spearman.txt", quote=F, sep="\t")



## Now Tumor Data
tumor.targets <- tumor$targets
tumor.targets.legacy <- tumor.legacy$targets

head(tumor.targets)
head(tumor.targets.legacy)

tumor.targets <- tumor.targets[order(tumor.targets$Case), ]
tumor.targets.legacy <- tumor.targets.legacy[order(tumor.targets.legacy$Case), ]
## Do we have them all?
size <- unique(rbind(dim(tumor.targets), dim(tumor.targets.legacy))[1])
cases <- cbind(as.character(tumor.targets$Case), as.character(tumor.targets.legacy$Case))
cases <- t(unique(t(cases)))
stopifnot(dim(cases) == c(size[1], 1))
tumor.df <- data.frame(tumor$Counts)
tumor.df.legacy <- data.frame(tumor.legacy$Counts)

head(tumor.df)
head(tumor.df.legacy)

results.tumor.df <- data.frame()
for (i in 1:nrow(tumor.targets)) {
  results.tumor.df <- rbind(results.tumor.df, 
                      spearmanCorr(tumor.df[, tumor.targets[i, "ID"]], 
                                   tumor.df.legacy[, tumor.targets.legacy[i, "ID"]]))
}
colnames(results.tumor.df) <- c("rho", "S", "p")
df <- data.frame(tumor.targets$Case, results.tumor.df$rho)
colnames(df) <- c("Case", "Rho")
df <- df[order(df$Rho), ]
df
ggplot(data=df, aes(x=Case, y=Rho, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point() 
mean(results.tumor.df$rho)

write.table(df, file="cancer_spearman.txt", quote=F, sep="\t")
## Listo, con 18547 genes