
setwd("~/Workspace/rnapipeline/tcga-rawcounts-comparison")
load(file = "cancerRawCleanLegacy.RData")
load(file = "healthyRawCleanLegacy.RData")
load(file = "annotEntrez.RData")
normal.legacy <- normal
tumor.legacy <- tumor
load(file = "cancerRawNew.RData")
load(file = "healthyRawNew.RData")

normal.df <- data.frame(normal$Counts)
tumor.df <- data.frame(tumor$Counts)

## Check if genes are in the same position
size <- unique(rbind(dim(normal$Counts), dim(tumor$Counts))[1])
genes <- cbind(tumor$Annot, normal$Annot)
genes <- t(unique(t(genes)))
stopifnot(dim(genes) == c(size[1], 1))

rownames(normal.df) <- normal$Annot
rownames(tumor.df) <- tumor$Annot
head(normal.df)
head(tumor.df)

dim(normal.df)
## [1] 60488   101
dim(tumor.df)
## [1] 60488   773

all.ensembl <- unique(merged.annot$EnsemblId)
length(all.ensembl)
## [1] 18563
out.ensembl.ids <- which(!(rownames(normal.df) %in% all.ensembl))
length(out.ensembl.ids)
## [1] 41941
normal.df <- normal.df[-out.ensembl.ids, ]
dim(normal.df)
##[1] 18547   101
tumor.df <- tumor.df[-out.ensembl.ids, ]
dim(tumor.df)
##[1] 18547   773

merged.annot <- merged.annot[merged.annot$EnsemblId %in% rownames(normal.df), ]
dim(merged.annot)
##[1] 18547     8
normal.df.legacy <- data.frame(normal.legacy$Counts)
rownames(normal.df.legacy) <- normal.legacy$Annot
tumor.df.legacy <- data.frame(tumor.legacy$Counts)
rownames(tumor.df.legacy) <- tumor.legacy$Annot

dim(normal.df.legacy)
## [1] 18563   101
dim(tumor.df.legacy)
## [1] 18563   773

normal.df.legacy <- normal.df.legacy[rownames(normal.df.legacy) %in% merged.annot$EnsemblId, ]
tumor.df.legacy <- tumor.df.legacy[rownames(tumor.df.legacy) %in% merged.annot$EnsemblId, ]

dim(normal.df.legacy)
## [1] 18547   101
dim(tumor.df.legacy)
## [1] 18547   101

## Order rows by EnsemblId

normal.df <- normal.df[order(rownames(normal.df)), ]
tumor.df <- tumor.df[order(rownames(tumor.df)), ]
normal.df.legacy <- normal.df.legacy[order(rownames(normal.df.legacy)), ]
tumor.df.legacy <- tumor.df.legacy[order(rownames(tumor.df.legacy)), ]
merged.annot <- merged.annot[order(merged.annot$EnsemblId), ]

head(normal.df)
head(tumor.df)
head(normal.df.legacy)
head(tumor.df.legacy)
head(merged.annot)

##save all
normal$Counts <- data.matrix(normal.df, rownames.force = F)
tumor$Counts <- data.matrix(tumor.df, rownames.force = F)
normal.legacy$Counts <- data.matrix(normal.df.legacy, rownames.force = F)
tumor.legacy$Counts <- data.matrix(tumor.df.legacy, rownames.force = F)

normal$Annot <- as.character(merged.annot$EnsemblId)
tumor$Annot <-  as.character(merged.annot$EnsemblId)
normal.legacy$Annot <- as.character(merged.annot$EnsemblId)
tumor.legacy$Annot <-  as.character(merged.annot$EnsemblId)

save(merged.annot, file="annotEntrez.RData", compress="xz")
save(normal, file="healthyRawClean.RData", compress="xz")
save(tumor, file="cancerRawClean.RData", compress="xz")
save(normal.legacy, file="healthyRawCleanLegacy.RData", compress="xz")
save(tumor.legacy, file="cancerRawCleanLegacy.RData", compress="xz")
