setwd("~/Workspace/raw-counts-comparison")

### Map new and old files for normal
healthy.legacy <- read.table("fileid_caseid_filename_legacy_healthy.txt", sep=",", header=TRUE, stringsAsFactors = FALSE)
colnames(healthy.legacy) <- c("filename", "submitterid", "caseid", "fileid")
healthy.new <- read.table("fileid_caseid_filename_new_healthy.txt", sep=",", header=TRUE, stringsAsFactors = FALSE)
colnames(healthy.new) <- c("filename", "submitterid", "caseid", "fileid")
### Merge
healthy.all <- merge(healthy.legacy, healthy.new, by="caseid",  all = TRUE, suffixes = c(".legacy",".new"))
## check if there are cases with different submitterid
which(healthy.all$submitterid.legacy != healthy.all$submitterid.new)
healthy.all$filename.new <- sub(".gz", "", healthy.all$filename.new)

dim(healthy.all)
## [1] 101   6
dim(healthy.legacy)
## [1] 101   4
dim(healthy.new)
## [1] 101   4
save(healthy.all, file="healthyFiles.RData", compress="xz")
###############################################################################
### Load New healthy Matrix
library("BiocParallel")
load(file="healthyFiles.RData")
normalFiles <- healthy.all$filename.new 
normalFiles <- sort(normalFiles)
is.unsorted(normalFiles)
normal<-bplapply(paste("new/healthy", normalFiles, sep="/"), read.delim, sep="\t", header=F, col.names=c("EnsemblID", "raw_counts"))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(normal, dim)))
stopifnot(nrow(size)==1)
cat('Normal samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(normal, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in normal samples match positions \n')

##Let's keep only the raw counts
normal<-do.call(cbind, lapply(normal, function(x)x[,"raw_counts"]))

##Order healthy.all by filename, check if files are in same order and add Case Ids
healthy.all <- healthy.all[order(healthy.all$filename.new), ]
normalFileNames <- cbind(normalFiles, healthy.all$filename.new)
normalFileNames <- t(unique(t(normalFileNames)))
stopifnot(dim(normalFileNames)==c(length(normalFiles), 1))

targets<-data.frame(File=normalFiles, Case=healthy.all$caseid, ID=paste("N", 1:length(normalFiles), sep=""))
colnames(normal) <- targets$ID

##Let's change the annotation. Remove Ensembl version
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))

##Save clean data
normal<-list(Counts=normal, Annot=genes[,1], targets=targets)
save(normal, file="healthyRawNew.RData", compress="xz")
cat('healthyRawNew.RData saved \n')
#########################################################################
library("BiocParallel")
load(file="healthyFiles.RData")
normalFiles <- healthy.all$filename.legacy 
normalFiles <- sort(normalFiles)
is.unsorted(normalFiles)
normal<-bplapply(paste("legacy/healthy", normalFiles, sep="/"), read.delim, sep="\t", header=T)

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(normal, dim)))
stopifnot(nrow(size)==1)
cat('Normal samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(normal, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in normal samples match positions \n')

##Let's keep only the raw counts
normal<-do.call(cbind, lapply(normal, function(x)x[,"raw_counts"]))

##Order healthy.all by filename, check if files are in same order and add Case Ids
healthy.all <- healthy.all[order(healthy.all$filename.legacy), ]
normalFileNames <- cbind(normalFiles, healthy.all$filename.legacy)
normalFileNames <- t(unique(t(normalFileNames)))
stopifnot(dim(normalFileNames)==c(length(normalFiles), 1))

targets<-data.frame(File=normalFiles, Case=healthy.all$caseid, ID=paste("N", 1:length(normalFiles), sep=""))
colnames(normal) <- targets$ID

##Let's change the annotation. Remove Ensembl version
genes<-do.call(rbind,sapply(genes[,1], strsplit, split="|", fixed=TRUE))
colnames(genes)<-c("Symbol", "EntrezId")
##Save clean data
normal<-list(Counts=normal, Annot=genes, targets=targets)
save(normal, file="healthyRawLegacy.RData", compress="xz")
cat('healthyRawLegacy.RData saved \n')

##########################################################################

