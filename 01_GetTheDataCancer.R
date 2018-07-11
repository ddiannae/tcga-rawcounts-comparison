setwd("~/Workspace/raw-counts-comparison")

### Map new and old files for tumor
cancer.legacy <- read.table("fileid_caseid_filename_legacy_cancer.txt", sep=",", header=TRUE, stringsAsFactors = FALSE)
colnames(cancer.legacy) <- c("filename", "submitterid", "caseid", "fileid")
cancer.new <- read.table("fileid_caseid_filename_new_cancer.txt", sep=",", header=TRUE, stringsAsFactors = FALSE)
colnames(cancer.new) <- c("filename", "submitterid", "caseid", "fileid")

####################### Remove duplicated CaseIds
dup.caseids <- unique(cancer.new[which(duplicated(cancer.new$caseid)), "caseid"])
cancer.new.dup.caseid.all <- cancer.new[cancer.new$caseid %in% dup.caseids, ] 
cancer.new.dup.caseid.all <- cancer.new.dup.caseid.all[order(cancer.new.dup.caseid.all$caseid, cancer.new.dup.caseid.all$submitterid), ]
cancer.legacy.dup.caseid.all <- cancer.legacy[cancer.legacy$caseid %in% dup.caseids, ]
cancer.legacy.dup.caseid.all <- cancer.legacy.dup.caseid.all[order(cancer.legacy.dup.caseid.all$caseid, cancer.legacy.dup.caseid.all$submitterid), ]

############################### Should remove submitters ids that are not in the cancer.legacy from cancer.new
fileids.to.remove <- cancer.new.dup.caseid.all[which(cancer.new.dup.caseid.all$submitterid %in% setdiff(cancer.new.dup.caseid.all$submitterid, cancer.legacy.dup.caseid.all$submitterid)), "fileid"]
cancer.new <- cancer.new[!(cancer.new$fileid %in% fileids.to.remove),]
dim(cancer.new)
#### [1] 779   4

### We still have duplicates
dup.caseids <- unique(cancer.new[which(duplicated(cancer.new$caseid)), "caseid"])
cancer.new.dup.caseid.all <- cancer.new[cancer.new$caseid %in% dup.caseids, ] 
cancer.new.dup.caseid.all <- cancer.new.dup.caseid.all[order(cancer.new.dup.caseid.all$caseid, cancer.new.dup.caseid.all$submitterid), ]
cancer.legacy.dup.caseid.all <- cancer.legacy[cancer.legacy$caseid %in% dup.caseids, ]
cancer.legacy.dup.caseid.all <- cancer.legacy.dup.caseid.all[order(cancer.legacy.dup.caseid.all$caseid, cancer.legacy.dup.caseid.all$submitterid), ]
##############################################
#                                                 filename      submitterid                               caseid                               fileid
# 361 c22879b5-4183-4c3c-804b-f6cef1df617e.htseq.counts.gz TCGA-A7-A26E-01A 011b9b2d-ebe5-42bf-9662-d922faccc7a1 2932cdbb-71af-4404-8b13-7ea73e8d48c5
# 631 6aa69f7f-28c6-4943-b907-62d8085e9c6c.htseq.counts.gz TCGA-A7-A26E-01A 011b9b2d-ebe5-42bf-9662-d922faccc7a1 83acd71c-c12a-4394-ba8d-05e9c5ad2cf1
# 316 4e7c8fb8-81c8-477a-b012-b6e4d89c88f2.htseq.counts.gz TCGA-A7-A26J-01A 6cb6f179-defd-4661-af0a-c353b74c0c49 4ee13efd-d94c-4a2f-a8db-3f94ed939902
# 775 6577de41-48a3-4933-94de-64553f13bf05.htseq.counts.gz TCGA-A7-A26J-01A 6cb6f179-defd-4661-af0a-c353b74c0c49 89d21293-106d-4f52-a519-88c365fa1471
# 219 11b92965-3650-4f3e-a6da-65728487d7ab.htseq.counts.gz TCGA-A7-A13D-01A 8785012f-f73e-4d68-87cf-1d804af32782 a26a5dc8-1f9d-4f6a-917a-7457f8f8c056
# 724 01cdf894-fd1d-42c0-866b-10c68cde7c54.htseq.counts.gz TCGA-A7-A13D-01A 8785012f-f73e-4d68-87cf-1d804af32782 d382ff76-879d-4071-9e12-442cef2fd973
# 103 aa42dd76-efb1-4fb7-9922-cf54540844e7.htseq.counts.gz TCGA-A7-A13E-01A 8c7e74e0-71ef-49b8-9217-94b8ef740ef9 75ee2909-b949-45f8-b3cb-87806d223efa
# 778 96833a09-14ce-4483-b17d-f4328e0b68f5.htseq.counts.gz TCGA-A7-A13E-01A 8c7e74e0-71ef-49b8-9217-94b8ef740ef9 c804e8f4-bfe7-4ef6-b1ae-ca8cefa53aa7
# 66  7af2075c-0386-4971-ae25-375330ef6cec.htseq.counts.gz TCGA-A7-A0DB-01A f130f376-5801-40f9-975d-a7e2f7b5670d c66bfdc3-873c-413b-b6d6-495f3b43c0e9
# 620 0c7d119e-22c0-47ed-b736-3404ab2c7065.htseq.counts.gz TCGA-A7-A0DB-01A f130f376-5801-40f9-975d-a7e2f7b5670d b1d84ba0-bcf9-43be-8413-a20078a866ef
###############################################
### Get AliquotesId from TCGA
cancer.new.dup.aliquotes <- read.csv("fileid_aliquotesid.csv", quote="")
colnames(cancer.new.dup.aliquotes) <- c("aliq.submitterid", "sample.type", "filename", "submitterid", 
                                "fileid", "data.category", "cases.submitterid", "caseid", "id")
cancer.legacy.dup.caseid.all$aliq.submitterid <- substr(cancer.legacy.dup.caseid.all$filename, 14, 41)
############################### Should remove  aliquotes submitters ids that are not in the cancer.legacy from cancer.new
fileids.to.remove <- cancer.new.dup.aliquotes [!(cancer.new.dup.aliquotes$aliq.submitterid %in% cancer.legacy.dup.caseid.all$aliq.submitterid), "fileid"]
cancer.new <- cancer.new[!(cancer.new$fileid %in% fileids.to.remove),]
dim(cancer.new)
dup.caseids <- unique(cancer.new[which(duplicated(cancer.new$caseid)), "caseid"])
stopifnot(length(dup.caseids)==0)
#####################################################################
## Now we have more files in the legacy set that in the new set
dim(cancer.new)
## 774   4
dim(cancer.legacy)
## [1] 777   4
####################################################################
## Remove cases that are not shared
setdiff(cancer.new$caseid, cancer.legacy$caseid)
setdiff(cancer.legacy$caseid, cancer.new$caseid)
cancer.legacy <- cancer.legacy[!(cancer.legacy$caseid %in% setdiff(cancer.legacy$caseid, cancer.new$caseid)), ]
setdiff(cancer.new$caseid, cancer.legacy$caseid)
setdiff(cancer.legacy$caseid, cancer.new$caseid)
dim(cancer.new)
## 774   4
dim(cancer.legacy)
## [1] 774   4

#####################################################################
## Merge
cancer.all <- merge(cancer.legacy, cancer.new, by="caseid", all = TRUE, suffixes = c(".legacy",".new"))
dim(cancer.all)
##########################################################
## Are there different submitterids? YES
which(cancer.all$submitterid.legacy != cancer.all$submitterid.new)
cancer.all.diff.submitterids <- cancer.all[which(cancer.all$submitterid.legacy != cancer.all$submitterid.new), ]

###################################################
# caseid
# 199 3b7b9c1e-a84c-47ed-983c-9e4b00cbf01a      ## Old submitter (aliquote) got removed
# 331 67c5f371-3fa9-47c5-8b15-c2dd9acc8519      ## Old submitter (aliquote) got removed
# 346 6cd9baf5-bbe0-4c1e-a87f-c53b3af22890      ## This caseId has two counts.file but only one got added in the list fileid_caseid_filename_new_cancer.txt. ¿TCGA error?
                                                ## Should have to be downloadad manually or removed. Remove it
# 554 b2ecbc0f-2c30-4200-8d5e-7b95424bcadb      ## Old submitter (aliquote) got removed
# filename.legacy
# 199    UNCID_591641.TCGA-A7-A26F-01A-21R-A169-07.110919_UNC12-SN629_0141_AB0104ABXX.6.trimmed.annotated.gene.quantification.txt
# 331 UNCID_1111853.TCGA-AC-A2QH-01A-11R-A18M-07.111101_UNC13-SN749_0132_AD0ETKABXX.3_2.trimmed.annotated.gene.quantification.txt
# 346    UNCID_450071.TCGA-A7-A13G-01A-11R-A13Q-07.110805_UNC11-SN627_0134_AC042TABXX.3.trimmed.annotated.gene.quantification.txt
# 554    UNCID_591312.TCGA-A7-A26I-01A-11R-A169-07.110919_UNC12-SN629_0142_BB00V1ABXX.1.trimmed.annotated.gene.quantification.txt
# submitterid.legacy                        fileid.legacy                                         filename.new  submitterid.new
# 199   TCGA-A7-A26F-01A dac427ae-ee20-4253-9948-02bb10f8d96d 28dfbaea-e781-4f64-ae7e-d12d49f33675.htseq.counts.gz TCGA-A7-A26F-01B
# 331   TCGA-AC-A2QH-01A bfb8bcde-8681-4c16-9d36-6fef577847d7 2656413d-a5d2-4812-9ce9-a02c15ab04bd.htseq.counts.gz TCGA-AC-A2QH-01B
# 346   TCGA-A7-A13G-01A aa72e2d1-fcf2-4c2d-962a-e332f52f03d4 eaee0d82-cd7e-47a5-859f-59b7bf5ad6a0.htseq.counts.gz TCGA-A7-A13G-01B
# 554   TCGA-A7-A26I-01A e2f4857d-d772-4209-b46a-541b5fb3f671 afe5f01a-b3c1-4d24-94cf-27ec0894e308.htseq.counts.gz TCGA-A7-A26I-01B
# fileid.new
# 199 c26de1c6-96c2-45da-8e2a-9868f21a9757
# 331 238b412c-dc67-4bed-8092-95052029d46f
# 346 26562f3a-6491-477c-888b-8cd9e994233a
# 554 e72ac810-c468-4e78-a10b-3e402996bb25
cancer.all <- cancer.all[cancer.all$caseid != "6cd9baf5-bbe0-4c1e-a87f-c53b3af22890", ]
dim(cancer.all)
cancer.all$filename.new <- sub(".gz", "", cancer.all$filename.new)
cancer.all <- subset(cancer.all, select=-c(filename.new))
colnames(cancer.all) <- c("submitterid", "filename.legacy", "caseid", "fileid.legacy",  "fileid.new", "filename.new")
cancer.all <- cancer.all[order(cancer.all$filename.new),]
save(cancer.all, file="cancerFiles.RData", compress="xz")

############################################################
library("BiocParallel")
load(file="cancerFiles.RData")
tumorFiles <- cancer.all$filename.new 
tumorFiles <- sort(tumorFiles)
is.unsorted(tumorFiles)
tumor <- bplapply(paste("new/cancer", tumorFiles, sep="/"), read.delim, sep="\t", header=F, col.names=c("EnsemblID", "raw_counts"))
length(intersect(tumorFiles, cancer.all$filename.new))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(tumor, dim)))
stopifnot(nrow(size)==1)
cat('Tumor samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(tumor, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in tumor samples match positions \n')

##Lets keep only the raw counts
tumor<-do.call(cbind, lapply(tumor, function(x)x[,"raw_counts"]))

##Order cancer.all by newfilenames, check if files are in same order and add Case Ids
cancer.all <- cancer.all[order(cancer.all$filename.new), ]
tumorFileNames <- cbind(tumorFiles, cancer.all$filename.new)
tumorFileNames <- t(unique(t(tumorFileNames)))
stopifnot(dim(tumorFileNames)==c(length(tumorFiles), 1))

targets<-data.frame(File=tumorFiles, ID=paste("T", 1:length(tumorFiles), sep=""), Case=cancer.all$caseid)
colnames(tumor)<-targets$ID

##Let's change the annotation
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))
colnames(genes)<-c("EnsemblId", "version")

tumor<-list(Counts=tumor, Annot=genes[,1], targets=targets)
save(tumor, file="cancerRawNew.RData", compress="xz")
cat('cancerRawNew.RData saved \n')
######################################################################3
library("BiocParallel")
load(file="cancerFiles.RData")
tumorFiles <- cancer.all$filename.legacy 
tumorFiles <- sort(tumorFiles)
is.unsorted(tumorFiles)
tumor <- bplapply(paste("legacy/cancer", tumorFiles, sep="/"), read.delim, sep="\t", header=T)
length(intersect(tumorFiles, cancer.all$filename.legacy))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(tumor, dim)))
stopifnot(nrow(size)==1)
cat('Tumor samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(tumor, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in tumor samples match positions \n')
dim(genes)
##Lets keep only the raw counts
tumor<-do.call(cbind, lapply(tumor, function(x)x[,"raw_counts"]))

##Order cancer.all by newfilenames, check if files are in same order and add Case Ids
cancer.all <- cancer.all[order(cancer.all$filename.legacy), ]
tumorFileNames <- cbind(tumorFiles, cancer.all$filename.legacy)
tumorFileNames <- t(unique(t(tumorFileNames)))
stopifnot(dim(tumorFileNames)==c(length(tumorFiles), 1))
targets<-data.frame(File=tumorFiles, ID=paste("T", 1:length(tumorFiles), sep=""), Case=cancer.all$caseid)
colnames(tumor)<-targets$ID

##Let's change the annotation
genes<-do.call(rbind,sapply(genes[,1], strsplit, split="|", fixed=TRUE))
colnames(genes)<-c("Symbol", "EntrezId")
tumor<-list(Counts=tumor, Annot=genes, targets=targets)
save(tumor, file="cancerRawLegacy.RData", compress="xz")
cat('cancerRawNew.RData saved \n')
