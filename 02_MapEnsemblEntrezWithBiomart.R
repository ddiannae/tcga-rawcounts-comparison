
setwd("~/Workspace/raw-counts-comparison")
load(file="healthyRawLegacy.RData")
load(file="cancerRawLegacy.RData")

annot <- read.delim(file="Biomart_EnsemblG92_GRCh38_p12_NCBI.txt", sep="\t")
names(annot) <- c("EnsemblId", "Chr", "Start", "End", "GC", "Type", "Symbol", "EntrezId")
dim(annot)
## [1] 64927     8
levels(annot$Chr)
annot<-annot[annot$Chr%in%c(as.character(1:22), "X", "Y"),]
annot$Chr<-droplevels(annot$Chr)
table(annot$Chr)
dim(annot)
## [1] 58519     8

### Genes with no Symbol
annot.nosymbol <- annot[annot$Symbol=="", ]
annot.nosymbol
dim(annot.nosymbol)
##[1] 20925     8
annot.entrezna <- annot[is.na(annot$EntrezId),]
annot.entrezna
dim(annot.entrezna)
## [1] 39054     8

###### Healthy Data Frame
normal.df <- data.frame(normal$Counts)
tumor.df <- data.frame(tumor$Counts)
head(normal.df)
dim(normal.df)
## [1] 20532   101
rownames(normal.df) <- as.numeric(normal$Annot[,2])
head(normal.df)
##[1] 20532     2
head(normal$Annot[,2])
setdiff(normal$Annot[,2], tumor$Annot[,2])
length(normal$Annot)
length(tumor$Annot)
all.entrez <- unique(annot$EntrezId)
length(all.entrez)
#[1] 19338

## Should remove entrez ids that are not in all.entrez
out.entrez.ids <- which(!(rownames(normal.df) %in% all.entrez))
normal.df <- normal.df[-out.entrez.ids, ]
dim(normal.df)
#[1] 18578   101
normal$Annot <- normal$Annot[-out.entrez.ids,]
dim(normal$Annot)
#[1] 18578     2

################ Merge and remove duplicates
annot.normal <- data.frame(normal$Annot, stringsAsFactors = FALSE)
##Add Biomart data
merged.annot<-merge(x=annot.normal, y=annot, by="EntrezId", all=F, sort=FALSE)
dim(merged.annot)
#[1] 18639     9
merged.annot[duplicated(merged.annot[,c("EntrezId", "EnsemblId")]), ]
duplicated.entrez.ensembl <- merged.annot[duplicated(merged.annot[,c("EntrezId", "EnsemblId")]), c("EntrezId", "EnsemblId")]
duplicated.entrez.ensembl.rows <- merged.annot[merged.annot$EntrezId %in% duplicated.entrez.ensembl$EntrezId |
                                                 merged.annot$EnsemblId %in% duplicated.entrez.ensembl$EntrezId, ]

duplicated.entrez.ensembl.rows 
# EntrezId Symbol.x       EnsemblId Chr    Start      End   GC           Type Symbol.y
# 2977     6349   CCL3L1 ENSG00000276085  17 36194869 36196758 52.7 protein_coding   CCL3L1 https://www.ncbi.nlm.nih.gov/gene/6349 No Ensembl Id
# 2978     6349   CCL3L1 ENSG00000276085  17 36194869 36196758 52.7 protein_coding   CCL3L3 https://www.ncbi.nlm.nih.gov/gene/6349 No Ensembl Id
# 2979   414062   CCL3L3 ENSG00000276085  17 36194869 36196758 52.7 protein_coding   CCL3L3 OK
# 2980   414062   CCL3L3 ENSG00000276085  17 36194869 36196758 52.7 protein_coding   CCL3L1 https://www.ncbi.nlm.nih.gov/gene/414062 Symbol: CCL3L3

ids.to.remove = c("2977", "2978", "2980")
merged.annot[ids.to.remove, ]
merged.annot <- merged.annot[!(rownames(merged.annot) %in% ids.to.remove), ]
dim(merged.annot)
## [1] 18636     9

##### Duplicates. EnsemblId
dup.ensembl.ids <- which(duplicated(merged.annot$EnsemblId))
dup.ensembl <-unique(merged.annot[dup.ensembl.ids, "EnsemblId"])
dup.ensembl.rows <- merged.annot[which(merged.annot$EnsemblId %in% dup.ensembl), ]
dim(dup.ensembl.rows)
###[1] 75  9
## We need to remove the duplicated EnsemblIds. Lets keep them but only one Entrez counts should count
## probaby take max
length(unique(dup.ensembl.rows$EntrezId))
## [1] 70
dup.entrez.in.dup.ensembl <- dup.ensembl.rows[which(duplicated(dup.ensembl.rows$EntrezId)), c("EntrezId", "EnsemblId")]
dup.entrez.in.dup.ensembl <- merged.annot[merged.annot$EntrezId %in% dup.entrez.in.dup.ensembl$EntrezId 
                                          | merged.annot$EnsemblId %in% dup.entrez.in.dup.ensembl$EnsemblId , ]

dup.entrez.in.dup.ensembl <- dup.entrez.in.dup.ensembl[order(dup.entrez.in.dup.ensembl$EnsemblId, dup.entrez.in.dup.ensembl$EntrezId), ]
dup.entrez.in.dup.ensembl
#       EntrezId Symbol.x       EnsemblId    Start      End Chr    GC           Type Symbol.y
# 2786    220869    CBWD5 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5
# 2784    445571    CBWD3 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5
# 2788    644019    CBWD6 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5
# 6720 100132565  GOLGA8F ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F
# 6722    283768  GOLGA8G ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F
# 8946      3963   LGALS7 ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B
# 8944    653499  LGALS7B ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B
# 6719 100132565  GOLGA8F ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G
# 6721    283768  GOLGA8G ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G
# 8947      3963   LGALS7 ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7
# 8945    653499  LGALS7B ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7
# 2787    644019    CBWD6 ENSG00000215126 41131306 41199261   9 36.14 protein_coding    CBWD6

### Duplicates. Entrez.Ids
dup.entrez.ids <-which(duplicated(merged.annot$EntrezId))
dup.entrez <- unique(merged.annot[dup.entrez.ids, "EntrezId"])
dup.entrez.rows <- merged.annot[which(merged.annot$EntrezId %in% dup.entrez), ]
dim(dup.entrez.rows)
###[1] 117   9
length(unique(dup.ensembl.rows$EnsemblId))
## [1] 35
dup.ensembl.in.dup.entrez <- dup.entrez.rows[which(duplicated(dup.entrez.rows$EnsemblId)), c("EntrezId", "EnsemblId")]
dup.ensembl.in.dup.entrez <- merged.annot[merged.annot$EntrezId %in% dup.ensembl.in.dup.entrez$EntrezId 
                                          | merged.annot$EnsemblId %in% dup.ensembl.in.dup.entrez$EnsemblId , ]
dup.ensembl.in.dup.entrez <- dup.ensembl.in.dup.entrez[order(dup.ensembl.in.dup.entrez$EnsemblId, dup.ensembl.in.dup.entrez$EntrezId), ]

identical(dup.ensembl.in.dup.entrez, dup.entrez.in.dup.ensembl)
######################################
## dup.esenbl.in.dup.entrez has 150472 - ENSG00000215126 extra. 
## check one by one
dup.ensembl.in.dup.entrez
#       EntrezId Symbol.x       EnsemblId    Start      End Chr    GC           Type Symbol.y
# 2786    220869    CBWD5 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5 https://www.ncbi.nlm.nih.gov/gene/220869 OK
# 2784    445571    CBWD3 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5 https://www.ncbi.nlm.nih.gov/gene/445571 CBWD3 ENSG00000196873 
# 2788    644019    CBWD6 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5 https://www.ncbi.nlm.nih.gov/gene/644019 CBWD6 ENSG00000215126
# 6720 100132565  GOLGA8F ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F https://www.ncbi.nlm.nih.gov/gene/100132565 OK
# 6722    283768  GOLGA8G ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F https://www.ncbi.nlm.nih.gov/gene/283768 GOLGA8G ENSG00000183629
# 8946      3963   LGALS7 ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B https://www.ncbi.nlm.nih.gov/gene/3963 LGALS7 ENSG00000205076
# 8944    653499  LGALS7B ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B https://www.ncbi.nlm.nih.gov/gene/653499 OK
# 6719 100132565  GOLGA8F ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G https://www.ncbi.nlm.nih.gov/gene/100132565 GOLGA8F ENSG00000153684
# 6721    283768  GOLGA8G ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G https://www.ncbi.nlm.nih.gov/gene/283768 OK
# 8947      3963   LGALS7 ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7 https://www.ncbi.nlm.nih.gov/gene/3963 OK
# 8945    653499  LGALS7B ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7 https://www.ncbi.nlm.nih.gov/gene/653499 LGALS7B ENSG00000178934
# 2783    150472    CBWD2 ENSG00000215126 41131306 41199261   9 36.14 protein_coding    CBWD6 https://www.ncbi.nlm.nih.gov/gene/150472 CBWD2 ENSG00000136682
# 2787    644019    CBWD6 ENSG00000215126 41131306 41199261   9 36.14 protein_coding    CBWD6 https://www.ncbi.nlm.nih.gov/gene/644019 OK

setdiff(rownames(dup.ensembl.in.dup.entrez), rownames(dup.entrez.in.dup.ensembl))
## [1] "2783"

merged.annot[merged.annot$EnsemblId == "ENSG00000196873", ]
merged.annot["2784", ] <- merged.annot["2785", ]
merged.annot[merged.annot$EnsemblId == "ENSG00000196873", ]
merged.annot[merged.annot$EntrezId == merged.annot["2784", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000215126", ]
merged.annot["2788", ] <- merged.annot["2787", ]
merged.annot[merged.annot$EnsemblId == "ENSG00000215126", ]
merged.annot[merged.annot$EntrezId == merged.annot["2788", "EntrezId"] , ]

merged.annot[merged.annot$EnsemblId == "ENSG00000183629", ]
merged.annot["6722", ] <- merged.annot["6721", ]
merged.annot[merged.annot$EnsemblId == "ENSG00000183629", ]
merged.annot[merged.annot$EntrezId == merged.annot["6722", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000205076", ]
merged.annot["8946", ] <- merged.annot["8947", ]
merged.annot[merged.annot$EnsemblId == "ENSG00000205076", ]
merged.annot[merged.annot$EntrezId == merged.annot["8946", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000153684", ]
merged.annot["6719", ] <- merged.annot["6720", ]
merged.annot[merged.annot$EntrezId == merged.annot["6719", "EntrezId"], ]
merged.annot[merged.annot$EnsemblId == "ENSG00000153684", ]

merged.annot[merged.annot$EnsemblId == "ENSG00000178934", ]
merged.annot["8945", ] <- merged.annot["8944", ]
merged.annot[merged.annot$EnsemblId == "ENSG00000178934", ]
merged.annot[merged.annot$EntrezId == merged.annot["8945", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000136682", ]
merged.annot["2783", ] <- merged.annot["2782", ]
merged.annot[merged.annot$EnsemblId == "ENSG00000136682", ]
merged.annot[merged.annot$EntrezId == merged.annot["2783", "EntrezId"], ]

### We can remove Entrezid - Ensemblid duplicates 
merged.annot[duplicated(merged.annot), ]
merged.annot <- merged.annot[!duplicated(merged.annot), ]
dim(merged.annot)
## [1] 18629     9

##### Duplicates. EnsemblId - Take2
dup.ensembl.ids <- which(duplicated(merged.annot$EnsemblId))
dup.ensembl <-unique(merged.annot[dup.ensembl.ids, "EnsemblId"])
dup.ensembl.rows <- merged.annot[which(merged.annot$EnsemblId %in% dup.ensembl), ]
dim(dup.ensembl.rows)
###[1] 62  9
length(unique(dup.ensembl.rows$EntrezId))
##[1] 62 YEI!!!
dup.ensembl.rows$EntrezId <- as.numeric(dup.ensembl.rows$EntrezId)
dup.ensembl.rows <- dup.ensembl.rows[order(dup.ensembl.rows$EntrezId), ]
dim(dup.ensembl.rows)
## [1] 62  9
dup.ensembl.rows
write.table(dup.ensembl.rows, file = "dup_ensembl_rows.tsv", quote=F, sep="\t")
### Check one by one
### Remove the ones with no EnsemblId
ids.to.remove = c("6338", "6341", "6340", "6073", "6337", "6339", "6343")
merged.annot[ids.to.remove, ]
merged.annot <- merged.annot[!(rownames(merged.annot) %in% ids.to.remove), ]
dim(merged.annot)
##[1] 18622     9

merged.annot[merged.annot$EnsemblId == "ENSG00000160014" & merged.annot$EntrezId == "808", ]
merged.annot["2652",] = merged.annot["2653",]
merged.annot["2652",]
merged.annot[merged.annot$EnsemblId == "ENSG00000160014", ]
merged.annot[merged.annot$EntrezId == merged.annot["2652", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000196735" & merged.annot$EntrezId == "3117", ]
merged.annot["7374",] = merged.annot["7373",]
merged.annot["7374",]
merged.annot[merged.annot$EnsemblId == "ENSG00000196735", ]
merged.annot[merged.annot$EntrezId == merged.annot["7374", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000231887" & merged.annot$EntrezId == "5554", ]
merged.annot["12754",] = merged.annot["12755", ]
merged.annot["12755", ]
merged.annot[merged.annot$EnsemblId == "ENSG00000231887", ]
merged.annot[merged.annot$EntrezId == merged.annot["12755", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000177954" & merged.annot$EntrezId == "6232", ]
merged.annot["13942",] = merged.annot["13941",]
merged.annot[merged.annot$EnsemblId == "ENSG00000177954", ]
merged.annot[merged.annot$EntrezId == merged.annot["13941", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000172062" & merged.annot$EntrezId == "6606", ]
merged.annot["15115",] = merged.annot["15114",]
merged.annot[merged.annot$EnsemblId == "ENSG00000172062", ]
merged.annot[merged.annot$EntrezId == merged.annot["15115", "EntrezId"] , ]

merged.annot[merged.annot$EnsemblId == "ENSG00000244682" & merged.annot$EntrezId == "9103", ]
merged.annot["5920", "Start"] <- 161581339
merged.annot["5920", "End"] <- 	161605662
merged.annot["5920", "GC"] <- 44.62
merged.annot["5920", "Type"] <- "polymorphic_pseudogene"
merged.annot["5920", "EnsemblId"] <- "ENSG00000244682"
merged.annot["5920", "Symbol.y"] <- "FCGR2C"
merged.annot[merged.annot$EnsemblId == "ENSG00000244682", ]
merged.annot[merged.annot$EntrezId == "9103", ]

merged.annot[merged.annot$EnsemblId == "ENSG00000256060" & merged.annot$EntrezId == "10597", ]
merged.annot["16795",] = merged.annot["16796",]
merged.annot[merged.annot$EnsemblId == "ENSG00000256060", ]
merged.annot[merged.annot$EntrezId == merged.annot["16795", "EntrezId"] , ]

merged.annot[merged.annot$EnsemblId == "ENSG00000169618" & merged.annot$EntrezId == "10887", ]
merged.annot["12825",]
merged.annot["12825", "EnsemblId"] <- "ENSG00000169618"
merged.annot["12825", "Symbol.y"] <- "PROKR1"
merged.annot["12825", "Start"] <- 68643589
merged.annot["12825", "End"] <- 	68658247
merged.annot["12825", "GC"] <- 46.47
merged.annot[merged.annot$EnsemblId == "ENSG00000169618", ]
merged.annot[merged.annot$EntrezId == "10887", ]

merged.annot[merged.annot$EnsemblId == "ENSG00000144134" & merged.annot$EntrezId == "11159", ]
merged.annot["13253",] <- merged.annot["13254",]
merged.annot[merged.annot$EnsemblId == "ENSG00000144134", ]
merged.annot[merged.annot$EntrezId == merged.annot["13253", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000188092" & merged.annot$EntrezId == "51463", ]
merged.annot["6879",] <- merged.annot["6880",]
merged.annot[merged.annot$EnsemblId == "ENSG00000188092", ]
merged.annot[merged.annot$EntrezId == merged.annot["6879", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000187272" & merged.annot$EntrezId == "83901", ]
merged.annot["8772",] = merged.annot["8771",]
merged.annot[merged.annot$EntrezId == merged.annot["8772", "EntrezId"] , ]
merged.annot[merged.annot$EnsemblId == "ENSG00000187272", ]

merged.annot[merged.annot$EnsemblId == "ENSG00000182986" & merged.annot$EntrezId == "162967", ]
merged.annot["18208",] <- merged.annot["18207",]
merged.annot[merged.annot$EnsemblId == "ENSG00000182986", ]
merged.annot[merged.annot$EntrezId == merged.annot["18208", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000158427" & merged.annot$EntrezId == "286527", ]
merged.annot["16591",] <- merged.annot["16592",]
merged.annot[merged.annot$EntrezId == merged.annot["16591", "EntrezId"], ]
merged.annot["16593",] <- merged.annot["16592",]
merged.annot[merged.annot$EnsemblId == "ENSG00000158427", ]

merged.annot[merged.annot$EnsemblId == "ENSG00000238269" & merged.annot$EntrezId == "389860", ]
merged.annot["11643",] <- merged.annot["11644",]
merged.annot[merged.annot$EnsemblId == "ENSG00000238269", ]
merged.annot[merged.annot$EntrezId == merged.annot["11643", "EntrezId"] , ]

merged.annot[merged.annot$EnsemblId == "ENSG00000278299" & merged.annot$EntrezId == "414060", ]
merged.annot["15942",] <- merged.annot["15943",]
merged.annot[merged.annot$EnsemblId == "ENSG00000278299", ]
merged.annot[merged.annot$EntrezId == merged.annot["15942", "EntrezId"] , ]

merged.annot[merged.annot$EnsemblId == "ENSG00000268940" & merged.annot$EntrezId == "541466", ]
merged.annot["4078",] <- merged.annot["4077",]
merged.annot[merged.annot$EnsemblId == "ENSG00000268940", ]
merged.annot[merged.annot$EntrezId == merged.annot["4078", "EntrezId"], ]

merged.annot[merged.annot$EnsemblId == "ENSG00000183474" & merged.annot$EntrezId == "728340", ]
merged.annot["7026",] <- merged.annot["7025",]
merged.annot[merged.annot$EntrezId == merged.annot["7026", "EntrezId"] , ]
merged.annot[merged.annot$EnsemblId == "ENSG00000183474", ]

merged.annot[merged.annot$EnsemblId == "ENSG00000099974" & merged.annot$EntrezId == "100037417", ]
merged.annot["4448",] <- merged.annot["4449",]
merged.annot[merged.annot$EnsemblId == "ENSG00000099974", ]
merged.annot[merged.annot$EntrezId == merged.annot["4448", "EntrezId"], ]

all(lapply(as.character(merged.annot$EnsemblId), nchar) == 15)
#[1] TRUE

### Duplicates Ensambl Take 3
merged.annot[duplicated(merged.annot), ]
merged.annot <- merged.annot[!duplicated(merged.annot), ]
dim(merged.annot)
## [1] 18605     9

dup.ensembl.ids <- which(duplicated(merged.annot$EnsemblId))
dup.ensembl <-unique(merged.annot[dup.ensembl.ids, "EnsemblId"])
dup.ensembl.rows <- merged.annot[which(merged.annot$EnsemblId %in% dup.ensembl), ]
dim(dup.ensembl.rows)
### [1] 16  9
length(unique(dup.ensembl.rows$EntrezId))
dup.ensembl.rows[order(dup.ensembl.rows$EnsemblId, dup.ensembl.rows$EntrezId), ]

# EntrezId Symbol.x       EnsemblId Chr     Start       End    GC           Type Symbol.y
# 13248    23637  RABGAP1 ENSG00000011454   9 122940833 123104866 38.49 protein_coding  RABGAP1 
# 6841      2844    GPR21 ENSG00000011454   9 122940833 123104866 38.49 protein_coding  RABGAP1 GPR21 is a small gene that overlaps with RABGAP1. Sum the counts
# 9381     23499    MACF1 ENSG00000127603   1  39081316  39487177 42.57 protein_coding    MACF1 
# 8377    643314 KIAA0754 ENSG00000127603   1  39081316  39487177 42.57 protein_coding    MACF1 KIAA0754 is a small gene that overlaps with MACF1. Sum the counts
# 6367     26290   GALNT8 ENSG00000130035  12   4720341   4851927 43.83 protein_coding   GALNT8
# 8191      3742    KCNA6 ENSG00000130035  12   4720341   4851927 43.83 protein_coding   GALNT8 GALNT8 and KCNA6 are contiguous. Start - End covers both. Sum the counts
# 9816     84953  MICALCL ENSG00000133816  11  12094008  12359144 45.98 protein_coding   MICAL2
# 9814      9645   MICAL2 ENSG00000133816  11  12094008  12359144 45.98 protein_coding   MICAL2 MICALCL and MICAL2 are contiguous. Start - End covers both. Sum the counts
# 3964    285464   CRIPAK ENSG00000163945   4   1347266   1395992 56.55 protein_coding    UVSSA
# 8426     57654 KIAA1530 ENSG00000163945   4   1347266   1395992 56.55 protein_coding    UVSSA CRIPAK and KIAA1530 are contiguous. Start - End covers both. Sum the counts
# 8393     23285 KIAA1107 ENSG00000189195   1  92080305  92184725 38.18 protein_coding    BTBD8 
# 1570    284697    BTBD8 ENSG00000189195   1  92080305  92184725 38.18 protein_coding    BTBD8 KIAA1107 and BTBD8 are contiguous. Start - End covers both. Sum the counts
# 8834     10715    LASS1 ENSG00000223802  19  18868545  18896727 57.34 protein_coding    CERS1 
# 6465      2657     GDF1 ENSG00000223802  19  18868545  18896727 57.34 protein_coding    CERS1 LASS1 and GDF1 overlap. Start - End covers both. Sum the counts
# 5338      2074    ERCC6 ENSG00000225830  10  49455368  49539538 40.54 protein_coding    ERCC6
# 12046   267004    PGBD3 ENSG00000225830  10  49455368  49539538 40.54 protein_coding    ERCC6 PGBD3 is a small gene that overlaps with ERCC6. Sum the counts

#### Get the Counts
dup.ensembl.ids <- which(rownames(normal.df) %in% dup.ensembl.rows$EntrezId)
dup.ensembl.rows$EE <- paste(dup.ensembl.rows$EnsemblId, dup.ensembl.rows$EntrezId, sep=".")
length(dup.ensembl.ids)
dup.ensembl.counts <- normal.df[dup.ensembl.ids, ] 
rownames(dup.ensembl.counts) <- dup.ensembl.rows$EE
dup.ensembl.counts <- dup.ensembl.counts[order(rownames(dup.ensembl.counts)), ]

dup.entrez.ids <-which(duplicated(merged.annot$EntrezId))
dup.entrez <- unique(merged.annot[dup.entrez.ids, "EntrezId"])
dup.entrez.rows <- merged.annot[which(merged.annot$EntrezId %in% dup.entrez), ]
dim(dup.entrez.rows)
