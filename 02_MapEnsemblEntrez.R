setwd("~/Workspace/raw-counts-comparison")
library("biomaRt")

load(file="healthyRawNew.RData")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

entrezIds <- normal$Annot[,2]
ensembl.entrez = getBM(attributes = c('entrezgene', 'ensembl_gene_id'), 
                              filters = 'entrezgene', 
                              values = entrezIds, 
                              mart = ensembl)
names(ensembl.entrez) <- c("EntrezId", "EnsemblId")
dim(ensembl.entrez)
##[1] 21030     2
length(entrezIds)
##[1] 20532
########################################################################
## Bring the other annotation to remove genes from non conventional chromosomes
## NOTE: NCBI Gene in Biomart = to EntrezId
## I should've used Biomart annotations and no biomaRt
annot <- read.delim(file="Biomart_EnsemblG92_GRCh38_p12.txt", sep="\t")
names(annot) <- c("EnsemblId", "Start", "End", "Chr", "GC", "Type", "Symbol")
dim(annot)
## [1] 63967     7
annot.entrez <- merge(ensembl.entrez, annot, by = "EnsemblId", all=FALSE)
dim(annot.entrez)
##[1] 21038     8
levels(annot.entrez$Chr)
annot.entrez<-annot.entrez[annot.entrez$Chr%in%c(as.character(1:22), "X", "Y"),]
annot.entrez$Chr<-droplevels(annot.entrez$Chr)
table(annot.entrez$Chr)
dim(annot.entrez)
## [1] 18639     8

### Genes with no Symbol
annot.nosymbol <- annot.entrez[annot.entrez$Symbol=="", ]
annot.nosymbol
# EnsemblId  EntrezId     Start       End Chr    GC                 Type Symbol
# 529   ENSG00000036549     26009  77562416  77683419   1 37.17       protein_coding       
# 5760  ENSG00000126749     10436   6970893   6979941  12 48.30       protein_coding       
# 11834 ENSG00000167355    282763   5303444   5505652  11 38.82       protein_coding       
# 11864 ENSG00000167524    124923  28607964  28614200  17 45.98       protein_coding       
# 16508 ENSG00000197991     64881  61409858  61427849  13 38.48       protein_coding       
# 18561 ENSG00000235994 100128124 167975924 167979776   6 58.16 processed_transcript       
# 19212 ENSG00000260548      8471 108722070 108732687   X 36.11       protein_coding       
# 19388 ENSG00000268861     23370   7382834   7472477  19 51.37       protein_coding       
# 19396 ENSG00000269226    286527 104063871 104076212   X 41.27       protein_coding       
# 19439 ENSG00000272617     84342  69326913  69339667  16 51.52       protein_coding       
# 19451 ENSG00000272916      8509  73796514  73811651  10 53.82       protein_coding       
# 20377 ENSG00000278500      3233 176151085 176173097   2 51.67       protein_coding       
# 20523 ENSG00000280987      9782 139273752 139331671   5 40.63       protein_coding       
# 20528 ENSG00000281028     55300  25160663  25277306   4 41.33       protein_coding       
# 20734 ENSG00000283597    283777  98437162  98547728  15 43.49       protein_coding       
# 20791 ENSG00000284024     51182  14838160  14847018  10 44.63       protein_coding       
# 20814 ENSG00000284194      9997  50523568  50526439  22 67.13       protein_coding       
# 20870 ENSG00000284741     50940 177628069 178108339   2 39.13       protein_coding       
# 20875 ENSG00000284762      8622  77086732  77427124   5 41.54       protein_coding       
# 20880 ENSG00000284779      3481   2129112   2158658  11 60.39       protein_coding       
# 20892 ENSG00000284844    220074  72105924  72109329  11 53.11       protein_coding       
# 20893 ENSG00000284849    132332 121758930 121765185   4 41.75       protein_coding       
# 20895 ENSG00000284862    339829 180602858 180679500   3 35.65       protein_coding       
# 20927 ENSG00000285043       226  30053090  30070420  16 54.42       protein_coding       
# 20937 ENSG00000285077     89839  30624494  30638805  15 40.07       protein_coding       
# 20949 ENSG00000285130      9414  69035750  69255187   9 44.06       protein_coding       
# 20960 ENSG00000285188      5143  18207961  18255419  19 54.12       protein_coding       
# 20981 ENSG00000285253     54816  56652488  56918571  15 37.87       protein_coding       
# 21014 ENSG00000285396      5888  40695129  40732158  15 43.33       protein_coding       
# 21022 ENSG00000285441      6648 159679119 159762529   6 42.38       protein_coding       
# 21025 ENSG00000285447    169834 112997120 113050043   9 37.57       protein_coding   
# They do have EnsemblId EntrezId, check qwith Biomart
extra.symbols = getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'), 
                      filters = 'ensembl_gene_id', 
                      values = annot.nosymbol, 
                      mart = ensembl)
extra.symbols 
### Todos NA, still keep them
save(annot.entrez, file="annotEntrez.RData", compress = "xz")
#################################################################
## Load Data
load(file="annotEntrez.RData")
load(file="cancerRawLegacy.RData")
load(file="healthyRawLegacy.RData")
######## Pasar a data frames para mapear ensemblids

normal.df <- data.frame(normal$Counts)
head(normal.df)
dim(normal.df)
## [1] 20532   101
rownames(normal.df) <- as.numeric(normal$Annot[,2])
head(normal.df)
##[1] 20532     2
head(normal$Annot[,2])

all.entrez <- unique(annot.entrez$EntrezId)
length(all.entrez)
#[1] 18578

## Should remove entrez ids that are not in all.entrez
out.entrez.ids <- which(!(rownames(normal.df) %in% all.entrez))
normal.df <- normal.df[-out.entrez.ids, ]
dim(normal.df)
#[1] 18578   101
normal$Annot <- normal$Annot[-out.entrez.ids,]
dim(normal$Annot)
#[1] 18578     2

########################################################
## Merge and remove duplicates
annot.normal <- data.frame(normal$Annot, stringsAsFactors = FALSE)
##Add Biomart data
merged.annot<-merge(x=annot.normal, y=annot.entrez, by="EntrezId", all=F, sort=FALSE)
dim(merged.annot)
#[1] 18639     9

merged.annot[duplicated(merged.annot[,c("EntrezId", "EnsemblId")]), ]
#         EntrezId Symbol.x       EnsemblId    Start      End Chr   GC           Type Symbol.y
# 2978     6349   CCL3L1 ENSG00000276085 36194869 36196758  17 52.7 protein_coding   CCL3L3
# 2980   414062   CCL3L3 ENSG00000276085 36194869 36196758  17 52.7 protein_coding   CCL3L3
## Two EntrezIds with different Symbols map to the same EnsemblId 
merged.annot <- merged.annot[!duplicated(merged.annot[,c("EntrezId", "EnsemblId")]), ]
dim(merged.annot)
#[1] 18637     9

##### Duplicates. EnsemblId
dup.ensembl.ids <- which(duplicated(merged.annot$EnsemblId))
dup.ensembl <-unique(merged.annot[dup.ensembl.ids, "EnsemblId"])
dup.ensembl.rows <- merged.annot[which(merged.annot$EnsemblId %in% dup.ensembl), ]
dim(dup.ensembl.rows)
###[1] 77  9
## We need to remove the duplicated EnsemblIds. Lets keep them but only one Entrez counts should count
## probaby take max
length(unique(dup.ensembl.rows$EntrezId))
## [1] 72
dup.entrez.in.dup.ensembl <- dup.ensembl.rows[which(duplicated(dup.ensembl.rows$EntrezId)), c("EntrezId", "EnsemblId")]
dup.entrez.in.dup.ensembl <- merged.annot[merged.annot$EntrezId %in% dup.entrez.in.dup.ensembl$EntrezId 
                                              | merged.annot$EnsemblId %in% dup.entrez.in.dup.ensembl$EnsemblId , ]

dup.entrez.in.dup.ensembl <- dup.entrez.in.dup.ensembl[order(dup.entrez.in.dup.ensembl$EnsemblId, dup.entrez.in.dup.ensembl$EntrezId), ]
dup.ensembl.entrez.counts.ids <- which(rownames(normal.df) %in% dup.entrez.in.dup.ensembl$EntrezId)

#       EntrezId Symbol.x       EnsemblId    Start      End Chr    GC           Type Symbol.y
# 2786    220869    CBWD5 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5
# 2784    445571    CBWD3 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5
# 2788    644019    CBWD6 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5
# 6720 100132565  GOLGA8F ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F
# 6721    283768  GOLGA8G ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F
# 8947      3963   LGALS7 ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B
# 8945    653499  LGALS7B ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B
# 6719 100132565  GOLGA8F ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G
# 6722    283768  GOLGA8G ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G
# 8946      3963   LGALS7 ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7
# 8944    653499  LGALS7B ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7
# 2787    644019    CBWD6 ENSG00000215126 41131306 41199261   9 36.14 protein_coding    CBWD6

### Duplicates. Entrez.Ids
dup.entrez.ids <-which(duplicated(merged.annot$EntrezId))
dup.entrez <- unique(merged.annot[dup.entrez.ids, "EntrezId"])
dup.entrez.rows <- merged.annot[which(merged.annot$EntrezId %in% dup.entrez), ]
dim(dup.entrez.rows)
###[1] 117   9
length(unique(dup.ensembl.rows$EnsemblId))
## [1] 36
dup.ensembl.in.dup.entrez <- dup.entrez.rows[which(duplicated(dup.entrez.rows$EnsemblId)), c("EntrezId", "EnsemblId")]
dup.ensembl.in.dup.entrez <- merged.annot[merged.annot$EntrezId %in% dup.ensembl.in.dup.entrez$EntrezId 
                                              | merged.annot$EnsemblId %in% dup.ensembl.in.dup.entrez$EnsemblId , ]
dup.ensembl.in.dup.entrez <- dup.ensembl.in.dup.entrez[order(dup.ensembl.in.dup.entrez$EnsemblId, dup.ensembl.in.dup.entrez$EntrezId), ]

identical(dup.ensembl.in.dup.entrez, dup.entrez.in.dup.ensembl)
######################################
## dup.esenbl.in.dup.entrez has 150472 - ENSG00000215126 extra. 
## check one by one

#       EntrezId Symbol.x       EnsemblId    Start      End Chr    GC           Type Symbol.y
# 2786    220869    CBWD5 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5 https://www.ncbi.nlm.nih.gov/gene/220869 OK
# 2784    445571    CBWD3 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5 https://www.ncbi.nlm.nih.gov/gene/445571 CBWD3 ENSG00000196873 
# 2788    644019    CBWD6 ENSG00000147996 65668805 65734041   9 36.05 protein_coding    CBWD5 https://www.ncbi.nlm.nih.gov/gene/644019 CBWD6 ENSG00000215126
# 6720 100132565  GOLGA8F ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F https://www.ncbi.nlm.nih.gov/gene/100132565 OK
# 6721    283768  GOLGA8G ENSG00000153684 28378621 28392018  15 49.12 protein_coding  GOLGA8F https://www.ncbi.nlm.nih.gov/gene/283768 GOLGA8G ENSG00000183629
# 8947      3963   LGALS7 ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B https://www.ncbi.nlm.nih.gov/gene/3963 LGALS7 ENSG00000205076
# 8945    653499  LGALS7B ENSG00000178934 38789211 38791749  19 58.80 protein_coding  LGALS7B https://www.ncbi.nlm.nih.gov/gene/653499 OK
# 6719 100132565  GOLGA8F ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G https://www.ncbi.nlm.nih.gov/gene/100132565 GOLGA8F ENSG00000153684
# 6722    283768  GOLGA8G ENSG00000183629 28519611 28533014  15 49.13 protein_coding  GOLGA8G https://www.ncbi.nlm.nih.gov/gene/283768 OK
# 8946      3963   LGALS7 ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7 https://www.ncbi.nlm.nih.gov/gene/3963 OK
# 8944    653499  LGALS7B ENSG00000205076 38770971 38773492  19 58.52 protein_coding   LGALS7 https://www.ncbi.nlm.nih.gov/gene/653499 LGALS7B ENSG00000178934
# 2783    150472    CBWD2 ENSG00000215126 41131306 41199261   9 36.14 protein_coding    CBWD6 https://www.ncbi.nlm.nih.gov/gene/150472 CBWD2 ENSG00000136682
# 2787    644019    CBWD6 ENSG00000215126 41131306 41199261   9 36.14 protein_coding    CBWD6 https://www.ncbi.nlm.nih.gov/gene/644019 OK

merged.annot["2784", "Symbol.y"] <- "CBWD3"
merged.annot["2784", "EnsemblId"] <- "ENSG00000196873"
merged.annot["2784", ]

merged.annot["2788", "Symbol.y"] <- "CBWD6"
merged.annot["2788", "EnsemblId"] <- "ENSG00000215126"
merged.annot["2788", ]

merged.annot["6721", "Symbol.y"] <- "GOLGA8G"
merged.annot["6721", "EnsemblId"] <- "ENSG00000183629"
merged.annot["6721", ]

merged.annot["8947", "Symbol.y"] <- "LGALS7"
merged.annot["8947", "EnsemblId"] <- "ENSG00000205076"
merged.annot["8947", ]

merged.annot["6719", "Symbol.y"] <- "GOLGA8F"
merged.annot["6719", "EnsemblId"] <- "ENSG00000153684"
merged.annot["6719", ]

merged.annot["8944", "Symbol.y"] <- "LGALS7B"
merged.annot["8944", "EnsemblId"] <- "ENSG00000178934"
merged.annot["8944", ]

merged.annot["2783", "Symbol.y"] <- "CBWD2"
merged.annot["2783", "EnsemblId"] <- "ENSG00000136682"
merged.annot["2783", ]

### We can remove Entrezid - Ensemblid duplicates 
merged.annot[duplicated(merged.annot[,c("EntrezId", "EnsemblId")]), ]
merged.annot <- merged.annot[!duplicated(merged.annot[,c("EntrezId", "EnsemblId")]), ]
dim(merged.annot)
## [1] 18630     9

##### Duplicates. EnsemblId - Take2
dup.ensembl.ids <- which(duplicated(merged.annot$EnsemblId))
dup.ensembl <-unique(merged.annot[dup.ensembl.ids, "EnsemblId"])
dup.ensembl.rows <- merged.annot[which(merged.annot$EnsemblId %in% dup.ensembl), ]
dim(dup.ensembl.rows)
###[1] 64  9
length(unique(dup.ensembl.rows$EntrezId))
##[1] 64 YEI!!!
write.table(dup.ensembl.rows, file="dupEnsembl.txt", quote=FALSE, sep="\t")
dup.ensembl.rows$EntrezId <- as.numeric(dup.ensembl.rows$EntrezId)
dup.ensembl.rows <- dup.ensembl.rows[order(dup.ensembl.rows$EntrezId), ]

#        EntrezId  Symbol.x       EnsemblId     Start       End Chr    GC           Type Symbol.y  CorrectEnsemblAccordingToNCBI
# 2650        801     CALM1 ENSG00000198668  90396502  90408261  14 47.09 protein_coding    CALM1 OK
# 2653        808     CALM3 ENSG00000198668  90396502  90408261  14 47.09 protein_coding    CALM1 ENSG00000160014 CALM3
# 4450       1652       DDT ENSG00000099977  23971365  23980469  22 54.22 protein_coding      DDT OK
# 5338       2074     ERCC6 ENSG00000225830  49455368  49539538  10 40.54 protein_coding    ERCC6 OK
# 5918       2212    FCGR2A ENSG00000143226 161505430 161524013   1 43.82 protein_coding   FCGR2A OK
# 6338       2574    GAGE2C ENSG00000236362  49551278  49568218   X 44.54 protein_coding  GAGE12F No Ensembl
# 6341       2576     GAGE4 ENSG00000236362  49551278  49568218   X 44.54 protein_coding  GAGE12F No Ensembl
# 6465       2657      GDF1 ENSG00000223802  18868545  18896727  19 57.34 protein_coding    CERS1 OK
# 6841       2844     GPR21 ENSG00000011454 122940833 123104866   9 38.49 protein_coding  RABGAP1 OK. GPR21
# 7027       2966    GTF2H2 ENSG00000145736  71032670  71067689   5 38.44 protein_coding   GTF2H2 OK.
# 7373       3117  HLA-DQA1 ENSG00000237541  32741342  32747215   6 43.14 protein_coding HLA-DQA2 ENSG00000196735 HLA-DQA1
# 7375       3118  HLA-DQA2 ENSG00000237541  32741342  32747215   6 43.14 protein_coding HLA-DQA2 OK
# 8191       3742     KCNA6 ENSG00000130035   4720341   4851927  12 43.83 protein_coding   GALNT8 OK. KCNA6
# 12754      5554      PRH1 ENSG00000111215  10845849  10849475  12 40.20 protein_coding     PRR4 ENSG00000231887 PRH1
# 13941      6232     RPS27 ENSG00000185088  63125872  63158021  15 39.85 protein_coding   RPS27L ENSG00000177954 RPS27
# 2977       6349    CCL3L1 ENSG00000276085  36194869  36196758  17 52.70 protein_coding   CCL3L1 OK
# 16797      6399   TRAPPC2 ENSG00000196459  13712244  13734635   X 42.18 protein_coding  TRAPPC2 OK
# 15115      6606      SMN1 ENSG00000205571  70049612  70078522   5 42.22 protein_coding     SMN2 ENSG00000172062 SMN1
# 15116      6607      SMN2 ENSG00000205571  70049612  70078522   5 42.22 protein_coding     SMN2 OK
# 5920       9103    FCGR2C ENSG00000143226 161505430 161524013   1 43.82 protein_coding   FCGR2A ENSG00000244682 FCGR2C
# 9814       9645    MICAL2 ENSG00000133816  12094008  12359144  11 45.98 protein_coding   MICAL2 OK
# 16795     10597 TRAPPC2P1 ENSG00000196459  13712244  13734635   X 42.18 protein_coding  TRAPPC2 ENSG00000256060 TRAPPC2B
# 8834      10715     LASS1 ENSG00000223802  18868545  18896727  19 57.34 protein_coding    CERS1 OK CERS1
# 12825     10887    PROKR1 ENSG00000169621  68467561  68655862   2 37.87 protein_coding     APLF ENSG00000169618 PROKR1
# 16590     11013   TMSB15A ENSG00000158164 102513676 102516784   X 43.49 protein_coding  TMSB15A OK
# 13255     11158    RABL2B ENSG00000079974  50767501  50783663  22 50.05 protein_coding   RABL2B OK
# 13253     11159    RABL2A ENSG00000079974  50767501  50783663  22 50.05 protein_coding   RABL2B ENSG00000144134 RABL2A
# 12876     11272      PRR4 ENSG00000111215  10845849  10849475  12 40.20 protein_coding     PRR4 OK
# 8393      23285  KIAA1107 ENSG00000189195  92080305  92184725   1 38.18 protein_coding    BTBD8 OK BTBD8
# 9381      23499     MACF1 ENSG00000127603  39081316  39487177   1 42.57 protein_coding    MACF1 OK
# 13248     23637   RABGAP1 ENSG00000011454 122940833 123104866   9 38.49 protein_coding  RABGAP1 OK
# 6367      26290    GALNT8 ENSG00000130035   4720341   4851927  12 43.83 protein_coding   GALNT8 OK
# 6340      26749    GAGE2E ENSG00000236362  49551278  49568218   X 44.54 protein_coding  GAGE12F No Ensembl
# 13940     51065    RPS27L ENSG00000185088  63125872  63158021  15 39.85 protein_coding   RPS27L OK
# 6879      51463    GPR89B ENSG00000117262 145607990 145670648   1 38.91 protein_coding   GPR89A ENSG00000188092 GPR89B
# 8426      57654  KIAA1530 ENSG00000163945   1347266   1395992   4 56.55 protein_coding    UVSSA OK UVSSA
# 8769      83899  KRTAP9-2 ENSG00000239886  41226648  41227652  17 51.74 protein_coding KRTAP9-2 OK
# 8771      83901  KRTAP9-8 ENSG00000239886  41226648  41227652  17 51.74 protein_coding KRTAP9-2 ENSG00000187272 KRTAP9-8
# 9816      84953   MICALCL ENSG00000133816  12094008  12359144  11 45.98 protein_coding   MICAL2 OK
# 18297     90333    ZNF468 ENSG00000204604  52838008  52890375  19 43.20 protein_coding   ZNF468 OK
# 18208    162967    ZNF320 ENSG00000204604  52838008  52890375  19 43.20 protein_coding   ZNF468 ENSG00000182986 ZNF320
# 782      200558      APLF ENSG00000169621  68467561  68655862   2 37.87 protein_coding     APLF OK
# 11645    203569     PAGE2 ENSG00000234068  55089008  55092842   X 39.14 protein_coding    PAGE2 OK
# 12046    267004     PGBD3 ENSG00000225830  49455368  49539538  10 40.54 protein_coding    ERCC6 OK
# 1570     284697     BTBD8 ENSG00000189195  92080305  92184725   1 38.18 protein_coding    BTBD8 OK
# 3964     285464    CRIPAK ENSG00000163945   1347266   1395992   4 56.55 protein_coding    UVSSA OK
# 16593    286527   TMSB15B ENSG00000158164 102513676 102516784   X 43.49 protein_coding  TMSB15A ENSG00000158427 TMSB15B
# 10664    340527     NHSL2 ENSG00000204131  71910818  72161750   X 45.33 protein_coding    NHSL2 OK
# 11643    389860    PAGE2B ENSG00000234068  55089008  55092842   X 39.14 protein_coding    PAGE2 ENSG00000238269 PAGE2B
# 6073     392490  FLJ44635 ENSG00000204131  71910818  72161750   X 45.33 protein_coding    NHSL2 No Ensembl
# 15941    414059   TBC1D3B ENSG00000274808  36165681  36176636  17 59.93 protein_coding  TBC1D3B OK
# 15943    414060   TBC1D3C ENSG00000274808  36165681  36176636  17 59.93 protein_coding  TBC1D3B ENSG00000278299 TBC1D3C
# 2979     414062    CCL3L3 ENSG00000276085  36194869  36196758  17 52.70 protein_coding   CCL3L1 OK
# 4080     441519    CT45A3 ENSG00000269096 135759846 135768191   X 42.77 protein_coding   CT45A3 OK
# 4078     541466    CT45A1 ENSG00000269096 135759846 135768191   X 42.77 protein_coding   CT45A3 ENSG00000268940 CT45A1
# 8377     643314  KIAA0754 ENSG00000127603  39081316  39487177   1 42.57 protein_coding    MACF1 OK
# 6337     645037    GAGE2B ENSG00000189064  49589529  49596827   X 44.17 protein_coding   GAGE2A No Ensembl
# 6878     653519    GPR89A ENSG00000117262 145607990 145670648   1 38.91 protein_coding   GPR89A OK
# 7025     728340   GTF2H2C ENSG00000145736  71032670  71067689   5 38.44 protein_coding   GTF2H2 ENSG00000183474 GTF2H2C
# 6339     729408    GAGE2D ENSG00000236362  49551278  49568218   X 44.54 protein_coding  GAGE12F No Ensembl
# 6336     729447    GAGE2A ENSG00000189064  49589529  49596827   X 44.17 protein_coding   GAGE2A OK
# 6332  100008586   GAGE12F ENSG00000236362  49551278  49568218   X 44.54 protein_coding  GAGE12F OK
# 4448  100037417      DDTL ENSG00000099977  23971365  23980469  22 54.22 protein_coding      DDT ENSG00000099974 DDTL
# 6343  100101629     GAGE8 ENSG00000189064  49589529  49596827   X 44.17 protein_coding   GAGE2A No Ensembl


### Remove the entries that don't have Enembl Id
rownames.to.remove <-c("6338", "6341", "6340", "6073", "6337", "6339", "6343") 
merged.annot[rownames(merged.annot) %in% rownames.to.remove, ]
merged.annot <- merged.annot[!rownames(merged.annot) %in% rownames.to.remove, ]
dim(merged.annot)
##[1] 18623     9

#### Get the Counts
dup.ensembl.ids <- which(rownames(normal.df) %in% dup.ensembl.rows$EntrezId)
dup.ensembl.rows$EE <- paste(dup.ensembl.rows$EnsemblId, dup.ensembl.rows$EntrezId, sep=".")
length(dup.ensembl.ids)
dup.ensembl.counts <- normal.df[dup.ensembl.ids, ] 
rownames(dup.ensembl.counts) <- dup.ensembl.rows$EE
dup.ensembl.counts <- dup.ensembl.counts[order(rownames(dup.ensembl.counts)), ]



##### Dup Entrez. Ids
dup.entrez.ids <-which(duplicated(merged.annot$EntrezId))
dup.entrez <- unique(merged.annot[dup.entrez.ids, "EntrezId"])
dup.entrez.rows <- merged.annot[which(merged.annot$EntrezId %in% dup.entrez), ]
dim(dup.entrez.rows)
###[1] 103   9
length(unique(dup.ensembl.rows$EnsemblId))
## [1] 30
dup.ensembl.in.dup.entrez <- dup.entrez.rows[which(duplicated(dup.entrez.rows$EnsemblId)), c("EntrezId", "EnsemblId")]
dup.ensembl.in.dup.entrez
## <0 rows> (or 0-length row.names) YEI!!!
dup.entrez.rows <- dup.entrez.rows[order(dup.entrez.rows$EntrezId, dup.entrez.rows$EnsemblId), ]

# head(dup.entrez.rows)
# EntrezId Symbol.x       EnsemblId     Start       End Chr    GC           Type Symbol.y
# 71     10061    ABCF2 ENSG00000285292 151207837 151227166   7 49.13 protein_coding    ABCF2 //Stable ID ENSG00000285292 not present in GRCh37.
# 72     10061    ABCF2 ENSG00000033050 151212490 151227230   7 48.51 protein_coding    ABCF2 OK 
# 512      226    ALDOA ENSG00000149925  30064164  30070457  16 59.26 protein_coding    ALDOA
# 513      226    ALDOA ENSG00000285043  30053090  30070420  16 54.42 protein_coding          //Stable ID ENSG00000285043 not present in GRCh37.
# 909    23370 ARHGEF18 ENSG00000268861   7382834   7472477  19 51.37 protein_coding         
# 910    23370 ARHGEF18 ENSG00000104880   7395113   7472477  19 51.44 protein_coding ARHGEF18








