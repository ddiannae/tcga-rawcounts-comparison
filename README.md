# tcga-rawcounts-comparison

The goal of this project is to compare raw RNASeq counts from the TCGA legacy database with the new database. 
Experimentally, the main change between them is that all samples were re-aligned to the new genome assembly version GRCh38/hg38 (Dec 2013). Legacy database was aligned using GRCh27/hg19 (Feb. 2009). 
The change is important because we want to make sure we'll catch the same features with the new database regarding the loss of trans-regulation project. 

We've tried to compare the same legacy/new raw counts file. In order to do that, case ids and aliquote submitter ids (when possible) were matched. We kept 773 tumor samples and 101 healthy samples. 

Gene IDs used by TCGA for legacy and new files are different. The legacy database uses Entrez ID (NCBI) while the new database uses Ensembl ID.  Mapping was done using BioMart and manual curation. We kept 18547 features. 

Spearman's correlation between each legacy/new sample was calculated. 


