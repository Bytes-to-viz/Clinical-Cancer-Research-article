#R version 3.6.1 (2019-07-05) -- "Action of the Toes"
#Copyright (C) 2019 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)



#TITLE: ANALYSIS OF THE TCGA-HNSC DATASET



# PART 0: PREPARE THE ENVIRONMENT --------------------------------

#1. Load (or install) all necessary packages.
library("readxl")
library("writexl")
library("readr")
library("ComplexHeatmap")
library("circlize")
library("ggplot2")
library("gridExtra")
library("gplots")
library("mclust")
library("GSEABase")
library("methods")
library("edgeR")
library("geneplotter")
library("genefilter")
library("BiocGenerics")
library("Biobase")
library("graph")
library("XML")
library("lattice")
library("limma")
library("shinythemes")
library("shiny")
library("RColorBrewer")
library("parallel")
library("cluster")
library("Matrix")
library("locfit")
library("snow")
library("GSVA") 
library("dplyr")
library("ggplot2")
library("data.table")
library ("remotes")
library("plyr")
library("magrittr")
library("OIsurv")
library("survival")
library("KMsurv")
library("splines")
library("survminer")
library("ggpubr")
library("survutils")
library("scales")
library("ggpubr")
library("tidyverse")
library("corrr")
library("igraph")
library("ggraph")
library("tidygraph")
library("CoxBoost")
library("glmnet") 
library("randomForest") 
library("class") 
library("dml") 
library("MASS")
library("readr")
library("Rtsne") 
library("stats") 
library("ggridges")
library("gdata")
library("ggrepel")
library("corrplot")
library("ggExtra")
library("gridExtra")
library("rstatix")   
library("ggpubr")
library("viridis")
library("corrplot")
library("xlsx")
library("plotly")
library("plot3D")
library("tidyr")
library("dplyr")



# PART 1: DATA IMPORTATION-----------------------------------------

#1. Import the TGFB and ALTEJ genelists. 
TGFBgeneset <- read.table("Input/TGFB_list.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
ALTEJgeneset <- read.table("Input/ALTEJ_list.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
TGFBlist <- as.list(TGFBgeneset[,1:1]) #55 genes (50 + synonyms)
ALTEJlist <- as.list(ALTEJgeneset[,1:1]) #45 genes (36 + synonyms)
Bothgenelists <- list(TGFBUPgeneset=TGFBlist, ALTEJgeneset=ALTEJlist)

#2. Import the files with genes' weights. 
TGFBimportance <- read_excel("Input/Relative importance of BAlt genes.xlsx", sheet = "TGFB")
ALTEJimportance <- read_excel("Input/Relative importance of BAlt genes.xlsx", sheet = "ALTEJ")
##These files are the ones generated based on the results from the "inhouse HNSC dataset", and can also be found in "1. Inhouse-HNSC/Output".

#3. Import the file with gene expression.
Expression <- read.delim("Input/ExpressionHiSeqV2_HNSC_UCSCXena_Nov2020.txt") 
##This file was downloaded from the dataset "TCGA.HNSC.sampleMap/HiSeqV2" of the UCSC Xena platform using the R package UCSCXenaTools in November 2020. 
##Contains gene expression measured by RNAseq with the platform	IlluminaHiSeq_RNASeqV2. 
##Values (Log2(norm_count+1)) had been RSEM normalized and log2(x+1) transformed.

#4. Import the clinical information.
Survival <- read_excel("Input/TCGA-CDR-SupplementalTableS1.xlsx")
Survival <- Survival[,c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "gender", "race", 
                        "ajcc_pathologic_tumor_stage", "clinical_stage", "histological_grade", 
                        "treatment_outcome_first_course", "OS", "OS.time", "PFI", "PFI.time")] #34->13 variables.
##This file was downloaded the dataset "TCGA-CDR-SupplementalTableS1.xl" from the GDC in November 2020. 
##This file is the main dataset with clinical information, such as OS/PFS, tumor stage and age.
Clin_additional1 <- read.delim("Input/hnsc_tcga_pan_can_atlas_2018_clinical_data.tsv")
Clin_additional1 <- Clin_additional1[,c("Subtype", "Patient.ID", "Sample.ID", 
                                        "Fraction.Genome.Altered", "Mutation.Count")]
##This file was downloaded from cBioPortal.  
##This file contains some additional clinical information (fraction of the genome altered, HPV status).
Clin_additional2 <- read_delim("Input/Phenotype_HNSC_UCSCXena_Nov2020.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Clin_additional2 <- Clin_additional2[,c("sampleID", "_PATIENT", "anatomic_neoplasm_subdivision", 
                                        "method_of_curative_tx", "additional_pharmaceutical_therapy", 
                                        "additional_radiation_therapy", "radiation_therapy",
                                        "history_of_neoadjuvant_treatment", "postoperative_rx_tx", 
                                        "targeted_molecular_therapy")]
##This file was downloaded from the UCSC Xena platform using the R package UCSCXenaTools in November 2020. 
##This file contains some additional clinical information (treatment).



# PART 2: ELIMINATE DUPLICATED, RECURRENT AND NORMAL TISSUE SAMPLES -----------------------------------------

#1. Change points for slashes in the column names of the expression dataset.
names(Expression) <- gsub(x = names(Expression), pattern = "\\.", replacement = "-") 

#2. Transpose the dataset with gene expression and create a variable that shows the sample types. 
rownames(Expression) <- Expression$sample
Expression$sample <-NULL
ExpressionT <- as.data.frame(t(Expression))
ExpressionT$Sampletype <- rownames(ExpressionT)
ExpressionT$Sampletype <- substr(ExpressionT$Sampletype, 13, 15) #Keep only sample type numbers.
table(ExpressionT$Sampletype) #520 primary tumors, 44 normal tissue samples, 2 metastasic tumor samples.

#3. Slect only primary solid tumor samples. 
ExpressionT <- ExpressionT[which(ExpressionT$Sampletype=="-01"),] #566->520 samples.
ExpressionT$Sampletype <- NULL

#4. Check if some patients have more than one sample.
ExpressionT$Patient <- rownames(ExpressionT)
ExpressionT$Patient <- substr(ExpressionT$Patient, 1, 12) #Keep only patient ID numbers.
Duplicated <- ExpressionT[duplicated(ExpressionT$Patient), ] #There are no duplicated patients.
ExpressionT$Patient <- NULL

#5. Transpose back the dataframe with gene expression.
Expression <- as.data.frame(t(ExpressionT))



# PART 3.A: SCORES CALCULATION - WEIGHTED BALT -----------------------------------------

#1. Turn the dataframe with gene expression into a matrix.
Expression <- as.matrix(Expression)

#2.  Mean normalize all values per gene by transforming them into Z-scores ((x - mean(Xcolumn)) / sd(Xcolumn)).
genesT <- t(Expression) #Transpose the matrix so that each gene is a column. 
genesT<-scale(genesT, center = TRUE, scale = TRUE) #calculate the Z-score of each gene.
genes<-t(genesT)

#3. Create a dataset with only the expression of the TGFB and ALTEJ genes. 
genes <- genes[which(rownames(genes) %in% TGFBlist | rownames(genes) %in% ALTEJlist), ] #86 genes.

#4. Edit gene names so that they match the ones from the genes' weight dataframe. 
genes <- as.data.frame(genes)
genes$Gene <- rownames(genes)
genes$Gene <- ifelse(genes$Gene=="APEX2", "APE2", 
                     ifelse(genes$Gene=="HRAS", "HRAS1",
                     ifelse(genes$Gene=="PRPF19", "PRP19",
                     ifelse(genes$Gene=="RAD51L3", "RAD51D",
                     ifelse(genes$Gene=="KAT5", "TIP60",
                     ifelse(genes$Gene=="ARHGAP32", "RICS",
                     ifelse(genes$Gene=="C19orf40", "FAAP24",
                     ifelse(genes$Gene=="CYTH1", "PSCD1",
                     ifelse(genes$Gene=="NUDT1", "MTH1",
                     ifelse(genes$Gene=="OBFC2B", "NABP2", 
                     ifelse(genes$Gene=="STAG1", "TMEPAI", genes$Gene)))))))))))

#5. Merge the genes and their weights in one dataframe. 
TGFBimportance$Gene <- TGFBimportance$Gene
genesimportance <- rbind(TGFBimportance[,c("Gene", "mean weight")], ALTEJimportance[,c("Gene", "mean weight")])
genes2 <- merge(genesimportance, genes, 
               by.x="Gene", by.y="Gene", all.x=TRUE, all.y=TRUE) #86 genes.

#6. Calculate the TGFB and ALTEJ weighted scores in each sample, by multiplying each gene by its factor.
##Factor = 1+weight if weight>0; or 1 if weight<0. 
genes2$factor <- ifelse(genes2$`mean weight`< 0, 1, genes2$`mean weight`)
genes2$factor <- ifelse(genes2$factor==1, 1, 1+genes2$factor)
rownames(genes2) <- genes2$Gene
scores <- sapply(genes2[,-which(colnames(genes2) %in% c("Gene", "mean weight", "factor"))], '*', (genes2$factor))
rownames(scores) <- rownames(genes2)
scores <- as.data.frame(t(scores))
scores$TGFBUPgeneset <- rowSums(scores[ , which(colnames(scores) %in% TGFBlist)], na.rm=TRUE)
scores$ALTEJgeneset <- rowSums(scores[ , which(colnames(scores) %in% ALTEJlist)], na.rm=TRUE)
scores <- scores[,c("TGFBUPgeneset", "ALTEJgeneset")]
cor.test(scores$TGFBUPgeneset, 
         scores$ALTEJgeneset, 
         method = "pearson") #Test if TGFB and ALTEJ are anticorrelated: PCC= -0.296.

#7. Create a new variable that is the weighted Balt score.
scores$balt <- sqrt((max(scores$ALTEJgeneset)-scores$ALTEJgeneset)^2+
                      (min(scores$TGFBUPgeneset)-scores$TGFBUPgeneset)^2) - 
               sqrt((min(scores$ALTEJgeneset)-scores$ALTEJgeneset)^2+
                      (max(scores$TGFBUPgeneset)-scores$TGFBUPgeneset)^2)
scores$balt <- scores$balt * -1 

#8. Create a new variable that are the weighted BAlt score tertiles. 
scores$tertile <- cut(scores$balt, quantile(scores$balt, c(0, 1/3, 2/3, 1)), na.rm=TRUE, include.lowest=TRUE, 
                      labels = c("Low", "Middle", "High"))
scores$Tertile <- ifelse(scores$tertile == "Low", "high TGFβ and low ALTEJ", 
                         ifelse(scores$tertile == "High", "low TGFβ and high ALTEJ", NA))
scores$Tertile <- as.character(scores$Tertile)

#9. Create a clean dataframe with the results. 
scores$balt_weighted <- scores$balt
scores$TGFBUPgeneset_weighted <- scores$TGFBUPgeneset
scores$ALTEJgeneset_weighted <- scores$ALTEJgeneset
scores$tertile_weighted <- scores$tertile
scores$Tertile_weighted <- scores$Tertile
scores_weighted <- scores[,c("balt_weighted", "TGFBUPgeneset_weighted", "ALTEJgeneset_weighted", 
                             "tertile_weighted", "Tertile_weighted")]



# PART 3.B: SCORES CALCULATION - ORIGINAL BALT AND OTHER VARIATIONS -----------------------------------------

#1. Create some additional genelists. 
##Short (list of all TGFB and altEJ genes with a positive weight).
TGFBshortlist <- genesimportance[which((genesimportance$Gene %in% TGFBlist) & genesimportance$`mean weight`>0),] 
TGFBshortlist <- as.list(TGFBshortlist$Gene) #39 genes.
ALTEJshortlist <- genesimportance[which((genesimportance$Gene %in% ALTEJlist) & genesimportance$`mean weight`>0),] 
ALTEJshortlist <- as.list(ALTEJshortlist$Gene) #32 genes.
Bothshortlists <- list(TGFBUPgeneset=TGFBshortlist, ALTEJgeneset=ALTEJshortlist)
##Top15 (list of the top 15 TGFB and altEJ genes with the highest weight).
TGFBtop15list <- genesimportance[which(genesimportance$Gene %in% TGFBlist),]
TGFBtop15list <- top_n(TGFBtop15list, n=15, TGFBtop15list$`mean weight`)
TGFBtop15list <- as.list(TGFBtop15list$Gene)
ALTEJtop15list <- genesimportance[which(genesimportance$Gene %in% ALTEJlist),]
ALTEJtop15list <- top_n(ALTEJtop15list, n=15, ALTEJtop15list$`mean weight`)
ALTEJtop15list <- as.list(ALTEJtop15list$Gene)
Bothtoplists <- list(TGFBUPgeneset=TGFBtop15list, ALTEJgeneset=ALTEJtop15list)
##Bottom15 (list of the bottom 15 TGFB and altEJ genes with the lowest weight).
TGFBbottom15list <- genesimportance[which(genesimportance$Gene %in% TGFBlist),]
TGFBbottom15list <- top_n(TGFBbottom15list, n=15, TGFBbottom15list$`mean weight`*-1)
TGFBbottom15list <- as.list(TGFBbottom15list$Gene)
ALTEJbottom15list <- genesimportance[which(genesimportance$Gene %in% ALTEJlist),]
ALTEJbottom15list <- top_n(ALTEJbottom15list, n=15, ALTEJbottom15list$`mean weight`*-1)
ALTEJbottom15list <- as.list(ALTEJbottom15list$Gene)
Bothbottomlists <- list(TGFBUPgeneset=TGFBbottom15list, ALTEJgeneset=ALTEJbottom15list)

#2. Create a matrix with the normalized (z-score transformed) expression of all genes. 
Expression <- as.matrix(Expression)
genesT <- t(Expression) #Transpose the matrix so that each gene is a column. 
genesT<-scale(genesT, center = TRUE, scale = TRUE) #Z-score transformation.
genes<-t(genesT)

#3. Calculate the Original Balt. 
ssgsea <- gsva(genes, Bothgenelists, method=c("ssgsea"))
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
cor.test(ssgsea$TGFBUPgeneset, ssgsea$ALTEJgeneset) #Test if TGFB and ALTEJ are anticorrelated: PCC -0.35.
ssgsea$balt <- sqrt((max(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (min(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2) - 
               sqrt((min(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (max(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2)
ssgsea$balt <- ssgsea$balt * -1 
ssgsea$balt_original <- ssgsea$balt
ssgsea$tertile <- cut(ssgsea$balt, quantile(ssgsea$balt, c(0, 1/3, 2/3, 1)), 
                      na.rm=TRUE, include.lowest=TRUE, labels = c("Low", "Middle", "High"))
ssgsea$Tertile_original <- ifelse(ssgsea$tertile == "Low", "high TGFβ and low ALTEJ", 
                         ifelse(ssgsea$tertile == "High", "low TGFβ and high ALTEJ", NA))
ssgsea_original <- ssgsea
ssgsea_original <- ssgsea_original[,c("balt_original", "Tertile_original")]

#4. Calculate the Short Balt (based on all TGFB and altEJ genes with a positive weight).
ssgsea <- gsva(genes, Bothshortlists, method=c("ssgsea"))
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
cor.test(ssgsea$TGFBUPgeneset, ssgsea$ALTEJgeneset) #Test if TGFB and ALTEJ are anticorrelated: PCC -0.37.
ssgsea$balt <- sqrt((max(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (min(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2) - 
               sqrt((min(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (max(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2)
ssgsea$balt <- ssgsea$balt * -1 
ssgsea$balt_short <- ssgsea$balt
ssgsea$tertile <- cut(ssgsea$balt, quantile(ssgsea$balt, c(0, 1/3, 2/3, 1)), 
                      na.rm=TRUE, include.lowest=TRUE, labels = c("Low", "Middle", "High"))
ssgsea$Tertile_short <- ifelse(ssgsea$tertile == "Low", "high TGFβ and low ALTEJ", 
                              ifelse(ssgsea$tertile == "High", "low TGFβ and high ALTEJ", NA))
ssgsea_short <- ssgsea
ssgsea_short <- ssgsea_short[,c("balt_short", "Tertile_short")]

#5. Calculate the Top15 Balt (based on the top 15 TGFB and altEJ genes with the highest weight).
ssgsea <- gsva(genes, Bothtoplists, method=c("ssgsea"))
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
cor.test(ssgsea$TGFBUPgeneset, ssgsea$ALTEJgeneset) #Test if TGFB and ALTEJ are anticorrelated: PCC -0.34.
ssgsea$balt <- sqrt((max(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (min(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2) - 
               sqrt((min(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (max(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2)
ssgsea$balt <- ssgsea$balt * -1 
ssgsea$balt_top <- ssgsea$balt
ssgsea$tertile <- cut(ssgsea$balt, quantile(ssgsea$balt, c(0, 1/3, 2/3, 1)), 
                      na.rm=TRUE, include.lowest=TRUE, labels = c("Low", "Middle", "High"))
ssgsea$Tertile_top <- ifelse(ssgsea$tertile == "Low", "high TGFβ and low ALTEJ", 
                            ifelse(ssgsea$tertile == "High", "low TGFβ and high ALTEJ", NA))
ssgsea_top <- ssgsea
ssgsea_top <- ssgsea_top[,c("balt_top", "Tertile_top")]

#6. Calculate the Bottom15 Balt (based on the bottom 15 TGFB and altEJ genes with the lowest weight).
ssgsea <- gsva(genes, Bothbottomlists, method=c("ssgsea"))
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
cor.test(ssgsea$TGFBUPgeneset, ssgsea$ALTEJgeneset) #Test if TGFB and ALTEJ are anticorrelated: PCC -0.32.
ssgsea$balt <- sqrt((max(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (min(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2) - 
               sqrt((min(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (max(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2)
ssgsea$balt <- ssgsea$balt * -1 
ssgsea$balt_bottom <- ssgsea$balt
ssgsea$tertile <- cut(ssgsea$balt, quantile(ssgsea$balt, c(0, 1/3, 2/3, 1)), 
                      na.rm=TRUE, include.lowest=TRUE, labels = c("Low", "Middle", "High"))
ssgsea$Tertile_bottom <- ifelse(ssgsea$tertile == "Low", "high TGFβ and low ALTEJ", 
                               ifelse(ssgsea$tertile == "High", "low TGFβ and high ALTEJ", NA))
ssgsea_bottom <- ssgsea
ssgsea_bottom <- ssgsea_bottom[,c("balt_bottom", "Tertile_bottom")]



# PART 4: DATA WRANGLING -----------------------------------------

#1.Merge all the scores together in one dataframe. 
ssgsea_original$ID <- rownames(ssgsea_original)
ssgsea_short$ID <- rownames(ssgsea_short)
ssgsea_top$ID <- rownames(ssgsea_top)
ssgsea_bottom$ID <- rownames(ssgsea_bottom)
scores_weighted$ID <- rownames(scores_weighted)
scores <- merge(scores_weighted, ssgsea_original, 
                by.x="ID", by.y="ID", 
                all.x=TRUE, all.y=TRUE) #520+520->520 samples.
scores <- merge(scores, ssgsea_short, 
                by.x="ID", by.y="ID", 
                all.x=TRUE, all.y=TRUE) #520+520->520 samples.
scores <- merge(scores, ssgsea_top, 
                by.x="ID", by.y="ID", 
                all.x=TRUE, all.y=TRUE) #520+520->520 samples.
scores <- merge(scores, ssgsea_bottom, 
                by.x="ID", by.y="ID", 
                all.x=TRUE, all.y=TRUE) #520+520->520 samples.
scores$Patient <- substr(scores$ID, 1, 12) 

#2. Add the clinical information into the dataframe with the scores. 
ALL <- merge(scores, Survival, 
             by.x = "Patient", by.y = "bcr_patient_barcode",
             all.x=TRUE, all.y=FALSE) #520 patients with gene expression, all of them with clinical data.
ALL <- merge(ALL, Clin_additional1,
             by.x = "Patient", by.y = "Patient.ID",
             all.x=TRUE, all.y=FALSE) #520 patients.
ALL <- merge(ALL, Clin_additional2,
             by.x = "ID", by.y = "sampleID",
             all.x=TRUE, all.y=FALSE) #520 patients.

#3. Group the tumor stages into fewer groups. 
table(ALL$ajcc_pathologic_tumor_stage)
ALL$Stage <- ifelse(ALL$ajcc_pathologic_tumor_stage=="Stage I", "I",
                    ifelse(ALL$ajcc_pathologic_tumor_stage=="Stage II", "II",
                    ifelse(ALL$ajcc_pathologic_tumor_stage=="Stage III", "III",
                    ifelse(ALL$ajcc_pathologic_tumor_stage=="Stage IV", "IV",
                    ifelse(ALL$ajcc_pathologic_tumor_stage=="Stage IVA", "IV",
                    ifelse(ALL$ajcc_pathologic_tumor_stage=="Stage IVB", "IV",
                    ifelse(ALL$ajcc_pathologic_tumor_stage=="Stage IVC", "IV", NA)))))))
ALL$stage <- ifelse(ALL$Stage=="I", "I-II",
                    ifelse(ALL$Stage=="II", "I-II",
                    ifelse(ALL$Stage=="III", "III-IV",
                    ifelse(ALL$Stage=="IV", "III-IV", NA))))

#4.Prepare the survival variables so that they are in the appropriate format. 
table(ALL$OS)
class(ALL$OS.time) #[1] "numeric"
ALL$OS.months <- ALL$OS.time/30.5
table(ALL$PFI)
class(ALL$PFI.time) #[1] "numeric"
ALL$PFI.months <- ALL$PFI.time/30.5

#5. Tidy other variables. 
ALL$age <- ALL$age_at_initial_pathologic_diagnosis
table(ALL$Subtype)
ALL$Subtype <- ifelse(ALL$Subtype == "HNSC_HPV-", "HPV-", 
                      ifelse(ALL$Subtype == "HNSC_HPV+", "HPV+", NA))

#6. Check if some patients have more than one sample. 
Duplicated <- ALL[duplicated(ALL$Patient), ] #There are no patients with more than one sample. 

#7. Export the "ALL" dataframe. 
write_xlsx(ALL, "Output/balt_scores_TCGAHNSC.xlsx")



# PART 5: DESCRIPTIVE DATA ANALYSIS -----------------------------------------

#1. Check HPV status. 
table(ALL$Subtype) #72 HPV+. 

#2. Check stage. 
table(ALL$Stage)

#3. Check tumor location. 
table(ALL$anatomic_neoplasm_subdivision)

#4. Explore the treatment variables.
table(ALL$method_of_curative_tx) #101 IQ, 36 RT, 44 RT-ChT, 7 ChT, 332 NA. 
##ChT
table(ALL$additional_pharmaceutical_therapy) #52Y, 64N. 
table(ALL$history_of_neoadjuvant_treatment) #10Y, 510N.
table(ALL$targeted_molecular_therapy) #148Y, 260N.
##RT
table(ALL$additional_radiation_therapy) #29Y, 87N.
table(ALL$postoperative_rx_tx) #64Y, 115N.
table(ALL$radiation_therapy) #291Y, 159N.



# PART 6: SCATTERPLOT OF THE CORRELATION BETWEEN TGFβ AND ALTEJ -----------------------------------------

#1. Scatterplot of TGFB versus ALTEJ weighted scores, colored by weighted Balt tertiles. 
p <- ggplot(ALL, aes(x=TGFBUPgeneset_weighted, y=ALTEJgeneset_weighted)) + 
  geom_point(size=5, alpha=0.7,shape=21, stroke=0.1, color="white", aes(fill=tertile_weighted)) +
  scale_fill_manual(values=c("indianred1", "bisque3", "steelblue4")) + #Or "grey50"
  labs(x = "TGFβ score", y = "alt-EJ score", col="βalt tertile") +
  geom_smooth(method=loess, fullrange=TRUE) + 
  #geom_smooth(method="glm", fullrange=TRUE) +
  stat_cor(method="pearson", size=5) +
  theme(text = element_text(size=15)) +
  theme(legend.title=element_text(size=14))
p <- p + rremove("legend")
p1 <- ggMarginal(p,  size=10, fill = "grey65"); p1 #PDF 5x6.5

#2. Calculate Pearson and Spearman correlation coefficients.
cor.test(ALL$TGFBUPgeneset, ALL$ALTEJgeneset, method = "pearson")
cor.test(ALL$TGFBUPgeneset, ALL$ALTEJgeneset, method = "spearman")



# PART 7: VIOLIN PLOTS OF THE ASSOCIATION OF THE WEIGHTED βALT WITH OTHER VARIABLES -----------------------------

#1. Create violinplots showing the relation between the weighted Balt score and HPV status.
ggplot(data=subset(ALL, !is.na(Subtype)), aes(x=Subtype, y=balt_weighted, fill=Subtype)) + 
  geom_violin() + 
  scale_fill_manual(values=c("grey65", "slateblue3")) + #Or slateblue.
  labs(x = "HPV status", y = "weighted βalt score") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, fill="black") +
  theme(text = element_text(size=15)) + theme(legend.title=element_text(size=14)) +
  rremove("legend")

#2. Do a U Mann Whitney test to compare the weighted Balt score between HPV+ and HPV- patients.
wilcox.test(ALL$balt_weighted ~ ALL$Subtype, data=ALL) 

#3. Create violinplots showing the relation between "fraction of the genome altered" and weighted Balt tertiles.
ggplot(data=subset(ALL, !(tertile_weighted=="Middle")), 
       aes(x=tertile_weighted, y=Fraction.Genome.Altered, fill=tertile_weighted)) + 
  geom_violin() + 
  scale_fill_manual(values=c("indianred1", "steelblue4")) + 
  labs(x = "weighted βalt group", y = "Fraction of the genome altered") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, fill="black") +
  theme(text = element_text(size=15)) + theme(legend.title=element_text(size=14)) +
  rremove("legend")

#4. Do a U Mann Whitney Test to compare genomic alterations between the weighted BAlt tertiles.
wilcox.test(ALL$Fraction.Genome.Altered ~ ALL$Tertile_weighted, data=ALL) 

#5. Create violinplots showing the relation between "mutation count" and weighted Balt tertiles.
ggplot(data=subset(ALL, !(tertile_weighted=="Middle")), 
       aes(x=tertile_weighted, y=Mutation.Count, fill=tertile_weighted)) + 
  geom_violin() + 
  scale_fill_manual(values=c("indianred1", "steelblue4")) + 
  labs(x = "weighted βalt group", y = "Mutation count") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, fill="black") +
  theme(text = element_text(size=15)) + theme(legend.title=element_text(size=14)) +
  rremove("legend") +
  scale_y_continuous(trans='log2')

#6. Do a U Mann Whitney Test to compare the mutation count between the weighted BAlt tertiles.
wilcox.test(ALL$Mutation.Count ~ ALL$Tertile_weighted, data=ALL) 



# PART 8.1: SURVIVAL CURVES - WEIGHTED βALT -----------------------------

#1. Create a dataframe of patients treated with RT.
ALL_RT <- ALL[which(ALL$radiation_therapy=="YES" | 
                 ALL$additional_radiation_therapy=="YES" | 
                 ALL$postoperative_rx_tx=="YES"),] #520->312 patients.

#2. Create a dataframe excluding patients whose primary curative treatment was surgery.
#The remaining patients probably received RT and/or genotoxic ChT, based on the standard of care. 
ALL_RTChT <- ALL[-which(ALL$method_of_curative_tx=="Surgery"),] #520->419 patients.

#3. Create separate datasets separating patients by their tumor stage. 
Stage1 <- ALL[which(ALL$stage_4groups=="I"),] #27 patients.
Stage2 <- ALL[which(ALL$stage_4groups=="II"),] #71 patients.
Stage3 <- ALL[which(ALL$stage_4groups=="III"),] #81 patients.
Stage4 <- ALL[which(ALL$stage_4groups=="IV"),] #266 patients.
Stage34 <- ALL[which(ALL$stage_4groups=="III" | ALL$stage_4groups=="IV"),] #347 patients.

#4. Select the dataset with patients presumably treated with RT and/or genotoxic ChT. 
X <- ALL_RTChT

#5. Plot and compare the overall survival curves between the weighted Balt score top and bottom tertiles.

##Without the intermediate tertile.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
my.fit<-survfit(my.surv.object~X$Tertile_weighted)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="TCGA-HNSC", subtitle = "All except those with 'curative.treatment=surgery' (n=419)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

##With the intermediate tertile.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
my.fit<-survfit(my.surv.object~X$tertile_weighted)
my.fit
p<-ggsurvplot(my.fit, data=X, pval = TRUE, size = 1.5, censor.size=7, conf.int = FALSE, break.time.by = 12,
              xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),
              legend.labs=c("low βalt", "high βalt", "mid βalt"), legend.title=" ", 
              font.main=15, palette=c("steelblue4", "indianred1", "goldenrod3"),
              risk.table = TRUE)
p$plot <- p$plot + labs(title="TCGA-HNSC", subtitle = "All except those with 'curative.treatment=surgery' (n=419)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#6. Calculate the overall survival hazard ratio between the weighted Balt score top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
cox<-coxph(my.surv.object ~  X$Tertile_weighted)
summary(cox)

#7. Plot the overall survival curves of patients treated with RT.
X <- ALL_RT
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
my.fit<-survfit(my.surv.object~X$Tertile_weighted)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="TCGA-HNSC", subtitle = "All those treated with RT (n=312)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#8. Plot the overall survival curves of patients from each stage. 
X <- ALL
fit <- survfit( Surv(time=OS.months, event=OS) ~ Tertile_weighted, data = X)
ggsurvplot_facet(fit, X, facet.by = "stage",
                 palette=c("indianred1", "steelblue4"), pval = TRUE, break.time.by = 12,
                 xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60))
ggsurvplot_facet(fit, X, facet.by = "Stage",
                 palette=c("indianred1", "steelblue4"), pval = TRUE, break.time.by = 12,
                 xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60))
table(X$Stage) #I: 27, II: 71, III: 81, IV: 266.



# PART 8.2: SURVIVAL CURVES - ORIGINAL βALT AND OTHER VARIATIONS -----------------------------

#1. Select the dataset with patients presumably treated with RT and/or genotoxic ChT. 
X <- ALL_RTChT

#2. Plot and compare the overall survival curves between the original Balt score top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
my.fit<-survfit(my.surv.object~X$Tertile_original)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="TCGA-HNSC: Original  βalt", subtitle = "All except those with 'curative.treatment=surgery' (n=419)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#3. Plot and compare the overall survival curves between the short Balt score top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
my.fit<-survfit(my.surv.object~X$Tertile_short)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="TCGA-HNSC: Short  βalt", subtitle = "All except those with 'curative.treatment=surgery' (n=419)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#4. Plot and compare the overall survival curves between the top Balt score top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
my.fit<-survfit(my.surv.object~X$Tertile_top)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="TCGA-HNSC: Top 15  βalt", subtitle = "All except those with 'curative.treatment=surgery' (n=419)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#5. Plot and compare the overall survival curves between the bottom Balt score top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
my.fit<-survfit(my.surv.object~X$Tertile_bottom)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="TCGA-HNSC: Bottom 15  βalt", subtitle = "All except those with 'curative.treatment=surgery' (n=419)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#6. Calculate the overall survival hazard ratios between the scores top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
##Original Balt. 
cox<-coxph(my.surv.object ~  X$Tertile_original)
summary(cox)
##Short Balt. 
cox<-coxph(my.surv.object ~  X$Tertile_short)
summary(cox)
##Top15 Balt. 
cox<-coxph(my.surv.object ~  X$Tertile_top)
summary(cox)
##Bottom15 Balt. 
cox<-coxph(my.surv.object ~  X$Tertile_bottom)
summary(cox)



# PART 9: COX REGRESSIONS - WEIGHTED BALT -----------------------------------------

#1. Select the dataset with patients presumably treated with RT and/or genotoxic ChT. 
X <- ALL_RTChT

#2. Calculate the univariate Cox regressions of the weighted Balt score. 
##As a categoric variable (tertiles 1 vs 3). 
my.surv.object <- Surv(time=X$OS.months, event=X$OS)
cox<-coxph(my.surv.object ~  X$Tertile_weighted)
summary(cox)
##As a continuous variable. 
cox<-coxph(my.surv.object ~  X$balt_weighted)
summary(cox)

#2. Calculate the multivariate Cox regressions of the weighted Balt score. 
##As a categoric variable (tertiles 1 vs 3).
cox<-coxph(my.surv.object ~  X$Tertile_weighted + X$age + X$stage + X$Subtype)
summary(cox)
##As a continuous variable. 
cox<-coxph(my.surv.object ~  X$balt_weighted + X$age + X$stage + X$Subtype)
summary(cox)




