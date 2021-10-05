#R version 3.6.1 (2019-07-05) -- "Action of the Toes"
#Copyright (C) 2019 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)



#TITLE: ANALYSIS OF THE NCI60 PANCANCER CELL LINES DATASET



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



# PART 1: DATA IMPORTATION -----------------------------------------

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

#3. Import cell lines files. 
Exp <- read_delim("Input/data_mRNA_median_Zscores.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Clin <- read_delim("Input/cellline_nci60_clinical_data.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
##These files were downloaded from cBioPortal in April 2021.  
SF2 <- read_excel("Input/SF2.xlsx")
##This file contains the surviving fraction after 2 Gy (SF2) of the NCI60 cell lines based on clonogenic assays.
##It was obtained from the article with doi: 10.1016/j.ijrobp.2009.05.056. 



# PART 2: SCORES CALCULATION - WEIGHTED BALT -----------------------------------------

#1. Check how many TGFB and ALTEJ genes are available. 
A <- Exp[which(Exp$Hugo_Symbol %in% TGFBlist), ] #50/50 genes.
A <- Exp[which(Exp$Hugo_Symbol %in% ALTEJlist), ] #36/36 genes.

#2. Select only the TGFB and ALTEJ genes. 
genes <- Exp[which(Exp$Hugo_Symbol %in% TGFBlist | Exp$Hugo_Symbol %in% ALTEJlist), ] #86 genes.

#3. Edit gene names so that they match the ones from the genes'weight dataset. 
genes <- as.data.frame(genes)
genes$Gene <- genes$Hugo_Symbol
genes$Gene <- ifelse(genes$Gene=="APEX2", "APE2", 
                     ifelse(genes$Gene=="HRAS", "HRAS1",
                     ifelse(genes$Gene=="MRE11", "MRE11A", 
                     ifelse(genes$Gene=="SPDL1", "CCDC99",
                     ifelse(genes$Gene=="PRPF19", "PRP19",
                     ifelse(genes$Gene=="RAD51L3", "RAD51D",
                     ifelse(genes$Gene=="KAT5", "TIP60",
                     ifelse(genes$Gene=="ARHGAP32", "RICS",
                     ifelse(genes$Gene=="C19orf40", "FAAP24",
                     ifelse(genes$Gene=="CYTH1", "PSCD1",
                     ifelse(genes$Gene=="NUDT1", "MTH1",
                     ifelse(genes$Gene=="OBFC2B", "NABP2", 
                     ifelse(genes$Gene=="PMEPAI", "TMEPAI",      
                     ifelse(genes$Gene=="STAG1", "TMEPAI", genes$Gene))))))))))))))

#4. Merge the genes and their weights in one dataset. 
genesimportance <- rbind(TGFBimportance[,c("Gene", "mean weight")], ALTEJimportance[,c("Gene", "mean weight")])
genes2 <- merge(genesimportance, genes, 
                by.x="Gene", by.y="Gene", all.x=TRUE, all.y=TRUE) #86 genes.
rownames(genes2) <- genes2$Gene
genes2$Hugo_Symbol <- NULL
genes2$Entrez_Gene_Id <- NULL

#5. Calculate the TGFB and ALTEJ weighted scores in each sample, by multiplying each gene by its factor.
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
         method = "pearson") #PCC=-0.53.

#6. Create a new variable that is the weighted Balt score.
scores$balt <- sqrt((max(scores$ALTEJgeneset)-scores$ALTEJgeneset)^2+
                      (min(scores$TGFBUPgeneset)-scores$TGFBUPgeneset)^2) - 
               sqrt((min(scores$ALTEJgeneset)-scores$ALTEJgeneset)^2+
                      (max(scores$TGFBUPgeneset)-scores$TGFBUPgeneset)^2)
scores$balt <- scores$balt * -1 

#7. Create a new variable that are the weighted BAlt score tertiles. 
scores$tertile <- cut(scores$balt, quantile(scores$balt, c(0, 1/3, 2/3, 1)), na.rm=TRUE, include.lowest=TRUE, 
                      labels = c("Low", "Middle", "High"))
scores$Tertile <- ifelse(scores$tertile == "Low", "high TGFβ and low ALTEJ", 
                         ifelse(scores$tertile == "High", "low TGFβ and high ALTEJ", NA))
scores$Tertile <- as.character(scores$Tertile)



# PART 3: DATA WRANGLING -----------------------------------

#1. Merge the "scores" and the clinical dataframes to put everything in one dataset. 
table(Clin$`Patient ID`)
SF2$ID <- ifelse(SF2$`Cell Line`=="BREAST_MCF7ADRr", "NCI_ADR_RES", SF2$ID)
Clin2 <- merge(SF2, Clin,
               by.x="ID", by.y="Patient ID", 
               all.x=FALSE, all.y=FALSE) #60 cell lines.
table(Clin2$`Cancer Type`) #There are 6 cell lines from blood cancers.
scores$ID <- rownames(scores)
scores$ID <- ifelse(scores$ID=="CSF_268", "SF_268", 
                    ifelse(scores$ID=="CSF_295", "SF_295",
                    ifelse(scores$ID=="CSF_539", "SF_539",
                    ifelse(scores$ID=="CSNB_19", "SNB_19",
                    ifelse(scores$ID=="CSNB_75", "SNB_75",
                    ifelse(scores$ID=="CU251", "U251", scores$ID))))))
ALL <- merge(scores, Clin2, 
                 by.x="ID", by.y="ID", 
                 all.x=FALSE, all.y=FALSE) #60+60->60 cell lines.

#2. Tidy some vriable names. 
ALL$TGFBUPgeneset_weighted <- ALL$TGFBUPgeneset
ALL$ALTEJgeneset_weighted <- ALL$ALTEJgeneset
ALL$balt_weighted <- ALL$balt
ALL$tertile_weighted <- ALL$tertile
ALL$Tertile_weighted <- ALL$Tertile
ALL$TGFBUPgeneset <- NULL
ALL$ALTEJgeneset <- NULL
ALL$balt <- NULL
ALL$tertile <- NULL
ALL$Tertile <- NULL

#3. Export the "ALL" dataset. 
write_xlsx(ALL, "Output/balt_scores_NCI60.xlsx")



# PART 4: SCATTERPLOTS OF THE CORRELATION BETWEEN βALT AND SF2 -----------------------

#1. Import the "ALL" dataset. Instead of running parts 1-3, can upload this and start here directly. 
ALL <- read_excel("Output/balt_scores_NCI60.xlsx")

#2. Scatterplot of TGFB versus ALTEJ weighted scores.
##Colored by cancer type.
p <- ggplot(ALL, aes(x=TGFBUPgeneset_weighted, y=ALTEJgeneset_weighted)) + 
  geom_point(size=5, alpha=0.7,  shape=21, stroke=0.1, color="white", aes(fill=ALL$`Cancer Type`)) +
  labs(x = "TGFβ score", y = "alt-EJ score", col="Cancer type") +
  geom_smooth(method="glm", fullrange=TRUE) +
  stat_cor(method="spearman", size=5) +
  theme(text = element_text(size=15)) +
  theme(legend.title=element_text(size=14)); p
p <- p + rremove("legend")
p1 <- ggMarginal(p,  size=10, fill = "grey65"); p1 #PDF 5x6.5.
##Not colored.
p <- ggplot(ALL, aes(x=TGFBUPgeneset_weighted, y=ALTEJgeneset_weighted)) + 
  geom_point(size=5, color="grey40") +
  labs(x = "TGFβ score", y = "alt-EJ score") +
  geom_smooth(method="glm", fullrange=TRUE) +
  stat_cor(method="spearman", size=5) +
  theme(text = element_text(size=15)); p
p1 <- ggMarginal(p,  size=10, fill = "grey65"); p1 #PDF 5x6.5.

#3. Calculate Pearson and Spearman correlation coefficients.
cor.test(ALL$TGFBUPgeneset_weighted, ALL$ALTEJgeneset_weighted, method = "pearson")
cor.test(ALL$TGFBUPgeneset_weighted, ALL$ALTEJgeneset_weighted, method = "spearman")

#4. Scatterplot of Balt versus SF2.
##Colored by cancer type.
p <- ggplot(ALL, aes(x=balt_weighted, y=`Recorded SF2`)) + 
  geom_point(size=5, alpha=0.7,  shape=21, stroke=0.1, color="white", aes(fill=ALL$`Cancer Type`)) +
  labs(x = c(expression(βalt[w]~score)), y = "SF2") +
  stat_cor(method="spearman", size=5) +
  theme(text = element_text(size=15)) +
  theme(legend.title=element_text(size=14)); p
p <- p + rremove("legend")
p1 <- ggMarginal(p,  size=10, fill = "grey65"); p1 #PDF 5x6.5.
##Not colored.
p <- ggplot(ALL, aes(x=balt_weighted, y=`Recorded SF2`)) + 
  geom_point(size=5, color="grey40") +
  labs(x = c(expression(βalt[w]~score)), y = "SF2") +
  geom_smooth(method="glm", fullrange=TRUE) +
  stat_cor(method="spearman", size=5) +
  theme(text = element_text(size=15)); p
p1 <- ggMarginal(p,  size=10, fill = "grey65"); p1 #PDF 5x6.5.

#5. Calculate Pearson and Spearman correlation coefficients.
cor.test(ALL$balt_weighted, ALL$`Recorded SF2`, method = "pearson")
cor.test(ALL$balt_weighted, ALL$`Recorded SF2`, method = "spearman")



