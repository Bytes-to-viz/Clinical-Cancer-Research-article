#R version 3.6.1 (2019-07-05) -- "Action of the Toes"
#Copyright (C) 2019 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)



#TITLE: ANALYSIS OF THE GSE41613-HNSC DATASET



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



# PART 1: DATA IMPORTATION (CAN SKIP THIS AND JUMP DIRECTLY TO PART 7) -----------------------------------------

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



# PART 2: DOWNLOAD DATA FROM GEO (CAN SKIP THIS AND JUMP DIRECTLY TO PART 7) -----------------------------------------

#1.Download the "GSE41613" series. 
library("GEOquery")
?getGEO()
gse <- getGEO("GSE41613")

#2. Check how many expression objects there are.
length(gse) #1 
gse <- gse[[1]] #Assign it to the same variable.

#3. Characterize briefly what we have.
show(gse)
head(exprs(gse)) #exprs(gse) is a matrix.
length(exprs(gse)) 

#4. Download a matrix with gene expression.
GEOexpression <- exprs(gse)
dim(GEOexpression) #97 patients x 54.613 genes.
colnames(GEOexpression)
GEOexpression[1:5,] #The row names in the expression data are specific IDs not gene names.
featureNames(gse)[1:10]

#5. Change the names of the genes' IDs to the gene names.
annotation(gse) #This is the platform used. GPL570.
platf <- getGEO(annotation(gse), AnnotGPL=TRUE) 
show(platf) #This is the platform metadata. Contains extensive information for each gene (ID, symbol, description, GO code...).
##The attr() function allows us to extract specific data from this table.
attr(dataTable(platf), "table")[1:100,c("ID", "Gene symbol")] #For example, we can get the IDs and symbols from the first 100 genes.
##Store a matrix of IDs.
IDs <- attr(dataTable(platf), "table")[,c("ID", "Gene symbol")]
##Merge the IDs with the expression file. 
GEOexpression2 <- as.data.frame(GEOexpression)
GEOexpression2$OldID <- rownames(GEOexpression2)
GEOexpression2 <- merge(IDs, GEOexpression2, 
                        by.x="ID", by.y="OldID", 
                        all.x=FALSE, all.y=TRUE) 
#Replace the rownames with the symbols.
GEOexpression3 <- as.matrix(GEOexpression2[,!(names(GEOexpression2) %in% c("ID", "Gene symbol"))])
rownames(GEOexpression3) <- GEOexpression2$`Gene symbol`
tail(rownames(GEOexpression3), 10)

#6. Download phenotypic data.
summary(pData(gse))
phenotype <- pData(gse) 



# PART 3: DEAL WITH DUPLICATED AND EMPTY GENES (CAN SKIP THIS AND JUMP DIRECTLY TO PART 7) -----------------------------------------

#1.Create a submatrix that includes only non-empty row names.
##Create a vector showing FALSE or TRUE on whether probes (genes) have IDs.
sel <- apply(as.matrix(IDs[,2], ncol=1), 1,
             function(x){if (x == "") return(FALSE) else return(TRUE)})
##Then select from the matrix.
selExpr <- c()
selSymb <- c()
for (i in 1:nrow(GEOexpression3)) {
  if (sel[i]) {
    selExpr <- rbind(selExpr, GEOexpression3[i,])
    selSymb <- rbind(selSymb, IDs[i,2])
  }
} #WARNING: Quite slow. 15 mins aprox.
rownames(selExpr) <- selSymb
dim(selExpr)  #97 patients.

#2. Check how many genes are duplicated. 
Duplicated <- selExpr[duplicated(rownames(selExpr)),] 
dim(Duplicated) #22.919 duplicated genes.

#3. Calculate the means of probes (rows) with the same gene IDs. 
selExpr2 <- as.data.frame(selExpr)
head(selExpr2)
selExpr2$geneID <- rownames(selExpr)
selExpr3 <- aggregate(.~geneID, selExpr2, FUN=mean) #45.108->22.189 genes. #WARNING: Quite slow.
rownames(selExpr3) <- selExpr3$geneID
selExpr3$geneID <- NULL



# PART 4: NORMALIZE GENE EXPRESSION (CAN SKIP THIS AND JUMP DIRECTLY TO PART 7) -----------------------------------------

#In GEO, when clicking in one of the samples one can see how the dataset was processed: 
##Processing:	We used gcRMA algorithm from Bioconductor to extract gene expression values and perform normalization. 
##Values: log2 gcRMA signal.

#1. Plot the distribution of the expression of a few random genes accross all samples, to see if log2 transformation is needed. 
GEODF <- as.data.frame(t(selExpr3[c("A1BG"),]))
ggplot(GEODF, aes(x=GEODF$A1BG)) + geom_density()
GEODF <- as.data.frame(t(selExpr3[c("TGFB1"),]))
ggplot(GEODF, aes(x=GEODF$TGFB1)) + geom_density()
GEODF <- as.data.frame(t(selExpr3[c("PARP1"),]))
ggplot(GEODF, aes(x=GEODF$PARP1)) + geom_density()
GEODF <- as.data.frame(t(selExpr3[c("LIG1"),]))
ggplot(GEODF, aes(x=GEODF$LIG1)) + geom_density()
GEODF <- as.data.frame(t(selExpr3[c("BRCA1"),]))
ggplot(GEODF, aes(x=GEODF$BRCA1)) + geom_density()
##CONCLUSION: The values have not been scaled but have been log2 transformed.

#2. Convert gene expression values into gene-centered Z-scores ((x - mean(Xcolumn)) / sd(Xcolumn)).
selExpr5 <- as.matrix(selExpr3) #Turn the dataset with gene expression into a matrix.
genesT <- t(selExpr5) #Transpose the matrix so that each gene is a column.
genesT<-scale(genesT, center = TRUE, scale = TRUE) ##Calculate the Z-scores.



# PART 5: SCORES CALCULATION - WEIGHTED AND ORIGINAL BALT (CAN SKIP THIS AND JUMP DIRECTLY TO PART 7) -----------------------------------------

#1. Check how many TGFB and ALTEJ genes are available. 
A <- selExpr5[which(rownames(selExpr5) %in% TGFBlist), ] #49/50 genes.
A <- selExpr5[which(rownames(selExpr5) %in% ALTEJlist), ] #36/36 genes.

#2. Select only the TGFB and ALTEJ genes. 
genes <- as.data.frame(t(genesT))
genes <- genes[which(rownames(genes) %in% TGFBlist | rownames(genes) %in% ALTEJlist), ] #85 genes.

#3. Edit gene names so that they match the ones from the genes'weight dataset. 
genes <- as.data.frame(genes)
genes$Gene <- rownames(genes)
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
                by.x="Gene", by.y="Gene", all.x=TRUE, all.y=TRUE) #86 genes, but RUNX1 with null values.
rownames(genes2) <- genes2$Gene

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
         method = "pearson") 

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

#8. Calculate also the TGFB and altEJ ssgsea scores and the original Balt. 
genes<-t(genesT)
genes <- as.data.frame(genes)
genes <- as.matrix(genes)
ssgsea <- gsva(genes, Bothgenelists, method=c("ssgsea"))
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
cor.test(ssgsea$TGFBUPgeneset, ssgsea$ALTEJgeneset) 
ssgsea$balt <- sqrt((max(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (min(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2) - 
               sqrt((min(ssgsea$ALTEJgeneset)-ssgsea$ALTEJgeneset)^2+
                      (max(ssgsea$TGFBUPgeneset)-ssgsea$TGFBUPgeneset)^2)
ssgsea$balt <- ssgsea$balt * -1 
ssgsea$balt_original <- ssgsea$balt
ssgsea_original <- ssgsea
ssgsea_original$tertile_original <- cut(ssgsea_original$balt_original, 
                                        quantile(ssgsea_original$balt_original, c(0, 1/3, 2/3, 1)), 
                                        na.rm=TRUE, include.lowest=TRUE, labels = c("Low", "Middle", "High"))
ssgsea_original$Tertile_original <- ifelse(ssgsea_original$tertile_original == "Low", "high TGFβ and low ALTEJ", 
                                           ifelse(ssgsea_original$tertile_original == "High", "low TGFβ and high ALTEJ", NA))
ssgsea_original$Tertile_original <- as.character(ssgsea_original$Tertile_original)



# PART 6: DATA WRANGLING (CAN SKIP THIS AND JUMP DIRECTLY TO PART 7) -----------------------------------

#1. Merge the "scores" and "phenotype" dataframes to put everything in one dataset.
scores$ID <- rownames(scores) 
phenotype$ID <- rownames(phenotype)
GEOALL <- merge(scores, phenotype, 
                by.x = "ID", by.y = "ID",
                all.x=FALSE, all.y=FALSE)
dim(GEOALL) #97 patients x 50 variables.
ssgsea_original$ID <- rownames(ssgsea_original) 
GEOALL <- merge(GEOALL, ssgsea_original[,c("ID", "balt_original", "Tertile_original", "tertile_original")], 
                by.x = "ID", by.y = "ID",
                all.x=FALSE, all.y=FALSE)
dim(GEOALL) #97 patients x 53 variables.

#2. Explore what clinical information is available. 
names(GEOALL)

#3. Prepare the survival variables so that they are in the appropriate format.  Note:FU time means follow-up time.
GEOALL$OS.time <- GEOALL$`fu time:ch1`
table(GEOALL$OS.time)
GEOALL$OS.months <- as.numeric(GEOALL$OS.time)
table(GEOALL$`vital:ch1`)
GEOALL$OS <- ifelse(GEOALL$`vital:ch1`=="Alive", 0,
                           ifelse(GEOALL$`vital:ch1`=="Dead-non OC", 1, 
                           ifelse(GEOALL$`vital:ch1`=="Dead-oral ca", 1,    
                           ifelse(GEOALL$`vital:ch1`=="Dead -unk cause", 1, NA))))
GEOALL$specific.OS <- ifelse(GEOALL$`vital:ch1`=="Alive", 0,
                             ifelse(GEOALL$`vital:ch1`=="Dead-non OC", 0, 
                             ifelse(GEOALL$`vital:ch1`=="Dead-oral ca", 1,    
                             ifelse(GEOALL$`vital:ch1`=="Dead -unk cause", 0, NA))))

#4. Check stage. 
table(GEOALL$`tumor stage:ch1`) #I-II: 41, III-IV: 56.

#5. Check treatment. 
table(GEOALL$`treatment:ch1`) #53 multimodal, 43 unimudal, unknown 1. 
table(GEOALL$`treatment:ch1`, GEOALL$`tumor stage:ch1`)

#6. Tidy other variables. 
GEOALL$balt_weighted <- GEOALL$balt
GEOALL$TGFBUPgeneset_weighted <- GEOALL$TGFBUPgeneset
GEOALL$ALTEJgeneset_weighted <- GEOALL$ALTEJgeneset
GEOALL$tertile_weighted <- GEOALL$tertile
GEOALL$Tertile_weighted <- GEOALL$Tertile
GEOALL$balt <- NULL
GEOALL$TGFBUPgeneset <- NULL
GEOALL$ALTEJgeneset <- NULL
GEOALL$tertile <- NULL
GEOALL$Tertile <- NULL
GEOALL$stage <- GEOALL$`tumor stage:ch1`

#7. Export the "GEOALL" dataset. 
write_xlsx(GEOALL, "Output/balt_scores_GSE41613HNSC.xlsx")



# PART 7: SCATTERPLOT OF THE CORRELATION BETWEEN TGFβ AND ALTEJ (CAN START HERE DIRECTLY) -----------------------

#1. Import the "GEOALL" dataset. Instead of running parts 1-6, can upload this and start here directly. 
GEOALL <- read_excel("Output/balt_scores_GSE41613HNSC.xlsx")

#2. Scatterplot of TGFB versus ALTEJ weighted scores, colored by weighted Balt tertiles. 
p <- ggplot(GEOALL, aes(x=TGFBUPgeneset_weighted, y=ALTEJgeneset_weighted)) + 
  geom_point(size=5, alpha=0.7,shape=21, stroke=0.1, color="white", aes(fill=tertile_weighted)) +
  scale_fill_manual(values=c("steelblue4", "indianred1", "bisque3")) + #Or "grey50"
  labs(x = "TGFβ score", y = "alt-EJ score", col="βalt tertile") +
  geom_smooth(method=loess, fullrange=TRUE) + 
  #geom_smooth(method="glm", fullrange=TRUE) +
  stat_cor(method="pearson", size=5) +
  theme(text = element_text(size=15)) +
  theme(legend.title=element_text(size=14))
p <- p + rremove("legend")
p1 <- ggMarginal(p,  size=10, fill = "grey65"); p1 #PDF 5x6.5

#3. Calculate Pearson and Spearman correlation coefficients.
cor.test(GEOALL$TGFBUPgeneset_weighted, GEOALL$ALTEJgeneset_weighted, method = "pearson")
cor.test(GEOALL$TGFBUPgeneset_weighted, GEOALL$ALTEJgeneset_weighted, method = "spearman")



# PART 8.1: SURVIVAL CURVES - WEIGHTED βALT -----------------------------------------------

#1. Plot and compare the specific overall survival curves between the weighted Balt score top and bottom tertiles.
X <- GEOALL

##Without the intermediate tertile.
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
my.fit<-survfit(my.surv.object~X$Tertile_weighted)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Specific Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="GSE41613-HNSC", subtitle = "All (n=97)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

##With the intermediate tertile.
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
my.fit<-survfit(my.surv.object~X$tertile_weighted)
my.fit
p<-ggsurvplot(my.fit, data=X, pval = TRUE, size = 1.5, censor.size=7, conf.int = FALSE, break.time.by = 12,
              xlab="Time (months)", ylab="Specific Surviving Fraction", xlim=c(0, 60),
              legend.labs=c("low βalt", "high βalt", "mid βalt"), legend.title=" ", 
              font.main=15, palette=c("steelblue4", "indianred1", "goldenrod3"),
              risk.table = TRUE)
p$plot <- p$plot + labs(title="GSE41613-HNSC", subtitle = "All (n=97)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#2. Calculate the specific overall survival hazard ratio between the weighted Balt score top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
cox<-coxph(my.surv.object ~  X$Tertile_weighted)
summary(cox)

#3. Plot the specific overall survival curves of patients from each stage. 
fit <- survfit( Surv(time=OS.months, event=specific.OS) ~ Tertile_weighted, data = X)
ggsurvplot_facet(fit,data=X, facet.by = "stage",  
                 palette=c("indianred1", "steelblue4"), break.time.by = 12,
                 #pval = TRUE, #When adding "pval = TRUE" it stops working
                 xlab="Time (months)", ylab="Specific Surviving Fraction", xlim=c(0, 60))
table(X$stage) #I-II: 41, III-IV: 56. 

#4. Since p-values could not be added to the former plot, calculate them. 
X <- GEOALL[which(GEOALL$stage=="I/II"),]
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
cox<-coxph(my.surv.object ~  X$Tertile_weighted)
summary(cox) #p=0.1.
X <- GEOALL[which(GEOALL$stage=="III/IV"),]
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
cox<-coxph(my.surv.object ~  X$Tertile_weighted)
summary(cox) #p=0.2.



# PART 8.2: SURVIVAL CURVES - ORIGINAL βALT -----------------------------------------------

#1. Plot and compare the specific overall survival curves between the original Balt score top and bottom tertiles.
X <- GEOALL

##Without the intermediate tertile.
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
my.fit<-survfit(my.surv.object~X$Tertile_original)
my.fit
p<- ggsurvplot(my.fit, data=X, pval = TRUE, conf.int = TRUE, break.time.by = 12,
               xlab="Time (months)", ylab="Specific Surviving Fraction", xlim=c(0, 60),
               legend.labs=c("low βalt", "high βalt"), legend.title=" ", 
               font.main=15, palette=c("indianred1", "steelblue4"),
               #surv.median.line = "hv",
               risk.table = TRUE) 
p$plot <- p$plot + labs(title="GSE41613-HNSC", subtitle = "All (n=97)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

##With the intermediate tertile.
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
my.fit<-survfit(my.surv.object~X$tertile_original)
my.fit
p<-ggsurvplot(my.fit, data=X, pval = TRUE, size = 1.5, censor.size=7, conf.int = FALSE, break.time.by = 12,
              xlab="Time (months)", ylab="Specific Surviving Fraction", xlim=c(0, 60),
              legend.labs=c("low βalt", "high βalt", "mid βalt"), legend.title=" ", 
              font.main=15, palette=c("steelblue4", "indianred1", "goldenrod3"),
              risk.table = TRUE)
p$plot <- p$plot + labs(title="GSE41613-HNSC", subtitle = "All (n=97)") + 
  theme(plot.title=element_text(hjust = 0.5), plot.subtitle=element_text(hjust = 0.5))
p1<-p;p1 #PDF 8x7

#2. Calculate the specific overall survival hazard ratio between the original Balt score top and bottom tertiles.
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
cox<-coxph(my.surv.object ~  X$Tertile_original)
summary(cox)



# PART 9: COX REGRESSIONS - WEIGHTED BALT -----------------------------------------------

#1. Calculate the univariate Cox regressions of the weighted Balt score.
X <- GEOALL
##As a categoric variable (tertiles 1 vs 3). 
my.surv.object <- Surv(time=X$OS.months, event=X$specific.OS)
cox<-coxph(my.surv.object ~  X$Tertile_weighted)
summary(cox)
##As a continuous variable. 
cox<-coxph(my.surv.object ~  X$balt_weighted)
summary(cox)

#2. Calculate the multivariate Cox regressions of the weighted Balt score. 
##As a categoric variable (tertiles 1 vs 3).
cox<-coxph(my.surv.object ~  X$Tertile_weighted + X$`age:ch1` + X$stage)
summary(cox)
##As a continuous variable. 
cox<-coxph(my.surv.object ~  X$balt_weighted + X$`age:ch1` + X$stage)
summary(cox)



