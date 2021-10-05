#R version 3.6.1 (2019-07-05) -- "Action of the Toes"
#Copyright (C) 2019 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)



#TITLE: ANALYSIS OF THE "INHOUSE OV DATASET"



# PART 0: PREPARE THE ENVIRONMENT ---------------------------------------------

#1. Load all necessary packages.
library("readxl")
library("writexl")
library("xlsx")
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



# PART 1: IMPORT THE FILES FOR THE ANALYSIS ---------------------------------------------

#1. Import the genelists. 
TGFBgeneset <- read.table("Input/TGFB_list.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
ALTEJgeneset <- read.table("Input/ALTEJ_list.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
TGFB_list <- as.list(TGFBgeneset[,1:1]) #55 genes (50 + synonyms)
ALTEJ_list <- as.list(ALTEJgeneset[,1:1]) #45 genes (36 + synonyms)
Bothgenelists <- list(TGFBUPgeneset=TGFB_list, ALTEJgeneset=ALTEJ_list)

#2. Import the file with gene expression.
Exp_original <- read.delim("Input/Expression_Nsolver_normalized.txt")
##This file contatins normalized gene expression of 48 OV cancer samples (23 FFPE + 25 snap-frozen) from 41 patients. 
##Normalization (as described in Methods): Processed in the Nsolver (background threshold + positive control normalization + normalization to housekeeping genes).
##This data still needs to be Log2 transformed -> Z-scored transformed grouping FFPE and snap-frozen samples separately. 



# PART 2: DATA WRANGLING ---------------------------------------------

#1. Eliminate unnecessary rows and columns from the dataframe.
Exp <- Exp_original[-2,]
Exp <- Exp[which(Exp$Class.Name  == "Endogenous" | Exp$Class.Name == "Housekeeping" | Exp$Class.Name == ""),]
Exp$Class.Name <- NULL
Exp$Annotation <- NULL
Exp$Accession.. <- NULL
Exp$NS.Probe.ID <- NULL
Exp$Analyte.Type <- NULL
Exp$X..Samples.above.Threshold <- NULL
Exp$Positive.Flag <- NULL
Exp$Avg.Count <- NULL
Exp$Species.Name <- NULL
rownames(Exp) <- Exp$Probe.Name
Exp$Probe.Name <- NULL

#2. Extract the sample IDs from the column names. 
Exp <- as.data.frame(t(Exp))
Exp$sample_type <- Exp$V1
Exp$V1 <- NULL
Exp$sample_type <- ifelse(Exp$sample_type=="UCSF", "frozen", "FFPE")
Exp$SampleID <- rownames(Exp)
Exp$SampleID_FFPE <- ifelse(Exp$sample_type == "FFPE", substring(Exp$SampleID,1,nchar(Exp$SampleID)-7), NA) #Exclude the last 7 characters.
Exp$SampleID_FFPE <-substring(Exp$SampleID_FFPE,nchar(Exp$SampleID_FFPE)-5, nchar(Exp$SampleID_FFPE)) #Select the last 5 characters.
Exp$SampleID_FFPE <- chartr(".", "-", Exp$SampleID_FFPE) #Turn the "."s into "-"s.
Exp$SampleID_FFPE <- sub("_", "", Exp$SampleID_FFPE) #Delete the "_"s characters.
Exp$SampleID_frozen <- ifelse(Exp$sample_type == "frozen", substring(Exp$SampleID,1,nchar(Exp$SampleID)-7), NA) #Exclude the last 7 characters.
Exp$SampleID_frozen <-substring(Exp$SampleID_frozen,21, nchar(Exp$SampleID_frozen)) #Remove the first 20 characters.
Exp$SampleID_frozen <- chartr(".", "-", Exp$SampleID_frozen) #Turn the "."s into "-"s.
Exp$SampleID_frozen <- ifelse(Exp$SampleID_frozen=="OV-5", "OV-2-2", Exp$SampleID_frozen) #Relabel this sample as explained in email Alexander Cheung from the NanoString core 20 July 2020, in which he specified that sample OV-5 was empty and he run OV-2 twice, once within its own tube and again by substituting it in place of OV-5. 
Exp$SampleID_frozen <- gsub("OV-", "OV1-", Exp$SampleID_frozen) #Turn "OV-" into "OV1-".
Exp$SampleID <- ifelse(Exp$sample_type == "FFPE", Exp$SampleID_FFPE, Exp$SampleID_frozen)
rownames(Exp) <- Exp$SampleID

#3. Create a list with the FFPE and snap-frozen samples. 
FFPE_list <- Exp[which(Exp$sample_type=="FFPE"),]
FFPE_list <- as.list(as.character(FFPE_list$SampleID))
frozen_list <- Exp[which(Exp$sample_type=="frozen"),]
frozen_list <- as.list(as.character(frozen_list$SampleID))

#4. Eliminate unnecessary rows and columns.
Exp$SampleID <- NULL
Exp$SampleID_FFPE <- NULL
Exp$SampleID_frozen <- NULL    
Exp$sample_type <- NULL 
Exp <- as.data.frame(t(Exp))
dim(Exp)  #150 genes, 48 samples.

#5. Turn the variables into numeric. 
ExpB <- as.data.frame(sapply(Exp, function(x) as.numeric(chartr(",", ".", as.character(x))))) #Turn values into character --> Turn the ","s into "."s --> Turn values into numeric.
rownames(ExpB) <- rownames(Exp)
Exp <- ExpB

#6. Exclusion of the sample "OV2-14" because it had a normalization flag in the nSolver, indicating that it did not pass the quality control.
Exp <- Exp[,-which(colnames(Exp)=="OV2-14")]

#7. Create a list with the housekeeping genes. 
HK_list <- Exp_original[which(Exp_original$Class.Name=="Housekeeping"),]
HK_list <- as.list(as.character(HK_list$Probe.Name))



# PART 3: RIDGELINE PLOTS ---------------------------------------------

#1. Turn the dataframe of gene expression into a long formatted dataframe of few columns.
Exp1 <- as.matrix(Exp) 
Exp2 <- data.frame(Gene=rownames(Exp1)[row(Exp1)],
                   Nanostring_ID=colnames(Exp1)[col(Exp1)],
                   Expression=c(Exp1))

#2. Add a "signature" column to the dataframe.
Signature <- as.data.frame(Exp2[,c("Gene")])
colnames(Signature) <- "Gene"
Signature$Signature <- ifelse(Signature$Gene %in% TGFB_list, "TGFB", 
                              ifelse(Signature$Gene %in% ALTEJ_list, "ALTEJ", 
                              ifelse(Signature$Gene %in% HK_list, "Housekeeping", "Other")))
Signature <- Signature[!duplicated(Signature$Gene), ] #Remove duplicated genes.
Exp3 <- merge(Exp2, Signature, by.x="Gene", by.y="Gene", all.x=TRUE, all.y=FALSE)
Exp3$Expression <- as.numeric(as.character(Exp3$Expression)) #Make sure that this variable is numeric.

#3. Add a "sample type" column to the dataframe.
sample_type <- as.data.frame(Exp3[,c("Nanostring_ID")])
colnames(sample_type) <- "Nanostring_ID"
sample_type$sample_type <- ifelse(sample_type$Nanostring_ID %in% FFPE_list, "FFPE", 
                                 ifelse(sample_type$Nanostring_ID %in% frozen_list, "frozen", NA))
sample_type <- sample_type[!duplicated(sample_type$Nanostring_ID), ] #Remove duplicated samples.
Exp4 <- merge(Exp3, sample_type, by.x="Nanostring_ID", by.y="Nanostring_ID", all.x=TRUE, all.y=FALSE)

#4. Create ridgeline plots of the genes from each signature. 
TGFB <- Exp4[which(Exp4$Signature == "TGFB"),] 
p1<-ggplot(TGFB, aes(x = Expression, y=Gene, color=sample_type, fill=sample_type)) + ggtitle("TGFβ genes") +
  geom_density_ridges(alpha=0.3) +rremove("legend"); p1
ALTEJ <- Exp4[which(Exp4$Signature == "ALTEJ"),] 
p2<-ggplot(ALTEJ, aes(x = Expression, y=Gene, color=sample_type, fill=sample_type)) + ggtitle("Alt-EJ genes") +
  geom_density_ridges(alpha=0.3) +rremove("legend"); p2
HK <- Exp4[which(Exp4$Signature == "Housekeeping"),] 
p3<-ggplot(HK, aes(x = Expression, y=Gene, color=sample_type, fill=sample_type)) + ggtitle("Housekeeping genes") +
  geom_density_ridges(alpha=0.3) +rremove("legend"); p3
Other <- Exp4[which(Exp4$Signature == "Other"),] 
p4<-ggplot(Other, aes(x = Expression, y=Gene, color=sample_type, fill=sample_type)) + ggtitle("Other genes") +
  geom_density_ridges(alpha=0.3) +rremove("legend"); p4
##Together
grid.arrange(p1, p2, p3, p4, nrow = 1)



# PART 4: DATA NORMALIZATION ---------------------------------------------

#1. Log2 transform the data. 
Exp_log2 <- log2(Exp)

#2. Repeat all the steps of "PART 3" but using "Exp_log2" instead of "Exp" in the first line of code.

#3. Transform the values into gene-centered z-scores ((x - mean(Xcolumn)) / sd(Xcolumn)) grouping FFPE and snap-frozen samples separately. 
genesT <- t(Exp_log2) #transpose the dataframe so that EACH GENE IS A COLUMN.
genesT <- as.data.frame(genesT)
genesT$sample_type <- ifelse(rownames(genesT) %in% FFPE_list, "FFPE", 
                            ifelse(rownames(genesT) %in% frozen_list, "frozen", NA)) #add a "sample_type" column.
genesT2 <- genesT %>% 
  group_by(sample_type) %>% 
  mutate_at(vars(1:150), scale)
rownames(genesT2) <- rownames(genesT)
genesT2$sample_type <- NULL
rownames(genesT2) <- rownames(genesT)
Exp_Normalized <- as.data.frame(t(genesT2))

#4. Repeat all the steps of "3.A" but using "Exp_Normalized" instead of "Exp" in the first line of code.

#5. Export the results. 
write.xlsx(Exp_Normalized, file="Output/Expression_OV_specimens.xlsx", col.names = TRUE, row.names = TRUE)

#6. Transpose the dataframe.
Exp_Normalized <- as.data.frame(t(Exp_Normalized))



# PART 5: CLUSTERING HEATMAP ---------------------------------------------

#1. Create a matrix with TGFB and ALTEJ genes' expression values. 
Miniexp <- Exp_Normalized[,which(colnames(Exp_Normalized) %in% TGFB_list | colnames(Exp_Normalized) %in% ALTEJ_list)] 
Miniexp <- as.data.frame(Miniexp)
Miniexp <- as.matrix(Miniexp)
Miniexp <- t(Miniexp) #Transpose it so that genes are rows. 

#2. Create row annotations: Signature. 
Signature <- data.frame(Gene=rownames(Miniexp))
Signature$signature <- ifelse(Signature$Gene %in% TGFB_list, "TGFB", "ALTEJ")
rownames(Signature) <- Signature$Gene
Signature$Gene <- NULL
SignatureColors=list(signature=c("TGFB"="hotpink2", "ALTEJ"="seagreen3"))
Signature_Annotation = HeatmapAnnotation(df  = Signature, col = SignatureColors, which = "row") #Create the heatmap annotation.

#3. Create column annotations: Sample type + Replicated.
Sampleinfo <- data.frame(Sample=colnames(Miniexp))
rownames(Sampleinfo) <- Sampleinfo$Sample
Sampleinfo$sample_type <- ifelse(Sampleinfo$Sample %in% FFPE_list, "FFPE", 
                             ifelse(Sampleinfo$Sample %in% frozen_list, "snap-frozen", NA))
Replicate_list <- as.list(c("OV1-2-2", "OV1-2", "OV1-10-2", "OV1-10", "OV1-5-2", "OV1-5", 
                            "OV1-6-2", "OV1-6", "OV1-7-2", "OV1-7", "OV1-8-2", "OV1-8", 
                            "OV1-9-2", "OV1-9"))
Sampleinfo$Replicated_samples <- ifelse(Sampleinfo$Sample %in% Replicate_list, "Replicated samples", "Non-replicated samples")
Sampleinfo$Sample <- NULL
Sampleinfo <- Sampleinfo[match(colnames(Miniexp), rownames(Sampleinfo)), ] #Order samples as in the expression file.
summary(Sampleinfo)
SampleColors=list(Replicated_samples=c("Non-replicated samples"="grey60", "Replicated samples"="red"), 
                   sample_type=c("FFPE"="gold", "snap-frozen"="deepskyblue")) #Choose colors for the annotations.
Sample_Annotation = HeatmapAnnotation(df = Sampleinfo, col = SampleColors, which = "column") #Create the heatmap annotation.

#4.  Create a hierarchical clustering heatmap. 
Heatmap = Heatmap(Miniexp,
                  clustering_distance_rows = "euclidean", 
                  clustering_method_rows = "ward.D2",
                  row_dend_width = unit(1.5, "cm"),
                  clustering_distance_columns = "euclidean", 
                  clustering_method_columns = "ward.D2",
                  column_dend_height = unit(1.5, "cm"),
                  row_names_gp = gpar(fontsize = 6.2),
                  column_names_gp = gpar(fontsize = 10), 
                  column_km = 2,
                  row_km = 2,
                  bottom_annotation = Sample_Annotation,
                  colorRamp2(c(4, 1.2, 0.5, 0, -0.5, -1.2, -4), brewer.pal(7,"YlGnBu"))); Heatmap
                  #colorRamp2(c(4, 1.1, 0.5, 0, -0.5, -1.1, -4), brewer.pal(7,"RdBu"))); Heatmap
Heatmap + Signature_Annotation + rowAnnotation(rn = anno_text(rownames(Miniexp), gp=gpar(fontsize=6.2))) #PDF9x11.1.



# PART 6: WEIGHTED GENE COEXPRESSION NETWORK ---------------------------------------------

#1. Create a dataset with only the expression values of the TGFB/ALTEJ genes. Rows=samples, columns=genes.
Exp1 <- Exp_Normalized #150 genes.
Exp2 <- data.frame(sapply(Exp1, function(x) as.numeric(as.character(x)))) #Convert the values into numeric.
rownames(Exp2) <- rownames(Exp1)
Exp2 <- Exp2[,which(colnames(Exp2) %in% TGFB_list | colnames(Exp2) %in% ALTEJ_list)] #NOTE: CAN CHANGE TO ONLY TGFB OR ALTEJ GENES.

#2. Exclude the duplicated samples that are replicates from other ones.
Duplicated_list <- as.list(c("OV1-2-2", "OV1-10-2", "OV1-5-2", "OV1-6-2", "OV1-7-2", "OV1-8-2", "OV1-9-2"))
Exp2 <- Exp2[-which(rownames(Exp1) %in% Duplicated_list),] #Remove duplicated samples.

#3. Create a matrix with the Pearson correlation coefficient between the expression of each pair of genes.
GeneCorr <- Exp2 %>% correlate() %>% stretch()
cor.test(Exp2$ABCG1, Exp2$AMIGO2, method="pearson") #Check that the results are correct.

#4.Create a tbl_graph object.
##Create a dataframe with the Edges.
Edges <- GeneCorr %>% filter(r > 0.007) 
Edges$class <- "Mixed edge"
Edges <- within(Edges, class[Edges$x %in% TGFB_list & Edges$y %in% TGFB_list] <- "TGFB edge")
Edges <- within(Edges, class[Edges$x %in% ALTEJ_list & Edges$y %in% ALTEJ_list] <- "ALTEJ edge")
##Create a dataframe with the Nodes.
Nodes <- data.frame(Gene=colnames(Exp2))
Nodes$signature <- ifelse(Nodes$Gene %in% TGFB_list, "TGFB signature", 
                          ifelse(Nodes$Gene %in% ALTEJ_list, "ALTEJ signature", NA))
rownames(Nodes) <- Nodes$Gene
##Create the tbl_graph object.
Tbl_graph1 <- tbl_graph(nodes = Nodes, edges = Edges, directed = FALSE)

#5. Plot the weighted gene coexpression network.
par(pty="s")
ggraph(Tbl_graph1, layout = "fr", weights = r) +
  geom_edge_link2(aes(colour=class), edge_alpha=0.2, edge_width=0.4) +
  geom_node_text(aes(label = Gene), colour="black", size=3.5, repel = TRUE) +
  geom_node_point(aes(fill=signature, size=centrality_degree(weights=r)), shape=21) +
  scale_edge_color_manual(values = c("springgreen4", "grey68", "hotpink2")) + 
  scale_fill_manual(values = c("seagreen4", "hotpink2")) + 
  scale_size(range = c(0.1, 10)) +
  theme_graph() + rremove("legend") #PDF7x8




