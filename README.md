# Under-review-article

This repository is part of a pre-publication release containing all the code used in the analyses from the project called "Validation of TGFβ signaling and alternative end‐joining DNA repair signatures that predict response to genotoxic cancer therapy". The results are currently under review for publication. 

These analyses report a functional relationship between TGFβ signaling and DNA damage repair that can be used to predict the response of cancer patients to genotoxic treatment. 

The gene signatures used in here were reported in our recent Science Translational Medicine publication from earlier this year (link: https://stm.sciencemag.org/content/13/580/eabc4465). Here we report their biological validation in primary head and neck squamous carcinoma specimens, to demonstrate that the signatures accurately report biological processes. Next, we used functional readouts to weight each signature gene by its association with biological responses to define a modified score -- termed βaltw-- that was retested for association with response to genotoxic therapies in independent datasets.

Each folder (Inhouse-HNSC, Inhouse-OV, TCGA-pancancer, TCGA-HNSC, GSE41613-HNSC, etc.) contains the code that was used to analyze each dataset. Among many others, some of the analytical processes that can be found in here are: 
1. Analysis of the association between the expression of multiple genes and biological readouts of TGFβ activity and DNA damage repair proficiency. 
2. Construction of a weighted gene coexpression network and calculation of genes' centrality degree within it. 
3. Assignment of a weight to each gene based on its centrlity degree and the strength of its association with biological readouts. 
4. Definition and calculation of the βaltw score. 
5. Comparision of the predictive performance of the βaltw score and the former version of it in 500 datasets generated with bootstrapping. 
6. Analysis of the association of the βaltw score with cancer patients' outcome after genotoxic treatment. 
7. Analysis of the association of the βaltw score with cell lines' radiosensitivity. 
