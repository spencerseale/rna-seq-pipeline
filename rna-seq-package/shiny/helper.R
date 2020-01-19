#Install 

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require(biomaRt)) BiocManager::install("biomaRt")
if(!require(edgeR)) BiocManager::install("edgeR")
if(!require(EnhancedVolcano)) BiocManager::install('EnhancedVolcano')
if(!require(DESeq2)) BiocManager::install("DESeq2")
if(!require(GO.db)) BiocManager::install("GO.db")
if(!require(org.Ce.eg.db)) BiocManager::install("org.Ce.eg.db")
if(!require(VennDiagram)) install.packages("VennDiagram", repos = "http://cran.us.r-project.org")
if(!require(RColorBrewer)) install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")

library(edgeR)
library(biomaRt)
library(org.Ce.eg.db)
library(GO.db)
library(EnhancedVolcano)
library(DESeq2)
library(RColorBrewer)
library(VennDiagram)
library(gridExtra)
library(grid)

setwd('.')

#reading in counts data, removing last 5 lines that don't contain gene counts
GenewiseCounts <- read.table("../merged_counts.tsv", header = TRUE, row.names = 1, sep="\t", check.names = FALSE)
GenewiseCounts <- GenewiseCounts[1:(nrow(GenewiseCounts)-5),]

#setting treatments
conditions <- colnames(GenewiseCounts)
group <- vector()
for (cond in conditions) {
  if (grepl("_A", cond)) {
    group <- append(group, "group_A")
  }
  else if (grepl("_B", cond)) {
    group <- append(group, "group_B")
  }
  else if (grepl("_C", cond)) {
    group <- append(group, "group_C")
  }
  else if (grepl("_D", cond)) {
    group <- append(group, "group_D")
  }
}

#determining number of pairwise comparisons
if ("group_D" %in% group) {
  treatment_num <- 4
} else if ("group_C" %in% group) {
  treatment_num <- 3
} else if ("group_B" %in% group) {
  treatment_num <- 2
}

#getting gene name from biomart using ensembl id and storing as var symbol
mart <- useDataset("celegans_gene_ensembl", useMart("ensembl"))
symbol <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "external_gene_name"),
  values=rownames(GenewiseCounts),
  mart=mart)

#creating DGE object, getting dispersion, fitting generalized linear model, and identifying gene names for each ensembl id

#setting dif gene expression list and normalizing
dge.er <- DGEList(counts=GenewiseCounts)
dge.er <- calcNormFactors(dge.er)

#setting design variable
design.er <- model.matrix(~0 + group)

#running all dispersion options
dge.er <- estimateGLMCommonDisp(dge.er, design.er)
dge.er <- estimateGLMTrendedDisp(dge.er, design.er)
dge.er <- estimateGLMTagwiseDisp(dge.er, design.er)

#fitting a negative binomial generalized log-linear model
fit.er <- glmFit(dge.er, design.er)

#3 treatments
#specifying contrasts, normalizing data, liklihood ratio tests, adding gene names, and producing venn diagram

if (treatment_num == 3){
  
  #setting all contrasts for DGE (including GO)
  contrasts <- makeContrasts(A_vs_B=groupgroup_A-groupgroup_B, A_vs_C=groupgroup_A-groupgroup_C, B_vs_C=groupgroup_B-groupgroup_C, levels=design.er)
  
  #conducting likelihood ratio tests for each pairwise comparison
  A_vs_B_lrt <- glmLRT(fit.er, contrast=contrasts[,"A_vs_B"])
  A_vs_C_lrt <- glmLRT(fit.er, contrast=contrasts[,"A_vs_C"])
  B_vs_C_lrt <- glmLRT(fit.er, contrast=contrasts[,"B_vs_C"])
  
  #swapping gene symbol in for ensembl id
  row.names(A_vs_B_lrt$table) <- symbol$external_gene_name
  row.names(A_vs_C_lrt$table) <- symbol$external_gene_name
  row.names(B_vs_C_lrt$table) <- symbol$external_gene_name
  
  #identifying significantly expressed genes based off of trusted edgeR p-values and stats
  A_vs_B_de <- decideTestsDGE(A_vs_B_lrt, adjust.method = "fdr")
  A_vs_C_de <- decideTestsDGE(A_vs_C_lrt, adjust.method = "fdr")
  B_vs_C_de <- decideTestsDGE(B_vs_C_lrt, adjust.method = "fdr")
  
  #heatmap normalization, gene name fetching, and pairwise fetching
  #creating logCPM table containing all samples to be used for heat maps
  logCPM <- cpm(dge.er, prior.count=2, log=TRUE)
  #rawCPM <- cpm(dge.er, prior.count=2, log=FALSE)
  
  #swapping in gene symbol for ensembl ID
  row.names(logCPM) <- symbol$external_gene_name
  #row.names(rawCPM)<- symbol$external_gene_name
  #write.table(logCPM[1:46904,1:6], file="logCPM.tsv", sep="\t")
  #write.table(rawCPM[1:46904,1:6],file="rawCPM.tsv", sep="\t")
  
  #subsetting each pairwise comparison
  logCPM_A_vs_B <- logCPM[,grepl("_A|_B", colnames(logCPM))]
  logCPM_A_vs_C <- logCPM[,grepl("_A|_C", colnames(logCPM))]
  logCPM_B_vs_C <- logCPM[,grepl("_B|_C", colnames(logCPM))]
}



#2 treatments
#specifying contrasts, normalizing data, liklihood ratio tests, adding gene names, and producing venn diagram

if (treatment_num == 2) {
  
  #print("In a 2-way comparison, a venn diagram is not outputted")
  
  #setting all contrasts for DGE (including GO)
  contrasts <- makeContrasts(A_vs_B=groupgroup_A-groupgroup_B,levels=design.er)
  
  #conducting likelihood ratio tests for each pairwise comparison
  A_vs_B_lrt <- glmLRT(fit.er, contrast=contrasts[,"A_vs_B"])
  
  #swapping gene symbol in for ensembl id
  row.names(A_vs_B_lrt$table) <- symbol$external_gene_name
  
  #identifying significantly expressed genes based off of trusted edgeR p-values and stats
  A_vs_B_de <- decideTestsDGE(A_vs_B_lrt, adjust.method = "fdr")
  
  #heatmap normalization, gene name fetching, and pairwise fetching
  #creating logCPM table containing all samples to be used for heat maps
  logCPM <- cpm(dge.er, prior.count=2, log=TRUE)
  #swapping in gene symbol for ensembl ID
  row.names(logCPM) <- symbol$external_gene_name
  
  #subsetting each pairwise comparison
  logCPM_A_vs_B <- logCPM[,grepl("_A|_B", colnames(logCPM))]
  
}

#4 treatments
#specifying contrasts, normalizing data, liklihood ratio tests, adding gene names, and producing venn diagram

if (treatment_num == 4) {
  
  
  #setting all contrasts for DGE (including GO)
  contrasts <- makeContrasts(A_vs_B=groupgroup_A-groupgroup_B, A_vs_C=groupgroup_A-groupgroup_C, B_vs_C=groupgroup_B-groupgroup_C, A_vs_D=groupgroup_A-groupgroup_D, B_vs_D=groupgroup_B-groupgroup_D, C_vs_D=groupgroup_C-groupgroup_D,levels=design.er)
  
  #conducting likelihood ratio tests for each pairwise comparison
  A_vs_B_lrt <- glmLRT(fit.er, contrast=contrasts[,"A_vs_B"])
  A_vs_C_lrt <- glmLRT(fit.er, contrast=contrasts[,"A_vs_C"])
  B_vs_C_lrt <- glmLRT(fit.er, contrast=contrasts[,"B_vs_C"])
  A_vs_D_lrt <- glmLRT(fit.er, contrast=contrasts[,"A_vs_D"])
  B_vs_D_lrt <- glmLRT(fit.er, contrast=contrasts[,"B_vs_D"])
  C_vs_D_lrt <- glmLRT(fit.er, contrast=contrasts[,"C_vs_D"])
  
  #swapping gene symbol in for ensembl id
  row.names(A_vs_B_lrt$table) <- symbol$external_gene_name
  row.names(A_vs_C_lrt$table) <- symbol$external_gene_name
  row.names(B_vs_C_lrt$table) <- symbol$external_gene_name
  row.names(A_vs_D_lrt$table) <- symbol$external_gene_name
  row.names(B_vs_D_lrt$table) <- symbol$external_gene_name
  row.names(C_vs_D_lrt$table) <- symbol$external_gene_name
  
  #identifying significantly expressed genes based off of trusted edgeR p-values and stats
  A_vs_B_de <- decideTestsDGE(A_vs_B_lrt, adjust.method = "fdr")
  A_vs_C_de <- decideTestsDGE(A_vs_C_lrt, adjust.method = "fdr")
  B_vs_C_de <- decideTestsDGE(B_vs_C_lrt, adjust.method = "fdr")
  A_vs_D_de <- decideTestsDGE(A_vs_D_lrt, adjust.method = "fdr")
  B_vs_D_de <- decideTestsDGE(B_vs_D_lrt, adjust.method = "fdr")
  C_vs_D_de <- decideTestsDGE(C_vs_D_lrt, adjust.method = "fdr")
  
  #heatmap normalization, gene name fetching, and pairwise fetching
  #creating logCPM table containing all samples to be used for heat maps
  logCPM <- cpm(dge.er, prior.count=2, log=TRUE)
  #swapping in gene symbol for ensembl ID
  row.names(logCPM) <- symbol$external_gene_name
  
  #subsetting each pairwise comparison
  logCPM_A_vs_B <- logCPM[,grepl("_A|_B", colnames(logCPM))]
  logCPM_A_vs_C <- logCPM[,grepl("_A|_C", colnames(logCPM))]
  logCPM_B_vs_C <- logCPM[,grepl("_B|_C", colnames(logCPM))]
  logCPM_A_vs_D <- logCPM[,grepl("_A|_D", colnames(logCPM))]
  logCPM_B_vs_D <- logCPM[,grepl("_B|_D", colnames(logCPM))]
  logCPM_C_vs_D <- logCPM[,grepl("_C|_D", colnames(logCPM))]
}


### Volcano Plots and Heatmaps

if (treatment_num == 3) {
  o1 <- order(A_vs_B_lrt$table$PValue)
  o2 <- order(A_vs_C_lrt$table$PValue)
  o3 <- order(B_vs_C_lrt$table$PValue)
}

if (treatment_num == 2) {
  o1 <- order(A_vs_B_lrt$table$PValue)
 }

if (treatment_num == 4) {
  o1 <- order(A_vs_B_lrt$table$PValue)
  o2 <- order(A_vs_C_lrt$table$PValue)
  o3 <- order(B_vs_C_lrt$table$PValue)
  o4 <- order(A_vs_D_lrt$table$PValue)
  o5 <- order(B_vs_D_lrt$table$PValue)
  o6 <- order(C_vs_D_lrt$table$PValue)
}

## Gene Ontology Analysis

#fitting glm and preparing data for GO
#getting entrez id from biomart using ensembl id, necessary for running GO
#mart var assigned earlier in script 
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=rownames(GenewiseCounts),
  mart=mart)

#swapping entrez id in for ensembl id
EntGenewiseCounts <- cbind(GenewiseCounts, genes$entrezgene_id)
#removing rows where entrez ID is NA, cannot run GO unless an entrez id exists
EntGenewiseCounts <- na.omit(EntGenewiseCounts)
#removing duplicate entrez IDs
entrezCountsUnique <- EntGenewiseCounts[!duplicated(EntGenewiseCounts[,ncol(EntGenewiseCounts)]),]
#setting row names as entrez ID rather than ensembl ID
row.names(entrezCountsUnique) <- entrezCountsUnique[,ncol(entrezCountsUnique)]
#removing entrez ID from the last column as it is now the row name
entrezCountsUnique <- entrezCountsUnique[,1:ncol(entrezCountsUnique)-1]

#####Setting DGE list for GO
GO_dge.er <- DGEList(counts=entrezCountsUnique)

#looking at number of genes in dge variable
GO_dge.er <- calcNormFactors(GO_dge.er)

#getting dispersions
#using same design.er as created earlier
GO_dge.er <- estimateGLMCommonDisp(GO_dge.er, design.er)
GO_dge.er <- estimateGLMTrendedDisp(GO_dge.er, design.er)
GO_dge.er <- estimateGLMTagwiseDisp(GO_dge.er, design.er)

#fitting linear model
GO_fit.er <- glmFit(GO_dge.er, design.er)

if (treatment_num == 3) {
  
  #conducting likelihood ratio tests for each pairwise comparison
  GO_A_vs_B_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"A_vs_B"])
  GO_A_vs_C_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"A_vs_C"])
  GO_B_vs_C_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"B_vs_C"])
}

if (treatment_num == 2) {
  
  #conducting likelihood ratio tests for each pairwise comparison
  GO_A_vs_B_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"A_vs_B"])
}

if (treatment_num == 4) {
  
  #conducting likelihood ratio tests for each pairwise comparison
  GO_A_vs_B_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"A_vs_B"])
  GO_A_vs_C_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"A_vs_C"])
  GO_B_vs_C_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"B_vs_C"])
  GO_A_vs_D_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"A_vs_D"])
  GO_B_vs_D_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"B_vs_D"])
  GO_C_vs_D_lrt <- glmLRT(GO_fit.er, contrast=contrasts[,"C_vs_D"])
}

if (treatment_num == 3) {
  
  #print("Top 15 up: Treatment A vs Treatment B")
  go_A_vs_B <- goana(GO_A_vs_B_lrt, FDR = 0.05, species="Ce")
  go_A_vs_C <- goana(GO_A_vs_C_lrt, FDR = 0.05, species="Ce")
  go_B_vs_C <- goana(GO_B_vs_C_lrt, FDR = 0.05, species="Ce")
}

if (treatment_num == 2) {
  
  #print("Top 15 up: Treatment A vs Treatment B")
  go_A_vs_B <- goana(GO_A_vs_B_lrt, FDR = 0.05, species="Ce")
}

if (treatment_num == 4) {
  
  #print("Top 15 up: Treatment A vs Treatment B")
  go_A_vs_B <- goana(GO_A_vs_B_lrt, FDR = 0.05, species="Ce")
  go_A_vs_C <- goana(GO_A_vs_C_lrt, FDR = 0.05, species="Ce")
  go_B_vs_C <- goana(GO_B_vs_C_lrt, FDR = 0.05, species="Ce")
  go_A_vs_D <- goana(GO_A_vs_D_lrt, FDR = 0.05, species="Ce")
  go_B_vs_D <- goana(GO_B_vs_D_lrt, FDR = 0.05, species="Ce")
  go_C_vs_D <- goana(GO_C_vs_D_lrt, FDR = 0.05, species="Ce")
}
