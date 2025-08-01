#### Differential Expression Analysis ####
## Date: 09/11/2023
## Author: Epitácio Farias 

#### Loading the Packages ####
pkg.bioconductor <- c("limma","SummarizedExperiment")
pkg.cran <- c("tidyverse","rtracklayer")

#check if each package is on the local machine
#if it is not installed, it will be installed and loaded
pkg.check <- lapply(pkg.bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)}
})
pkg.check <- lapply(pkg.cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)}
})

rm(pkg.cran, pkg.bioconductor, pkg.check)

#### Data from ICGC - RECA  ####
load("~/Tesis_Project/Kidney/0_data/ICGC_DEA_Data.RData")

uniques <- donor_data_final %>% group_by(submitted_donor_id) %>% tally() %>% filter(n==1) %>% pull(submitted_donor_id)
donor_data_sub <- donor_data_final[!(donor_data_final$submitted_donor_id %in% uniques),]

#### Count Data from ICGC - RECA ####
exp_data <- count_final_ICGC

sum(duplicated(exp_data$gene_id))

rownames(exp_data) <- exp_data$gene_id

#### Clinical Data from TCGA - KIRC ####
## Normal Clinical Data
rownames(donor_data_sub) <- stringr::str_c(donor_data_sub$submitted_donor_id,"",donor_data_sub$tissue_type)

#### Selecting pairs for DEA Analysis ####
exp_data <- exp_data[,colnames(exp_data) %in% rownames(donor_data_sub)]

rm(donor_data_final, count_final_ICGC, uniques)

#### DEA with limma ####

## Organazing the data for each DEA 
# Normal x Primary 

sapply(donor_data_sub, class)
donor_data_sub$submitted_donor_id <- as.factor(donor_data_sub$submitted_donor_id)
donor_data_sub$pathologic_t <- as.factor(donor_data_sub$pathologic_t)
donor_data_sub$pathologic_n <- as.factor(donor_data_sub$pathologic_n)
donor_data_sub$pathologic_m <- as.factor(donor_data_sub$pathologic_m)
donor_data_sub$tissue_type <- as.factor(donor_data_sub$tissue_type)
donor_data_sub$gender <- as.factor(donor_data_sub$gender)

donor_data_sub %>%
  dplyr::count(tissue_type) 

donor_data_sub %>%
  dplyr::count(pathologic_m) 

exp_data <- exp_data %>% select(order(colnames(exp_data)))

index_order <- order(rownames(donor_data_sub))
donor_data_sub <- donor_data_sub[index_order,]

rownames(donor_data_sub) == colnames(exp_data)

## Desing Matrix
subject <- donor_data_sub$submitted_donor_id
treat <- donor_data_sub$tissue_type
design <- model.matrix(~treat)

## Reading and Normalizating the Data
y <- edgeR::DGEList(count = exp_data, genes = rownames(exp_data))
keep <- edgeR::filterByExpr(y, design)

y <- y[keep,,keep.lib.sizes=FALSE]

y <- edgeR::calcNormFactors(y)

## Differential Expression: voom
v <- voom(y, design, plot = T)

# Fitting the genewise glms:
fit <-  lmFit(v, design)

fit <-  eBayes(fit)

res_DGE_RAxRT_limma <- topTable(fit, coef = ncol(design), number = "inf")

## Differentially Expressed Genes
#name of differentially expressed genes
res_DGE_RAxRT_limma <- res_DGE_RAxRT_limma %>% 
  mutate(
    adj.P.Val = ifelse(is.na(adj.P.Val), 1, adj.P.Val),
    signif = 
      case_when(
        adj.P.Val <= 0.01 & logFC < 0 ~ "Downregulated",
        adj.P.Val <= 0.01 & logFC > 0 ~ "Upregulated",
        adj.P.Val > 0.01 ~ "Not significant"
      )
  )

ggplot(res_DGE_RAxRT_limma, aes(x = logFC, y = -log10(adj.P.Val), color = signif)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "blue", "Not significant" = "grey", "Upregulated" = "red"))


DGE_RAxRT_limma <- res_DGE_RAxRT_limma |> filter(adj.P.Val < 0.01) |> arrange(desc(logFC))

## Correcting gene names
DGE_RAxRT_limma <- DGE_RAxRT_limma[!duplicated(DGE_RAxRT_limma$genes),]

rm(design, donor_data_sub,exp_data,fit,v,y,index_order, keep,subject,treat)

#### Saving the DE ####
save(DGE_RAxRT_limma,res_DGE_RAxRT_limma, file = "~/Tesis_Project/Kidney/0_data/DGE_limma.RData")
