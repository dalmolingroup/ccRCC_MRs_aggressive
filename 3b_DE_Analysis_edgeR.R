#### Differential Expression Analysis ####
## Date: 09/11/2023
## Author: Epitácio Farias 

#### Loading the Packages ####
pkg.bioconductor <- c("edgeR","SummarizedExperiment")
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

#### DEA with DESeq2 ####

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

## Reading and Normalizating the Data
y <- DGEList(count = exp_data, genes = rownames(exp_data))
y <- normLibSizes(y)

## Data Exploration
plotMDS(y)

## Desing Matrix
subject <- donor_data_sub$submitted_donor_id
treat <- donor_data_sub$tissue_type
design <- model.matrix(~treat)

## Dispersion Estimation
y2_RTRA <- estimateDisp(y,design)
plotBCV(y2_RTRA)

## Differential Expression
# Fitting the genewise glms:
fit_RTRA <-  glmQLFit(y2_RTRA, design)

# Likelihood ratio test (RA vs RT) and the top genes:
lrt_RTRA <- glmLRT(fit_RTRA)
topTags(lrt_RTRA)

#  Total of DE at 5% of FDR:
summary(decideTests(lrt_RTRA))

# Plotting logFC vs logCPM:
plotMD(lrt_RTRA)
abline(h=c(-1, 1), col="blue")

## Differentially Expressed Genes
res_DGE_RAxRT_edgeR_unpaired <- lrt_RTRA[["table"]]
res_DGE_RAxRT_edgeR_unpaired <- rownames_to_column(res_DGE_RAxRT_edgeR_unpaired, var = "ENSEMBL_ID")
res_DGE_RAxRT_edgeR_unpaired$Padjust <- p.adjust(res_DGE_RAxRT_edgeR_unpaired$PValue, method = "fdr")

res_DGE_RAxRT_edgeR_unpaired <- res_DGE_RAxRT_edgeR_unpaired %>% 
  mutate(
    Padjust = ifelse(is.na(Padjust), 1, Padjust),
    signif = 
      case_when(
        Padjust <= 0.01 & logFC < 0 ~ "Downregulated",
        Padjust <= 0.01 & logFC > 0 ~ "Upregulated",
        Padjust > 0.01 ~ "Not significant"
      )
  )

ggplot(res_DGE_RAxRT_edgeR_unpaired, aes(x = logFC, y = -log10(Padjust), color = signif)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "blue", "Not significant" = "grey", "Upregulated" = "red"))

DGE_RAxRT_edgeR_unpaired <- res_DGE_RAxRT_edgeR_unpaired |> filter(Padjust < 0.01) |> arrange(desc(logFC))

## Correcting gene names
DGE_RAxRT_edgeR_unpaired <- DGE_RAxRT_edgeR_unpaired[!duplicated(DGE_RAxRT_edgeR_unpaired$ENSEMBL_ID),]

#### Saving the DE ####
save(res_DGE_RAxRT_edgeR_unpaired, DGE_RAxRT_edgeR_unpaired, file = "~/Tesis_Project/Kidney/0_data/DGE_edgeR.RData")
