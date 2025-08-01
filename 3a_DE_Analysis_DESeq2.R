#### Differential Expression Analysis ####
## Date: 09/11/2023
## Author: Epitácio Farias 

#### Loading the Packages ####
pkg.bioconductor <- c("DESeq2","SummarizedExperiment")
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

## DEA for paired data with DESeq2
dds_unpaired <- DESeqDataSetFromMatrix(countData = as.matrix(exp_data),
                                       colData = donor_data_sub,
                                       design = ~ tissue_type + submitted_donor_id)

dds_unpaired <- DESeq(dds_unpaired)

resultsNames(dds_unpaired)

plotPCA(vst(dds_unpaired), intgroup="tissue_type")

stand_unpaired <- vst(dds_unpaired)

# Computar PCA
tt_unpaired <- PCAtools::pca(assay(stand_unpaired), metadata = as.matrix(colData(dds_unpaired)))

# Plotar PC1 e PC2
PCAtools::biplot(tt_unpaired, colby = "tissue_type", shape = "tissue_type", legendPosition = "left")

# Proporção de variabilidade de cada componente
PCAtools::screeplot(tt_unpaired)

# PCAtools::eigencorplot(tt_unpaired, components = PCAtools::getComponents(tt_unpaired, c(1:10)), metavars = c("sizeFactor", "tissue_type"))

#Primary x normal paired
res_unpaired <- results(dds_unpaired, contrast = c("tissue_type","RT","RA"))

#PRIMARY X NORMAL PAIRED
res_unpaired <- res_unpaired %>% as.data.frame() %>% tidyr::drop_na()

res_unpaired <- tibble::rownames_to_column(res_unpaired,"gene_name")

#duplicated genes
res_unpaired[duplicated(res_unpaired$gene_name), "padj"]

#name of differentially expressed genes
DEG_unpaired <- res_unpaired %>% 
  mutate(
    padj = ifelse(is.na(padj), 1, padj),
    signif = 
      case_when(
        padj <= 0.01 & log2FoldChange < 0 ~ "Downregulated",
        padj <= 0.01 & log2FoldChange > 0 ~ "Upregulated",
        padj > 0.01 ~ "Not significant"
      )
  )

ggplot(DEG_unpaired, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point() +
  scale_color_manual(values = c("Downregulated" = "blue", "Not significant" = "grey", "Upregulated" = "red"))

#foldchange
fold_unpaired <- res_unpaired$log2FoldChange
names(fold_unpaired) <- res_unpaired$gene_name

DEG_unpaired <- DEG_unpaired %>%
  filter(padj <= 0.01) %>%
  arrange(desc(log2FoldChange))

symbol_paired <- DEG_unpaired$gene_name

save(DEG_unpaired, fold_unpaired, symbol_paired, res_unpaired, tt_unpaired, stand_unpaired,dds_unpaired,  file = "~/Tesis_Project/Kidney/0_data/DEG_DESeq2.RData")
