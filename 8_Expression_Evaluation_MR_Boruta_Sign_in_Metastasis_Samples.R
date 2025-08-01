###########################################
##### Expression Evaluation of the MR #####
##### Author: Farias, Epitácio  ###########
##### Date: 08/05/2024          ###########

## Install packages
pkg.bioconductor <- c("TCGAbiolinks", "TCGAWorkflow", "clusterProfiler","org.Hs.eg.db" )
pkg.cran <- c("tidyverse", "skimr", "tableone", "survminer", "survival","vroom",
              "finalfit", "DT", "data.table")
# optional pkgs: "finalfit", "DT", "data.table"

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
rm(pkg.bioconductor, pkg.check, pkg.cran)

#### TCGA-KIRC Data set ######
#TCGA: primary and normal samples
query_kirc <- GDCquery(project = "TCGA-KIRC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       sample.type = c("Primary Tumor", "Solid Tissue Normal"))

results_kirc <- getResults(query_kirc)

GDCdownload(query = query_kirc,
            method = "api")

KIRC <- GDCprepare(query = query_kirc,
                   summarizedExperiment = T,
                   directory = "./GDCdata/")

#### Clinical TCGA-KIRC Data Adjustment ####
## Retrieving and Taming the Clinical Data
clin_data <- as.data.frame(KIRC@colData)

glimpse(clin_data)

#### Transcription  TCGA-KIRC Data Adjustment ####
## Retrieving the sequencing Gene Library 
gene_library <- as.data.frame(KIRC@rowRanges@elementMetadata@listData)

gene_library <- gene_library[,c(5:7)]

## Retrieving the Raw Count Data
count_data_kirc <- SummarizedExperiment::assay(KIRC) |> as.data.frame() |> rownames_to_column("ENSEMBL_ID")

count_data_kirc <- merge(count_data_kirc, gene_library, by.x = "ENSEMBL_ID", by.y = "gene_id", all.x = T)

count_data_kirc <- count_data_kirc[,c(1,615,616,2:614)]

sum(duplicated(count_data_kirc$ENSEMBL_ID))

sum(duplicated(count_data_kirc$gene_name))

count_data_kirc$ENSEMBL_ID <- sub("\\..*", "", count_data_kirc$ENSEMBL_ID)

sum(duplicated(count_data_kirc$ENSEMBL_ID))

count_data_kirc <- count_data_kirc[!duplicated(count_data_kirc$ENSEMBL_ID),]

rownames(count_data_kirc) <- NULL

count_data_kirc <- count_data_kirc |> column_to_rownames("ENSEMBL_ID")

## Adjusting the gene library 
gene_library$gene_id <- sub("\\..*", "", gene_library$gene_id)

sum(!duplicated(gene_library$gene_id))

gene_library <- gene_library[!duplicated(gene_library$gene_id),]

all(rownames(count_data_kirc) %in% gene_library$gene_id)

## Selecting only the protein_codings 
gene_library2 <- gene_library[gene_library$gene_type == "protein_coding",]
count_data_kirc2 <- count_data_kirc[rownames(count_data_kirc) %in% gene_library2$gene_id, ]

#### Saving the Data ####
save(count_data_kirc2, gene_library2, clin_data, file = "~/Kidney/0_data/KIRC_MR_expression_eval_data.RData" )
rm(list = ls())

#### Evaluating the MR Expression TCGA-KIRC ####
selected_features <- read.table("~/Kidney/0_data/selected_features_boruta_balanced_2.txt","\n", header = F)
colnames(selected_features) <- "ENSEMBL_ID"

selected_features <- bitr(selected_features$ENSEMBL_ID, fromType= "ENSEMBL", toType = c("SYMBOL","ENTREZID","ENSEMBL"),
                          OrgDb = org.Hs.eg.db, drop = T)

selected_features[11,] <- c("ENSG00000215271", "HOMEZ", "20164")

load("~/Kidney/0_data/KIRC_MR_expression_eval_data.RData")

inteserct_pati <- intersect(colnames(count_data_kirc2), rownames(clin_data))

count_df <- count_data_kirc2[,colnames(count_data_kirc2) %in% inteserct_pati]

clin_data <- clin_data |> rownames_to_column("patients") |> as.data.frame()

clin_data <- clin_data[clin_data$patients %in% inteserct_pati,]

clin_data$ajcc_pathologic_m <- as.factor(clin_data$ajcc_pathologic_m)

clin_data <-  clin_data[clin_data$ajcc_pathologic_m != "MX",]

a <- clin_data[!is.na(clin_data$ajcc_pathologic_m),]

rownames(a) <- a$patients

b <- count_df[, colnames(count_df) %in% a$patients]

dds_kirc <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(b),
                                           colData = a,
                                           design = ~ ajcc_pathologic_m)
dds_kirc <- DESeq2::DESeq(dds_kirc)

res <- DESeq2::results(dds_kirc, tidy=TRUE, contrast=c("ajcc_pathologic_m", "M1", "M0")) %>%
  arrange(padj, pvalue) %>%
  tibble::as.tibble()
res

stopifnot(all(selected_features$ENSEMBL %in% names(dds_kirc)))
selected_features

tcounts_kirc <- t(log2((counts(dds_kirc[selected_features$ENSEMBL, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(SummarizedExperiment::colData(dds_kirc), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(selected_features$ENSEMBL)+1):ncol(.))

tcounts_kirc %>% 
  dplyr::select(Row.names, sample_type, gene, expression) %>% 
  head %>% 
  knitr::kable()

tcounts_kirc <- merge(tcounts_kirc, gene_library2, by.x = "gene", by.y = "gene_id", all.x = T)

c1 <- c("HOMEZ","PBX1","ZFP28","ZNF136","ZNF658")
c2 <- c("E2F2", "FOXM1", "MYBL2", "NFE2L3", "PAX6","TFAP2A")
tcounts_kirc$Cluster <- ifelse(tcounts_kirc$gene_name %in% c1, "Cluster 1","Cluster 2")

library(rstatix)
stat_test_kirc <- tcounts_kirc %>%
  dplyr::group_by(gene_name,Cluster) %>%
  rstatix::wilcox_test(expression ~ ajcc_pathologic_m) %>%
  rstatix::adjust_pvalue(method = "BY") %>%
  rstatix::add_significance("p.adj")

stat_test_kirc <- stat_test_kirc %>%
  add_xy_position(x = "gene_name", dodge = 0.8)

save(tcounts_kirc,stat_test_kirc, file = "tcounts_kirc_M0M1_Boruta.RData")

#Ordenar
tcounts_kirc_2 <- tcounts_kirc %>% dplyr::arrange(factor(Cluster, levels = unique(stat_test_kirc$Cluster)))

#Boxplot
bxp <- ggboxplot(tcounts_kirc_2, x = "gene_name", y = "expression", 
                 color = "ajcc_pathologic_m")

tiff(file = "~/Kidney/Figures/Expression_Pattern_Behaviour_Boruta_Sign_Metastatic_Samples.tiff", res = 300, width = 2500, height = 1250)
bxp + 
  stat_pvalue_manual(stat_test_kirc, label = "p.adj.signif", 
                     tip.length = 0.01,bracket.nudge.y = 6, vjust = 0) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1, size = 10)) +
  labs(x="", y="Expression (log normalized counts)") +
  scale_color_discrete(name = "Groups",
                       labels = c("M0", "M1"))
dev.off()

#### Evaluating the MR Expression ICGC_RECA ####
selected_features <- read.table("~/Kidney/0_data/selected_features_boruta_balanced_2.txt","\n", header = F)
colnames(selected_features) <- "ENSEMBL_ID"

selected_features <- bitr(selected_features$ENSEMBL_ID, fromType= "ENSEMBL", toType = c("SYMBOL","ENTREZID","ENSEMBL"),
                          OrgDb = org.Hs.eg.db, drop = T)
library(DESeq2)

load("~/Kidney/0_data/ICGC_DEA_Data.RData")

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

donor_data_sub <- donor_data_sub[donor_data_sub$pathologic_m != "MX", ]

exp_data <- exp_data %>% dplyr::select(order(colnames(exp_data)))

index_order <- order(rownames(donor_data_sub))
donor_data_sub <- donor_data_sub[index_order,]

rownames(donor_data_sub) == colnames(exp_data)

exp_data <- exp_data[,colnames(exp_data) %in% rownames(donor_data_sub)]

dds_icgc <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(exp_data),
                                           colData = donor_data_sub,
                                           design = ~ pathologic_m)
dds_icgc <- DESeq2::DESeq(dds_icgc)

res <- DESeq2::results(dds_icgc, tidy=TRUE, contrast=c("pathologic_m", "M1", "M0")) %>%
  arrange(padj, pvalue) %>%
  tbl_df()
res

stopifnot(all(selected_features$ENSEMBL %in% names(dds_icgc)))
selected_features

selected_features[!(selected_features$ENSEMBL %in% names(dds_icgc)),]

test <- selected_features[selected_features$ENSEMBL %in% names(dds_icgc),]

tcounts_icgc <- t(log2((counts(dds_icgc[test$ENSEMBL, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(SummarizedExperiment::colData(dds_icgc), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-9+1):ncol(.))

tcounts_icgc %>% 
  dplyr::select(Row.names, tissue_type, gene, expression) %>% 
  head %>% 
  knitr::kable()

tcounts_icgc <- merge(tcounts_icgc, gene_library2, by.x = "gene", by.y = "gene_id", all.x = T)

c1 <- c("PBX1","ZFP28","ZNF136")
c2 <- c("E2F2", "FOXM1", "MYBL2", "NFE2L3", "PAX6","TFAP2A")
tcounts_icgc$Cluster <- ifelse(tcounts_icgc$gene_name %in% c1, "Cluster 1","Cluster 2")

library(rstatix)
stat_test_icgc <- tcounts_icgc %>%
  dplyr::group_by(gene_name, Cluster) %>%
  rstatix::wilcox_test(expression ~ pathologic_m) %>%
  rstatix::adjust_pvalue(method = "BY") %>%
  rstatix::add_significance("p.adj")

stat_test_icgc <- stat_test_icgc %>%
  add_xy_position(x = "gene_name", dodge = 0.8)

save(tcounts_icgc,stat_test_icgc, file = "tcounts_ICGC_M0M1_Boruta.RData")

#Ordenar
tcounts_icgc_2 <- tcounts_icgc %>% dplyr::arrange(factor(Cluster, levels = unique(stat_test_icgc$Cluster)))

#Boxplot
bxp <- ggboxplot(tcounts_icgc_2, x = "gene_name", y = "expression", 
                 color = "pathologic_m")

tiff(file = "~/Kidney/Figures/Expression_Pattern_Behaviour_Boruta_Sign_ICGC_Metastatic_Samples.tiff", res = 300, width = 2500, height = 1250)
bxp + 
  stat_pvalue_manual(stat_test_icgc, label = "p.adj.signif", 
                     tip.length = 0.01,bracket.nudge.y = 6, vjust = 0) +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1, size = 10)) +
  labs(x="", y="Expression (log normalized counts)") +
  scale_color_discrete(name = "Groups",
                       labels = c("M0", "M1"))
dev.off()
