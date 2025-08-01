#### UpSetPlot for Intersection ####
## Author: Epitácio Farias
## Date: 11/11/2023

library("tidyverse")
if(!require("ComplexHeatmap")){BiocManager::install("ComplexHeatmap")}
# install.packages("devtools")
if(!require("patchwork")){devtools::install_github("thomasp85/patchwork")}
library("ggplot2")
library("ggrepel")

## Loading the Data 
load("~/Tesis_Project/Kidney/0_data/DEG_DESeq2.RData")
load("~/Tesis_Project/Kidney/0_data/DGE_limma.RData")
load("~/Tesis_Project/Kidney/0_data/DGE_edgeR.RData")

#### Removing the Non-Significant Genes 
DEG_unpaired <-  DEG_unpaired[DEG_unpaired$signif %in% c("Upregulated", "Downregulated"),]
DGE_RAxRT_edgeR_unpaired <-  DGE_RAxRT_edgeR_unpaired[DGE_RAxRT_edgeR_unpaired$signif %in% c("Upregulated", "Downregulated"),]
DGE_RAxRT_limma <-  DGE_RAxRT_limma[DGE_RAxRT_limma$signif %in% c("Upregulated", "Downregulated"),]

## UpSet Analysis
inter_list <- list(edgeR = DGE_RAxRT_edgeR_unpaired$ENSEMBL_ID,
                   limma = DGE_RAxRT_limma$genes,
                   DESeq2 = DEG_unpaired$gene_name)

inter_mx <- make_comb_mat(inter_list, mode = "intersect")

tiff(filename = "~/Tesis_Project/Kidney/Figures/UpsetPlot_DEA_Genes.tiff", res = 300, width = 5000, height = 2000)
UpSet(inter_mx, comb_order = order(comb_size(inter_mx)), 
      top_annotation = upset_top_annotation(inter_mx, add_numbers = T), right_annotation = upset_right_annotation(inter_mx, add_numbers = T))
dev.off()

distin_mx <- make_comb_mat(inter_list, mode = "distinct")

UpSet(distin_mx, comb_order = order(comb_size(distin_mx)), 
      top_annotation = upset_top_annotation(distin_mx, add_numbers = T), right_annotation = upset_right_annotation(distin_mx, add_numbers = T))

intersect_DE <- extract_comb(inter_mx,"111")

View(intersect_DE)

#### Selecting the genes on the intersection of the DEA methods ###
load("~/Tesis_Project/Kidney/0_data/KIRC_data_updated.RData")

DEG_genes_final <- DEG_unpaired[DEG_unpaired$gene_name %in% intersect_DE,]

DEG_genes_final <- DEG_genes_final[DEG_genes_final$padj != 0,]

DEG_genes_final <- DEG_genes_final %>% 
  mutate(
    padj = ifelse(is.na(padj), 1, padj),
    signif = 
      case_when(
        padj <= 0.01 & log2FoldChange < -1.5 ~ "Downregulated",
        padj <= 0.01 & log2FoldChange > 1.5 ~ "Upregulated",
        padj > 0.01  ~ "Not significant"
      )
  )

DEG_genes_final <- merge(DEG_genes_final, gene_library, by.x = "gene_name",by.y = "gene_id", all.x = T)

cutoff_pc <- sort(DEG_genes_final$padj)[10]

DEG_genes_final <- DEG_genes_final |> mutate(TopGeneLabel=ifelse(padj<=cutoff_pc, gene_name.y, ""))

p1 <- ggplot(DEG_genes_final, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
      geom_point() + labs(title = "Differentially Expressed Protein Coding Genes" ) +
      scale_color_manual(values = c("Downregulated" = "blue", "Not significant" = "grey", "Upregulated" = "red")) + guides(fill = "none") +
      geom_label_repel(aes(label=TopGeneLabel), 
                   seed = 123,
                   max.time = 3,
                   max.iter = Inf,
                   size = 5,
                   box.padding = 2, 
                   max.overlaps = Inf) + coord_cartesian(ylim = c(0,350))

tiff(filename = "~/Tesis_Project/Kidney/Figures/DEA_Genes.tiff", res = 300, width = 5000, height = 2000)
p1
dev.off()

save(DEG_genes_final, file = "~/Tesis_Project/Kidney/0_data/DEG_intersec.RData")
