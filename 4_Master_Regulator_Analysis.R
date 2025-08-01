#### TF Signature Construction ####
## Date: 09/14/2023
## Author: Epitácio Farias

#### Loading packages and data ####
library("RTN")

load("~/Tesis_Project/Kidney/0_data/DEG_DESeq2.RData")
load("~/Tesis_Project/Kidney/0_data/DEG_intersec.RData")
load("~/Tesis_Project/Kidney/0_data/tf_list.RData")
load("~/Tesis_Project/Kidney/0_data/rtn_TF_data.RData")

#### Correcting some issues ####
fold_paired <- DEG_genes_final$log2FoldChange
names(fold_paired) <- DEG_genes_final$gene_name

symbol_paired <- DEG_genes_final$gene_name

#### Creating RTN Analysis Object #### 
rtna_kirc_prim_norm <- tni2tna.preprocess(object = rtni_TF,
                                          phenotype = fold_paired,
                                          hits = symbol_paired)
#### Master Regulator Analysis ####
rtna_kirc_prim_norm <- tna.mra(rtna_kirc_prim_norm, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

rtna_kirc_prim_norm <- tna.gsea1(rtna_kirc_prim_norm, nPermutations=1000)
gsea1 <- tna.get(rtna_kirc_prim_norm, what="gsea1", ntop = -1)

rtna_kirc_prim_norm <- tna.gsea2(rtna_kirc_prim_norm, nPermutations=1000)
gsea2 <- tna.get(rtna_kirc_prim_norm, what="gsea2", ntop = -1)
gsea2_diff <- gsea2$differential

mra_kirc_prim_norm <- tna.get(rtna_kirc_prim_norm, what = "mra")

gsea1_MR <- gsea1[gsea1$Regulon %in% mra_kirc_prim_norm$Regulon,]
gsea2_diff_MR <- gsea2_diff[gsea2_diff$Regulon %in% mra_kirc_prim_norm$Regulon,]

mra_kirc_prim_norm[!(mra_kirc_prim_norm$Regulon %in% gsea1_MR$Regulon),]
mra_kirc_prim_norm[!(mra_kirc_prim_norm$Regulon %in% gsea2_diff_MR$Regulon),]

#### Saving Files ####
save(rtna_kirc_prim_norm, mra_kirc_prim_norm, file = "~/Kidney/0_data/MRA.RData")
write.csv(mra_kirc_prim_norm, file = "~/Kidney/0_data/MRA_TFs.csv")