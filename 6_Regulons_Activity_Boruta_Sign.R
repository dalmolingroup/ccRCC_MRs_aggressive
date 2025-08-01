#### Regulons Activity Analysis ####
## Date: 04/08/2024
## Author: Epitácio Farias 

#### Loading Packages ####
if(!require("RTN")){BiocManager::install("RTN")}
if(!require("RTNduals")){BiocManager::install("RTNduals")}
if(!require("RTNsurvival")){BiocManager::install("RTNsurvival")}
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("finalfit")){install.packages("finalfit")}
if(!require("survminer")){install.packages("survminer")}
if(!require("clusterProfiler")){BiocManager::install("clusterProfiler")}
if(!require("org.Hs.eg.db")){BiocManager::install("org.Hs.eg.db")}

#### Loading the data ####
selected_features <- read.table("~/Kidney/0_data/selected_features_boruta_balanced_2.txt","\n", header = F)
colnames(selected_features) <- "ENSEMBL_ID"
load("~/Kidney/0_data/rtn_TF_data.RData")
load("~/Kidney/0_data/MRA.RData")

selected_features <- bitr(selected_features$ENSEMBL_ID, fromType= "ENSEMBL", toType = c("SYMBOL","ENTREZID","ENSEMBL"),
                          OrgDb = org.Hs.eg.db, drop = T)
selected_features[11,] <- c("ENSG00000215271", "HOMEZ", "20164")

#### Retrieving the Regulons ####
targets_regs <- regulons_TF[names(regulons_TF) %in% selected_features$ENSEMBL]

extract_names <- function(my_list, vector_name) {
  if (vector_name %in% names(my_list)) {
    names(my_list[[vector_name]])
  } else {
    NULL
  }
}

for (i in 1:length(selected_features$ENSEMBL)){
  targets_regs[[i]] <- paste(extract_names(targets_regs, names(targets_regs[i])))
  names(targets_regs)[i] <- names(targets_regs)[i]
}

rm(regs, regulons_TF, extract_names, i)

#### Regulon Interactions ####
g <- tni.graph(rtni_TF, regulatoryElements = names(targets_regs))

g_all <- tni.graph(rtni_TF, regulatoryElements = regs$GENE_ID,tnet = "dpi", 
                   minRegulonSize = 15, gtype = "amapDend")

g_MR <- tni.graph(rtni_TF, tnet = "dpi", minRegulonSize = 15, gtype = "amapDend",
                  regulatoryElements = mra_kirc_prim_norm$Regulon)

g_r <- tni.graph(rtni_TF, tnet = "dpi", minRegulonSize = 15,
                 gtype = "amap",
                 regulatoryElements = names(targets_regs), amapFilter = "phyper",
                 amapCutoff = 1)

g_tal <- tni.graph(rtni_TF, tnet = "dpi", minRegulonSize = 15, gtype = "amapDend")

save(g, g_r, g_tal, g_MR, file = "~/Kidney/0_data/rtn_graph_Boruta_sign.RData")

#### Co-regulation analysis and inference of dual regulons ####
rmbr_TF <- tni2mbrPreprocess(rtni_TF)

rmbr_TF <- mbrAssociation(rmbr_TF)

mbrGet(rmbr_TF, what="summary")

## Filtering by the TF from the best model constructed
## Overlap Situation
overlap <- mbrGet(rmbr_TF, what="dualsOverlap")

overlap2 <- overlap[overlap$Regulon1 %in% selected_features$SYMBOL | overlap$Regulon2 %in% selected_features$SYMBOL, ]

## Correlation
correlation <- mbrGet(rmbr_TF, what="dualsCorrelation")
correlation2 <- correlation[correlation$Regulon1 %in% selected_features$SYMBOL | correlation$Regulon2 %in% selected_features$SYMBOL, ]

## Duals Relation
all(order(rownames(correlation2)) %in% order(rownames(overlap2)))
all(order(rownames(overlap2)) %in% order(rownames(correlation2)))

duals <- mbrGet(rmbr_TF, what="dualRegulons")

duals2 <- duals[duals %in% rownames(overlap2)]

save(rmbr_TF,correlation,correlation2, duals, duals2, file = "~/Kidney/0_data/coregulation_data_Boruta_sign.RData")

#### Survival Analysis ####
## Taming the survival data 
clin_data <- rtni_TF@colAnnotation

clin_data <- clin_data[,c(1,2,9,11,15:17,21,22)]

kirc_surv <- vroom::vroom("~/Kidney/0_data/TCGA-KIRC.survival.tsv")
kirc_surv <- kirc_surv[,c(1,3,2,4)]

kirc_pheno <- vroom::vroom("~/Kidney/0_data/TCGA-KIRC.GDC_phenotype.tsv.gz")
kirc_pheno <- kirc_pheno[,c(1,11,34,39:41,70,74:76,94,117)]

clin_GDC <- merge(kirc_pheno, kirc_surv, by.x = "submitter_id.samples", by.y = "sample", all.y = T)

clin_GDC <- clin_GDC[!duplicated(clin_GDC$submitter_id),]

clin_data <- merge(clin_data, clin_GDC, by.x = "sample", by.y = "submitter_id.samples", all.x = T)

clin_data[clin_data$sample == "TCGA-BP-4337-01A", 4] <- "2"

clin_data[clin_data$sample == "TCGA-B2-5635-01B", "OS"] <- 0
clin_data[clin_data$sample == "TCGA-B2-5633-01B", "OS"] <- 0
clin_data[clin_data$sample == "TCGA-B2-3923-01B", "OS"] <- 0
clin_data[clin_data$sample == "TCGA-B2-3924-01B", "OS"] <- 0
clin_data[clin_data$sample == "TCGA-B0-5104-01A", "OS"] <- 1
clin_data[clin_data$sample == "TCGA-CJ-5689-01A", "OS"] <- 1

clin_data[clin_data$sample == "TCGA-B2-5635-01B", "OS.time"] <- 755
clin_data[clin_data$sample == "TCGA-B2-5633-01B", "OS.time"] <- 963
clin_data[clin_data$sample == "TCGA-B2-3923-01B", "OS.time"] <- 992
clin_data[clin_data$sample == "TCGA-B2-3924-01B", "OS.time"] <- 1092

surv_data <- clin_data[,c(1:3,5:7,9,22,23)]
surv_data$ajcc_pathologic_m <- ifelse(surv_data$ajcc_pathologic_m == "0", "M0","M1")

rownames(surv_data) <- NULL
surv_data <- surv_data |> column_to_rownames("patient_ID")

table(surv_data$ajcc_pathologic_t)
table(surv_data$ajcc_pathologic_stage)

surv_data <- surv_data %>%
  mutate(ajcc_pathologic_stage = fct_collapse(ajcc_pathologic_stage, "G1" = "Stage I",
                                              "G2" = "Stage II", "G3" = "Stage III", "G4" = "Stage IV"))

surv_data <- surv_data %>%
  dplyr::rename("tumor_grade" = ajcc_pathologic_stage ,  "age" = age_at_index,
                "time" = OS.time, "event" = OS)

surv_data <- surv_data %>%
  mutate(ajcc_pathologic_t = as.character(ajcc_pathologic_t),
         ajcc_pathologic_n = as.character(ajcc_pathologic_n),
         ajcc_pathologic_m = as.character(ajcc_pathologic_m),
         tumor_grade = as.character(tumor_grade),
         event = as.character(event))

surv_data <- surv_data |> 
  pivot_wider(names_from = ajcc_pathologic_m, 
              values_from = ajcc_pathologic_m, values_fill = "0")

surv_data$M1 <- ifelse(surv_data$M1 == "M1", 1, 0)
surv_data$M0 <- ifelse(surv_data$M0 == "M0", 1, 0)

surv_data <- surv_data |> 
  pivot_wider(names_from = ajcc_pathologic_n, 
              values_from = ajcc_pathologic_n, values_fill = "0")

surv_data$NX <- ifelse(surv_data$NX == "NX", 1, 0)
surv_data$N1 <- ifelse(surv_data$N1 == "N1", 1, 0)
surv_data$N0 <- ifelse(surv_data$N0 == "N0", 1, 0)

surv_data <- surv_data |> 
  pivot_wider(names_from = ajcc_pathologic_t, 
              values_from = ajcc_pathologic_t, values_fill = "0")

surv_data$T1 <- ifelse(surv_data$T1 == "T1", 1, 0)
surv_data$T2 <- ifelse(surv_data$T2 == "T2", 1, 0)
surv_data$T3 <- ifelse(surv_data$T3 == "T3", 1, 0)
surv_data$T4 <- ifelse(surv_data$T4 == "T4", 1, 0)

surv_data <- surv_data |> 
  pivot_wider(names_from = tumor_grade, 
              values_from = tumor_grade, values_fill = "0")

surv_data$G1 <- ifelse(surv_data$G1 == "G1", 1, 0)
surv_data$G2 <- ifelse(surv_data$G2 == "G2", 1, 0)
surv_data$G3 <- ifelse(surv_data$G3 == "G3", 1, 0)
surv_data$G4 <- ifelse(surv_data$G4 == "G4", 1, 0)
surv_data$`NA` <- NULL

surv_data <- surv_data %>%
  mutate_if(is.character, as.integer)

surv_data$sample <- clin_data$patient_ID
surv_data <- surv_data |> column_to_rownames("sample")

surv_data <- surv_data[,c(3,2,1,4:16)] |> as.data.frame()

surv_data <- surv_data %>% 
  mutate(time = round(time/30.417, digit=0))

surv_data$group <- ifelse(surv_data$'1' == "0", "Primary", "Metastatic")

## Contructing the RTNsurvival object
rtns_TF <- tni2tnsPreprocess(rtni_TF, survivalData = surv_data, time = 3, event = 2)

rtns_TF <- tnsGSEA2(rtns_TF)

## Cox Analysis
rtns_TF <- tnsCox(rtns_TF)

regs <- tnsGet(rtns_TF,what = "regulatoryElements") |> as.data.frame()
colnames(regs)[1] <- "Regulon"

coxTab <- tnsGet(rtns_TF,what = "coxTable")

tiff(filename = "Kidney/Figures/cox_analysis_regs_boruta.tiff", height = 1000, width = 1000, res = 300)
tnsPlotCox(rtns_TF, regs = selected_features$SYMBOL)
dev.off()

## Kaplan-Meier Analysis
rtns_TF <- tnsKM(rtns_TF, sections = 1,)

attributes <- list(c("M0","M1"))

#,c("N0","N1","NX"),c("T1","T2","T3","T4"),c("G1","G2","G3","G4")

for (i in 1:length(selected_features$SYMBOL)) {
  tiff(filename = paste0("~/Kidney/Figures/Boruta_Sign_Figures/KM_Boruta_Sign/KM_analysis_",selected_features$SYMBOL[i],".tiff"), height = 1000, width = 2000, res = 300)
  tnsPlotKM(rtns_TF, regs = selected_features$SYMBOL[i],attribs = attributes)
  dev.off()
}

save(rtns_TF, file = "~/Kidney/0_data/survival_data_Boruta_Sign.RData")

#### Odds Ratio to Metastasis ####
expr_data <- as.data.frame(rtni_TF@gexp)

expr_data <- expr_data |> t() |> as.data.frame() |> rownames_to_column("Samples")


data_OR <- merge(expr_data,clin_data, by.x = "Samples", by.y = "patient_ID", all.x = T)

data_OR <- data_OR[,c(1,19946:19967,2:19945)]
data_OR <- data_OR |> column_to_rownames("Samples")
data_OR <- data_OR[,-c(1,3,7:20)]

explanatory = selected_features$ENSEMBL
  
dependent = 'ajcc_pathologic_m'

tiff(filename = "~/Kidney/Figures/OR_analysis_regsbm_Boruta_Sign.tiff", height = 1000, width = 1000, res = 300)
data_OR |> or_plot(dependent, explanatory)
dev.off()

#### Computing Regulons Activity ####
## Figures 
library("tidyverse")
if(!require("GDCRNATools")){BiocManager::install("GDCRNATools")}
library("RTN")
library("pheatmap")
library("patchwork")
library("paletteer")
load("~/Kidney/0_data/DEG_DESeq2.RData")
load("~/Kidney/0_data/rtni_TF.RData")
load("~/Kidney/0_data/ICGC_DEA_Data.RData")

selected_features <- vroom::vroom("~/Kidney/0_data/selected_features_boruta_balanced_2.txt",delim = "\t",col_names = F)

colnames(selected_features)[1] <- "ENSEMBL_ID"

###### HEATMAP TCGA-KIRC DATASET #####
# Compute regulon activity for individual samples
rtni_TF <- tni.gsea2(rtni_TF, regulatoryElements = selected_features$ENSEMBL_ID)
tcga_regact <- tni.get(rtni_TF, what = "regulonActivity")

c2 <- c("HOMEZ","ZNF658","PBX1","ZFP28","ZNF136")
annotation_row = data.frame(
  MRs = rownames(t(tcga_regact$dif))
)

annotation_row$cluster <- ifelse(annotation_row$MRs %in% c2, "C2","C1")
annotation_row <- annotation_row %>%  column_to_rownames("MRs")
my_colors <- paletteer_d("MexBrewer::Alacena")

# Get sample attributes from the 'rtni1st' dataset
kirc_annot <- tni.get(rtni_TF, "colAnnotation")
kirc_annot$ajcc_pathologic_m <- ifelse(kirc_annot$ajcc_pathologic_m == "1", "M1","M0")

# Get ER+/- and PAM50 attributes for pheatmap
attribs <- c("ajcc_pathologic_m","vital_status")
kirc_annot <- kirc_annot[,attribs]

tiff(filename = "~/Kidney/Figures/Boruta_Sign_Figures/Heatmap_Act_TCGA.tiff",
     res = 300, height = 1000, width = 3000)
pheatmap(t(tcga_regact$dif), 
         main="TCGA cohort (n=505 samples)",
         #annotation_col = kirc_annot, 
         show_colnames = FALSE, annotation_legend = T, 
         clustering_method = "ward.D2", fontsize_row = 12,
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         annotation_row = annotation_row,
         annotation_colors = list (cluster = c("C1" = "#1a9850", "C2" = "#542788")),
         cutree_rows = 2,
         color = my_colors, annotation_names_row = F)
dev.off()

## ICGC Cohort 
norm_ICGC_count <- DESeq2::vst(dds_unpaired) |> SummarizedExperiment::assay()

donor_data_final$patient <- paste0(donor_data_final$submitted_donor_id, donor_data_final$tissue_type)

uniques <- donor_data_final %>% group_by(submitted_donor_id) %>% tally() %>% filter(n==1) %>% pull(submitted_donor_id)
donor_data_sub <- donor_data_final[!(donor_data_final$submitted_donor_id %in% uniques),]
rownames(donor_data_sub) <-  NULL
donor_data_sub <- donor_data_sub |> column_to_rownames("patient")

norm_ICGC_count <- norm_ICGC_count[,colnames(norm_ICGC_count) %in% rownames(donor_data_sub)]

# Replace samples
rtni1st_icgcsamples <- tni.replace.samples(rtni_TF, norm_ICGC_count)

# Compute regulon activity for the new samples
rtni1st_icgcsamples <- tni.gsea2(rtni1st_icgcsamples, regulatoryElements = selected_features$ENSEMBL_ID)

icgc_regact <- tni.get(rtni1st_icgcsamples, what = "regulonActivity")

# Get sample attributes from the 'rtni1st_icgcsamples' dataset
donor_data_sub <- donor_data_sub %>% rownames_to_column("patient")

icgc_annot <- tni.get(rtni1st_icgcsamples, "colAnnotation")
icgc_annot <- merge(icgc_annot, donor_data_sub, by.x = "ID", by.y = "patient", all.x = T)

# Adjust PAM50 attributes for pheatmap
icgc_annot2 <- icgc_annot |> 
  pivot_wider(names_from = pathologic_n, 
              values_from = pathologic_n, values_fill = "0")

icgc_annot2$NX <- ifelse(icgc_annot2$NX == "NX", 1, 0)
icgc_annot2$N1 <- ifelse(icgc_annot2$N1 == "N1", 1, 0)
icgc_annot2$N0 <- ifelse(icgc_annot2$N0 == "N0", 1, 0)

icgc_annot2 <- icgc_annot2 |> 
  pivot_wider(names_from = pathologic_t, 
              values_from = pathologic_t, values_fill = "0")

icgc_annot2$T1 <- ifelse(icgc_annot2$T1 == "T1", 1, 0)
icgc_annot2$T2 <- ifelse(icgc_annot2$T2 == "T2", 1, 0)
icgc_annot2$T3 <- ifelse(icgc_annot2$T3 == "T3", 1, 0)
icgc_annot2$T4 <- ifelse(icgc_annot2$T4 == "T4", 1, 0)


icgc_annot2 <- icgc_annot2 |> 
  pivot_wider(names_from = pathologic_m, 
              values_from = pathologic_m, values_fill = "0")

icgc_annot2$MX <- ifelse(icgc_annot2$MX == "MX", 1, 0)
icgc_annot2$M1 <- ifelse(icgc_annot2$M1 == "M1", 1, 0)
icgc_annot2$M0 <- ifelse(icgc_annot2$M0 == "M0", 1, 0)

icgc_annot2 <- icgc_annot2 |> 
  pivot_wider(names_from = tissue_type, 
              values_from = tissue_type, values_fill = "0")

icgc_annot2$RT <- ifelse(icgc_annot2$RT == "RT", 1, 0)
icgc_annot2$RA <- ifelse(icgc_annot2$RA == "RA", 1, 0)

icgc_annot2 <- icgc_annot2 |> column_to_rownames("ID")

attribs <- c("obs.time","age","status","M0","M1","NX","N0","N1","T1","T2","T3","T4","RT","RA")

icgc_annot2 <- icgc_annot2[,attribs]

icgc_annot <- icgc_annot |> column_to_rownames("ID")

# Plot regulon activity profiles
tiff(filename = "~/Kidney/Figures/Boruta_Sign_Figures/Heatmap_Act_ICGC.tiff",
     res = 300, height = 1000, width = 3000)
pheatmap(t(icgc_regact$dif), 
               main="ICGC cohort (n=90 samples)",
               #annotation_col = icgc_annot[,c("status","pathologic_m")], 
               show_colnames = FALSE, annotation_legend = T, 
               clustering_method = "ward.D2", fontsize_row = 12,
               clustering_distance_rows = "manhattan",
               clustering_distance_cols = "manhattan",
         annotation_row = annotation_row,
         annotation_colors = list (cluster = c("C1" = "#1a9850", "C2" = "#542788")),
         cutree_rows = 2,
         color = my_colors, border_color = NA, annotation_names_row = F)
dev.off()

## Contructing the RTNsurvival object
icgc_annot2$status <- ifelse(icgc_annot2$status == "deceased", 1, 0)
rtns_TF_ICGC <- tni2tnsPreprocess(rtni1st_icgcsamples, survivalData = icgc_annot2, time = 1, event = 3)

rtns_TF_ICGC <- tnsGSEA2(rtns_TF_ICGC)

## Cox Analysis
rtns_TF_ICGC <- tnsCox(rtns_TF_ICGC)

regs_ICGC <- tnsGet(rtns_TF_ICGC,what = "regulatoryElements") |> as.data.frame()
colnames(regs_ICGC)[1] <- "Regulon"

coxTab_ICGC <- tnsGet(rtns_TF_ICGC,what = "coxTable")

ls_mr <- selected_features[selected_features$SYMBOL != c("HOMEZ", "ZNF658"),]

tiff(filename = "Kidney/Figures/cox_analysis_regs_boruta.tiff", height = 1000, width = 1000, res = 300)
tnsPlotCox(rtns_TF_ICGC, regs = ls_mr$SYMBOL)
dev.off()

## Kaplan-Meier Analysis
rtns_TF_ICGC <- tnsKM(rtns_TF_ICGC, sections = 1)

attributes <- list(c("0","1"),c("N0","N1","NX"),c("T1","T2","T3","T4"),
                   c("G1","G2","G3","G4"))

for (i in 1:length(ls_mr$SYMBOL)) {
  tiff(filename = paste0("~/Kidney/Figures/Boruta_Sign_Figures/KM_Boruta_Sign/KM_analysis_ICGC_",ls_mr$SYMBOL[i],".tiff"), height = 1000, width = 2000, res = 300)
  tnsPlotKM(rtns_TF_ICGC, regs = ls_mr$SYMBOL[i])
  dev.off()
}

save(rtns_TF_ICGC, file = "~/Kidney/0_data/survival_data_Boruta_Sign_ICGC.RData")