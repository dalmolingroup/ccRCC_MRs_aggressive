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

surv_data_m <- surv_data[surv_data$ajcc_pathologic_m == "1",]
surv_data_wm <- surv_data[surv_data$ajcc_pathologic_m != "1",]

## Metastatic Samples

surv_data_m <- surv_data_m |> 
  pivot_wider(names_from = ajcc_pathologic_m, 
              values_from = ajcc_pathologic_m, values_fill = "0")

surv_data_m <- surv_data_m |> 
  pivot_wider(names_from = ajcc_pathologic_n, 
              values_from = ajcc_pathologic_n, values_fill = "0")

surv_data_m$NX <- ifelse(surv_data_m$NX == "NX", 1, 0)
surv_data_m$N1 <- ifelse(surv_data_m$N1 == "N1", 1, 0)
surv_data_m$N0 <- ifelse(surv_data_m$N0 == "N0", 1, 0)

surv_data_m <- surv_data_m |> 
  pivot_wider(names_from = ajcc_pathologic_t, 
              values_from = ajcc_pathologic_t, values_fill = "0")

surv_data_m$T1 <- ifelse(surv_data_m$T1 == "T1", 1, 0)
surv_data_m$T2 <- ifelse(surv_data_m$T2 == "T2", 1, 0)
surv_data_m$T3 <- ifelse(surv_data_m$T3 == "T3", 1, 0)
surv_data_m$T4 <- ifelse(surv_data_m$T4 == "T4", 1, 0)

surv_data_m <- surv_data_m |> 
  pivot_wider(names_from = tumor_grade, 
              values_from = tumor_grade, values_fill = "0")

# surv_data_m$G1 <- ifelse(surv_data_m$G1 == "G1", 1, 0)
# surv_data_m$G2 <- ifelse(surv_data_m$G2 == "G2", 1, 0)
# surv_data_m$G3 <- ifelse(surv_data_m$G3 == "G3", 1, 0)
surv_data_m$G4 <- ifelse(surv_data_m$G4 == "G4", 1, 0)
surv_data_m$`NA` <- NULL

clin_data_m <- clin_data[clin_data$sample %in% surv_data_m$sample,]

surv_data_m$sample <- clin_data_m$patient_ID

surv_data_m <- surv_data_m |> column_to_rownames("sample")

surv_data_m <- surv_data_m[,c(3,2,1,4:12)] |> as.data.frame()

surv_data_m <- surv_data_m %>% 
  mutate(time = round(time/30.417, digit=0))

surv_data_m <- surv_data_m %>%
  mutate_if(is.character, as.integer)

## Contructing the RTNsurvival object
rtns_TF_m <- tni2tnsPreprocess(rtni_TF, survivalData = surv_data_m, time = 1, event = 2)

rtns_TF_m <- tnsGSEA2(rtns_TF_m)

## Cox Analysis
rtns_TF_m <- tnsCox(rtns_TF_m)

regs <- tnsGet(rtns_TF_m,what = "regulatoryElements") |> as.data.frame()
colnames(regs)[1] <- "Regulon"

coxTab <- tnsGet(rtns_TF_m,what = "coxTable")

tiff(filename = "Kidney/Figures/cox_analysis_regs_boruta.tiff", height = 1000, width = 1000, res = 300)
tnsPlotCox(rtns_TF_m, regs = selected_features$SYMBOL)
dev.off()

## Kaplan-Meier Analysis
rtns_TF_m <- tnsKM(rtns_TF_m, sections = 1)

attributes <- list(c("0","1"),c("N0","N1","NX"),c("T1","T2","T3","T4"),
                   c("G1","G2","G3","G4"))

for (i in 1:length(selected_features$SYMBOL)) {
  tiff(filename = paste0("~/Kidney/Figures/Boruta_Sign_Figures/KM_Boruta_Sign/KM_analysis_MS_",selected_features$SYMBOL[i],".tiff"), height = 1000, width = 2000, res = 300)
  tnsPlotKM(rtns_TF_m, regs = selected_features$SYMBOL[i])
  dev.off()
}

save(rtns_TF_m, file = "~/Kidney/0_data/survival_data_Boruta_Sign_metastatic_samples.RData")

## Without Metastasis Samples
surv_data_wm <- surv_data_wm |> 
  pivot_wider(names_from = ajcc_pathologic_m, 
              values_from = ajcc_pathologic_m, values_fill = "0")

surv_data_wm <- surv_data_wm |> 
  pivot_wider(names_from = ajcc_pathologic_n, 
              values_from = ajcc_pathologic_n, values_fill = "0")

surv_data_wm$NX <- ifelse(surv_data_wm$NX == "NX", 1, 0)
surv_data_wm$N1 <- ifelse(surv_data_wm$N1 == "N1", 1, 0)
surv_data_wm$N0 <- ifelse(surv_data_wm$N0 == "N0", 1, 0)

surv_data_wm <- surv_data_wm |> 
  pivot_wider(names_from = ajcc_pathologic_t, 
              values_from = ajcc_pathologic_t, values_fill = "0")

surv_data_wm$T1 <- ifelse(surv_data_wm$T1 == "T1", 1, 0)
surv_data_wm$T2 <- ifelse(surv_data_wm$T2 == "T2", 1, 0)
surv_data_wm$T3 <- ifelse(surv_data_wm$T3 == "T3", 1, 0)
surv_data_wm$T4 <- ifelse(surv_data_wm$T4 == "T4", 1, 0)

surv_data_wm <- surv_data_wm |> 
  pivot_wider(names_from = tumor_grade, 
              values_from = tumor_grade, values_fill = "0")

surv_data_wm$G1 <- ifelse(surv_data_wm$G1 == "G1", 1, 0)
surv_data_wm$G2 <- ifelse(surv_data_wm$G2 == "G2", 1, 0)
surv_data_wm$G3 <- ifelse(surv_data_wm$G3 == "G3", 1, 0)
surv_data_wm$G4 <- ifelse(surv_data_wm$G4 == "G4", 1, 0)
surv_data_wm$`NA` <- NULL

clin_data_wm <- clin_data[clin_data$sample %in% surv_data_wm$sample,]

surv_data_wm$sample <- clin_data_wm$patient_ID

surv_data_wm <- surv_data_wm |> column_to_rownames("sample")

surv_data_wm <- surv_data_wm[,c(3,2,1,4:15)] |> as.data.frame()

surv_data_wm <- surv_data_wm %>% 
  mutate(time = round(time/30.417, digit=0))

surv_data_wm <- surv_data_wm %>%
  mutate_if(is.character, as.integer)

## Contructing the RTNsurvival object
rtns_TF_wm <- tni2tnsPreprocess(rtni_TF, survivalData = surv_data_wm, time = 1, event = 2)

rtns_TF_wm <- tnsGSEA2(rtns_TF_wm)

## Cox Analysis
rtns_TF_wm <- tnsCox(rtns_TF_wm)

regs <- tnsGet(rtns_TF_wm,what = "regulatoryElements") |> as.data.frame()
colnames(regs)[1] <- "Regulon"

coxTab <- tnsGet(rtns_TF_wm,what = "coxTable")

tiff(filename = "Kidney/Figures/cox_analysis_regs_boruta.tiff", height = 1000, width = 1000, res = 300)
tnsPlotCox(rtns_TF_wm, regs = selected_features$SYMBOL)
dev.off()

## Kaplan-Meier Analysis
rtns_TF_wm <- tnsKM(rtns_TF_wm, sections = 1)

attributes <- list(c("0","1"),c("N0","N1","NX"),c("T1","T2","T3","T4"),
                   c("G1","G2","G3","G4"))

for (i in 1:length(selected_features$SYMBOL)) {
  tiff(filename = paste0("~/Kidney/Figures/Boruta_Sign_Figures/KM_Boruta_Sign/KM_analysis_WMS_",selected_features$SYMBOL[i],".tiff"), height = 1000, width = 2000, res = 300)
  tnsPlotKM(rtns_TF_wm, regs = selected_features$SYMBOL[i])
  dev.off()
}

save(rtns_TF_wm, file = "~/Kidney/0_data/survival_data_Boruta_Sign_non_metastatic_samples.RData")
