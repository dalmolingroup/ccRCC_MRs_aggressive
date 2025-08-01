#####################################
##### Constructing the Datasets #####
##### Author: Farias, Epitácio  #####
##### Date: 07/12/2024          #####

## Install packages
pkg.bioconductor <- c("TCGAbiolinks", "TCGAWorkflow")
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

#### TCGA-KIRC Data set ######
#TCGA: primary and normal samples
query_kirc <- GDCquery(project = "TCGA-KIRC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts",
                       sample.type = c("Primary Tumor"))
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
view(clin_data)

## 2. Cleaning data 

# Cleaning variables 
clin_data %>% missing_glimpse()
skim(clin_data) 

logical.vars <- names(clin_data |> select_if(is.logical))
list.vars <- names(clin_data |> select_if(is.list))

clin_data <- clin_data %>% 
  select(-all_of(c(logical.vars,list.vars)))

## Remove duplicate observations (patient_id or other id variable)
clin_data <- clin_data %>%
  distinct_at('sample', .keep_all = TRUE)

## Remove numeric variables with unique observations (or <50% complete)
# need to maintain 'days_to_death' for survival anaysis?
clin_data %>%
  select_if(is.numeric) %>%
  skim()

clin_data <- clin_data  %>%
  select(!c('days_to_collection', 'initial_weight', 'cigarettes_per_day', 'pack_years_smoked',
            'days_to_death','years_smoked','year_of_death','days_to_death', 'paper_mRNA_cluster',
            'paper_microRNA_cluster'))

## Remove character variables with unique observations 
clin_data %>%
  select_if(is.character) %>%
  skim()

clin_data <- clin_data  %>%
  select(!c('alcohol_history', 'ajcc_clinical_m', 
            'ajcc_staging_system_edition', 
            'progression_or_recurrence', 'tumor_grade', 
            'site_of_resection_or_biopsy', 'icd_10_code', 
            'classification_of_tumor', 'morphology', 'tissue_or_organ_of_origin',
            'last_known_disease_status','specimen_type','oct_embedded','preservation_method',
            'state','composition',))

## Remove variables with similar information - check each one!
clin_data$submitter_id == clin_data$bcr_patient_barcode
clin_data$age_at_index == as.integer(clin_data$age_at_diagnosis/365)
clin_data$bcr_patient_barcode == clin_data$sample_submitter_id
clin_data$sample_submitter_id == clin_data$sample
clin_data$definition == clin_data$tumor_descriptor
clin_data$sample_type == clin_data$tissue_type

clin_data <- clin_data  %>%
  select(!c('submitter_id', 'age_at_diagnosis','bcr_patient_barcode', 
            'sample_submitter_id','shortLetterCode','definition', 
            'tumor_descriptor','tissue_type','barcode', 'paper_patient',
            'pathology_report_uuid','patient','primary_diagnosis')) 

## Check and remove other variables
names.id <- names(clin_data %>%
                    select(ends_with("_id")))
names.id <- names.id[-1]

clin_data <- clin_data %>% 
  select(-all_of(names.id))


## 4. Taming data 
clin_data <- clin_data %>%
  mutate_if(is.character, as.factor) %>%
  mutate(sample = as.character(sample))

## 5. Checking NA patterns 
# Check distincts types of NAs: MCAR, MAR, MNAR
clin_data  %>%
  missing_plot()


## 6. Checking numeric variables 
# Check data distribution, unplausible values, outliers.
# Never delete an unusual value if the value is a possible one. 
# Deleting unusual values will bias your results and cause you to underestimate the variability in the observations.
# If the median and mean are similar, the distribution is likely roughly symmetrical. 
# Otherwise, it will be skewed to the right or to the left.

clin_data %>% 
  select_if(is.numeric) %>%
  summary()

# Histograms or density plots
ggplot(clin_data, aes(age_at_index)) +
  geom_histogram(bins = 20, alpha = 0.5, color = "red")

ggplot(clin_data, aes(year_of_birth)) +
  geom_histogram(bins = 20, alpha = 0.5, color = "red")

ggplot(clin_data, aes(year_of_diagnosis)) +
  geom_bar(stat= "count", alpha = 0.5, color = "blue")

ggplot(clin_data, aes(days_to_last_follow_up)) +
  geom_density( alpha = 0.5, color = "purple")

ggplot(clin_data, aes(x ='', y=age_at_index)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(clin_data$age_at_index)

ggplot(clin_data, aes(x ='', y=days_to_last_follow_up)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
boxplot.stats(clin_data$days_to_last_follow_up)

## 7. Checking categorical variables
# Check frequency, lables and levels 
# Cancer staging: https://www.cancer.gov/about-cancer/diagnosis-staging/staging
clin_data %>%
  select_if(is.factor) %>%
  summary() 

# agregating levels

clin_data <- clin_data %>%
  mutate(ajcc_pathologic_t = fct_collapse(ajcc_pathologic_t, T1=c('T1', 'T1a', 'T1b'), T2=c('T2','T2a', 'T2b'), 
                                          T3=c('T3a', 'T3b','T3c')))

## Graphics

# categorical frequencies graphics
clin_data %>% 
  dplyr::count(ajcc_pathologic_m) %>% 
  knitr::kable()


clin_data %>% 
  ggplot(aes(x = ajcc_pathologic_m, y = age_at_index, fill = ajcc_pathologic_m)) + 
  geom_boxplot() +
  theme_bw(15) +
  labs(title = 'Tissue Site', x = "tissue site", y = "age") + 
  facet_wrap(~ gender) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')



## 8. Saving dataset 
skim(clin_data)
write_csv(clin_data, "clin_data.csv")


## 9. Summary table with tableone
clin_data2 <- clin_data %>% 
  select(!where((is.character))) %>% 
  # select(!(updated_datetime)) %>% 
  mutate(days_to_last_follow_up = as.integer(days_to_last_follow_up),
         year_of_diagnosis = as.integer(year_of_diagnosis),
         year_of_birth = as.integer(year_of_birth))

CreateTableOne(data=clin_data2)

myVars <- names(clin_data2)
catVars <- names(clin_data2 %>% 
                   select_if(is.factor))
biomarkers <- c("days_to_last_follow_up", "year_of_diagnosis", "year_of_birth")

clin_tab <- CreateTableOne(vars=myVars, data=clin_data2, factorVars=catVars)
print(clin_tab, nonnormal=biomarkers, showAllLevels=TRUE, formatOptions=list(big.mark = ","))


## 10. Survival analysis 
clin_data3 <- clin_data

clin_data3 <- clin_data3[clin_data3$ajcc_pathologic_m != "MX",]

clin_data3 <- clin_data3[na.omit(clin_data3$sample),]

# Survival pkg
# Dichotomize and change data labels
clin_data3$ajcc_pathologic_m <- as.numeric(ifelse(clin_data3$ajcc_pathologic_m == "M0", 0, 1))
head(clin_data3$ajcc_pathologic_m)

# Fit survival data using the Kaplan-Meier method
# Time records survival time. Status indicates if patient is death (status=1) or was censored (status=0).
# A “+” after the time in the print out of km indicates censoring.
surv_object <- Surv(time = as.numeric(clin_data3$days_to_last_follow_up), event = clin_data3$ajcc_pathologic_m)
surv_object 

fit <- survfit(surv_object ~ race, data = clin_data3)
fit
summary(fit)
# summary(fit1, times = c(1,30,60,90*(1:10)))
summary(fit)$table

ggsurvplot(fit, data = clin_data3, pval = TRUE)
ggsurvplot(fit, data = clin_data3,
           conf.int = TRUE,
           pval = TRUE,
           fun = "pct",
           risk.table = TRUE,
           size = 1,
           linetype = "strata",
           palette = c("#E7B800", "#2E9FDF","pink","red"),
           legend = "bottom",
           legend.title = "Race",
           legend.labs = c("Not Reported", "White","Asian","Black of African American"))

# Fit a Cox proportional hazards model # ERROR
fit.coxph <- coxph(surv_object ~ gender + prior_malignancy, data = clin_data3)
ggforest(fit.coxph, data = clin_data3)


#### Transcription  TCGA-KIRC Data Adjustment ####
## Retrieving the sequencing Gene Library 
gene_library <- as.data.frame(KIRC@rowRanges@elementMetadata@listData)

gene_library <- gene_library[,c(5:7)]

## Retrieving the Raw Count Data
count_data_kirc <- assay(KIRC) |> as.data.frame() |> rownames_to_column("ENSEMBL_ID")

count_data_kirc <- merge(count_data_kirc, gene_library, by.x = "ENSEMBL_ID", by.y = "gene_id", all.x = T)

count_data_kirc <- count_data_kirc[,c(1,543,544,2:542)]

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
count_data_kirc2 <- count_data_kirc[count_data_kirc$gene_type == "protein_coding", ]

#### Saving the Data ####
save(clin_data,clin_data2,clin_data3, count_data_kirc,count_data_kirc2,gene_library, file = "~/Kidney/0_data/KIRC_data_updated.RData")
rm(list = ls())

#### ICGC-RECA Data set ####

url <- "https://dcc.icgc.org/api/v1/download?fn=/current/Projects/RECA-EU/exp_seq.RECA-EU.tsv.gz"
destfile <- "~/Kidney/0_data/exp_seq.RECA-EU.tsv.gz"
download.file(url, destfile)

url <- "https://dcc.icgc.org/api/v1/download?fn=/current/Projects/RECA-EU/donor.RECA-EU.tsv.gz"
destfile <- "~/Kidney/0_data/donor.RECA-EU.tsv.gz"
download.file(url, destfile)

donor_data <- vroom("~/Kidney/0_data/donor.RECA-EU.tsv.gz")
count_data <- vroom("~/Kidney/0_data/exp_seq.RECA-EU.tsv.gz")

#### Transcription ICGC-RECA Dataset Adjustment ####
count_data2 <- count_data %>% 
  dplyr::select(c("submitted_sample_id", "gene_id", "raw_read_count")) %>%
  tidyr::pivot_wider(names_from = "submitted_sample_id", values_from = "raw_read_count") %>% as.data.frame()

count_data3 <- merge(count_data2, gene_library, by.x = "gene_id", by.y = "gene_id")

count_data3 <- count_data3[,c(1,138:139,2:137)]

sum(duplicated(count_data3$gene_id))

count_data4 <- count_data3[!duplicated(count_data3$gene_id),]

## Final Count Dataset
count_final_ICGC <- count_data4[count_data4$gene_type =="protein_coding",]

#### Clinical ICGC-RECA Dataset Adjustment ####
exp_patients <- substr(colnames(count_final_ICGC), 1,5) %>% as.data.frame()
colnames(exp_patients)[1] <- "submitted_donor_id"
exp_patients$tissue_type <- substr(colnames(count_final_ICGC), 6,7)
exp_patients <- exp_patients[-c(1:3),]

donor_data2 <- donor_data %>%
  dplyr::select(c("submitted_donor_id", "donor_survival_time", "donor_vital_status", "donor_age_at_diagnosis", "donor_sex",  "donor_tumour_stage_at_diagnosis" )) %>% 
  dplyr::rename(obs.time = 'donor_survival_time',
                status = 'donor_vital_status',
                age = 'donor_age_at_diagnosis', 
                gender = 'donor_sex')

donor_data2$pathologic_t <- substr(donor_data2$donor_tumour_stage_at_diagnosis, 1, 2)
donor_data2$pathologic_n <- substr(donor_data2$donor_tumour_stage_at_diagnosis, 3, 4)
donor_data2$pathologic_m <- substr(donor_data2$donor_tumour_stage_at_diagnosis, 5, 6)
donor_data2$donor_tumour_stage_at_diagnosis <- NULL

donor_data3 <- merge(donor_data2, exp_patients, by.x = "submitted_donor_id", by.y = "submitted_donor_id")

donor_data_final <- donor_data3[,c(1,9,2:8)]

## Saving the data
save(count_data, count_data2,count_data3,count_data4,donor_data,donor_data2,donor_data3, file = "~/Kidney/0_data/ICGC_data.RData")
save(count_final_ICGC, donor_data_final, file = "~/Kidney/0_data/ICGC_DEA_Data.RData")