#library(devtools)
devtools::install_github("dalpozz/unbalanced")
suppressMessages(library(unbalanced))
library("tidyverse")

mra_tf = vroom::vroom("Kidney/0_data/MRA_TFs.csv")

metadata = vroom::vroom("Kidney/0_data/metadata_KIRC.csv")

load("Kidney/0_data/rtni_constructor_TF.RData")

count_data = rtni_constructor_TF@gexp %>% as.data.frame()

count_data <- count_data %>%  t() %>%  as.data.frame()

count_data <- count_data %>%  rownames_to_column("Sample")

metadata_icgc = vroom::vroom("Kidney/0_data/meta_icgc_data.csv")

count_data_icgc = vroom::vroom("Kidney/0_data/count_icgc_data.csv")

count_data_icgc <- GDCRNATools::gdcVoomNormalization(count_data_icgc, filter = F)

## Balancing TCGA dataset
predictor_variables <- metadata[,-17] # Select everything except response
response_variable <- metadata$ajcc_pathologic_m |> as.factor()

levels(response_variable) <- c('0', '1') 

undersampled_data <- ubBalance(predictor_variables, 
                               response_variable, 
                               type='ubUnder',         # Option for undersampling
                               verbose = TRUE)

metadata_balanced <- cbind(undersampled_data$X,    # combine output
                           undersampled_data$Y)

names(metadata_balanced)[names(metadata_balanced) == "undersampled_data$Y"] <- "ajcc_pathologic_m" # change name to class
levels(metadata_balanced$ajcc_pathologic_m) <- c('M0', 'M1')

count_data_balanced <- count_data[count_data$Sample %in% rownames(metadata_balanced),]
rownames(count_data_balanced) <-  NULL

count_data_balanced <- count_data_balanced |> column_to_rownames("Sample")

count_data_balanced <- count_data_balanced[,colnames(count_data_balanced) %in% mra_tf$...1]

metadata_balanced <- metadata_balanced |> column_to_rownames("...1")

write.csv(metadata_balanced, file = "metadata_KIRC_balanced.csv")
write.csv(count_data_balanced, file = "count_data_KIRC_pc_balanced.csv")

## Balancing ICGC dataset
predictor_variables <- metadata_icgc[,-9] # Select everything except response
response_variable <- metadata_icgc$pathologic_m |> as.factor()

levels(response_variable) <- c('0', '1') 

undersampled_data <- ubBalance(predictor_variables, 
                               response_variable, 
                               type='ubUnder',         # Option for undersampling
                               verbose = TRUE)

metadata_icgc_balanced <- cbind(undersampled_data$X,    # combine output
                                undersampled_data$Y)

names(metadata_icgc_balanced)[names(metadata_icgc_balanced) == "undersampled_data$Y"] <- "pathologic_m" # change name to class
levels(metadata_icgc_balanced$pathologic_m) <- c('M0', 'M1')

library("tidyverse")

count_data_icgc <- count_data_icgc |> column_to_rownames("...1")

count_data_icgc_balanced <- count_data_icgc[,colnames(count_data_icgc) %in% metadata_icgc_balanced$patient]

count_data_icgc_balanced <- count_data_icgc_balanced[rownames(count_data_icgc_balanced) %in% mra_tf$...1,]

metadata_icgc_balanced <- metadata_icgc_balanced[,-1] |> column_to_rownames("patient")

write.csv(metadata_icgc_balanced, file = "metada_icgc_balanced.csv")
write.csv(count_data_icgc_balanced, file = "count_data_icgc_balanced.csv")
