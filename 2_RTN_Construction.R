#### Regulatory Networks Reconstruction ####
## Author: Epitácio Farias
## Date: 09/08/2023 (MM/DD/YY)

BiocManager::install("RTN")
library("RTN")
library("tidyverse")

#### Genes in the cohort ####
load("~/Kidney/0_data/KIRC_data_updated.RData")

#### 3. RTN ####
ColAnnot <- clin_data3 |> rownames_to_column("patient_ID")
RowAnnot <- gene_library[gene_library$gene_id %in% rownames(count_data_kirc2),]
names(RowAnnot)[3] <- "SYMBOL"

#Regulatory elements
load("~/Kidney/0_data/tf_list.RData")
tf_list <- tf_list[!is.na(tf_list$Ensembl),]
tf_ensembl <- tf_list$Ensembl

#Normalizing the expression data 
exp_final <- GDCRNATools::gdcVoomNormalization(count_data_kirc2, filter = F)
exp_final <- exp_final[, colnames(exp_final) %in% ColAnnot$patient_ID]

#rtni object
rtni_constructor_TF <- tni.constructor(expData = exp_final,
                                       regulatoryElements = tf_ensembl,
                                       colAnnotation = ColAnnot,
                                       rowAnnotation = RowAnnot)

save(rtni_constructor_TF, file = "~/Kidney/0_data/rtni_constructor_TF.RData")

library("parallel")
options(cluster=snow::makeCluster(spec=32, "SOCK"))

#permutation
rtni_TF <- tni.permutation(rtni_constructor_TF,  nPermutations = 1000)
save(rtni_TF, file = "rtni_TF_permu.RData")

#bootstrap
rtni_TF <- tni.bootstrap(rtni_TF)
save(rtni_TF, file = "rtni_TF_boot.RData")

#ARACNe
rtni_TF <- tni.dpi.filter(rtni_TF)
save(rtni_TF, file = "rtni_TF_aracne.RData")

stopCluster(getOption("cluster"))

save(rtni_TF, file = "~/Kidney/0_data/rtni_TF.RData")

######### Summary
load(file = "~/Kidney/0_data/rtni_TF.RData")

tni.regulon.summary(rtni_TF)
regulons_TF <- tni.get(rtni_TF, what = "regulons.and.mode")
regs <- tni.get(rtni_TF, what = "regulonSize", idkey = "GENE_ID")

#count target
count <- 0
for (i in 1:length(regulons_TF)){
  if (length(regulons_TF[[i]]) >= 15){
    count <- count + 1
  }
}

save(rtni_TF, regulons_TF,regs, file = "~/Kidney/0_data/rtn_TF_data.RData")
