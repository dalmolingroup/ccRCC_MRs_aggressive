#### Checking the DE Profile ####
## Author: Epitácio Farias
## Date: 04/11/2024

#### Load packages ####
if(!require("tidyverse")){install.packages("tidyverse")}

#### Load the data ####
selected_features <- vroom::vroom("~/Kidney/0_data/selected_features_boruta_balanced_2.txt",delim = "\t",col_names = "ENSEMBL_ID")
load("~/Kidney/0_data/rtn_TF_data.RData")
load("~/Kidney/0_data/DEG_intersec.RData")

#### Retrieving the Regulons ####

targets_regs <- regulons_TF[names(regulons_TF) %in% selected_features$ENSEMBL_ID]

extract_names <- function(my_list, vector_name) {
  if (vector_name %in% names(my_list)) {
    names(my_list[[vector_name]])
  } else {
    NULL
  }
}

for (i in 1:length(selected_features$ENSEMBL_ID)){
  targets_regs[[i]] <- paste(extract_names(targets_regs, names(targets_regs[i])))
  names(targets_regs)[i] <- names(targets_regs)[i]
}

rm(regs, regulons_TF,extract_names, i)

#### Tamimg the DE Data ####
DEG_genes_final <- DEG_genes_final %>% 
  mutate(
    padj = ifelse(is.na(padj), 1, padj),
    signif = 
      case_when(
        padj <= 0.05 & log2FoldChange < 0 ~ "Downregulated",
        padj <= 0.05 & log2FoldChange > 0 ~ "Upregulated",
        padj > 0.05 ~ "Not significant"
      )
  )

#### Separating the expression Patern from each regulon ####

table_fig <- data.frame(matrix(0, nrow=9, ncol=3))

colnames(table_fig) <- c("Regulon","Upregulated", "Downregulated")

for (i in 1:length(targets_regs)) {
  print(paste0("Quantity of target Up and Downregulated for the regulon of the TF:",names(targets_regs)[i]))
  reg <- DEG_genes_final[DEG_genes_final$gene_name %in% targets_regs[[i]],]
  table_fig$Regulon[i] <- names(targets_regs)[i]
  table_fig$Upregulated[i] <- sum(reg$signif == "Upregulated")
  table_fig$Downregulated[i] <- sum(reg$signif == "Downregulated")
  print(table(reg$signif))
}

table_fig <- table_fig |> pivot_longer(cols = c("Upregulated", "Downregulated"), names_to = "Significance" )

tiff(filename = "~/Kidney/Figures/targ_genes_up_down_Boruta_sign.tiff", res = 300, width = 2500, height = 1000)
ggplot(table_fig, aes(x = Regulon, y = value, fill = Significance)) + 
  geom_col(position = "dodge") +
  theme(axis.ticks=element_blank()) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))+
  ggtitle("Downregulated and Upregulated Target Genes in Each Regulon")+
  ylab("Number") +
  xlab("Regulons")
dev.off()