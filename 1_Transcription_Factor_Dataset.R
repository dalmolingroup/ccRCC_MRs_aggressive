#### Transcription Factor List ####
## Date: 09/14/2023
## Author: Epitácio Farias

#### Loading Librarys ####
library("tidyverse")

#### Dowloading the data ####
url <- "http://bioinfo.life.hust.edu.cn/AnimalTFDB4/static/download/TF_list_final/Homo_sapiens_TF"
destfile <- "Homo_sapiens_TF.txt"
download.file(url, destfile)

#### Reading the data ####
tf_list <- read_delim("Homo_sapiens_TF.txt", delim = "\t")

save(tf_list, file = "tf_list.RData")