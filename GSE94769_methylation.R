library(ChAMP)
library(tidyverse)
library(ggplot2)
library(stringr)

# methylation data

Methy_THYM_Kajiura<-read_tsv(file="C:/Users/yuhon/Documents/Kajiura_thy_met/GSE94769_Average_Beta_Pval.txt",)
head(Methy_THYM_Kajiura)
Methy_THYM_Kajiura_beta<- Methy_THYM_Kajiura %>% dplyr::select("TargetID", str_subset(string = colnames(.),pattern = "AVG_Beta"))
Methy_THYM_Kajiura_det_p<- Methy_THYM_Kajiura %>% dplyr::select("TargetID", str_subset(string = colnames(.),pattern = "Detection.Pval"))
