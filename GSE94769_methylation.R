library(ChAMP)
library(tidyverse)
library(ggplot2)
library(stringr)

# methylation data

Methy_THYM_Kajiura<-read_tsv(file="C:/Users/yuhon/Documents/Kajiura_thy_met/GSE94769_Average_Beta_Pval.txt",) 

Methy_THYM_Kajiura_det_p<- Methy_THYM_Kajiura %>% 
  dplyr::select("TargetID", str_subset(string = colnames(.),pattern = "Detection.Pval"))

head(Methy_THYM_Kajiura) 
Methy_THYM_Kajiura_beta<- Methy_THYM_Kajiura %>% 
  dplyr::select("TargetID", str_subset(string = colnames(.),pattern = "AVG_Beta"))


colnames(Methy_THYM_Kajiura_beta) =  
  str_remove(string = colnames(Methy_THYM_Kajiura_beta),pattern=".AVG_Beta")

thymoma_ID<-Methy_THYM_Kajiura_beta

Methy_Thymoma_Kajiura_beta<- select ()

# meta data
meta_data_thy_kajiura<-read_tsv(file="C:/Users/yuhon/Documents/Kajiura_thy_met/meta_data.txt")

meta_data_thy_kajiura$Masaoka_classification = str_remove(string = meta_data_thy_kajiura$Masaoka_classification,pattern="Stage: ")
meta_data_thy_kajiura$histology_WHO_classification = str_remove(string = meta_data_thy_kajiura$histology_WHO_classification,pattern="tissue: ")
meta_data_thy_kajiura$tissue = str_remove(string = meta_data_thy_kajiura$tissue, pattern="tissue: ")

# thymoma ID list
thymoma_ID<-meta_data_thy_kajiura$Sample[meta_data_thy_kajiura$histology_WHO_classification == "B3 thymoma" |meta_data_thy_kajiura$histology_WHO_classification == "type A thymoma"]

Methy_Thymoma_Kajiura_beta<-Methy_THYM_Kajiura_beta %>% dplyr::select("TargetID",thymoma_ID)
meta_data_thymoma_kajiura<-meta_data_thy_kajiura %>% filter(Sample %in% thymoma_ID)



## methy that has data
met_THYM_gse<- Methy_Thymoma_Kajiura_beta %>% dplyr::select("TargetID",meta_data_thymoma_kajiura$Sample)


met_THYM_m <- met_THYM_gse %>% dplyr::select(-TargetID) %>% as.matrix()

row.names(met_THYM_m)<-met_THYM_gse$TargetID

met_THYM_m<-na.omit(met_THYM_m)




## methy that has

##filter
myload_met_gse <- champ.filter(beta=met_THYM_m,pd=meta_data_thymoma_kajiura,fixOutlier = T,filterBeads = F,filterDetP = F,autoimpute = F)

##quality control
#champ.QC(beta=myload_met$beta,pheno =myload_met$pd$patient_id)
## normalization
myNorm_met <- champ.norm(beta = myload_met_gse$beta,arraytype="450K", cores=8)

write.csv(myNorm_met,"gse94769_met_normalization_update_MG_status_20220805.csv")

# myNorm_met<-read.csv("tcga_met_normalization_update_MG_status_20220716.csv",check.names = F)
row.names(myNorm_met)<-myNorm_met$gene_id
myNorm_met<-myNorm_met[,-1]

myDMP <- champ.DMP(beta = myNorm_met,pheno=meta_data_thymoma_kajiura$Myasthenia_gravis,adjPVal = 1)

head(myload_met$beta)
head(myNorm_met)


myDMP$no_to_yes %>% filter(gene == "TTN")
myDMP$no_to_yes %>% filter(gene == "CHRNA1")
myDMP$no_to_yes %>% filter(gene == "NEFM")
myDMP$no_to_yes %>% filter(gene == "RYR1")
myDMP$no_to_yes %>% filter(gene == "RYR2")
myDMP$no_to_yes %>% filter(gene == "RYR3")
myDMP$no_to_yes %>% filter(gene == "CHRND")
myDMP$no_to_yes %>% filter(gene == "CHRNG")
myDMP$no_to_yes %>% filter(gene == "CHRNE")
myDMP$no_to_yes %>% filter(gene == "CHRNB1")
#DMPs 分析
myDMP <- champ.DMP(beta = myNorm_met,pheno=met_meta$MG,adjPVal = 1)
