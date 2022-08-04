library("ChAMP")
library(tidyverse)
#methylation
Methy_THYM<-read_tsv(file="C:/Users/yuhon/Documents/TCGApaper/GDCdata/TCGA-THYM/legacy/TCGA.THYM.sampleMap_HumanMethylation450/HumanMethylation450",
)
colnames(Methy_THYM)[2:length(Methy_THYM)]<-substr(colnames(Methy_THYM)[2:length(Methy_THYM)],1,12)
colnames(Methy_THYM)

head(Methy_THYM)


Methy_annotation<-read_csv(file="C:/Users/yuhon/Documents/TCGApaper/TCGA database/数据/humanmethylation450_15017482_v1_2_trimmed.csv")
head(Methy_annotation)

View(Methy_annotation[1:100,])


##clinical data
clinical_THYM<-read_tsv(file="./TCGApaper/GDCdata/TCGA-THYM/harmonized/Clinical/Clinical_Supplement/7cca5722-26cf-4ac8-a4a6-b803459f1861/nationwidechildrens.org_clinical_patient_thym_2.txt",
)
clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-3G-AB14"]<-"NO"
clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-X7-A8DF"]<-"NO"
#View(clinical_THYM)
meta_data<-clinical_THYM %>% dplyr::select(patient_id=bcr_patient_barcode,
                                    gender,height,weight,
                                    histo_type = histological_type,
                                    MG = history_myasthenia_gravis,
                                    section_MG = section_myasthenia_gravis,
                                    age_patho_diagnosis = 
                                      age_at_initial_pathologic_diagnosis) %>%
  filter(MG == "YES"|MG =="NO" ) %>%
  arrange(patient_id)

met_meta<-meta_data %>% filter(meta_data$patient_id %in% (colnames(Methy_THYM)[-1]))

## methy that has data
met_THYM<- Methy_THYM %>% dplyr::select("sample",met_meta$patient_id)

met_THYM_m <- met_THYM %>% dplyr::select(-sample) %>% as.matrix()
row.names(met_THYM_m)<-met_THYM$sample

met_THYM_m<-na.omit(met_THYM_m)

head(met_THYM_m)

View(met_THYM_m)





## methy that has

  ##filter
#myload_met <- champ.filter(beta=met_THYM_m,pd=met_meta,fixOutlier = T,filterBeads = F,filterDetP = F,autoimpute = F)

  ##quality control
#champ.QC(beta=myload_met$beta,pheno =myload_met$pd$patient_id)
  ## normalization
#myNorm_met <- champ.norm(beta = myload_met$beta,arraytype="450K", cores=8)

#write.csv(myNorm_met,"tcga_met_normalization_update_MG_status_20220716.csv")

myNorm_met<-read.csv("tcga_met_normalization_update_MG_status_20220716.csv",check.names = F)
row.names(myNorm_met)<-myNorm_met$gene_id
myNorm_met<-myNorm_met[,-1]

head(myload_met$beta)
head(myNorm_met)

#DMPs 分析
myDMP <- champ.DMP(beta = myNorm_met,pheno=met_meta$MG,adjPVal = 1)

myDMP$NO_to_YES %>% filter(gene == "TTN")
myDMP$NO_to_YES %>% filter(gene == "CHRNA1")
myDMP$NO_to_YES %>% filter(gene == "NEFM")
myDMP$NO_to_YES %>% filter(gene == "RYR1")
myDMP$NO_to_YES %>% filter(gene == "RYR2")
myDMP$NO_to_YES %>% filter(gene == "RYR3")
myDMP$NO_to_YES %>% filter(gene == "CHRND")
myDMP$NO_to_YES %>% filter(gene == "CHRNG")
myDMP$NO_to_YES %>% filter(gene == "CHRNE")
myDMP$NO_to_YES %>% filter(gene == "CHRNB1")

myDMP$NO_to_YES %>% filter(gene == "KCNA4")


?champ.DMP

#DMR 分析
myDMR <- champ.DMR(beta=myNorm_met,pheno=met_meta$MG,method="Bumphunter",cores = 14)

DMR.GUI(DMR=myDMR,
         beta=myNorm_met,
         pheno=myload_met$pd$MG,
         runDMP=TRUE,
         compare.group=NULL,
         arraytype="450K")

#GSEA 分析
myGSEA <- champ.GSEA(beta=myNorm_met,DMP=myDMP[[1]], DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")

View(myGSEA$DMP)


## DMP gene and RNA expression gene

DMP_gene_list<-unique(myDMP$NO_to_YES$gene)

RNA_gene_list<-unique(result_tmg$gene_name[result_tmg$padj<0.1])

both_gene<-intersect(DMP_gene_list,RNA_gene_list)


both_gene_entr<-en2ENSE %>% filter(SYMBOL %in% both_gene) %>% dplyr::select(ENTREZID) %>% pull()

both_gene<-enrichGO(gene          = both_gene_entr,
         universe      = na.omit(result_tmg$ENTREZID),
         OrgDb         = org.Hs.eg.db,
         ont           = "BP",
         pAdjustMethod = "none",
         minGSSize = 1,
         maxGSSize = 150,
         pvalueCutoff  = 0.05,
         qvalueCutoff  = 0.05,
         readable      = TRUE)

View(both_gene@result)
barplot(both_gene,showCategory = 47)



## gene exon position

gene_exon<-read.csv("gene_exon.csv",header = T,sep = "\t")

gene_uscs_mRNA_id<-read.csv("uscs_mRNA_id.txt",header = T,sep = "\t")

gene_exon_ucsc_gene_name<- gene_exon %>% left_join(gene_uscs_mRNA_id,by = (c("name"="kgID"))) 



CHRNA1_from<-gene_exon_ucsc_gene_name %>% filter(geneSymbol == "CHRNA1") %>% 
  select(exonStarts,exonEnds) %>% .[2,1] %>% str_split(",") %>% unlist() %>% as.numeric() %>% na.omit

CHRNA1_to<-gene_exon_ucsc_gene_name %>% filter(geneSymbol == "CHRNA1") %>% 

    select(exonStarts,exonEnds) %>% .[2,2] %>% str_split(",") %>% unlist() %>% as.numeric() %>% na.omit()



TTN_from<-gene_exon_ucsc_gene_name %>% filter(geneSymbol == "TTN") %>% 
  select(exonStarts,exonEnds) %>% .[13,1] %>% str_split(",") %>% unlist() %>% as.numeric() %>% na.omit
TTN_to<-gene_exon_ucsc_gene_name %>% filter(geneSymbol == "TTN") %>% 
  select(exonStarts,exonEnds) %>% .[13,2] %>% str_split(",") %>% unlist() %>% as.numeric() %>% na.omit()


NEFM_from<-gene_exon_ucsc_gene_name %>% filter(geneSymbol == "NEFM") %>% 
  select(exonStarts,exonEnds) %>% .[1,1] %>% str_split(",") %>% unlist() %>% as.numeric() %>% na.omit
NEFM_to<-gene_exon_ucsc_gene_name %>% filter(geneSymbol == "NEFM") %>% 
  select(exonStarts,exonEnds) %>% .[1,2] %>% str_split(",") %>% unlist() %>% as.numeric() %>% na.omit


gene_exon_ucsc_gene_name %>% filter(geneSymbol == "TTN")


CHRNA1_met<-myDMP$NO_to_YES %>% filter(gene == "CHRNA1")
TTN_met<-myDMP$NO_to_YES %>% filter(gene == "TTN")
NEFM_met<-myDMP$NO_to_YES %>% filter(gene == "NEFM")

sp<-ggplot(CHRNA1_met,aes(x = MAPINFO , y = YES_AVG)) + geom_point()
sp_TTN<-ggplot(TTN_met,aes(x = MAPINFO , y = YES_AVG)) + geom_point()
sp_NEFM<-ggplot(NEFM_met,aes(x = MAPINFO , y = YES_AVG)) + geom_point()

# Change x and y axis labels, and limits
sp + scale_x_continuous(breaks=c(CHRNA1_from,CHRNA1_to))
sp_TTN + scale_x_continuous(breaks=c(TTN_from,TTN_to))
sp_NEFM + scale_x_continuous(breaks=c(NEFM_from,NEFM_to)) +geom_vline(xintercept = c(NEFM_from,NEFM_to), size = 0.5,
                                                                      color = "firebrick", linetype = "dashed")

?geom_ribbon

?geom_hline


?scale_x_continuous
scale_x_discrete
