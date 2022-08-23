library(tidyverse)
library(ggplot2)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)




getwd()
# methylation data

setwd("C:/Users/yuhon/Documents/TMG_methylation_paper")


Methy_THYM_Kajiura<-read_tsv(file="./gse94769_meth_data/GSE94769_Average_Beta_Pval.txt",) 

Methy_THYM_Kajiura_det_p<- Methy_THYM_Kajiura %>% 
  dplyr::select("TargetID", str_subset(string = colnames(.),pattern = "Detection.Pval"))

head(Methy_THYM_Kajiura)
Methy_THYM_Kajiura_beta<- Methy_THYM_Kajiura %>% 
  dplyr::select("TargetID", str_subset(string = colnames(.),pattern = "AVG_Beta"))


colnames(Methy_THYM_Kajiura_beta) =  
  str_remove(string = colnames(Methy_THYM_Kajiura_beta),pattern=".AVG_Beta")

thymoma_ID<-Methy_THYM_Kajiura_beta


# meta data
meta_data_thy_kajiura<-read_tsv(file="./gse94769_meth_data//meta_data.txt")

meta_data_thy_kajiura$Masaoka_classification = str_remove(string = meta_data_thy_kajiura$Masaoka_classification,pattern="Stage: ")
meta_data_thy_kajiura$histology_WHO_classification = str_remove(string = meta_data_thy_kajiura$histology_WHO_classification,pattern="tissue: ")
meta_data_thy_kajiura$tissue = str_remove(string = meta_data_thy_kajiura$tissue, pattern="tissue: ")

   # thymoma( classification belong A and B) ID list
thymoma_ID<-meta_data_thy_kajiura$Sample[meta_data_thy_kajiura$histology_WHO_classification == "B3 thymoma" |meta_data_thy_kajiura$histology_WHO_classification == "type A thymoma"]

Methy_Thymoma_Kajiura_beta<-Methy_THYM_Kajiura_beta %>% dplyr::select("TargetID",thymoma_ID)
meta_data_thymoma_kajiura<-meta_data_thy_kajiura %>% filter(Sample %in% thymoma_ID)



## thymoma (A and B) methylation beta that has data
met_THYM_gse<- Methy_Thymoma_Kajiura_beta %>% dplyr::select("TargetID",meta_data_thymoma_kajiura$Sample)


met_THYM_m_gse <- met_THYM_gse %>% dplyr::select(-TargetID) %>% as.matrix()

row.names(met_THYM_m_gse)<-met_THYM_gse$TargetID

met_THYM_m_gse<-na.omit(met_THYM_m_gse)




##filter
myload_met_gse <- champ.filter(beta=met_THYM_m_gse,pd=meta_data_thymoma_kajiura,
                               fixOutlier = T,filterBeads = F,filterDetP = F,autoimpute = F)

##quality control
#champ.QC(beta=myload_met$beta,pheno =myload_met$pd$patient_id)


## normalization (beta value) and get M value
## myNorm_met_gse <- champ.norm(beta = myload_met_gse$beta,arraytype="450K", cores=8)

## write.csv(myNorm_met_gse,"./analy_data/gse94769_met_normalization_20220815.csv")

myNorm_met_gse<-read.csv("./analy_data/gse94769_met_normalization_20220815.csv",check.names = F)
colnames(myNorm_met_gse)[1]<-"geneID"
row.names(myNorm_met_gse)<-myNorm_met_gse$geneID
myNorm_met_gse<-myNorm_met_gse[,-1] %>% as.matrix

myNorm_met_gse<-na.omit(myNorm_met_gse)

  # get M value

M_myNorm_met_gse<-logit2(myNorm_met_gse)

##  MDS analysis

pal <- brewer.pal(8,"Dark2")   #creat a color panel
par(mfrow=c(1,2))
 # thymoma subtype
plotMDS(M_myNorm_met_gse, top=1000, gene.selection="common", 
        col=pal[factor(meta_data_thymoma_kajiura$histology_WHO_classification)])
legend("top", legend=levels(factor(meta_data_thymoma_kajiura$histology_WHO_classification)), text.col=pal,
       bg="white", cex=0.7)
 # MG status
plotMDS(M_myNorm_met_gse, top=1000, gene.selection="common", 
        col=pal[factor(meta_data_thymoma_kajiura$Myasthenia_gravis)])
legend("top", legend=levels(factor(meta_data_thymoma_kajiura$Myasthenia_gravis)), text.col=pal,
       bg="white", cex=0.7)
   
 # Masaoka classification
plotMDS(M_myNorm_met_gse, top=1000, gene.selection="common", 
        col=pal[factor(meta_data_thymoma_kajiura$Masaoka_classification)],dim=c(3,4))
legend("top", legend=levels(factor(meta_data_thymoma_kajiura$Masaoka_classification)), text.col=pal,
       bg="white", cex=0.7)

 # Gender
plotMDS(M_myNorm_met_gse, top=1000, gene.selection="common", 
        col=pal[factor(meta_data_thymoma_kajiura$Gender)],dim=c(2,3))
legend("top", legend=levels(factor(meta_data_thymoma_kajiura$Gender)), text.col=pal,
       bg="white", cex=0.7)



# dmrcate CPG annotation
group_gse<-factor(meta_data_thymoma_kajiura$Myasthenia_gravis) 
design_gse<-model.matrix(~group_gse)

myannotation <- cpg.annotate(datatype = "array", 
                             object = M_myNorm_met_gse, design = design_gse, coef = ncol(design_gse), fdr = 1,pcutoff = 1,
                              analysis.type = "differential", arraytype="450K", 
                             annotation = c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"), 
                             what = "M")

# dmrcate DMR analysis and plot NEFM 


dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2,pcutoff = 1)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")

which(results.ranges$overlapping.genes=="NEFM") # 5488


groups <- c(yes="magenta", no="forestgreen")
cols <- groups[meta_data_thymoma_kajiura$Myasthenia_gravis]
plot<-DMR.plot(ranges=results.ranges, dmr=5488, CpGs=myNorm_met_gse, what="Beta",
               arraytype = "450K", genome="hg19",phen.col=cols)


##### champ DMP analysis

myNorm_met<-na.omit(myNorm_met)
library(ChAMP)
myDMP_gse <- champ.DMP(beta = myNorm_met_gse,pheno=meta_data_thymoma_kajiura$Myasthenia_gravis,adjPVal = 1)


myDMP$no_to_yes %>% filter(gene == "TTN")
myDMP$no_to_yes %>% filter(gene == "CHRNA1")
myDMP_gse$no_to_yes %>% dplyr::filter(gene == "NEFM")
myDMP$no_to_yes %>% filter(gene == "RYR1")
myDMP$no_to_yes %>% filter(gene == "RYR2")
myDMP_gse$no_to_yes %>% filter(gene == "RYR3")
myDMP$no_to_yes %>% filter(gene == "CHRND")
myDMP$no_to_yes %>% filter(gene == "CHRNG")
myDMP$no_to_yes %>% filter(gene == "CHRNE")
myDMP$no_to_yes %>% filter(gene == "CHRNB1")


## DMR figure NEFM

 ## data transform into long form
library(reshape2)
myDMP_gse_NEFM %>% 
  pivot_wider(names_from = variable, values_from = value) %>%
ggplot(aes(x = MAPINFO)) + 
  geom_ribbon(aes(ymin = if_else(no_AVG > yes_AVG, yes_AVG, no_AVG), 
                  ymax = if_else(no_AVG > yes_AVG, no_AVG, yes_AVG)),
              fill = "#decbe4", colour = "black")

 ## 非线性拟合
myDMP_gse_NEFM<-myDMP_gse$no_to_yes %>% dplyr::filter(gene == "NEFM") %>% 
  filter ((MAPINFO < 24774000) & (MAPINFO > 24771000)) %>% select(MAPINFO,no_AVG, yes_AVG) %>% 
  melt(id.vars= "MAPINFO",measure.vars = c("no_AVG","yes_AVG"))

g <- ggplot(myDMP_gse_NEFM, aes(x = MAPINFO, y = value,color=variable))+
  geom_point() + geom_smooth()


 ##曲线下面积
myDMP_gse_NEFM<-myDMP_gse$no_to_yes %>% dplyr::filter(gene == "NEFM")  %>% 
  select(MAPINFO,no_AVG, yes_AVG) %>% 
  melt(id.vars= "MAPINFO",measure.vars = c("no_AVG","yes_AVG"))
g


# CHAMP DMR analysisis

myDMR_gse <- champ.DMR(beta=myNorm_met,pheno=as.factor(myload_met_gse$pd$Myasthenia_gravis),
                   method="Bumphunter",cores = 8,arraytype = "450K")



write.csv(myDMR_gse$BumphunterDMR,"gse94769_met_normalization_DMR_20220809.csv")





