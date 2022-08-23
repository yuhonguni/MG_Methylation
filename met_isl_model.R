library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(stringr)
library(data.table)

data(Islands.UCSC)
data(Locations)
data(Other)
#methylation island data from ill450hg19 package

isl_UC<- data.frame(Islands.UCSC) %>% 
  mutate(me_id = row.names(.)) 

me_ge_info<- data.frame(Other) %>% 
  select("UCSC_RefGene_Name","UCSC_RefGene_Accession","UCSC_RefGene_Group","Enhancer") %>%
  mutate(me_id = row.names(.)) 

isl_info<-data.frame(Locations) %>% 
  mutate(me_id = row.names(.)) %>% 
  inner_join(isl_UC,by="me_id") %>% 
  inner_join(me_ge_info,by="me_id") 

# extract methylation island methylation loci data,extract the cpg island with at least 3 methylation loci
# TCGA
me_da_t<-myNorm_met %>% data.frame(check.names=F) %>% mutate(me_id=row.names(.))
me_da_isl<-isl_info %>% right_join(me_da_t, by="me_id") %>% 
  filter(Islands_Name!="" & Relation_to_Island =="Island")
me_num_pisl<-me_da_isl %>% group_by(Islands_Name) %>% summarise(count = length(Islands_Name)) %>% 
  filter(count>=3) %>% arrange(Islands_Name)
me_num_pisl

me_da_isl<-me_da_isl %>% filter(Islands_Name %in% me_num_pisl$Islands_Name) %>% arrange(Islands_Name)

me_da_isl$Islands_Name

## TCGA
me_list<-list()
Sum_table<-data.frame()
j<-1

for (i in c(1:length(me_num_pisl$Islands_Name))) {
  me_da_t<-me_da_isl[(j:((me_num_pisl$count[i])+j-1)),]
  
  me_da_t2<-me_da_t %>% melt(colnames(.)[1:10],colnames(.)[11:length(colnames(.))],
                             variable.name="patient_id",value.name="me_beta") %>% arrange(pos) %>% 
    left_join(meta_data_thymoma_kajiura, by=c("patient_id"="Sample") )%>% 
    mutate(me_isl_id=paste0(me_id,"-",patient_id))
  
  
  expr_isl_t<-me_da_t2 %>% select(me_isl_id,me_beta) %>% 
    spread(key=me_isl_id,value=me_beta) %>% 
    select(me_da_t2$me_isl_id)
  
  met_loci<-factor(me_da_t2$me_id)
  MG<-factor(me_da_t2$MG,levels=c("YES","NO"))
  design<-model.matrix(~met_loci+MG)
  
  fit <- lmFit(expr_isl_t, design)
  fit <- eBayes(fit)
  topTable<-topTable(fit, coef="MGNO") %>% 
    mutate(chr=me_da_t$chr[1], 
           Islands_Name= me_da_t$Islands_Name[1], 
           me_id = paste(me_da_t$me_id,collapse= ";"), 
           UCSC_RefGene_Name= paste(me_da_t$UCSC_RefGene_Name,collapse= ";"),
           me_count = me_num_pisl$count[i])
  
  Sum_table<-rbind(Sum_table,topTable)
  j<-j+me_num_pisl$count[i]
}

Sum_table_tcga<-Sum_table %>% select(Islands_Name,chr,me_id,me_count,UCSC_RefGene_Name,
                                     everything())
length(Sum_table_tcga$Islands_Name)

write.csv(Sum_table_tcga,"CPG_island_differential_linear_model_TCGA.csv")


## GSE

me_da_g<-myNorm_met_gse %>% data.frame(check.names=F) %>% mutate(me_id=row.names(.))
me_da_islg<-isl_info %>% right_join(me_da_g, by="me_id") %>% 
  filter(Islands_Name!="" & Relation_to_Island =="Island")
me_num_pislg<-me_da_islg %>% group_by(Islands_Name) %>% summarise(count = length(Islands_Name)) %>% 
  filter(count>=3) %>% arrange(Islands_Name)
me_num_pislg

me_da_islg<-me_da_islg %>% filter(Islands_Name %in% me_num_pislg$Islands_Name) %>% arrange(Islands_Name)
me_da_islg$Islands_Name


me_da_islg[me_da_islg$Islands_Name=="chr8:24770908-24772547",]

me_num_pislg[me_num_pislg$Islands_Name=="chr8:24770908-24772547",]





## GSE

Sum_table<-data.frame()
j<-1
len<-length(me_num_pislg$Islands_Name)

for (i in c(16256:len)) {
  me_da_t<-me_da_islg[(j:((me_num_pislg$count[i])+j-1)),]
  
  me_da_t2<-me_da_t %>% melt(colnames(.)[1:10],colnames(.)[11:length(colnames(.))],
                             variable.name="patient_id",value.name="me_beta") %>% arrange(pos) %>% 
    left_join(meta_data_thymoma_kajiura, by=c("patient_id"="Sample") )%>% 
    mutate(me_isl_id=paste0(me_id,"-",patient_id))
  
  expr_isl_t<-me_da_t2 %>% select(me_isl_id,me_beta) %>% 
    spread(key=me_isl_id,value=me_beta) %>% 
    select(me_da_t2$me_isl_id)
  
  met_loci<-factor(me_da_t2$me_id)
  MG<-factor(me_da_t2$Myasthenia_gravis,levels=c("yes","no"))
  design<-model.matrix(~met_loci+MG)
  
  fit <- lmFit(expr_isl_t, design)
  fit <- eBayes(fit)
  topTable<-topTable(fit, coef="MGno") %>% 
    mutate(chr=me_da_t$chr[1], 
           Islands_Name= me_da_t$Islands_Name[1], 
           me_id = paste(me_da_t$me_id,collapse= ";"), 
           UCSC_RefGene_Name= paste(me_da_t$UCSC_RefGene_Name,collapse= ";"),
           me_count = me_num_pislg$count[i])
  
  Sum_table<-rbind(Sum_table,topTable)
  j<-j+me_num_pislg$count[i]
}

Sum_table_gse<-Sum_table %>% select(Islands_Name,chr,me_id,me_count,UCSC_RefGene_Name,
                                    everything())
write.csv(Sum_table_gse,"CPG_island_differential_linear_model_gse.csv")


Sum_table_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFM"))
Sum_table_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "NEFM"))


Sum_table_gse %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))
Sum_table_tcga %>% filter(str_detect(.$UCSC_RefGene_Name, pattern= "RYR3"))

