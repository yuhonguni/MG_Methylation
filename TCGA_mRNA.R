library(tidyverse)
library(magrittr)
library(knitr)
library(kableExtra)
library(ggfortify)
library(ggrepel)

getwd()


setwd("C:/Users/yuhon/Documents") #laptop
#setwd("C:/Users/Administrator/Documents") # desktop
source('./PostDocProject/NSC_RNAseq/RNA_seq/functions.R')
R.utils::setOption("clusterProfiler.download.method",'auto')


#TCGA thymus RNA 
RNA_expr_THYM<-read_tsv(file="./TCGApaper/GDCdata/TCGA-THYM/legacy/Gene_expression/Gene_expression_quantification_UCSC_htseq/TCGA-THYM.htseq_counts.tsv"
)

RNA_expr_THYM<-RNA_expr_THYM[-c(60484:60488),]
   #log2(count+1) tranform to raw count
RNA_expr_THYM[,2:122]<-2^RNA_expr_THYM[,2:122]-1

RNA_expr_THYM$Ensembl_ID<-substr(RNA_expr_THYM$Ensembl_ID,1,15)

colnames(RNA_expr_THYM)[1]<-"gene_id"
RNA_expr_THYM <- RNA_expr_THYM %>% select(order(colnames(RNA_expr_THYM)))

 #remove normal tissue
RNA_expr_THYM<-RNA_expr_THYM %>% select(-c("TCGA-X7-A8D7-11A","TCGA-X7-A8D6-11A"))
colnames(RNA_expr_THYM)[2:length(colnames(RNA_expr_THYM))]<-substr(colnames(RNA_expr_THYM)[2:length(colnames(RNA_expr_THYM))],1,12)



#clinical

clinical_THYM<-read_tsv(file="./TCGApaper/GDCdata/TCGA-THYM/harmonized/Clinical/Clinical_Supplement/7cca5722-26cf-4ac8-a4a6-b803459f1861/nationwidechildrens.org_clinical_patient_thym_2.txt",
)

#update myasthenia gravis information

clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-3G-AB14"]<-"NO"
clinical_THYM$history_myasthenia_gravis[clinical_THYM$bcr_patient_barcode=="TCGA-X7-A8DF"]<-"NO"
clinical_THYM$history_myasthenia_gravis

meta_data<-clinical_THYM %>% select(patient_id=bcr_patient_barcode,
                                    gender,height,weight,
                                    histo_type = histological_type,
                                    MG = history_myasthenia_gravis,
                                    section_MG = section_myasthenia_gravis,
                                    age_patho_diagnosis = 
                                      age_at_initial_pathologic_diagnosis) %>%
                              filter(MG == "YES"|MG =="NO" ) %>%
                              arrange(patient_id)

meta_data_RNA<-meta_data %>% filter(meta_data$patient_id %in% (colnames(RNA_expr_THYM)[-1]))


#TCGA_paper_clinical_data

clinical_THYM_paper<-read_csv(file="./TCGApaper/TCGA database/TCGApaper/TCGA_paper.csv")




#看TCGA文章里面和数据库里面标本的临床资料是否一???

clinical_THYM_paper$patient_short_barcode %in% clinical_THYM$bcr_patient_barcode

clinical_THYM_paper$patient_short_barcode %in% clinical_THYM$bcr_patient_barcode

setdiff(clinical_THYM_paper$patient_short_barcode,meta_data$patient_id)
#"TCGA-3G-AB14" 1 MG NO dead 2 MG unknown alive  
#"TCGA-5K-AAAP" 1 MG unknown 2 MG unknown 一致
#"TCGA-X7-A8DF" 1 MG NO      2 MG not available
# 根据TCGA paper内容更新

setdiff(meta_data$patient_id,clinical_THYM_paper$patient_short_barcode)

intersect(clinical_THYM_paper$patient_short_barcode,
          clinical_THYM$bcr_patient_barcode)

intersect(clinical_THYM_paper$patient_short_barcode,
          meta_data$patient_id)

intersect(clinical_THYM_paper$patient_short_barcode[clinical_THYM_paper$history_myasthenia_gravis=="YES"],
          meta_data$patient_id[meta_data$MG == "YES"])

intersect(clinical_THYM_paper$patient_short_barcode[clinical_THYM_paper$history_myasthenia_gravis=="NO"],
          meta_data$patient_id[meta_data$MG == "NO"])



## RNA that has data
RNA_expr_THYM_meta<- RNA_expr_THYM %>% select("gene_id",meta_data_RNA$patient_id)

## geneName import
geneNames<-read.csv("./TCGApaper/geneNames.csv")%>% 
  mutate(gene_id = substr(.$gene_id,1,15) ) # remove ensemble gene_id version

#mitogenes
MitoGenes <- geneNames[geneNames$seqnames=="chrM",]
RiboGenes <- geneNames[geneNames$gene_type %in% c("rRNA", "Mt_rRNA"),]
ProtGenes <- geneNames[geneNames$gene_type=="protein_coding",]

# Now we remove the chr names from the geneNames because there are duplicate
# gene names in chrX and chrY 
geneNames %<>% dplyr::select(gene_id, gene_name, gene_type) %>%
  unique 


## 2 Quality control
## 2.1 Transcript biotypes

# Counts for mitochondrial and rRNA genes
MitoCountSum <- colSums(RNA_expr_THYM_meta %>% filter(gene_id %in% MitoGenes$gene_id) %>% dplyr::select(-gene_id))
MitoCountSum <- data.frame(patient_id=names(MitoCountSum), mito_counts=round(MitoCountSum), stringsAsFactors=FALSE)
RiboCountSum <- colSums(RNA_expr_THYM_meta %>% filter(gene_id %in% RiboGenes$gene_id) %>% dplyr::select(-gene_id))
RiboCountSum <- data.frame(patient_id=names(RiboCountSum), ribo_counts=round(RiboCountSum), stringsAsFactors=FALSE)
ProtCountSum <- colSums(RNA_expr_THYM_meta %>% filter(gene_id %in% (ProtGenes %>% filter(seqnames!="chrM"))$gene_id) %>% dplyr::select(-gene_id))
ProtCountSum <- data.frame(patient_id=names(ProtCountSum), prot_counts=round(ProtCountSum), stringsAsFactors=FALSE)
TotCountSum  <- colSums(RNA_expr_THYM_meta %>% dplyr::select(-gene_id))
TotCountSum <- data.frame(patient_id=names(TotCountSum), counts=TotCountSum, stringsAsFactors=FALSE)

ReadFractions <- inner_join(MitoCountSum, RiboCountSum, by="patient_id") %>%
  inner_join(ProtCountSum, by="patient_id") %>%
  inner_join(TotCountSum, by="patient_id") %>%
  mutate(mito_fraction=mito_counts/counts, ribo_fraction=ribo_counts/counts, prot_fraction=prot_counts/counts,
         other_fraction=1-(mito_fraction+ribo_fraction+prot_fraction), lib_size=as.integer(round(counts))) %>%
  dplyr::select(patient_id, mito_fraction, ribo_fraction, prot_fraction, other_fraction, lib_size=counts)

meta_data_RNA %<>% left_join(ReadFractions,by="patient_id")

ReadFractions %>% dplyr::select(-lib_size) %>% gather("Gene_biotype", "Proportion", -patient_id) %>%
  arrange(Proportion) %>% mutate(Gene_biotype=factor(Gene_biotype, levels=c("ribo_fraction","mito_fraction","other_fraction","prot_fraction"))) %>%
  ggplot(aes(x=patient_id, y=Proportion, fill=Gene_biotype)) +
  geom_bar(position="stack", stat="identity", colour="white") +
  labs(y="Proportion of reads", x="Sample") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank())



##2.2 Highly expressed genes
## highly expressed 

# Counts excluding mitochondrial genes and mitochondrial reads fraction over
# total lib. size
minProp <- 0.01
MitoCountFiltered <- RNA_expr_THYM_meta %>% filter(!gene_id %in% MitoGenes$gene_id)
CountFreqs <- prop.table(as.matrix(MitoCountFiltered[,-1]), 2)
rownames(CountFreqs) <- MitoCountFiltered$gene_id

#CountFreqs[CountFreqs < minProp] <- NA

HighExpGenes <- CountFreqs[apply(CountFreqs, 1, function(x) any(x >= minProp)),] %>%
  as.data.frame %>% mutate(gene_id=rownames(.)) %>%
  left_join(geneNames) %>% 
  gather("patient_id", "Proportion", -gene_id, -gene_name, -gene_type)

s_order <- HighExpGenes %>% group_by(patient_id) %>% summarise(sum=sum(Proportion)) %>%
  arrange(sum) %>% pull(patient_id) 

g_order <- HighExpGenes %>% group_by(gene_name) %>% summarise(sum=sum(Proportion)) %>%
  arrange(sum) %>% pull(gene_name) 

HighExpGenes$patient_id <- factor(HighExpGenes$patient_id, levels=s_order)
HighExpGenes$gene_name <- factor(HighExpGenes$gene_name, levels=g_order)

ggplot(HighExpGenes, aes(patient_id, Proportion, fill=gene_name)) +
  geom_bar(stat = "identity", colour="white") +
  labs(y="Proportion of reads", x="Patient") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank())


##Now that we have checked for highly abundant transcripts, we can transform the count data to log2(CPM) and do some more filtering on the non-mitochondrial genes. Genes are filtered out if:

##  no Gene Symbol (or duplicated)
##  expression is below the median expression for more than 20% of the samples

MitoCountFiltered<-RNA_expr_THYM_meta %>% filter(!gene_id %in% MitoGenes$gene_id)

countMatrixFiltered <- MitoCountFiltered

GenesKept <- rep(0, 3)
names(GenesKept) <- c("Total nuclear transcripts", "Trancripts above threshold", "With Gene Symbol")

# Create log2 CPM matrix after removal of mitochondria-encoded genes
cpmMatrixFiltered <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), log2(Count2CPM(countMatrixFiltered[,-1])+1),check.names=F)
GenesKept[1] <- nrow(cpmMatrixFiltered)

# Add gene symbols
ExpDataCPM <- left_join(cpmMatrixFiltered, geneNames %>% dplyr::select(gene_id, gene_name)) %>%
  dplyr::select(ensemblID=gene_id, GeneSymbol=gene_name, everything())

# Filter low-expressed genes similar to https://f1000research.com/articles/5-1438
# (although a bit stricter: keep genes that have at least 10 reads in 25% of
# the samples)
min_lib_size <- min(meta_data_RNA$lib_size)
min_cpm <- round(10/(min_lib_size/1e+6),1)
keep <- rowSums(cpmMatrixFiltered[,-1] > min_cpm) >= round(nrow(meta_data_RNA)*.25)
#keep %>% table

# Filter genes below noise level
ExpDataCPM <- ExpDataCPM [keep,]
GenesKept[2] <- nrow(ExpDataCPM)

# Filter out genes with no Gene Symbol or with duplicated Gene Symbol (highest
# expressed option of duplicates is kept)
ExpDataCPM %<>% filter(GeneSymbol != "", !is.na(GeneSymbol)) %>%
  mutate(Sum = do.call(pmax, select_if(., is.numeric))) %>%
  arrange(desc(Sum)) %>% 
  distinct(GeneSymbol, .keep_all=TRUE) %>%
  dplyr::select(-Sum) %>%
  arrange(ensemblID)

GenesKept[3] <- nrow(ExpDataCPM)
prot_coding_n <- c(geneNames %>% filter(gene_id %in% ExpDataCPM$ensemblID) %>%
                     pull(gene_type) %>% table %>% prop.table)["protein_coding"]


GenesKept %>% kable(caption="Number of genes passing filters") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

# Reorder samples in CPM matrix
colnames(ExpDataCPM)[3:length(colnames(ExpDataCPM))] == meta_data_RNA$patient_id


if ( ! all(colnames(ExpDataCPM)[-c(1,2)] == meta_data_RNA$patient_id) ) {
  warning("Reordering CPM matrix...")
  ExpDataCPM <- ExpDataCPM[,c(1,2,(match(meta_data_RNA$patient_id, colnames(ExpDataCPM[,-c(1,2)]))+2))]
} else {
  message("No need to reorder CPM matrix")
}

# Filter low-expressed genes from counts matrix and reorder if necessary
Counts <- data.frame(countMatrixFiltered,check.names=F) %>%
  filter(gene_id %in% ExpDataCPM$ensemblID) %>%
  mutate_if(is.numeric, function(x) as.integer(round(x)))

if ( ! all(colnames(Counts)[-1] == meta_data_RNA$patient_id) ) {
  warning("Reordering count matrix...")
  Counts <- Counts[,c(1,(match(meta_data_RNA$patient_id, colnames(Counts[,-1]))+1))]
} else {
  message("No need to reorder count matrix")
}


############################################################################################
require("DESeq2")


# DE analysis

Mt <- meta_data_RNA %>% filter(MG == "NO" | MG == "YES")
Mt %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)

table(Mt$MG)

Ct <- Counts[,colnames(Counts) %in% Mt$patient_id]
rownames(Ct) <- Counts$gene_id

Md <- as.formula(paste0("~ MG"))
dds <- DESeq2RUN(Ct, Mt, Md)
?DESeq()
?p.adjust

?results

res_tmg_t<-results(dds, alpha = 0.05, format = "DataFrame", 
                   independentFiltering = T,pAdjustMethod = "fdr") %>% 
  data.frame() %>% dplyr::arrange(pvalue)

result_tmg_t<-data.frame(res_tmg_t) %>% mutate(gene_id=rownames(.)) %>% right_join(geneNames,by="gene_id") %>% arrange(pvalue)


result_tmg_t[result_tmg_t$gene_name=='CHRNA1',]

#CHRNA1 foldchange

lg2fc_chrna1<-result_tmg_t[result_tmg_t$gene_name=='CHRNA1',]$log2FoldChange
2^lg2fc_chrna1

result_tmg_t[result_tmg_t$gene_name=='CHRNG',]

result_tmg_t[result_tmg_t$gene_name=='NEFM',]# "similiarity with CHRNA1"

result_tmg_t[result_tmg_t$gene_name=='RYR1',]
result_tmg_t[result_tmg_t$gene_name=='RYR2',]

result_tmg_t[result_tmg_t$gene_name=='RYR3',]

Counts[Counts$gene_id=="ENSG00000138435",]## CHRNA1 foldchange

chrna1_counts_mg<-Counts %>% filter(gene_id=="ENSG00000138435") %>% dplyr::select(meta_data_RNA$patient_id[meta_data_RNA$MG=="YES"])


chrna1_counts_nmg<-Counts %>% filter(gene_id=="ENSG00000138435") %>% dplyr::select(meta_data_RNA$patient_id[meta_data_RNA$MG=="NO"])

median(t(chrna1_counts_mg))/median(t(chrna1_counts_nmg))

meta_data$patient_id[meta_data$MG=="YES"]

meta_data$patient_id[meta_data$MG=="NO"]



View(meta_data)



result_tmg_t[result_tmg_t$gene_name=='TTN',]
Counts[Counts$gene_id=="ENSG00000155657",]

## pathway analysis 

library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)

# extract ensemble and entrze gene ID mapping files
k<-keys(org.Hs.eg.db,keytype='ENSEMBL')
en2ENSE<-AnnotationDbi::select(org.Hs.eg.db,keys=k,
                               columns=c("ENTREZID",'SYMBOL'),
                               keytype = 'ENSEMBL')

result_tmg <- result_tmg_t %>% left_join(en2ENSE[,1:2],by=(c('gene_id'='ENSEMBL'))) %>% 
  dplyr::select(gene_name,gene_id,ENTREZID,gene_type,log2FoldChange,pvalue,padj)

result_tmg_up<-result_tmg %>% filter(padj<0.1 & log2FoldChange>0)
result_tmg_down<-result_tmg %>% filter(padj<0.1 & log2FoldChange<0)

## KEGG enrichment analysis
KEGG_up_tmg<-enrichKEGG(na.omit(result_tmg_up$ENTREZID), organism = "hsa", keyType = "kegg",
           pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1,
           universe = na.omit(result_tmg$ENTREZID))
setReadable(KEGG_up_tmg, org.Hs.eg.db, 
            keyType = "ENTREZID")
barplot(KEGG_up_tmg,showCategory = 23)


KEGGM_up_tmg<-enrichMKEGG(na.omit(result_tmg_up$ENTREZID), organism = "hsa", minGSSize=1,
                          pvalueCutoff = 0.05, pAdjustMethod = "none",qvalueCutoff = 1,
            universe = na.omit(result_tmg$ENTREZID))
setReadable(KEGGM_up_tmg, org.Hs.eg.db, 
            keyType = "ENTREZID")
barplot(KEGG_up_tmg,showCategory = 8)

##  GO analysis

## bp
go_up_tmg_bp <- enrichGO(gene          = na.omit(result_tmg_up$ENTREZID),
                universe      = na.omit(result_tmg$ENTREZID),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "none",
                minGSSize = 10,
                maxGSSize = 150,
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)

?enrichGO
barplot(go_up_tmg_bp,showCategory = 15)

View(go_up_tmg_bp@result)

go_down_tmg_bp <- enrichGO(gene          = na.omit(result_tmg_down$ENTREZID),
                         universe      = na.omit(result_tmg$ENTREZID),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "BP",
                         pAdjustMethod = "none",
                         minGSSize = 5,
                         maxGSSize = 150,
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
View(go_up_tmg_bp@result)


barplot(go_down_tmg_bp,showCategory = 60)

View(go_down_tmg_bp@result)


## mf
go_up_tmg_mf <- enrichGO(gene          = na.omit(result_tmg_up$ENTREZID),
                         universe      = na.omit(result_tmg$ENTREZID),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "MF",
                         pAdjustMethod = "none",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
barplot(go_up_tmg_mf,showCategory = 15)
View(go_up_tmg_mf@result)

go_down_tmg_mf <- enrichGO(gene          = na.omit(result_tmg_down$ENTREZID),
                         universe      = na.omit(result_tmg$ENTREZID),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "MF",
                         pAdjustMethod = "none",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
barplot(go_down_tmg_mf,showCategory = 50)
View(go_down_tmg_mf@result)


## cc
go_up_tmg_cc <- enrichGO(gene          = na.omit(result_tmg_up$ENTREZID),
                         universe      = na.omit(result_tmg$ENTREZID),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "CC",
                         pAdjustMethod = "none",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
barplot(go_up_tmg_cc,showCategory = 15)

go_down_tmg_cc <- enrichGO(gene          = na.omit(result_tmg_down$ENTREZID),
                         universe      = na.omit(result_tmg$ENTREZID),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "CC",
                         pAdjustMethod = "none",
                         pvalueCutoff  = 0.05,
                         readable      = TRUE)
barplot(go_down_tmg_cc,showCategory = 50)




?enrichGO


## GO analysis
GO_ws5a_cp2a_enrich<-function(a) {
  ego <- enrichGO(gene          = de_ws5a_cp2a_list[[a]]$ENTREZID,
                  universe      = bg_genes,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  minGSSize = 20,
                  maxGSSize = 150,
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(ego)
}



## KEGG module analysis
MKEGG_ws5a_cp2a_enrich<-function(a) {
  enrich <- enrichMKEGG(de_ws5a_cp2a_list[[a]]$ENTREZID, organism = "hsa", minGSSize=1,
                        universe = bg_genes)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

heatplot(MKEGG_ws5a_cp2a_enrich('nsc_up'))
heatplot(MKEGG_ws5a_cp2a_enrich('nsc_down'))

heatplot(MKEGG_ws5a_cp2a_enrich('ipsc_up'))
heatplot(MKEGG_ws5a_cp2a_enrich('ipsc_down'))

##import wiki pathway annotation file 

wpgmtfile <- system.file("extdata/wikipathways-20180810-gmt-Homo_sapiens.gmt", package="clusterProfiler")
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%") %>% 
  filter(gene %in% Gene_ID_ENTREZ$ENTREZID)

wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


## wiki enrichment analysis
wiki_ws5a_cp2a_enrich<-function(a) {
  enrich <- enricher(de_ws5a_cp2a_list[[a]]$ENTREZID, pvalueCutoff = 0.05,pAdjustMethod = "BH",
                     qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 200,
                     TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
  enrich <- setReadable(enrich, org.Hs.eg.db, 
                        keyType = "ENTREZID")
  return(enrich)
}

heatplot(wiki_ws5a_cp2a_enrich('nsc_up'),showCategory = 23)
barplot(wiki_ws5a_cp2a_enrich('nsc_down'),showCategory = 12)

wiki_ws5a_cp2a_enrich('ipsc_up')
wiki_ws5a_cp2a_enrich('ipsc_down')





