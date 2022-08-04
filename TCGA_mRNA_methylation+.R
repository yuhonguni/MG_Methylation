library(readr)
library(dplyr)
library(tidyverse)
library(magrittr)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggfortify)
library(ggrepel)
setwd("C:/Users/yuhon/Documents/")
#TCGA thymus RNA 
RNA_expr_THYM<-read_tsv(file="C:/Users/yuhon/Documents/TCGApaper/GDCdata/TCGA-THYM/legacy/Gene_expression/Gene_expression_quantification_UCSC_htseq/TCGA-THYM.htseq_counts.tsv"
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

colnames(RNA_expr_THYM)[-1]



#clinical
clinical_THYM<-read_tsv(file="C:/Users/yuhon/Documents/TCGApaper/GDCdata/TCGA-THYM/harmonized/Clinical/Clinical_Supplement/7cca5722-26cf-4ac8-a4a6-b803459f1861/nationwidechildrens.org_clinical_patient_thym.txt",
)

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

## RNA that has data
RNA_expr_THYM_meta<- RNA_expr_THYM %>% select("gene_id",meta_data_RNA$patient_id)

## geneName import
geneNames<-read.csv("C:/Users/yuhon/Documents/geneNames.csv")%>% 
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
MitoCountSum <- colSums(RNA_expr_THYM %>% filter(gene_id %in% MitoGenes$gene_id) %>% dplyr::select(-gene_id))
MitoCountSum <- data.frame(patient_id=names(MitoCountSum), mito_counts=round(MitoCountSum), stringsAsFactors=FALSE)
RiboCountSum <- colSums(RNA_expr_THYM %>% filter(gene_id %in% RiboGenes$gene_id) %>% dplyr::select(-gene_id))
RiboCountSum <- data.frame(patient_id=names(RiboCountSum), ribo_counts=round(RiboCountSum), stringsAsFactors=FALSE)
ProtCountSum <- colSums(RNA_expr_THYM %>% filter(gene_id %in% (ProtGenes %>% filter(seqnames!="chrM"))$gene_id) %>% dplyr::select(-gene_id))
ProtCountSum <- data.frame(patient_id=names(ProtCountSum), prot_counts=round(ProtCountSum), stringsAsFactors=FALSE)
TotCountSum  <- colSums(RNA_expr_THYM %>% dplyr::select(-gene_id))
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
MitoCountFiltered <- RNA_expr_THYM %>% filter(!gene_id %in% MitoGenes$gene_id)
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
cpmMatrixFiltered <- data.frame(gene_id=as.character(countMatrixFiltered$gene_id), log2(Count2CPM(countMatrixFiltered[,-1])+1))
GenesKept[1] <- nrow(cpmMatrixFiltered)

# Add gene symbols
ExpDataCPM <- left_join(cpmMatrixFiltered, geneNames %>% dplyr::select(gene_id, gene_name)) %>%
  dplyr::select(ensemblID=gene_id, GeneSymbol=gene_name, everything())

# Filter low-expressed genes similar to https://f1000research.com/articles/5-1438
# (although a bit stricter: keep genes that have at least 10 reads in 25% of
# the samples)
min_lib_size <- min(Metadata$lib_size)
min_cpm <- round(10/(min_lib_size/1e+6),1)
keep <- rowSums(cpmMatrixFiltered[,-1] > min_cpm) >= round(nrow(Metadata)*.25)
#keep %>% table

# Filter genes below noise level
ExpDataCPM <- ExpDataCPM[keep,]
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
substring(colnames(ExpDataCPM)[3:length(colnames(ExpDataCPM))],9,12) == substring(meta_data_RNA$patient_ID,9,12)

meta_data_RNA$patient_id <-colnames(ExpDataCPM)[3:length(colnames(ExpDataCPM))]

if ( ! all(colnames(ExpDataCPM)[-c(1,2)] == meta_data_RNA$patient_id) ) {
  warning("Reordering CPM matrix...")
  ExpDataCPM <- ExpDataCPM[,c(1,2,(match(meta_data_RNA$patient_id, colnames(ExpDataCPM[,-c(1,2)]))+2))]
} else {
  message("No need to reorder CPM matrix")
}

# Filter low-expressed genes from counts matrix and reorder if necessary
Counts <- data.frame(countMatrixFiltered) %>%
  filter(gene_id %in% ExpDataCPM$ensemblID) %>%
  mutate_if(is.numeric, function(x) as.integer(round(x)))

if ( ! all(colnames(Counts)[-1] == meta_data_RNA$patient_id) ) {
  warning("Reordering count matrix...")
  Counts <- Counts[,c(1,(match(meta_data_RNA$patient_id, colnames(Counts[,-1]))+1))]
} else {
  message("No need to reorder count matrix")
}


##2.4 Sample outliers

## Now that we have filtered the genes, we can look at the samples and identify outliers. 
## To that end, we calculate the pairwise correlations in expression between samples. 
## We characterise each sample by its median correlation to the rest of the samples. 
## Based on this value, we mark as outliers the samples outside the interval Q1-1.5*IQR,R3+1.5*IQR 
## (Q1=first quartile, Q3=third quartile, IQR=interquartile range). 
## We restrict the correlation to protein-coding genes.

sc <- ExpDataCPM %>%
  filter(ensemblID %in% ProtGenes$gene_id) %>%
  dplyr::select(-ensemblID, -GeneSymbol) %>%
  as.matrix %>%
  cor
diag(sc) <- NA
medians <- matrixStats::rowMedians(sc, na.rm=TRUE)
names(medians) <- colnames(sc)

iqrange <- IQR(medians)
quartiles <- quantile(medians, c(0.25, 0.75))
lim_low <- quartiles[1]-(1.5*iqrange)
lim_high <- quartiles[2]+(1.5*iqrange)

sample.outliers <- (medians < lim_low) | (medians > lim_high)

# Update object
meta_data_RNA <- bind_cols(meta_data_RNA, sample_outlier=sample.outliers)

View(meta_data_RNA)






####################################################################################













## 3 Sample characterisation
## 3.1 Sample clustering
## We now plot all the samples and their correlation to further explore outliers.

library(ComplexHeatmap)


# Create temporal Metadata data.frame
dfMeta <- meta_data_RNA %>% dplyr::select(-patient_id, -height, -weight, -section_MG, -age_patho_diagnosis) %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at(vars(matches("fraction$")), scales::rescale) %>%
  mutate(lib_size=scales::rescale(lib_size)) %>%
  as.data.frame
rownames(dfMeta) <- meta_data_RNA$patient_id

col_fun_type <- circlize::colorRamp2(c(0, 1), c("black", "orange3"))

col_fun_libsize <- circlize::colorRamp2(c(0,1), c("red3", "springgreen3"))
library(RColorBrewer)


column_ha = HeatmapAnnotation(GeneTypeFraction=cbind(mtDNA_Genes=dfMeta$mito_fraction,
                                                     Ribo_Genes=dfMeta$ribo_fraction,
                                                     Prot_Genes=dfMeta$prot_fraction,
                                                     Other_Genes=dfMeta$other_fraction),
                              Library_Size=dfMeta$lib_size,
                              Pathology=dfMeta$histo_type,
                              MG=dfMeta$MG,
                              #CellType=dfMeta$cell_type,
                             # Mutation=dfMeta$mutation,
                              #Homo_Het=dfMeta$homo_het,
                              #Subname_Mutation=dfMeta$subname_mut,
                              #Clone=dfMeta$clone,
                              #Source_Clone=dfMeta$source_of_clone,
                              #Age=dfMeta$age,
                              #Sex=dfMeta$sex,
                              #Sample_Outlier=dfMeta$sample.outlier,
                              na_col="white",
                              col=list(#DV200=col_fun_dv200,
                                       GeneTypeFraction=col_fun_type,
                                       Library_Size=col_fun_libsize
                                       #CellType=cols_celltype,
                                       #Mutation=cols_mutation,
                                       #Age=cols_age,
                                       #Sex=cols_sex
                                       #Sample_Outlier=c("TRUE"="red", "FALSE"="grey")
                              ))
Heatmap(sc, column_title="Sample correlation in gene expression", col=viridis::viridis(10), top_annotation=column_ha,
        show_row_names=FALSE, show_column_names=TRUE, show_row_dend=FALSE, column_dend_height=unit(5,"cm"))

## 3.2 Correlation between variables

## We can plot the pair-wise Pearson’s correlation between the variables to have 
## a general picture of which ones are more associated.

library(corrplot)

dfDummy <- fastDummies::dummy_cols(dfMeta) %>%
  select_if(is.numeric)
rownames(dfDummy) <- rownames(dfMeta)

M <- cor(dfDummy)
res <- cor.mtest(dfDummy, conf.level=0.95)
diag(M) <- NA

corrplot.mixed(M, tl.pos="lt", tl.col="black", tl.srt=45, order="original", na.label=" ")

##3.3 PCA over all samples and transcripts

## Now we plot the PCA of all samples.

sample_pca <- prcomp(t(ExpDataCPM[,-c(1,2)]))

pcameta<-dfMeta 
#

autoplot(sample_pca, data=pcameta,colour="MG")


assocPCA(M=t(ExpDataCPM[,-c(1,2)]), variables=dfDummy, plot.return=TRUE, ncomp=5, plot.alpha=0.05, plot.colour="orangered") +
  ggtitle("Association of PC 1-5 with experimental variables")






############################################################################################
require("DESeq2")


# DE analysis

Mt <- meta_data_RNA %>% filter(MG == "NO" | MG == "YES") %>%
  filter(!sample_outlier)
Mt %>% kable(digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F)


Ct <- Counts[,colnames(Counts) %in% Mt$patient_id]
rownames(Ct) <- Counts$gene_id

Md <- as.formula(paste0("~ MG"))
dds <- DESeq2RUN(Ct, Mt, Md)

res_tmg_t<-results(dds, alpha = 0.05, format = "DataFrame", independentFiltering = T) %>% arrange(pvalue)


result_tmg_t<-data.frame(res_tmg_t) %>% mutate(gene_id=rownames(.)) %>% left_join(geneNames,by="gene_id") %>% arrange(pvalue)