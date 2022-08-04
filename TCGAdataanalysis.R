BiocManager::install("TCGAbiolinks",force=TRUE)
BiocManager::install("maftools",force=TRUE)
library("TCGAbiolinks")
library(dplyr)
library(data.table)
library(SummarizedExperiment)
library(tidyverse)
library(maftools)
# packageVersion("TCGAbiolinks")
getwd()

?GDCprepare()
#?̳?
# browseVignettes("TCGAbiolinks")
#????????????
TCGAbiolinks:::getProjectSummary("TCGA-THYM")

#?ٴ?????
clinical <- GDCquery_clinic(project = "TCGA-THYM", type = "clinical")
View(clinical)

##tcga THYMUS DATA
query <- GDCquery(
  project = "TCGA-THYM",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
?GDCquery()
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)

patient_info <- clinical.BCRtab.all$clinical_patient_thym

View(clinical.BCRtab.all)
View(patient_info)

#???ݷ???
# RNA sequence data
query_RNAseq <- GDCquery(
  project = "TCGA-THYM",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)

GDCdownload(query_RNAseq,method="client")

THYMRnaseqSE <- GDCprepare(query_RNAseq)

THYMMatrix <- assay(THYMRnaseqSE,"raw_count")

# DNA methylation data
query_DNAmethy <- GDCquery(
  project = "TCGA-THYM",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  legacy = F
)
GDCdownload(query_DNAmethy,method="client")

## miRNA data
query_miRNA <- GDCquery(
  project = "TCGA-THYM",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  legacy = FALSE
)


# RNA sequence data
query_RNAseq <- GDCquery(
  project = "TCGA-THYM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  legacy = F
)
THYMRnaseqSE <- GDCprepare(query_RNAseq)
getwd()
?GDCquery
GDCdownload(query_RNAseq,method="client")

#SNP data
query_SNV <- GDCquery(
  project = "TCGA-THYM",
  data.category = "Simple nucleotide variation",
  data.type = "Simple somatic mutation",
  legacy = T
)
GDCdownload(query_SNV)





