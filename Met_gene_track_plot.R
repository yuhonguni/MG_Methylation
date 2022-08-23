library(Gviz)

# gse dataset NEFM plot
# indicate which genome is being used
gen <- "hg19"
# the index of the DMR that we will plot 
which(results.ranges$overlapping.genes=="NEFM") # 5488
dmrIndex <- 5488


# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <-  24770000 #as.numeric(start(results.ranges[dmrIndex]))
end <-  24776000 #as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))


# CpG islands (download from Wu etal  http://www.haowulab.org/software/makeCGI/index.html )
islandHMM <- read.csv("reference_data/model-based-cpg-islands-hg19.txt",
                      sep="\t", stringsAsFactors=FALSE, header=T)
head(islandHMM)


islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
islandData

## ideogram track and genome axis track
iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")


##build gene region track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
rTrack <- GeneRegionTrack(txdb, chromosome = "chr8", 
                        start = 24770000,  end = 24778000)

##Ensure that the methylation data is ordered by chromosome and base position.

myDMP_gse_NEFM<-myDMP_gse$no_to_yes %>% dplyr::filter(gene == "NEFM") %>% 
  dplyr:: filter ((MAPINFO < 24774000) & (MAPINFO > 24771000))

myDMP_gse_NEFM <- myDMP_gse_NEFM[order(myDMP_gse_NEFM$CHR,myDMP_gse_NEFM$MAPINFO),]

myDMP_tcga_NEFM<-myDMP$NO_to_YES %>% dplyr::filter(gene == "NEFM") %>% 
  dplyr:: filter ((MAPINFO < 24774000) & (MAPINFO > 24771000))

myDMP_tcga_NEFM <- myDMP_gse_NEFM[order(myDMP_tcga_NEFM$CHR,myDMP_gse_NEFM$MAPINFO),]

#  methylation data
NEFM_beta_gse <- myNorm_met_gse[match(rownames(myDMP_gse_NEFM),rownames(myNorm_met_gse)),]
NEFM_beta_tcga <- myNorm_met[match(rownames(myDMP_tcga_NEFM),rownames(myNorm_met)),]


# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(rep("chr8",19)),  ##改一下
                   ranges=IRanges(start=myDMP_gse_NEFM$MAPINFO, end=myDMP_gse_NEFM$MAPINFO),
                   strand=Rle(rep("*",nrow(myDMP_gse_NEFM))),
                   betas=NEFM_beta_gse)


cpgData_tcga <- GRanges(seqnames=Rle(rep("chr8",19)),  ##改一下
                   ranges=IRanges(start=myDMP_tcga_NEFM$MAPINFO, end=myDMP_tcga_NEFM$MAPINFO),
                   strand=Rle(rep("*",nrow(myDMP_tcga_NEFM))),
                   betas=NEFM_beta_tcga)

# extract data on CpGs in DMR
#cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])


# methylation data track
methTrack_gse <- DataTrack(range=cpgData, groups=meta_data_thymoma_kajiura$Myasthenia_gravis,genome = gen,
                       chromosome=chrom, ylim=c(0,0.6), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)

methTrack_tcga <- DataTrack(range=cpgData_tcga, groups=met_meta$MG,genome = gen,
                       chromosome=chrom, ylim=c(0,0.8), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)

# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")


# DMR position data track
dmrTrack <- AnnotationTrack(start=as.numeric(start(results.ranges[dmrIndex])), end=as.numeric(end(results.ranges[dmrIndex])), genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")

#Set up the track list and indicate the relative sizes of the different tracks. Finally, draw the plot using the plotTracks function (Figure 11).

tracks <- list(iTrack, gTrack, rTrack, islandTrack, dmrTrack,methTrack_tcga, methTrack_gse)

sizes <- c(2,2,4,2,2,12,12) # set up the relative sizes of the tracks
plotTracks(tracks, from=24770000, to=24777000, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))




