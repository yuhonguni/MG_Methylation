### FUNCTIONS

readPicardCollapsed <- function( files=list() ){
  if(length(files)==0) stop("No ", pattern, " files found in ", path, call.=FALSE)

  out1 <- vector("list", length(files))

  for(i in seq_along( files)){
    x <- readr::read_lines( files[i] )
    x <- x[x!=""]

    n1 <- grep("^PF_BASES", x)
  # warn and skip instead?
    if( length(n1) == 0)  stop("Line starting with PF_BASES not found in ", n[i] )
    x1 <- readr::read_tsv( paste(x[n1:(n1+1)], collapse ="\n") )
    names(x1) <- tolower(names(x1))
    ## long format, convert to numeric, remove NAs and add sample name
    out1[[i]] <- tidyr::gather(x1, key="key", value="value", -c("sample", "library", "read_group"))  %>%
        dplyr::mutate(sample=basename(files[i])) %>%
        dplyr::mutate(value = as.numeric(value))
  }
  message("Read ", length(files), " files")
  dplyr::bind_rows(out1)
}


Count2CPM <- function(countData){
  sapply(countData, function(smp) (10^6)*smp/sum(smp))
}

CountEGL2TPM <- function(countData, geneLength){
    if (!all(dim(countData) == dim(geneLength))) stop("Matrices need to have the same dimensions")
    A <- 10^3 * countData / geneLength
    10^6*t(t(A)/colSums(A))
}



GetMarkers <- function(region){
    region <- as.character(region)
    CellType_genes <- mouseMarkerGenesCombined[[region]]

    for (i in 1:length(CellType_genes)) {
        CellType_genes[[i]] <- as.vector(mouse2human(CellType_genes[[i]])$humanGene)
    }
    
    names(CellType_genes) <- sapply(names(CellType_genes), function(x) paste(x, "Genes", sep="_"))
    
    # Exclude specific genes: 5HTR3A - also expressed in VIP negative 5HTR3A cells.
    if (region == "Cortex") {
        CellType_genes$GabaVIPReln_Genes <- CellType_genes$GabaVIPReln_Genes[!CellType_genes$GabaVIPReln_Genes %in% "HTR3A"]
    }
    
    names(CellType_genes) <- sapply(names(CellType_genes), function(x) gsub(" ", "_", x))
    return(CellType_genes)
} 

correct_sign <- function(Celltype, object){
  sign = sum(object$rotation[,1])  
  if(sign < 0) {
    object$rotation[,1] <- -1*object$rotation[,1]
    object$x[,1] <- -1*object$x[,1]
  }
  return(object)
}

PCA_genes_All_based <- function(dataset, group, CellType_genes, NoiseThreshold, contName="Control"){
    colnames(group) <- c("id", "group")
    group %<>% filter(id %in% colnames(dataset))
    groups <- unique(as.character(group[[2]]))
    groupList <- split(group$id, group$group)
    PCAresults <- list()
    PCAresults$ControlOnly <- list()
    PCAresults$All <- list()
    PCAresults$modified <- list()
    for(i in 1:length(CellType_genes)){
        #message(paste0("Cell type: ", names(CellType_genes)[i]))
        Data2 <- dataset[dataset$GeneSymbol %in% CellType_genes[[i]],]
        # Remove genes with expression level bellow the noise threshold in 95%
        # of control samples. Ths step is done to ensure that the marker genes
        # can be detected in human bulk tissue
        geneContExp <- apply(Data2 %>% select(as.character(groupList[[contName]])), 1, function(x) quantile(x, 0.05))
        Data2 <- Data2[geneContExp > NoiseThreshold,]
        if (nrow(Data2) > 2){
            PCAresults$ControlOnly[[i]] <- Data2 %>%
                select(as.character(groupList[[contName]])) %>% t %>% prcomp(scale=T)
            rownames(PCAresults$ControlOnly[[i]]$rotation) <- Data2$GeneSymbol
            PCAresults$ControlOnly[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$ControlOnly[[i]])
            
            PCAresults$All[[i]] <- Data2 %>% 
                select(as.character(group$id)) %>% t %>% prcomp(scale=TRUE) 
            PCAresults$All[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$All[[i]])
            rownames(PCAresults$All[[i]]$rotation) <- Data2$GeneSymbol
            while (sum(PCAresults$All[[i]]$rotation[,1] > 0) < nrow(PCAresults$All[[i]]$rotation)){
                if(sum(PCAresults$All[[i]]$rotation[,1]) < 0){
                  PCAresults$All[[i]]$rotation[,1] <- -1*PCAresults$All[[i]]$rotation[,1]
                }
                minorGene <- rownames(PCAresults$All[[i]]$rotation)[PCAresults$All[[i]]$rotation[,1] < 0]
                Data2 %<>% filter(!GeneSymbol %in% minorGene)
                if(nrow(Data2) > 2){
                  rownames(Data2) <- Data2$GeneSymbol
                  PCAresults$All[[i]] <- Data2 %>% select(as.character(group$id)) %>% t %>% prcomp(scale=TRUE) 
                  PCAresults$All[[i]] <- correct_sign(names(CellType_genes)[i], PCAresults$All[[i]])
                } else {
                  warning(paste("WARNING: ", names(CellType_genes)[i], "has less than 3 genes after sign exclusion"))
                  x <- matrix(nrow=length(group$id), ncol=1)
                  rownames(x) <- group$id
                  colnames(x)="x"
                  PCAresults$ControlOnly[[i]] <- list(x)
                  names(PCAresults$ControlOnly[[i]]) <- "x"
                  PCAresults$All[[i]] <- list(x)
                  names(PCAresults$All[[i]]) <- "x"
                  PCAresults$modified[[i]] <- list(x)
                  names(PCAresults$modified[[i]]) <- "x"
                  break
                }
            }
            PCAresults$modified[[i]] <- PCAresults$All[[i]]$x %>% list
            PCAresults$modified[[i]] <- apply(PCAresults$modified[[i]][[1]],2,
                                        function(x) rescale(x,c(0,1))) %>% list
            names(PCAresults$modified[[i]]) <- "x"
        } else {
            warning(paste("WARNING: No genes for", names(CellType_genes)[i]))
            x <- matrix(nrow=length(group$id), ncol=1)
            rownames(x) <- group$id
            colnames(x)="x"
            PCAresults$ControlOnly[[i]] <- list(x)
            names(PCAresults$ControlOnly[[i]]) <- "x"
            PCAresults$All[[i]] <- list(x)
            names(PCAresults$All[[i]]) <- "x"
            PCAresults$modified[[i]] <- list(x)
            names(PCAresults$modified[[i]]) <- "x"
        }
    }
    for(i in 1:3){
        names(PCAresults[[i]]) <- names(CellType_genes)
    }
    return(PCAresults)
}

DESeq2RUN <- function(data, Meta, model){
    DESeqDS <- DESeqDataSetFromMatrix(countData = data, colData = Meta, design = model)
    DESeqOut <- DESeq(DESeqDS)
    return(DESeqOut)
}

GetDESeq2Results <- function(DESeqOut, coef, alpha=0.05, indepFilter=TRUE, geneNames){
    DEresults <- results(DESeqOut, name = coef, alpha = alpha, format = "DataFrame", independentFiltering = indepFilter)
    DEresults$GeneSymbol <- geneNames$gene_name[match(rownames(DEresults), geneNames$gene_id)]
    DEresults$EnsemblID <- rownames(DEresults)
    DEresults %<>% data.frame %>% filter(GeneSymbol != "") %>% dplyr::select(GeneSymbol, EnsemblID, log2FoldChange, pvalue, padj, everything())
    return(DEresults)
}

GetOneSidedPval <- function(ResultsObj, adjust="BH", logFCcol="log2FoldChange", GeneCol="GeneSymbol", pvalCol="pvalue"){
    DESeqResultsDF <- data.frame(ResultsObj)
    #Just for now - remove the duplicated genes (5 at this point) 
    DESeqResultsDF <- DESeqResultsDF[!duplicated(DESeqResultsDF$GeneSymbol),]
    DESeqResultsDF$DownPval <- apply(DESeqResultsDF %>% select(logFCcol, pvalCol), 1, function(x){
        if(x[1] < 0){
            x[2]/2
        } else {
            1-x[2]/2
        }
    })
    DESeqResultsDF$DownPvalAdj <- p.adjust(DESeqResultsDF$DownPval, "BH")
    DESeqResultsDF$UpPval <- apply(DESeqResultsDF %>% select(logFCcol, pvalCol), 1, function(x){
        if(x[1] > 0){
            x[2]/2
        } else {
            1-x[2]/2
        }
    })
    DESeqResultsDF$UpPvalAdj <- p.adjust(DESeqResultsDF$UpPval, "BH")
    rownames(DESeqResultsDF) <- as.character(DESeqResultsDF[[GeneCol]])
    return(DESeqResultsDF)
}

GetAnnoFiles <- function(platform){
    platform <- as.character(platform)
    if(length(list.files(pattern=platform)) == 0){
        download.file(paste0("https://gemma.msl.ubc.ca/annots/", platform, "_noParents.an.txt.gz"), destfile=paste0(platform, ".gz"))
    }
    if(length(list.files(pattern=platform)) > 1){
        print("Multiple annotation files matching the platform exist")
    } else {
        warning("Using existing annotation file, consider updating")
    }
    Anno_file <- read.table(paste0(platform, ".gz"), comment="#", header=T, quote='"', sep="\t")
    return(Anno_file)
}

GetAdjVal <- function(model, preserve){
    Resid <- resid(model)
    Coef <- t(coef(model)) %>% as.vector
    data <- model.frame(model)
    if (!preserve %in% colnames(data)) stop("preserve arguments must be a column in the model matrix")
    mod <- as.matrix(model.matrix(model)) 
    if ("factor" %in% class(data[[preserve]])) preserve <- paste0(preserve, levels(data[[preserve]])[length(levels(data[[preserve]]))])
    if (!preserve %in% colnames(mod)) stop("Something unexpected went wrong, sorry")
    Columns <- !grepl(paste0("^",preserve,"$"), colnames(mod), perl=TRUE)
    #new_values <- colMeans(mod)
    for (i in 2:length(Columns)) {
        if (Columns[i]) mod[,i] <- mean(mod[,i])
    }
    AdjValue <- (mod %*% Coef) + Resid
    return(AdjValue)
}

assocPCA <- function(M, variables, center=TRUE, scale=TRUE, ncomp=4, plot.return=FALSE, plot.alpha=0.05, plot.font.size=3, plot.colour="#c70039") {
    M <- as.matrix(M)
    pca <- prcomp(M, center=center, scale.=scale)
    ncomp <- min(ncomp, dim(pca$x)[2])
    var.exp <- summary(pca)$importance[2,1:ncomp]
    var.exp.str <- paste0(formatC(round(var.exp*100)), '%')
    res <- sapply(1:ncomp, function(i) {
                  unlist(lapply(1:ncol(variables), function(j) {
                             kept <- !is.na(variables[,j])
                             vars <- model.matrix(~.-1, variables[kept,j, drop=FALSE])
                             if (ncol(vars) == 2){vars <- vars[,1,drop=FALSE]}
                             ret <- sapply(1:ncol(vars), function(jp) {
                                        summary(lm(pca$x[kept,i] ~ vars[,jp]))$coefficients[2,4]
                                    })
                             names(ret) <- colnames(vars)
                             ret
                  }))
          })
    colnames(res) <- paste0('PC', seq(1:ncomp), ' (', var.exp.str, ')')
    res <- rbind(res,var.exp)
    if (plot.return==TRUE) {
        Data <- res %>% as.data.frame %>% rownames_to_column("variable") %>%
            #select(-var.exp) %>%
            gather(key='PC', value='p.value', -variable) %>%
            mutate(PC=factor(as.character(PC))) %>%
            mutate(sig=ifelse(p.value<plot.alpha, "SIG", "NOTSIG"), log10p=-log10(p.value)) %>% 
            mutate(p.lab=ifelse(p.value<(plot.alpha*2), format(p.value, digits=2), NA))
        Data$variable <- factor(Data$variable, levels=rev(colnames(variables)))
        ggplot(data=Data %>% filter(variable!="var.exp"), aes(x=PC, y=variable)) +
            geom_tile(aes(alpha=log10p, fill=sig), colour="white") +
            geom_text(aes(label=p.lab), size=plot.font.size) +
            scale_fill_manual(values=c(SIG=plot.colour, NOTSIG="black"), guide=FALSE) +
            scale_alpha(guide=FALSE, name="-log10(p-value)") +
            scale_x_discrete(position="top", expand=c(0,0)) +
            scale_y_discrete(expand=c(0,0)) + 
            theme_classic() +
            theme(axis.title.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.ticks=element_blank(),
                  axis.line=element_blank()) -> Plot
            return(Plot)
    } else {
        return(res)
    }
}

draw_treemap <- function(dtf, title="", vp=NULL) {
    treemap::treemap(
        dtf=dtf,
        index=c("Name", "members"), 
        vSize="NumGenes", 
        vColor="log10P",
        sortID="color",
        title=title,
        type="value",
        palette=rev(brewer.pal(10, "RdBu")),
        border.col=c("white", "black"),
        border.lwds=c(5,1),
        overlap.labels=1,
        bg.labels="#FFFFFFCC",
        force.print.labels=FALSE,
        vp=vp)
}

gene_pathway_matrix <- function(input){
    m <- table(input[[2]], input[[1]]) == 1
    m[,match(unique(input[[1]]), colnames(m))]
}
sim_matrix <- function (gmmat){
    #require(parallel)
    Kappas <- sapply(1:ncol(gmmat), function(p1) {
                  sapply(1:ncol(gmmat), function(p2) irr::kappa2(gmmat[,c(p1,p2)])$value)
    })
    rownames(Kappas) <- colnames(gmmat)
    colnames(Kappas) <- colnames(gmmat)
    return(Kappas)
}
jaccard_matrix <- function (gmmat){
    as.matrix(1 - dist(t(gmmat), method="binary"))
}
overlap_matrix <- function (gmmat){
    #require(parallel)
    Overlap <- sapply(1:ncol(gmmat), function(p1) {
                   sapply(1:ncol(gmmat), function(p2) {
                       sum(gmmat[,p1] & gmmat[,p2]) / min(sum(gmmat[,p1]), sum(gmmat[,p2]))
                   })
               })
    rownames(Overlap) <- colnames(gmmat)
    colnames(Overlap) <- colnames(gmmat)
    return(Overlap)
}
update_simmat <- function(simmat, gmmat, failed, passed){
    diag(simmat) <- 1
    #gene_pathway_mat <- gene_pathway_matrix(pathways)
    genes_in_cluster <- rep(0, length(rownames(gmmat)))
    names(genes_in_cluster) <- rownames(gmmat)
    # genes in the union of the two pathways
    merged <- unique(c(rownames(gmmat)[gmmat[,failed]],
                       rownames(gmmat)[gmmat[,passed]]))
    genes_in_cluster <- ifelse(names(genes_in_cluster) %in% merged, 1, 0)
    names(genes_in_cluster) <- rownames(gmmat)
    new_kappa_col <- sapply(1:nrow(simmat), function(x){
            irr::kappa2(data.frame(gmmat[,x],genes_in_cluster))$value
    })
    simmat[,passed] <- new_kappa_col
    simmat <- simmat[-which(rownames(simmat)==failed), -which(rownames(simmat)==failed)]
    diag(simmat) <- NA
    return(simmat)
}
cluster_pathways <- function(pathways, threshold=0.4, subsetsize=length(unique(pathways[[1]])),
                            lowThreshold=50, highThreshold=1000, method=c("kappa", "jaccard", "overlap")) {
    # check pathways input is data.frame
    stopifnot("data.frame" %in% class(pathways))
    # max number of pathways allowed are 1000, safety measure
    totPaths <- length(unique(pathways[[1]]))
    if (subsetsize > totPaths) {
        warning(paste0("Subset size chosen (", subsetsize, ") greater than the number of pathways (", totPaths, "), setting subsetsize <- ", totPaths))
    	subsetsize <- totPaths
    }
    if (subsetsize > 1000){
        warning("Maximum subset size allowed is 1000, setting subsetsize <- 1000")
    	subsetsize <- 1000
    }

    method <- method[1]
    
    # Filter to only top pathways
    pathways[,1] <- as.character(pathways[,1])
    pathways[,2] <- as.character(pathways[,2])
    pathways <- pathways[pathways[[1]] %in% head(unique(pathways[[1]]), n=subsetsize),]

    # membership of genes to pathways and pathway sizes
    membership_matrix <- gene_pathway_matrix(pathways)
    pSizes <- colSums(membership_matrix)

    # similarity of pathways
    if (method == "kappa") {
        simmat <- sim_matrix(membership_matrix)
    } else if (method == "jaccard") {
        simmat <- jaccard_matrix(membership_matrix)
    } else if (method == "overlap") {
        simmat <- overlap_matrix(membership_matrix)
    } else {
        stop(paste0("Method \"", method, "\" unknown"))
    }
    diag(simmat) <- NA

    pathwayNames <- colnames(simmat)   # pathway names
    if (any(pathwayNames != unique(pathways[[1]]))) stop("HORRIBLE ERROR")

    cluster <- matrix(rep(0, subsetsize*subsetsize), nrow=subsetsize)
    diag(cluster) <- 1
    colnames(cluster) <- pathwayNames
    rownames(cluster) <- pathwayNames

    iter <- 1
    while(max(simmat[upper.tri(simmat, diag=FALSE)]) >= threshold){
        message(paste0("Iteration ", iter, "..."))
        index <- which(simmat == max(simmat[upper.tri(simmat, diag = FALSE)]), arr.ind=TRUE)
        t_i = index[1,1] # first row index
        t_j = index[1,2] # first col index
        t_both <- c(t_i, t_j)
        pair_size <- c(pSizes[rownames(simmat)[t_i]], pSizes[rownames(simmat)[t_j]])
        pair_names <- c(rownames(simmat)[t_i], rownames(simmat)[t_j])
        if ( (pair_size[1] <= lowThreshold || pair_size[1] >= highThreshold) && (pair_size[2] <= lowThreshold || pair_size[2] >= highThreshold) ) {
            failed <- rownames(simmat)[max(t_both)]
            passed <- rownames(simmat)[min(t_both)]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else if ( pair_size[1] <= lowThreshold || pair_size[1] >= highThreshold ) {
            failed <- rownames(simmat)[t_i]
            passed <- rownames(simmat)[t_j]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else if ( pair_size[2] <= lowThreshold || pair_size[2] >= highThreshold ) {
            failed <- rownames(simmat)[t_j]
            passed <- rownames(simmat)[t_i]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else {
            failed <- rownames(simmat)[max(t_both)]
            passed <- rownames(simmat)[min(t_both)]
            simmat <- update_simmat(simmat, membership_matrix, failed, passed)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        }
        iter <- iter + 1
    }

    # remodel "cluster" structure to extract the titles
    cluster_idx <- which(diag(cluster)==1)
    clusters <- apply(cluster[,cluster_idx], 2, function(x) names(which(x==1)))
    data.frame(title=rep(names(clusters), sapply(clusters, length)),
               member=unlist(clusters, , FALSE),
               stringsAsFactors=FALSE,
               row.names=NULL)
}

RunGSR <- function(scores, scoreColumn, annotation, bigIsBetter=FALSE, logTrans=TRUE, aspects=c("B", "M", "C"), iterations=100000) {
    gsr(scores=res_df, scoreColumn=scoreColumn, bigIsBetter=bigIsBetter, logTrans=logTrans,
        annotation=annotation, aspects=aspects, iterations=iterations)$results %>%
    arrange(Pval)
}

treePlotSummary <- function(gsea_up, gsea_down, pval=0.05, maxCluster=200, maxPlot=100, title="") {
    clusters <- list()
    # UP
    obj <- gsea_up %>% filter(CorrectedPvalue < pval) %>%
        arrange(CorrectedPvalue) %>% select(Name, GeneMembers, CorrectedPvalue, NumGenes)
    pathways <- apply(obj, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
    names(pathways) <- obj$Name
    pathways <- qdapTools::list2df(pathways, col1="gene", col2="Name")[,c(2,1)]
    
    clustersUP <- cluster_pathways(pathways, method="overlap", subsetsize=maxCluster)
    names(clustersUP) <- c("Name", "members")
    clustersUP <- left_join(clustersUP, obj %>% select(-GeneMembers, members=Name), by="members") %>%
        mutate(log10P=-log10(CorrectedPvalue)) %>%
        arrange(CorrectedPvalue)
    clusters[["UP"]] <- clustersUP
    
    # DOWN
    obj <- gsea_down %>% filter(CorrectedPvalue < pval) %>%
        arrange(CorrectedPvalue) %>%
        select(Name, GeneMembers, CorrectedPvalue, NumGenes)
    pathways <- apply(obj, 1, function(x) unlist(strsplit(x[[2]], '\\|')))
    names(pathways) <- obj$Name
    pathways <- qdapTools::list2df(pathways, col1="gene", col2="Name")[,c(2,1)]
    
    clustersDOWN <- cluster_pathways(pathways, method="overlap", subsetsize=maxCluster)
    names(clustersDOWN) <- c("Name", "members")
    clustersDOWN <- left_join(clustersDOWN, obj %>% select(-GeneMembers, members=Name), by="members") %>%
        mutate(log10P=-log10(CorrectedPvalue)) %>%
        arrange(CorrectedPvalue)
    clusters[["DOWN"]] <- clustersDOWN
    
    draw_treemap(dtf=bind_rows(clusters[["UP"]],
                               clusters[["DOWN"]] %>% mutate(log10P=-log10P)) %>%
                         arrange(desc(abs(log10P))) %>%
                         head(n=maxPlot),
                 title=title)

}
