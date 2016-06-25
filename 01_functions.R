unsupervisedPCA <- function(inMatrix, #specify columns of tlRPKMs matrix to be compared
                               outPDF='~/Desktop/uniquenessHeatmap.pdf', #name graph and define file path
                               forceGenes = NULL, #a spicific set of genes to use for the PCA
                               topN=100, #the top N most unique genes that will be selected per sample
                               trimThresh=10, #genes will be skipped if there isn't at least one cell per embryo with an RPKM above this value
                               repSets="byColName", #if the samples are not to be grouped by standard replicates, specify a list of vectors, each grouping the column numbers (ie. list(1:4, 5:8, 9:12))
                               buff=25, #this number will be added to all log2(rpkms) before to dampen the low-expression differences between samples
                               crossEmbryoCommon=3, #how many embryos need to have this gene as top ranked for the gene to get plotted
                               falseNeg=1, #how many replicates can not have a gene for that gene to make it past the first filtration step
                               plotWidth=5, #width of pdf
                               plotHeight=5, #height of pdf
                               returnGenes=F, #if true, spits out the genes used for PCA
                               labelSamples=F #should the tiny 'st###' labels be written over each dot?
){


    #make repSeps list, if it is byColName
    suppressWarnings(if(repSets=="byColName"){
        repSets<-makeRepSetList(inMatrix)
    })
    
    #subsetMatrix <- trim big huge matrix down to the genes that have at least one cell over the threshold in each replicate
    subsetMatrix <- multiFilterCutoff(inMatrix=inMatrix, cutValue=trimThresh, groups=repSets, falseNeg=falseNeg)
    #this matrix should be transfered into log2(n+2)
    subsetMatrix <- log2(subsetMatrix+4)
    #report the size of this new matrix
    print(paste("After rejecting genes without ", trimThresh, 
                "RPKM in at least one cell in at least ", (length(repSets)-falseNeg),
                " embryos, the remaining matrix has ", dim(subsetMatrix)[1], 
                " genes and ", dim(subsetMatrix)[2], " samples.",sep=""))
    
    
    sharedGeneList <- NULL
    if(length(forceGenes)>0){
        sharedGeneList <- intersect(forceGenes, rownames(subsetMatrix))
    } else {
        
        #create a matrix with uniqueness scores for each gene and sample 
        uniqueScoreMatrix <- RPKMmatr2uniqMatr(subsetMatrix, repSets=repSets, buff=buff)
        #turn that uniqueness matrix into ranked lists for each sample, as long as the specified cutoff
        topNlist <- uniqMatr2topNlist(uniqueScoreMatrix, topNcut=topN, repSets=repSets)
        sampleLabels <- makeLabelVector(inMatrix)
        
        #make a master vector of all the informative genes
        masterGeneList = NULL
        for(i in 1:length(topNlist)){
            masterGeneList <- union(masterGeneList, names(topNlist[[i]]))
        }
        print(paste("The union of the ", topN, " genes most unique to every cell adds up to ", length(masterGeneList), sep=""))
        #make a subset that only includes genes that come up in more than one embryo
        sharedGeneList <- selectCommons(topNlist, repSets, crossEmbryoCommon)   
        print(paste("Of those genes, ", length(sharedGeneList), " genes are shared by ", crossEmbryoCommon, " or more embryos", sep=""))
    }
    princ=NULL
    print(paste("Performing PCA using gene expression data from ", length(sharedGeneList), " genes", sep=""))
    if(returnGenes){
        return(sharedGeneList)
    }
    princ <- prcomp(subsetMatrix[sharedGeneList,])
    repCols<- makeRepColors(subsetMatrix) 
    #return(princ)
    library(scales)
    pdf(outPDF, width=plotWidth, height=plotHeight, family="Times")
    plot(princ$rotation[,1], 
         princ$rotation[,2],
         pch=16,
         cex=2, col=alpha(repCols,0.5), 
         xaxt='n', 
         yaxt='n',
         ann=FALSE)
    if(labelSamples){
        text(princ$rotation[,1], princ$rotation[,2],
             colnames(subsetMatrix), cex=.5, col=repCols)

        }
    dev.off()
}

makeRepSetList <- function(inMatrix){
    sampNums <- which(row.names(tlDesign) %in% colnames(inMatrix))
    repLevels <- levels(as.factor(tlDesign$replicate[sampNums]))
    repSets <- list()
    for(i in 1:length(repLevels)){
        colNums = which(tlDesign$replicate %in% repLevels[i])
        repSets[[i]]=which(sampNums %in% colNums)
    }
    return(repSets)
}

multiFilterCutoff <- function(inMatrix, cutValue, groups, falseNeg="none") {
    if(falseNeg=="none"){
        posEmbryos = length(groups)
    } else {
        posEmbryos = length(groups)-falseNeg
    }
    cutoff = NULL
    outMatrix = matrix()
    #for every gene in the matrix
    for (gene in 1:dim(inMatrix)[1]) {
        counter = 0
        for (grp in 1:length(groups)) {
            if (max(inMatrix[gene,groups[[grp]]]) > cutValue) {
                counter = counter+1
            }
        }
        cutoff[gene] = counter
    }
    outMatrix <- inMatrix[which(cutoff >= posEmbryos),]
    return(outMatrix)
}

RPKMmatr2uniqMatr <- function(inMatrix, buff=3, repSets="byColName"){
    #make a new matrix the same size as the inMatrix, for uniqueness scores
    uniqueScoreMatrix <- matrix(NA, nrow=dim(inMatrix)[1], ncol=dim(inMatrix)[2],
                                dimnames=dimnames(inMatrix))
    #if repSets is the default, generate a list representing the replicates
    #otherwise, just use the list entered in the command to represent the replicates
    suppressWarnings(if(repSets=="byColName"){
        repSets <- makeRepSetList(inMatrix)
    })
    #for each repSet
    for(rep in 1:length(repSets)){
        #for each gene
        for(gene in 1:dim(inMatrix)[1]){
            #run the pre-existing uniqueScore function on all cells
            uniqueScoreMatrix[gene,repSets[[rep]]] <- uniqueScore(inMatrix[gene,repSets[[rep]]], buff=buff)
        }       
    }
    uniqueScoreMatrix[is.na(uniqueScoreMatrix)] <- 0
    return(uniqueScoreMatrix)
}

uniqueScore <- function(vect, buff=3){
    #make temporary vector to find the haves/havenots cutoff
    vect <- as.numeric(vect+buff)
    tempVect = NULL
    for(permut in 1:(length(vect)-1)){
        havesMean <- mean(vect[order(vect, decreasing=TRUE)[1:permut]])
        havenotsMean <- mean(vect[order(vect, decreasing=TRUE)[(permut+1):length(vect)]])
        tempVect[permut] <- havesMean/havenotsMean
    }
    #lower score for genes that are shared by multiple cell types
    theLine <- order(tempVect, decreasing=TRUE)[1]
    theTotal <- length(vect)
    theNorm <- (1-((theLine/theTotal)^2))
    theScore<-(tempVect[theLine]*theNorm)
    #with that cutoff, assign a uniqueScore to the haves, and an NA to the havenots
    outVect = NULL
    outVect[order(vect, decreasing=TRUE)[1:theLine]]<- theScore
    outVect[order(vect, decreasing=TRUE)[(theLine+1):length(vect)]] <- NA
    return(outVect)
}

uniqMatr2topNlist <- function(uniqueScoreMatrix, topNcut, repSets){
    topNlist<-list()
    #for each column, return the top N genes, and their scores
    for(colmn in colnames(uniqueScoreMatrix)){
        topNlist[[colmn]]<-uniqueScoreMatrix[order(uniqueScoreMatrix[,colmn], decreasing=TRUE)[1:topNcut],colmn]
    }
    return(topNlist)
}

makeLabelVector <- function(inMatrix){
    sampNums = NULL
    for(name in colnames(inMatrix)){
        sampNums <- c(sampNums, which(row.names(tlDesign) %in% name))
    }
    labelVect <- paste(tlDesign$stage[sampNums], " (", tlDesign$replicate[sampNums], ")", sep="")
    return(labelVect)
}

selectCommons <- function(topNlist, repSets, crossEmbryoCommon){
    sharedGeneScore <- vector()
    for(embryo in 1:length(repSets)){
        embryoVect <- repSets[[embryo]]
        embryoGenes = vector()
        for(cell in embryoVect){
            embryoGenes <- union(embryoGenes, names(topNlist[[cell]]))
        }
        sharedGeneScore <- c(sharedGeneScore, embryoGenes)
    }
    overTwo = vector()
    for(indiv in levels(as.factor(sharedGeneScore))){
        if(length(which(sharedGeneScore==indiv))>=crossEmbryoCommon){
            overTwo = c(overTwo, indiv)
        }
    }
    return(overTwo)
}

makeRepColors <- function(subsetMatrix, repColors=c("dark blue", "red", "orange", "purple", "brown", 
                                                    "darkcyan", "darksalmon", "darkgreen", "violetred")){
    repCols=NULL
    for(cell in colnames(subsetMatrix)){
        thisCol <- repColors[tlDesign[cell,"replicate"]]
        repCols<- c(repCols, thisCol)
    }
    return(repCols)
}

ERdiffExV1 <- function(matr=tlRPKMs[,typeCells(c("AB", "P1"))], metaD = tlDesign,
                       mastMatr = tlRPKMs, # change if using erin or yanai's data
                       sig=.1, thresh=0, #thresh is in RPKM. genes to throw out before doing analysis
                       plus=0, postThresh=0, #genes to throw out when plotting, post-analysis.
                       group="default", #otherwise, provide (as factors) the conditions for comparison for each sample
                       plotWidth = 8, plotHeight = 8,
                       overlay = NULL, # this is a vector of gene names that will be circled in teal
                       label = NULL, # this is a vector of gene names that will be labeled in white
                       to="window",     #if 'to' is "list", it will output a ranked list,
                       #if it is 'plottable data', it will output tab, sig, plus, samps, overlay, label, and postThresh
                       # if it is anything else it will output that pdf.
                       visibleLabs=F, # T grays out plot and writes labels big
                       opac = .4,
                       forceY = F, # ylim for MA plot
                       by = "PValue") { #when exporting a list, should it be ranked 
    # by "PValue", "logCPM", or "logFC".
    
    
    suppressWarnings(if(group=="default"){
        group <- factor(metaD[colnames(matr), "ID"])
    })
    group <- factor(group)
    samps <- levels(group)
    
    cellsUsed <- colnames(matr)
    threshWhich <- which(apply(matr, 1, function(x){length(x[which(x>thresh)])>0}))
    #threshWhich <- which(apply(mastMatr[,cellsUsed], 1, function(x){summary(x)["Max."]})>thresh)
    matrThreshed <- matr[threshWhich,]
    
    postThreshGenes <- which(apply(matr, 1, function(x){length(x[which(x>postThresh)])>1}))
    postThresh <- names(postThreshGenes)
    
    overlay = as.character(overlay)
    label = as.character(label)
    
    #firstThreshed <- which(apply(matr[,group==samps[1]], 1, mean) > thresh+plus)
    #secondThreshed <- which(apply(matr[,group==samps[2]], 1, mean) > thresh+plus)
    #    matrThreshed <- matr[union(firstThreshed, secondThreshed),]
    
    x <- matrThreshed+plus
    library(edgeR)
    
    y <- DGEList(counts=x+plus, group=group)
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    tab <- et$table
    #return(et$table)
    
    if(to=="list"){
        outList = list()
        #make a FC down subset
        downSub <- tab[which(tab$PValue < sig & tab$logFC < 0 & rownames(tab) %in% postThresh),]
        outList[[samps[1]]] = downSub[order(downSub[,by]),]
        upSub <- tab[which(tab$PValue < sig & tab$logFC > 0 & rownames(tab) %in% postThresh),]
        outList[[samps[2]]] = upSub[order(upSub[,by]),]
        return(outList)
    }
    else if(to=="plottable data"){
        return(list(tab, sig, plus, samps, overlay, label, postThresh))
    }
    else if(to=="window"){
        outList = list(tab, sig, plus, samps, overlay, label, postThresh)
        #return(outList)
        plotTabSig(tab, sig, plus, samps, overlay, label, postThresh, opac, forceY=forceY)
    }else{
        pdf(to, plotWidth, plotHeight)
        plotTabSig(tab, sig, plus, samps, overlay, label, postThresh, opac, visibleLabs=visibleLabs, forceY=forceY)
        dev.off()
    }
}  


plotTabSig <- function(tab, sig, plus, samps, overlay, label, postThresh, opac, visibleLabs=F, cellHead="", genesHead="", forceY=F){
    library(scales)
    #Plot all genes in gray
    if(forceY){
        plot(tab$logCPM, tab$logFC, ylim=forceY,
             col=alpha("black", opac), pch=16, 
             xlab="Log ( Average Count )", ylab="Log ( Fold Change )", 
             main=cellHead)
    } else {
        plot(tab$logCPM, tab$logFC, 
             col=alpha("black", opac), pch=16, 
             xlab="Log ( Average Count )", ylab="Log ( Fold Change )", 
             main=cellHead)
    }
    mtext(genesHead)
    #Set a color scheme along that range of Pvalues
    colFunc <- colorRampPalette(c("red", "purple"))
    
    #Draw points, and color them by their Pvalue, using the color scheme
    binBreaks = vector()
    for(i in 1:11){binBreaks[i]=((i-1)*sig/10)}
    for(i in 1:10){
        lowerSig = binBreaks[i]
        upperSig = binBreaks[i+1]
        sigSubset = tab[which(tab$PValue >= lowerSig & tab$PValue<upperSig & rownames(tab) %in% postThresh),]
        points(sigSubset[,"logCPM"], sigSubset[,"logFC"],
               col=alpha(colFunc(10)[i], 2*opac), pch=16)
    }
    
    usr <- par( "usr" )
    tempLeft = usr[1]
    tempRight = usr[2]
    tempTop = usr[4]
    tempBottom = usr[3]
    
    text(tempLeft+(tempRight-tempLeft)/25, tempBottom+(tempTop-tempBottom)/25, adj = c( 0, 0 ),
         samps[1])
    text(tempLeft+(tempRight-tempLeft)/25, tempTop-(tempTop-tempBottom)/25, adj = c( 0, 1 ),
         samps[2])
    
    # overlay teal circles
    points(tab[overlay, "logCPM"], 
           tab[overlay,"logFC"], 
           col=alpha("#000000", 1), cex=1.2)
    points(tab[overlay, "logCPM"], 
           tab[overlay,"logFC"], 
           col=alpha("cyan", 1), cex=1.1)
    
    if(visibleLabs){
        rect(tempLeft, tempBottom, tempRight, tempTop, col = "#FFFFFF50")
    }
    
    # overlay white gene names
    if(length(label) > 0){
        if(visibleLabs){
            text(tab[label,"logCPM"], tab[label, "logFC"], labels=label, cex=.6, col="black")
        } else {
            text(tab[label, "logCPM"], 
                 tab[label, "logFC"], 
                 labels=label, cex =.1, col = "white")
        }
    }
}

ABxxxPCA <- function(geneList, outPDF, cells="all", 
                     width=5, height=5, labelSamps=F, splitReps=F){
    # define cells to use
    suppressWarnings(if(cells=="all"){
        cells <- row.names(tlDesign[which(tlDesign$stage=="16-cell" 
                                          & tlDesign$quality=="Use" 
                                          & !(tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2"))),])
    })
    princ <- prcomp(prop.table(as.matrix(tlRPKMs[geneList,cells]+1), margin=2))
    repSets<-makeRepSetList(tlRPKMs[geneList,cells])
    repColors<- makeRepColors(tlRPKMs[geneList,cells]) 
    #return(princ)
    
    # determine axis limits
    lims <- c(range(princ$rotation[,1]), range(princ$rotation[,2]))
    
    pdf(outPDF, width, height, family="Times")
    plot(princ$rotation[,1], 
         princ$rotation[,2],
         xlim=lims[1:2], ylim=lims[3:4],
         pch=16,
         cex=2, col=alpha(repColors,0.5), 
         xlab="Principle Component 1", 
         ylab="Principle Component 2")
    if(labelSamps){
        text(princ$rotation[,1], princ$rotation[,2],
             cells, cex=.5, col=repColors,
             adj = c(0,1))
    }
    dev.off()
    if(splitReps){
        for(repl in 1:length(repSets)){
            dots <- repSets[[repl]]
            outName <- strsplit(outPDF, split = ".pdf")[[1]]
            outName <- paste(outName, "_", repl, ".pdf", sep = "")
            pdf(outName, width, height, family="Times")
            plot(princ$rotation[dots,1],
                 princ$rotation[dots,2],
                 xlim=lims[1:2], ylim=lims[3:4],
                 pch=16,
                 cex=2, col=alpha(repColors[dots],0.5), 
                 xlab="Principle Component 1", 
                 ylab="Principle Component 2")
            dev.off()
        }
    }
}


mapExpression <- function(gene, # can be single gene or c() vector of genes
                          mapType = "cell", # other option = lineage
                          logTrans=F, # log10 transform RPKMs?
                          to="window", # anything else will be interpreted as file path
                          legend=T,
                          legendBottom=F,
                          stat = "avg", # can also be "med"
                          fillRange="standard", # otherwise, a c(min, max) vector
                          safeCell = T,
                          width="default", height="default"
){

    # average gene expression for each cell type
    cellNames = NULL
    if(safeCell){cellNames = unique(cellMapShapes[["nameRef"]]$safeDesignNames)}
    else{cellNames = unique(cellMapShapes[["nameRef"]]$designNames)}
    cellNames = cellNames[!is.na(cellNames)]
    exprFrame = as.matrix(data.frame(row.names = cellNames, avgExp = rep(0,length(cellNames))))
    # define statistical type
    matr = NULL
    if (stat=="avg" & safeCell==T) {matr = tlSafeCellAvg}
    if (stat=="avg" & safeCell==F) {matr = tlCellAvg}
    if (stat=="med" & safeCell==T) {matr=tlSafeCellMed}
    if (stat=="med" & safeCell==F) {matr=tlCellMed}
    for(cellType in cellNames){
        if(length(gene)==1){
            exprFrame[cellType, "avgExp"] = matr[gene,cellType]
        } else {
            exprFrame[cellType, "avgExp"] = mean(matr[gene, cellType])
        }
    }
    if(logTrans){
        exprFrame[,"avgExp"] <- log10(exprFrame[,"avgExp"]+1)
    }
    
    # make new mappable dataframe cell names that correspond to my IDs
    if(mapType=="cell"){
        origMatr <- cellMapShapes[["mapAllCells"]]
        tempMatr <- cellMapShapes[["mapAllCells"]]
        safeMatr <- cellMapShapes[["mapAllSafeCells"]]
    } else if (mapType=="lineage"){
        origMatr <- cellMapShapes[["mapLineage"]]
        tempMatr <- cellMapShapes[["mapLineage"]]
        safeMatr <- cellMapShapes[["mapLineage"]]
    } else if (mapType=="horz"){
        origMatr <- cellMapShapes[["mapCellsHorz"]]
        tempMatr <- cellMapShapes[["mapCellsHorz"]]
        safeMatr <- cellMapShapes[["mapSafeHorz"]]
    }
    tempMatr$id <- as.character(tempMatr$id)
    for(rowNum in 1:length(tempMatr$long)){
        # switch cell type name to ID's I have info on
        tempMapName = tempMatr[rowNum,"id"]
        tempDesignName = NULL
        if(safeCell){ tempDesignName = cellMapShapes[["nameRef"]][which(cellMapShapes[["nameRef"]]$mapNames==tempMapName),"safeDesignNames"]}
        else {tempDesignName = cellMapShapes[["nameRef"]][which(cellMapShapes[["nameRef"]]$mapNames==tempMapName),"designNames"]}
        tempMatr[rowNum,"id"] = tempDesignName
    }
    
    # add the average expression data to that dataframe
    tempMatr$expr = 0
    for(cellType in cellNames){
        exprValue = exprFrame[cellType, "avgExp"]
        tempMatr[which(tempMatr$id == cellType), "expr"] = exprValue
    }
    
    # clear polygons off the tempMatr that have no cell type
    tempMatr <- tempMatr[which(!is.na(tempMatr$id)),]
    
    # plot that matrix onto the cell map
    
    # set legend title to log or not
    legendTitle = NULL
    if(logTrans){ legendTitle <- "log10\n(RPKM)"} else { legendTitle <- "RPKM"}
    
    # set the heatmap range differently if desired
    suppressWarnings(
        if(fillRange=="standard"){
            fillRange = c(min(tempMatr[,"expr"]), max(tempMatr[,"expr"]))
        } else if(fillRange=="lowExpr"){
            if(logTrans){
                fillRange = c(min(tempMatr[,"expr"]), 1.7)
            } else {fillRange=c(min(tempMatr[,"expr"]), 50)}
        }
    )
    
    c <- ggplot()+
        geom_polygon(data=tempMatr, aes(x=long, y=lat, group=id, fill=expr)) +
        coord_equal(ratio=1)+
        theme_nothing(legend=legend)+
        scale_fill_continuous(name=legendTitle, low = "white", high = "#005c99", limits=fillRange,
                              guide=guide_legend(reverse=FALSE))
    if(legendBottom){
        c <- c+
            theme(legend.position="bottom")
    }
    
    if(safeCell){
        c <- c+
            geom_polygon(data=origMatr, aes(x=long, y=lat, group=id), fill=NA, color=alpha("black", .15))+
            geom_polygon(data=safeMatr, aes(x=long, y=lat, group=id), fill=NA, color="black")
    } else {
        c <- c+
            geom_polygon(data=origMatr, aes(x=long, y=lat, group=id), fill=NA, color="black")
    }
    
    if(to=="window"){
        c
    } else {
        if(mapType=="cell"){
            if(width=="default"){width<-3}
            if(height=="default"){height<-6}
            ggsave(c, file = to, width = width, height = height)
        } else if (mapType=="lineage"){
            if(width=="default"){width<-6}
            if(height=="default"){height<-3}
            ggsave(c, file = to, width = width, height = height)
        } else if (mapType=="horz"){
            if(width=="default"){width<-8}
            if(height=="default"){height<-2}
            ggsave(c, file = to, width = width, height = height)
        }
    }
}


mapStats <- function(dataIn, # one column name of stats matrix
                     mapType = "cell", # other option = "lineage"
                     logTrans=F, # log10 transform RPKMs?
                     to="window", # anything else will be interpreted as file path
                     legend=T, 
                     fillRange = "standard", # otherwise, a vector of 2 values: min and max
                     legendTitle = "stat",
                     safeCell = T
){
    
    if(legendTitle=="stat"){
        legendTitle <- dataIn
    }    
    
    dataIn <- tlCellStats[,dataIn]
    names(dataIn) <- rownames(tlCellStats)
    
    # transform data
    if(logTrans){
        dataIn <- log10(dataIn+1)
    }
    
    if(mapType=="cell"){
        origMatr <- cellMapShapes[["mapAllCells"]]
        tempMatr <- cellMapShapes[["mapAllCells"]]
        safeMatr <- cellMapShapes[["mapAllSafeCells"]]
    } else if (mapType=="lineage"){
        origMatr <- cellMapShapes[["mapLineage"]]
        tempMatr <- cellMapShapes[["mapLineage"]]
    } else if (mapType=="horz"){
        origMatr <- cellMapShapes[["mapCellsHorz"]]
        tempMatr <- cellMapShapes[["mapCellsHorz"]]
        safeMatr <- cellMapShapes[["mapSafeHorz"]]
    }
    
    # rename the map to match the cellTypes I have ID's
    tempMatr$id <- as.character(tempMatr$id)
    for(rowNum in 1:length(tempMatr$long)){
        tempMapName = tempMatr[rowNum,"id"]
        tempDesignName = NULL
        if(safeCell){tempDesignName = cellMapShapes[["nameRef"]][which(cellMapShapes[["nameRef"]]$mapNames==tempMapName),"safeDesignNames"]}
        else {tempDesignName = cellMapShapes[["nameRef"]][which(cellMapShapes[["nameRef"]]$mapNames==tempMapName),"designNames"]}
        tempMatr[rowNum,"id"] = tempDesignName
    }
    
    # add the average expression data to that dataframe
    tempMatr$stat = 0
    for(cellType in safeCellTypes){
        statValue = dataIn[cellType]
        tempMatr[which(tempMatr$id == cellType), "stat"] = statValue
    }
    
    # clear polygons off the tempMatr that have no cell type
    tempMatr <- tempMatr[which(!is.na(tempMatr$id)),]
    
    # plot that matrix onto the cell map
    
    suppressWarnings(
        if(fillRange=="standard"){
            fillRange = c(min(dataIn), max(dataIn)*1.1)
            
        } else if(fillRange=="skewed"){
            if(logTrans){
                trueMin = 10^(min(dataIn))
                fillRange = c(log10(.75*trueMin), max(dataIn)*1.1)
                
            } else {
                fillRange = c(.75*min(dataIn), max(dataIn)*1.1)
            }
        } else if(fillRange=="zeroed"){
            fillRange = c(0, max(dataIn)*1.1)
        }
    )
    
    c <- ggplot()+
        geom_polygon(data=tempMatr, aes(x=long, y=lat, group=id, fill=stat)) +
        scale_fill_continuous(name=legendTitle, low = "white", high = "#cc0000", limits=fillRange,
                              #scale_fill_distiller(name=legendTitle, palette = "Reds", limits=fillRange, 
                              breaks=(pretty_breaks(3)((10*min(dataIn)):(10*max(dataIn))))/10, guide=guide_legend(reverse=FALSE)) +
        coord_equal(ratio=1)+
        theme_nothing(legend=legend)+
        theme(legend.text=element_text(size=16),
              legend.title=element_text(size=20, face="plain"))
    
    if(mapType=="cell" & legend==T){
        c <- c+
            theme(legend.position="bottom", 
                  legend.direction="vertical")
    }
    
    if(safeCell){
        c <- c+
            geom_polygon(data=origMatr, aes(x=long, y=lat, group=id), fill=NA, color=alpha("black", .15)) +
            geom_polygon(data=safeMatr, aes(x=long, y=lat, group=id), fill=NA, color="black")
    }
    else {
        c <- c+
            geom_polygon(data=origMatr, aes(x=long, y=lat, group=id), fill=NA, color="black")
    }
    
    if(to=="window"){
        c
    } else {
        if(mapType=="cell"){
            ggsave(c, file = to, width = 4, height = 8)
        } else if (mapType=="lineage"){
            ggsave(c, file = to, width = 4, height = 2)
        } else if (mapType=="horz"){
            ggsave(c, file = to, width=8, height = 2)
        }
        
    }
}

typeCells <- function(cells){
    return(rownames(tlDesign[which((tlDesign$ID %in% cells | tlDesign$safeID %in% cells) & tlDesign$quality=="Use"),]))
}

makeEnrichTable <- function(constant=1){
    enrichMatr <- matrix(0, ncol=3, nrow=length(all2cellFCgenes), 
                         dimnames = list(all2cellFCgenes, c("Nishimura", "Hashimshony", "Tintori")))
    for(gene in all2cellFCgenes){
        for(study in c("Nishimura", "Hashimshony", "Tintori")){
            enrichMatr[gene, study] <- all2cellFCs[gene,study]*all2cellCPMs[gene,study]^constant
        }
    }
    return(enrichMatr)
}

ABcompPlot <- function(firstStudy, secondStudy, # Hashimshony, Nishimura, or Tintori
                       constant=1, axisLim="default"){
    # plot all genes
    library("scales")
    if(axisLim=="default"){
        axisLim <- max(abs(range(enrichMatr)))*1.05
    }
    plot(enrichMatr[,firstStudy], enrichMatr[,secondStudy],
         xlab = paste("AB/P1 enrichment from\n", firstStudy, " et al.", sep=""), 
         ylab = paste("AB/P1 enrichment from\n", secondStudy, " et al.", sep=""),
         pch=16, col= alpha("black", .3), xlim=c(-axisLim,axisLim), ylim=c(-axisLim, axisLim)
    )
    abline(h=0)
    abline(v=0)
    # add linear regression lines
    tempLinR <- lm(enrichMatr[,firstStudy] ~ enrichMatr[,secondStudy])
    abline(tempLinR[[1]][[1]], tempLinR[[1]][[2]], col=alpha("#000000", .3))
    
    # set up correlation calculation
    corVect <- c("Total"=0, "Tintori"=0, "Nishimura"=0, "Hashimshony"=0)
    corVect["Total"] <- round(cor(enrichMatr[,firstStudy], enrichMatr[,secondStudy]), digits = 2)
    colVect <- c("Total"="#000000", "Tintori"="#FF4100", "Nishimura"="#C000FF", "Hashimshony"="#00BEFF")
    sigGeneVect <- list("Total"=all2cellFCgenes, "Tintori"=my2cellIDgenes, "Nishimura"=erinIDgenes, "Hashimshony"=hashimIDgenes)
    
    # put linear regression lines and  circles and calculate corr coef 
    # for just the ones each study called significant
    for(tempStudy in c("Tintori", "Nishimura", "Hashimshony")){
        if(tempStudy %in% c(firstStudy, secondStudy)){
            subEnrichMatr <- enrichMatr[intersect(sigGeneVect[[tempStudy]], all2cellFCgenes),]
            
            tempLinR <- lm(subEnrichMatr[,firstStudy] ~ subEnrichMatr[,secondStudy])
            abline(tempLinR[[1]][[1]], tempLinR[[1]][[2]], col=alpha(colVect[tempStudy], .3))
            
            points(subEnrichMatr[,firstStudy], subEnrichMatr[,secondStudy],
                   pch=16, col=alpha(colVect[tempStudy], .3))
            
            corVect[tempStudy] = round(cor(subEnrichMatr[,firstStudy], subEnrichMatr[,secondStudy]), digits = 2)
        }
    }
    # plot correlation coefficients
    usr <- par( "usr" )
    tempLeft = usr[1]
    tempRight = usr[2]
    tempTop = usr[4]
    tempBottom = usr[3]
    txt <- c(paste("R = ", corVect["Total"], sep = ""),
             paste("R = ", corVect[firstStudy], sep = ""),
             paste("R = ", corVect[secondStudy], sep = ""))
    txtCol <- colVect[c("Total", firstStudy, secondStudy)]
    thisy <- tempTop-(tempTop-tempBottom)/10
    for(i in 1:3){
        text(tempLeft+(tempRight-tempLeft)/10, thisy, txt[i], col=txtCol[i], adj=c(0,1))
        thisy <- thisy-1.75*strheight(txt[i])
    }
}

rearrangeNameMatrix <- function(nameMatrix, priorityVector){
    # rearrange each row of the gene name matrix to have genes in the priority vector come first (and be what genes get renamed to)
    newMatrix <- matrix(NA, ncol = ncol(nameMatrix), nrow = nrow(nameMatrix))
    for(i in 1:dim(nameMatrix)[1]){
        leadGene <- NULL
        leadGene <- nameMatrix[i,which(nameMatrix[i,] %in% priorityVector)]
        tempLen <- length(unique(nameMatrix[i,]))
        if(length(leadGene)>0){
            newMatrix[i,1:tempLen] <- unique(c(leadGene, nameMatrix[i,]))
        } else { newMatrix[i,] <- nameMatrix[i,]}
    }
    return(newMatrix)
}

renameGeneVect <- function(oldGeneVect, nameMatrix = geneNameMatr){
    # rename genes in a vector so that they match the first column of the geneNameMatr
    newGeneVect <- oldGeneVect
    for(i in 1:length(oldGeneVect)){
        if(!newGeneVect[i] %in% nameMatrix[,1]){
            tempCoords <- which(nameMatrix == oldGeneVect[i], arr.ind = T)
            if(length(tempCoords)>0){
                newGeneVect[i] <- nameMatrix[tempCoords[1],1]
            }
        }
    }
    return(newGeneVect)
}

exportGeneExpMatrix <- function(outFile="04_analysesOut/SuppFig3/matrix.txt", 
                                geneSet=rownames(tlRPKMs), hm2=F, avgOn=F, log=T, plus1 = F){
    # filter cells
    if(avgOn){
        tempCellSubMatr = tlSafeCellAvg[geneSet,]
    } else {
        tempCellSubMatr = tlRPKMs[geneSet,rownames(tlDesign[which(tlDesign$quality=="Use"),])]
        
        # replace stNNN names with descriptive names
        newColNames = NULL
        for(colmn in 1:dim(tempCellSubMatr)[2]){
            sample = colnames(tempCellSubMatr)[colmn]
            newName = paste(tlDesign[sample, "stage"], tlDesign[sample, "replicate"], tlDesign[sample, "ID"],sep="_")
            newColNames[colmn] = newName
        }
        colnames(tempCellSubMatr) = newColNames
        
        # reorder columns by cell type
        colNameCellTypes = NULL
        for(colmn in 1:dim(tempCellSubMatr)[2]){
            cellType = strsplit(colnames(tempCellSubMatr)[colmn], "_")[[1]][3]
            colNameCellTypes[colmn] = cellType
        }
        
        colOrder = NULL
        for(i in 1:length(safeCellTypes)){
            colOrder = c(colOrder, which(colNameCellTypes==safeCellTypes[i]))
        }
        
        tempCellSubMatr <- tempCellSubMatr[,colOrder]
    }
    
    # log the matrix
    if(log){
        tempCellSubMatr <- log10(tempCellSubMatr+1)
    }
    
    if(plus1){
        tempCellSubMatr <- tempCellSubMatr+1
    }
    
    if(hm2){
        library('gplots')
        heatmap.2(tempCellSubMatr, trace = "none", col=colorRampPalette(c("white", "black")), Rowv=F, Colv=F)
        return(NA)
    }
    
    # export this matrix
    write.table(tempCellSubMatr, file = outFile, sep="\t", col.names = NA, quote=F)
    
    # open in cluster 3.0
    # cluster genes (don't necessarily have to cluster cells)
    # uncentered similarity metrics seemed to work for this dataset
    # average linkage
    # then view in java tree view:
    #     open data (must close other data first)
    # setting -> pixel setting 0=white, neg=blue, fill -> fill -,
    # global -> x+y = fill,
    # contrast -> make red/blue visible by messing w this
    # average or complete
    
    rm(tempCellSubMatr, newColNames, colmn, sample, newName, 
       colNameCellTypes, cellType, colOrder)
}
