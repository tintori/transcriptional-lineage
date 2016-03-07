#######
# 2015.12.31
#   All functions for transcriptional lineage paper


feedNewRPKMs <- function(directory, by="rpkm"){
    #this script makes a data matrix of either counts or rpkms from 
    #   every file in a directory of rpkm files.
    
    dat=NULL
    if(by=="rpkm"){dat <- 4}
    if(by=="count"){dat <- 2}
    rpkmStringMatrix <- matrix(NA)
    colmNames = vector()
    
    for(file in list.files(directory, recursive = TRUE)){
        #first extract column name to match the matrix
        splitFile <- as.vector(strsplit(file, split = "/"))
        splitFile <- splitFile[[1]][length(splitFile[[1]])]
        splitFile <- as.vector(strsplit(splitFile, split="_"))
        #splitFile <- as.vector(strsplit(file, split="_"))
        colm <- splitFile[[1]][1]
        
        #read in file
        tempTable <- read.table(paste(directory, file, sep="/"), header=FALSE, sep="\t")
        colnames(tempTable) <- c("gene", "count", "length", "rpkm")
        
        #if this is the first entry, make a new matrix
        if(dim(rpkmStringMatrix)[1]<10){
            #load the file into v2allRPKMs
            rpkmStringMatrix <- as.matrix(tempTable[,c(1,dat)])
        } else {
            rpkmStringMatrix <- merge(rpkmStringMatrix, tempTable[,c(1,dat)], by="gene", all=TRUE)         
        }
        colmNames <- c(colmNames, colm)
        colnames(rpkmStringMatrix) <- c("gene", colmNames)
        
    }
    rpkmMatrix <- as.matrix(rpkmStringMatrix[,2:dim(rpkmStringMatrix)[2]])
    class(rpkmMatrix) <- "numeric"
    colnames(rpkmMatrix) <- colmNames
    row.names(rpkmMatrix) <- rpkmStringMatrix[,1]
    return(rpkmMatrix)
}


typeCells <- function(cells){
    return(rownames(tlDesign[which((tlDesign$ID %in% cells | tlDesign$safeID %in% cells) & tlDesign$quality=="Use"),]))
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
                       by = "PValue") { #when exporting a list, should it be ranked 
    # by "PValue", "logCPM", or "logFC".
    
    
    suppressWarnings(if(group=="default"){
        group <- factor(metaD[colnames(matr), "ID"])
    })
    group <- factor(group)
    samps <- levels(group)
    
    cellsUsed <- colnames(matr)
    threshWhich <- which(apply(mastMatr[,cellsUsed], 1, function(x){length(x[which(x>thresh)])>0}))
    #threshWhich <- which(apply(mastMatr[,cellsUsed], 1, function(x){summary(x)["Max."]})>thresh)
    matrThreshed <- matr[threshWhich,]
    
    postThreshGenes <- which(apply(mastMatr[,cellsUsed], 1, function(x){length(x[which(x>postThresh)])>1}))
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
        plotTabSig(tab, sig, plus, samps, overlay, label, postThresh, opac)
    }else{
        pdf(to, plotWidth, plotHeight)
        plotTabSig(tab, sig, plus, samps, overlay, label, postThresh, opac, visibleLabs=visibleLabs)
        dev.off()
    }
}  


plotTabSig <- function(tab, sig, plus, samps, overlay, label, postThresh, opac, visibleLabs=F, cellHead="", genesHead=""){
    library(scales)
    #Plot all genes in gray
    plot(tab$logCPM, tab$logFC, 
         col=alpha("black", opac), pch=16, 
         xlab="Log ( Average Count )", ylab="Log ( Fold Change )", 
         main=cellHead)
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


singleCellMAplots <- function(geneList = "", #genes to circle in cyan
                              outDir="Desktop/whoopsgraduateschool archive/", to="pdf",
                              geneTypeName="", cellTypeList="all", #which cells are you making cell v embryo plots?
                              plotHeight=8, plotWidth=8, 
                              stage="self", # other options are "1cell", "2cell", etc
                              cellName="none", #use if youre combining cell types (eg Ex1 and Ex2), and comparing those to the rest of the embryo
                              thresh=0, postThresh=25, sig=.1,
                              safeCell=T, # Cx1 and Cx2, or Cx?
                              opac = .4 ){ # darkness of the black dots
    # creates single cell vs. embryo MA plots for all cell types, with a specified list of genes circled
    
    # identify all cell types to make plots for
    if("all" %in% cellTypeList){
        if(safeCell){cellTypeList <- unique(tlDesign$safeID)}
        else {cellTypeList <- unique(tlDesign$ID)}
        cellTypeList <- setdiff(cellTypeList, c("tossed", "P0", NA))
    }
    
    # for each cell type, make a plot
    for(cell in cellTypeList){
        if(stage=="self"){
            stage=unique(tlDesign[rownames(which(tlDesign[,c("ID", "safeID")]==cell, arr.ind = T)), "stage"])
        }
        stageCells = unique(c(rownames(tlDesign[which(tlDesign$stage==stage & tlDesign$quality=="Use"),]),
                              rownames(which(tlDesign[,c("ID", "safeID")] == cell & tlDesign[,"quality"]=="Use", arr.ind = T))))
        if(cellName=="none"){
            cellName <- cell
        }
        group = makeGroups(cell, stageCells, cellName = cellName, restName = paste(stage, "embryo", sep=" "))
        # make filenames for the two types of plots that will be generated
        outNameFewGenes = paste(outDir, cellName, "cells_", geneTypeName, ".pdf", sep="")
        outNameAllGenes = paste(outDir, cellName, "cells_", geneTypeName, "_allGenesLabelled.pdf", sep="")
        # identify genes to circle
        if(is.list(geneList)){
            geneList = geneList[[stage]]
        }
        
        if(to=="pdf"){
            # plot one with all genes labelled
            ERdiffExV1(matr=tlCounts[,stageCells], group=group, sig=sig, overlay=geneList, label=geneList,
                       to=outNameFewGenes, plotWidth=plotWidth, plotHeight=plotHeight, thresh=thresh, postThresh=postThresh, opac=opac)
            # plot one with only circled genes labelled
            ERdiffExV1(matr=tlCounts[,stageCells], group=group, sig=sig, overlay=geneList, label=rownames(v3uniqueCounts),
                       to=outNameAllGenes, plotWidth=plotWidth, plotHeight=plotHeight, thresh=thresh, postThresh=postThresh, opac=opac)
        }
        else if(to=="list"){
            sigGeneList = ERdiffExV1(matr=tlCounts[,stageCells], group=group, sig=sig, overlay=geneList, label=geneList,
                                     to="list", plotWidth=plotWidth, plotHeight=plotHeight, thresh=thresh, postThresh=postThresh)
            return(sigGeneList)
        }
        else if(to=="plottable data"){
            plottableSet = ERdiffExV1(matr=tlCounts[,stageCells], group=group, sig=sig, overlay=geneList, label=geneList,
                                      to="plottable data", plotWidth=plotWidth, plotHeight=plotHeight, thresh=thresh, postThresh=postThresh)
            return(plottableSet)
        }
    }    
}


makeGroups <- function(cellType, cellSet, cellName=cellType, restName="embryo"){
    outFact = rep(NA, length(cellSet))
    cells = which(cellSet %in% rownames(which(tlDesign[,c("ID", "safeID")]==cellType, arr.ind = T)))
    outFact[cells] = cellName
    outFact[setdiff(1:length(cellSet), cells)] = restName
    outFact = factor(outFact, levels = c(cellName, restName))
    return(outFact)
}


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
    vect<-(vect+buff)
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
        thisCol <- repColors[v3Design[cell,"replicate"]]
        repCols<- c(repCols, thisCol)
    }
    return(repCols)
}

exportGeneExpMatrix <- function(outFile="Desktop/whoopsgraduateschool archive/matrix.txt", 
                                geneSet=rownames(tlRPKMs), hm2=F, avgOn=F, log=T, plus1 = F){
    # filter cells
    if(avgOn){
        tempCellSubMatr = tlCellAvg[geneSet,]
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
        for(i in 1:length(cellTypes)){
            colOrder = c(colOrder, which(colNameCellTypes==cellTypes[i]))
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
}

ABxxxPCA <- function(geneList, outPDF, cells="all", width=8, height=8, labelSamps=F, splitReps=F){
    # define cells to use
    suppressWarnings(if(cells=="all"){
        cells <- row.names(tlDesign[which(tlDesign$stage=="16-cell" 
                                          & tlDesign$quality=="Use" 
                                          & !(tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2"))),])
    })
    princ <- prcomp(prop.table(tlRPKMs[geneList,cells]+1, margin=2))
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










cellSpecificHoms <- function(cell, #if compound, put in the form of a list,
                             cellName="none",
                             outDir="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/MAplotsWithHoms/", 
                             forceDate="none",
                             returnSets = F,
                             sig=.1, 
                             elim="equivalent" # for reducing redundant CHOMs
                             ){
    
    #generate a list of the genes that are enriched
    plotData = singleCellMAplots(geneList = c(), outDir= "", to="plottable data", geneTypeName = "null", 
                                 cellTypeList = cell, cellName=cellName, sig=sig)
    if(cellName=="none"){
        cellName <- cell
    }
    #take plot data and make a list of enriched genes
    singleCellSigGenes = list()
    #make a FC down subset
    downSub <- plotData[[1]][which(plotData[[1]]$PValue < plotData[[2]] & plotData[[1]]$logFC < 0 & rownames(plotData[[1]]) %in% plotData[[7]]),]
    singleCellSigGenes[[plotData[[4]][1]]] = downSub
    #    upSub <- plotData[[1]][which(plotData[[1]]$PValue < plotData[[2]] & plotData[[1]]$logFC > 0 & rownames(plotData[[1]]) %in% plotData[[7]]),]
    #    singleCellSigGenes[[plotData[[4]][2]]] = upSub
    
    cellEnriched <- rownames(singleCellSigGenes[[cellName]])
    
    #make a list of KOGs that show up in the gene set, 
    #   with the genes that match that KOG within the list
    CHOMlist = list()
    for(gene in cellEnriched){
        # identify the CHOM
        tempCHOMs <- as.character(ceHoms[which(ceHoms[,1]==gene),2])
        for(tempCHOM in tempCHOMs){
            # if something is there, push it with the gene to the list
            if(length(tempCHOM)>0){
                #add it if the listing already exists
                if(tempCHOM %in% names(CHOMlist)){
                    CHOMlist[[tempCHOM]] = c(CHOMlist[[tempCHOM]], gene)
                } 
                #otherwise, create a listing
                else{
                    CHOMlist[[tempCHOM]] = gene
                }
            }
        }
    }
    
    #print some stats
    library('Biobase')
    listLenKOG <- listLen(CHOMlist)
    if(length(listLenKOG) > 0){
        for(dups in 1:range(listLenKOG)[2]){
            print(paste("There are ", as.character(length(which(listLenKOG==dups))), " KOGs with ", as.character(dups), " genes represented", sep=""))
        }
    }
    
    
    #remove CHOM listings that have only one gene
    for(CHOM in names(CHOMlist)){
        if(length(CHOMlist[[CHOM]])<2){
            CHOMlist <- CHOMlist[-which(names(CHOMlist) == CHOM)]
        }
    }
    
    #remove CHOM listings where less than half the genes of that CHOM are enriched.
    #   if the CHOM still passes, add all the other genes to the group
    for(tempCHOM in unique(ceHoms[,2])){
        if(tempCHOM %in% names(CHOMlist)){
            tempShortLen <- length(CHOMlist[[tempCHOM]])
            tempLongLen <- length(ceHoms[which(ceHoms[,2]==tempCHOM), 1])
            if(tempShortLen > (.5*tempLongLen)){
                CHOMlist[[tempCHOM]] = union(CHOMlist[[tempCHOM]], ceHoms[which(ceHoms[,2]==tempCHOM), 1])
            } else {
                CHOMlist <- CHOMlist[-which(names(CHOMlist) == tempCHOM)]
            }
        }
    }
    
    # clear out the doubles
    tempList <- list()
    tempList[[cellName]] <- CHOMlist
    CHOMlist <- clearRedundCHOMs(tempList, elim = elim, enrichedGenes = cellEnriched)[[cellName]]
    
    # then turn that plot data into a series of plots... (and maybe a list of sets)
    returnList <- list()
    for(tempCHOM in names(CHOMlist)){
        cellGenes <- CHOMlist[[tempCHOM]]
        cellGenes <- intersect(cellGenes, cellEnriched)
        if(length(cellGenes)>=2){
            if(returnSets){
                returnList[[tempCHOM]] = cellGenes
            }
            overlay=CHOMlist[[tempCHOM]]
            label=CHOMlist[[tempCHOM]]
            if(forceDate=="none"){
                forceDate=as.character(format(Sys.Date(), format="%Y%m%d"))
            }
            outTo = paste(outDir, forceDate, "_", plotData[[4]][1], "-", tempCHOM, ".pdf", sep="")
            genesSubTi <- paste(cellGenes, collapse=", ")
            pdf(file=outTo, 5,5, family="Times")
            plotTabSig(plotData[[1]], plotData[[2]], 
                       plotData[[3]], plotData[[4]], 
                       overlay, label, plotData[[7]], visibleLabs=T, cellHead=cellName, genesHead=genesSubTi)
            dev.off()
        }
    }
    return(returnList)
}



batchHomplots <- function(cellList = "all", cellNameVect = "none", 
                          outDir="Desktop/whoopsgraduateschool\ archive/", 
                          forceDate="none", returnSets=F,
                          sig=.1, elim = "overlap" # for reducing redundant CHOMs (opt = "equivalent)
                          ){
    #make list of all cell types to be evaluated
    dir.create(outDir, showWarnings = FALSE)
    if(cellList=="all"){
        cellList = unique(tlDesign$ID)[-which(unique(tlDesign$ID) %in% c("tossed", "P0", NA))]
    }
    if(cellNameVect=="none"){ cellNameVect <- unlist(cellList) }
    #make set of plots for each cell type
    returnables = list()
    for(cellNum in 1:length(cellList)){
        print(paste("Generating plots for ", cellNameVect[cellNum], sep=""))
        cellDir = paste(outDir, "/", cellNameVect[cellNum], "/", sep="")
        dir.create(file.path(cellDir), showWarnings = FALSE)
        returnables[[cellNameVect[cellNum]]] <- cellSpecificHoms(cellList[cellNum], cellName = cellNameVect[cellNum], 
                                                                cellDir, forceDate, returnSets=returnSets, sig=sig, elim=elim)
    }
    if(returnSets){
        return(returnables)
    }
}


clearRedundCHOMs <- function(CHOMlist, elim = "overlap", # other option is "equivalent"
                             enrichedGenes=rownames(tlRPKMs)) {
    #for each cell or group that was analysed
    for(cell in names(CHOMlist)){    # cell = "Ex"
        #for each CHOM that was identified
        for(tempCHOM in names(CHOMlist[[cell]])){    # tempCHOM = "CHOM07551"
            #compare it to every other CHOM for that cell
            for(compCHOM in names(CHOMlist[[cell]][-which(names(CHOMlist[[cell]])==tempCHOM)])){
                # make some categories: total genes and enriched genes
                totTempGenes <- sort(CHOMlist[[cell]][[tempCHOM]])
                enrTempGenes <- intersect(totTempGenes, enrichedGenes)
                totCompGenes <- sort(CHOMlist[[cell]][[compCHOM]])
                enrCompGenes <- intersect(totCompGenes, enrichedGenes)
                #and if it fits inside of, or equal to, that other CHOM, throw it out
                if(elim=="overlap"){
                    if(length(setdiff(enrTempGenes, enrCompGenes))==0){
                        # if the enriched sets are exactly the same, get rid of the one with more total genes
                        if (identical(enrTempGenes, enrCompGenes)){
                            lenTemp <- length(totTempGenes)
                            lenComp <- length(totCompGenes)
                            if(lenTemp > lenComp){
                                CHOMlist[[cell]] <- CHOMlist[[cell]][-which(names(CHOMlist[[cell]])==tempCHOM)]
                            }else{
                                CHOMlist[[cell]] <- CHOMlist[[cell]][-which(names(CHOMlist[[cell]])==compCHOM)]
                            }
                        } else {
                        # otherwise, get rid of the one that fits inside the other
                            CHOMlist[[cell]] <- CHOMlist[[cell]][-which(names(CHOMlist[[cell]])==tempCHOM)]
                        }
                        break
                    }
                } else if(elim=="equivalent"){
                    if(identical(enrTempGenes, enrCompGenes)){
                        lenTemp <- length(totTempGenes)
                        lenComp <- length(totCompGenes)
                        if(lenTemp > lenComp){
                            CHOMlist[[cell]] <- CHOMlist[[cell]][-which(names(CHOMlist[[cell]])==tempCHOM)]
                        }else{
                            CHOMlist[[cell]] <- CHOMlist[[cell]][-which(names(CHOMlist[[cell]])==compCHOM)]
                        }
                        break
                    }
                    
                }
            }
        }
    }
    return(CHOMlist)
}









mapExpression <- function(gene, # can be single gene or c() vector of genes
         mapType = "cell", # other option = lineage
         logTrans=F, # log10 transform RPKMs?
         to="window", # anything else will be interpreted as file path
         legend=T,
         stat = "avg", # can also be "med"
         fillRange="standard", # otherwise, a c(min, max) vector
         safeCell = T,
         width="default", height="default"
){
    # load packages
    library("ggplot2")
    library("ggmap")
    library("scales")
    
    # average gene expression for each cell type
    cellNames = NULL
    if(safeCell){cellNames = unique(cellMapNameRef$safeDesignNames)}
    else{cellNames = unique(cellMapNameRef$designNames)}
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
        origMatr <- mapAllCells
        tempMatr <- mapAllCells
        safeMatr <- mapAllSafeCells
    } else if (mapType=="lineage"){
        origMatr <- mapLineage
        tempMatr <- mapLineage
    } else if (mapType=="horz"){
        origMatr <- mapCellsHorz
        tempMatr <- mapCellsHorz
        safeMatr <- mapSafeHorz
    }
    tempMatr$id <- as.character(tempMatr$id)
    for(rowNum in 1:length(tempMatr$long)){
        # switch cell type name to ID's I have info on
        tempMapName = tempMatr[rowNum,"id"]
        tempDesignName = NULL
        if(safeCell){ tempDesignName = cellMapNameRef[which(cellMapNameRef$mapNames==tempMapName),"safeDesignNames"]}
        else {tempDesignName = cellMapNameRef[which(cellMapNameRef$mapNames==tempMapName),"designNames"]}
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
    #scale_fill_distiller(name = legendTitle, palette = "Blues", limits=fillRange)
    
   
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
    
    # load packages
    library("ggplot2")
    library("ggmap")
    library("scales")

    if(legendTitle=="stat"){
        legendTitle <- dataIn
    }    
    
    dataIn <- tlCellStats[,dataIn]

    # transform data
    if(logTrans){
        dataIn <- log10(dataIn+1)
    }

    if(mapType=="cell"){
        origMatr <- mapAllCells
        tempMatr <- mapAllCells
        safeMatr <- mapAllSafeCells
    } else if (mapType=="lineage"){
        origMatr <- mapLineage
        tempMatr <- mapLineage
    } else if (mapType=="horz"){
        origMatr <- mapCellsHorz
        tempMatr <- mapCellsHorz
        safeMatr <- mapSafeHorz
    }
    
    # rename the map to match the cellTypes I have ID's
    tempMatr$id <- as.character(tempMatr$id)
    for(rowNum in 1:length(tempMatr$long)){
        tempMapName = tempMatr[rowNum,"id"]
        tempDesignName = NULL
        if(safeCell){tempDesignName = cellMapNameRef[which(cellMapNameRef$mapNames==tempMapName),"safeDesignNames"]}
        else {tempDesignName = cellMapNameRef[which(cellMapNameRef$mapNames==tempMapName),"designNames"]}
        tempMatr[rowNum,"id"] = tempDesignName
    }
    
    # add the average expression data to that dataframe
    tempMatr$stat = 0
    for(cellType in cellTypes){
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



findGenesOfPattern <- function(patternList, # 3 items in list: first group of cells, > or <, second group of cells
         thresh=25, fact=1, # fact is the number-fold by which higher expression must exceed lower
         compar = "bulk", # "bulk" will average expression between all cells in a group, 
         #and "singly" will consider each of them separately, in every pairwise comparison.
         wiggleRoom = 0 # how many comparisons are allowed to not match the pattern for a single gene
){
    # average reads, if they are to be considered in bulk
    subMatr <- matrix(nrow=dim(tlCellAvg)[1], ncol=2)
    rownames(subMatr) <- rownames(tlCellAvg)
    if(compar == "bulk"){
        if(length(patternList[[1]])>1){
            subMatr[,1] <- apply(tlCellAvg[,patternList[[1]]], 1, mean)
        } else {subMatr[,1] <- tlCellAvg[,patternList[[1]]]}
        if(length(patternList[[3]])>1){
            subMatr[,2] <- apply(tlCellAvg[,patternList[[3]]], 1, mean)
        } else {subMatr[,2] <- tlCellAvg[,patternList[[3]]]}
        colnames(subMatr) <-c("firstGroup", "secondGroup")
        # rename patternList groups
        patternList[[1]] <- "firstGroup"
        patternList[[3]] <- "secondGroup"
    } else if(compar =="singly"){
        subMatr <- cbind(tlCellAvg[,as.character(patternList[[1]])], tlCellAvg[,as.character(patternList[[3]])])
        colnames(subMatr) <- c(patternList[[1]], patternList[[3]])
    }
    subMatr <- subMatr[which(apply(subMatr, 1, max)>thresh),]+10
    tempGeneScores = rep(0, length(rownames(subMatr)))
    # for each unique combination of cell types from both cells, and for each gene...
    for(first in patternList[[1]]){
        for(second in patternList[[3]]){
            for(geneNum in 1:length(rownames(subMatr))){
                # test the logical operation, add a point to the gene if it doesn't match the settings
                if(patternList[[2]]==">"){
                    if(subMatr[geneNum,first]<=(fact*subMatr[geneNum,second])){
                        tempGeneScores[geneNum] <- tempGeneScores[geneNum]+1
                    }
                } else if (patternList[[2]]=="<"){
                    if((fact*subMatr[geneNum,first])>=subMatr[geneNum,second]){
                        tempGeneScores[geneNum] <- tempGeneScores[geneNum]+1
                    }
                }
            }
        }
    }
    outGenes <- rownames(subMatr)[which(tempGeneScores<(wiggleRoom+1))]
    return(outGenes)
}


makeEnrichTable <- function(constant=1){
    enrichMatr <- matrix(0, ncol=3, nrow=length(all2cellFCgenes), 
                         dimnames = list(all2cellFCgenes, c("Nishimura", "Hashimshony", "Tintori")))
    for(gene in all2cellFCs){
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

findGenesOfPattern <- function(patternList, # 3 items in list: first group of cells, > or <, second group of cells
         thresh=25, fact=1, # fact is the number-fold by which higher expression must exceed lower
         compar = "bulk", # "bulk" will average expression between all cells in a group, 
         #and "singly" will consider each of them separately, in every pairwise comparison.
         wiggleRoom = 0 # how many comparisons are allowed to not match the pattern for a single gene
){
    # average reads, if they are to be considered in bulk
    subMatr <- matrix(nrow=dim(tlCellAvg)[1], ncol=2)
    rownames(subMatr) <- rownames(tlCellAvg)
    if(compar == "bulk"){
        if(length(patternList[[1]])>1){
            subMatr[,1] <- apply(tlCellAvg[,patternList[[1]]], 1, mean)
        } else {subMatr[,1] <- tlCellAvg[,patternList[[1]]]}
        if(length(patternList[[3]])>1){
            subMatr[,2] <- apply(tlCellAvg[,patternList[[3]]], 1, mean)
        } else {subMatr[,2] <- tlCellAvg[,patternList[[3]]]}
        colnames(subMatr) <-c("firstGroup", "secondGroup")
        # rename patternList groups
        patternList[[1]] <- "firstGroup"
        patternList[[3]] <- "secondGroup"
    } else if(compar =="singly"){
        subMatr <- cbind(tlCellAvg[,as.character(patternList[[1]])], tlCellAvg[,as.character(patternList[[3]])])
        colnames(subMatr) <- c(patternList[[1]], patternList[[3]])
    }
    subMatr <- subMatr[which(apply(subMatr, 1, max)>thresh),]+10
    tempGeneScores = rep(0, length(rownames(subMatr)))
    # for each unique combination of cell types from both cells, and for each gene...
    for(first in patternList[[1]]){
        for(second in patternList[[3]]){
            for(geneNum in 1:length(rownames(subMatr))){
                # test the logical operation, add a point to the gene if it doesn't match the settings
                if(patternList[[2]]==">"){
                    if(subMatr[geneNum,first]<=(fact*subMatr[geneNum,second])){
                        tempGeneScores[geneNum] <- tempGeneScores[geneNum]+1
                    }
                } else if (patternList[[2]]=="<"){
                    if((fact*subMatr[geneNum,first])>=subMatr[geneNum,second]){
                        tempGeneScores[geneNum] <- tempGeneScores[geneNum]+1
                    }
                }
            }
        }
    }
    outGenes <- rownames(subMatr)[which(tempGeneScores<(wiggleRoom+1))]
    return(outGenes)
}

cellSpecifGenes <- function(cellType, oper = ">", fact = 5, expSettings=F, compTo = "all"){
    # make left handed side of operation
    relationsAvail <- intersect(names(fakeLinRelations[[cellType]]), c("daughters", "parent", "granddaughters", "grandparent"))
    firstSet = c(cellType)
    for(relation in relationsAvail){
        firstSet <- c(firstSet, fakeLinRelations[[cellType]][[relation]])
    }
    secondSet=NULL
    cellTypes = setdiff(unique(tlDesign$ID), c("tossed", NA))
    if(compTo=="all"){
        secondSet = cellTypes[which(!cellTypes%in%firstSet)]
    }
    if(compTo=="bro"){
        stage = unique(tlDesign[which(tlDesign$ID==cellType),"stage"])
        stageCells = unique(tlDesign[which(tlDesign$stage == stage & tlDesign$quality=="Use"), "ID"])
        secondSet = cellTypes[which(!cellTypes%in%firstSet & cellTypes %in% stageCells)]
    }
    
    # if testing settings, export them now
    if(expSettings){
        outList <- list(cellType, oper, secondSet, fact)
        return(outList)
    }

    #perform operation
    tempOut= NULL
    tempOut <- findGenesOfPattern(list(cellType, oper, secondSet), fact=fact, compar = "singly")
    return(tempOut)
}


ERdiffVsWhole <- function(cellType, reps="all", sig=.1, to="window", 
                          label="", overlay="", postThresh=0){
    #set stage
    stage = as.character(tlDesign[which(tlDesign$ID==cellType), "stage"][1])
    otherCellTypes = setdiff(unique(tlDesign[which(tlDesign$stage==stage),"ID"]), 
                             c("tossed", cellType, NA))
    #set reps
    suppressWarnings(if(reps=="all"){
        reps=as.character(unique(tlDesign[which(v3Design$ID==cellType), "replicate"]))
    })
    #make temp matrix
    tempCellCols <- NULL
    tempCellCols <- c(tempCellCols, row.names(tlDesign[which(tlDesign$ID==cellType & 
                                                                 tlDesign$replicate %in% reps &
                                                                 tlDesign$quality=="Use"),]))
    tempBroCols <- NULL
    tempBroCols <- c(tempBroCols, row.names(tlDesign[which(tlDesign$ID %in% otherCellTypes & 
                                                                tlDesign$replicate %in% reps &
                                                                tlDesign$quality=="Use"),]))
    
    tempMatr <- cbind(tlRPKMs[,tempCellCols], tlRPKMs[,tempBroCols])
    tempGroups <- factor(c(rep(cellType, length(tempCellCols)), 
                           rep("embryo", length(tempBroCols))), 
                         levels = c(cellType, "embryo"))
    ERdiffExV1(matr=tempMatr, group=tempGroups, sig=sig, to=to, 
               overlay=overlay, label=label, postThresh=postThresh)
}

whichNeighborsHaveGene <- function(centerCell, gene, thresh="default", stat = "avg", safeCells=F){
    # takes a cell and a gene, and returns the cell types in the vicinity that also have that gene
    avgMatr = NULL
    cellTypes = NULL
    if(safeCells){
        cellTypes <- unique(tlDesign$safeID)
        if(stat=="avg"){avgMatr <- tlSafeCellAvg}
        if(stat=="med"){avgMatr <- tlSafeCellMed}
    } else {
        cellTypes <- unique(tlDesign$ID)
        if(stat=="avg"){avgMatr <- tlCellAvg}
        if(stat=="med"){avgMatr <- tlCellMed}
    }
    cellsWithGene <- NULL
    relatCellTypes <- NULL
    for(relat in c("sister", "parent", "daughters", "cousins")){
        relatCellTypes <- c(relatCellTypes, fakeLinRelations[[centerCell]][[relat]])
    }
    for(relatCellType in intersect(relatCellTypes, cellTypes)){
        # if the relative has expression over the raw threshold
        if(thresh=="default"){
            if(avgMatr[gene, relatCellType] > 50){
                if(length(removeGenesUnderAvg(gene, relatCellType, removed="under"))==1){
                    cellsWithGene <- c(cellsWithGene, relatCellType)
                }
            }
        } else {
            if(avgMatr[gene, relatCellType] > thresh){
                cellsWithGene <- c(cellsWithGene, relatCellType)
            }
        }
    }
    cellsWithGene <- c(cellsWithGene, centerCell)
    cellsWithGene <- unique(cellsWithGene)
    return(cellsWithGene)
}

patternToPattName <- function(gene, origList = spatDynList){
    pattName <- origList[[gene]][1]
    for(i in 2:length(origList[[gene]])){
        pattName <- paste(pattName, origList[[gene]][i], sep="_")
    }
    return(pattName)
}

explSpatDyn <- function(rankNum, embryoNum="all", stat="med"){
    tempPatt <- names(spatDynPattScores[rankNum])
    tempGenes <- spatDynPattList[[tempPatt]]
    print(spatDynPattScores[rankNum])
    if(embryoNum=="all"){
        return(mapExpression(tempGenes, safeCell = T, stat=stat))
    } else {
        print(tempGenes[embryoNum])
        return(mapExpression(tempGenes[embryoNum], safeCell = T, stat=stat))
    }
}




# # # # FAKE LINEAGE LIST # # # # 


# make a database of lineal relationships for each cell:

fakeLinRelations <- list()
cellTypes <- factor(c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", 
                      "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                      "ABalx", "ABarx", "ABplx", "ABprx",
                      "MSx", "MSa", "MSp", "MSx1", "MSx2", "Ex", "Ea", "Ep", "Ex1", "Ex2",
                      "Cx", "Ca", "Cp", "Cx1", "Cx2", "D", "P4"), 
                    levels=c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", 
                             "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                             "ABalx", "ABarx", "ABplx", "ABprx",
                             "MSx", "MSa", "MSp", "MSx1", "MSx2", "Ex", "Ea", "Ep", "Ex1", "Ex2",
                             "Cx", "Ca", "Cp", "Cx1", "Cx2", "D", "P4"))
for(cell in cellTypes){
    fakeLinRelations[[cell]] <- list()
}

# DAUGHTERS
fakeLinRelations[["P0"]][["daughters"]] <- c("AB", "P1")

fakeLinRelations[["P1"]][["daughters"]] <- c("EMS", "P2")
fakeLinRelations[["AB"]][["daughters"]] <- c("ABa", "ABp")

fakeLinRelations[["ABa"]][["daughters"]] <- c("ABal", "ABar")
fakeLinRelations[["ABp"]][["daughters"]] <- c("ABpl", "ABpr")
fakeLinRelations[["EMS"]][["daughters"]] <- c("E", "MS")
fakeLinRelations[["P2"]][["daughters"]] <- c("C", "P3")

fakeLinRelations[["ABal"]][["daughters"]] <- c("ABalx")
fakeLinRelations[["ABar"]][["daughters"]] <- c("ABarx")
fakeLinRelations[["ABpl"]][["daughters"]] <- c("ABplx")
fakeLinRelations[["ABpr"]][["daughters"]] <- c("ABprx")
fakeLinRelations[["MS"]][["daughters"]] <- c("MSx1", "MSx2", "MSx")
fakeLinRelations[["E"]][["daughters"]] <- c("Ep", "Ea", "Ex")
fakeLinRelations[["C"]][["daughters"]] <- c("Cx1", "Cx2", "Cx")
fakeLinRelations[["P3"]][["daughters"]] <- c("D", "P4")

# PARENTS
fakeLinRelations[["P1"]][["parent"]] <- c("P0")
fakeLinRelations[["AB"]][["parent"]] <- c("P0")

fakeLinRelations[["ABa"]][["parent"]] <- c("AB")
fakeLinRelations[["ABp"]][["parent"]] <- c("AB")
fakeLinRelations[["EMS"]][["parent"]] <- c("P1")
fakeLinRelations[["P2"]][["parent"]] <- c("P1")

fakeLinRelations[["ABal"]][["parent"]] <- c("ABa")
fakeLinRelations[["ABar"]][["parent"]] <- c("ABa")
fakeLinRelations[["ABpl"]][["parent"]] <- c("ABp")
fakeLinRelations[["ABpr"]][["parent"]] <- c("ABp")
fakeLinRelations[["MS"]][["parent"]] <- c("EMS")
fakeLinRelations[["E"]][["parent"]] <- c("EMS")
fakeLinRelations[["C"]][["parent"]] <- c("P2")
fakeLinRelations[["P3"]][["parent"]] <- c("P2")

fakeLinRelations[["ABalx"]][["parent"]] <- c("ABal")
fakeLinRelations[["ABarx"]][["parent"]] <- c("ABar")
fakeLinRelations[["ABplx"]][["parent"]] <- c("ABpl")
fakeLinRelations[["ABprx"]][["parent"]] <- c("ABpr")
fakeLinRelations[["MSx1"]][["parent"]] <- c("MS")
fakeLinRelations[["MSx2"]][["parent"]] <- c("MS")
fakeLinRelations[["MSx"]][["parent"]] <- c("MS")
fakeLinRelations[["Ep"]][["parent"]] <- c("E")
fakeLinRelations[["Ea"]][["parent"]] <- c("E")
fakeLinRelations[["Ex"]][["parent"]] <- c("E")
fakeLinRelations[["Cx1"]][["parent"]] <- c("C")
fakeLinRelations[["Cx2"]][["parent"]] <- c("C")
fakeLinRelations[["Cx"]][["parent"]] <- c("C")
fakeLinRelations[["D"]][["parent"]] <- c("P3")
fakeLinRelations[["P4"]][["parent"]] <- c("P3")

# GRAND DAUGHTERS
fakeLinRelations[["P0"]][["granddaughters"]] <- c("AB", "P1", "ABa", "ABp", "EMS", "P2")

fakeLinRelations[["P1"]][["granddaughters"]] <- c("E", "MS", "C", "P3")
fakeLinRelations[["AB"]][["granddaughters"]] <- c("ABal", "ABar", "ABpl", "ABpr")

fakeLinRelations[["ABa"]][["granddaughters"]] <- c("ABalx", "ABarx")
fakeLinRelations[["ABp"]][["granddaughters"]] <- c("ABplx", "ABprx")
fakeLinRelations[["EMS"]][["granddaughters"]] <- c("Ep", "Ea", "Ex", "MSx1", "MSx2", "MSx")
fakeLinRelations[["P2"]][["granddaughters"]] <- c("Cx1", "Cx2", "Cx", "D", "P4")

# GRAND PARENTS
fakeLinRelations[["ABa"]][["grandparent"]] <- c("P0")
fakeLinRelations[["ABp"]][["grandparent"]] <- c("P0")
fakeLinRelations[["EMS"]][["grandparent"]] <- c("P0")
fakeLinRelations[["P2"]][["grandparent"]] <- c("P0")

fakeLinRelations[["ABal"]][["grandparent"]] <- c("AB")
fakeLinRelations[["ABar"]][["grandparent"]] <- c("AB")
fakeLinRelations[["ABpl"]][["grandparent"]] <- c("AB")
fakeLinRelations[["ABpr"]][["grandparent"]] <- c("AB")
fakeLinRelations[["MS"]][["grandparent"]] <- c("P1")
fakeLinRelations[["E"]][["grandparent"]] <- c("P1")
fakeLinRelations[["C"]][["grandparent"]] <- c("P1")
fakeLinRelations[["P3"]][["grandparent"]] <- c("P1")

fakeLinRelations[["ABalx"]][["grandparent"]] <- c("ABa")
fakeLinRelations[["ABarx"]][["grandparent"]] <- c("ABa")
fakeLinRelations[["ABplx"]][["grandparent"]] <- c("ABp")
fakeLinRelations[["ABprx"]][["grandparent"]] <- c("ABp")
fakeLinRelations[["MSx1"]][["grandparent"]] <- c("EMS")
fakeLinRelations[["MSx2"]][["grandparent"]] <- c("EMS")
fakeLinRelations[["MSx"]][["grandparent"]] <- c("EMS")
fakeLinRelations[["Ep"]][["grandparent"]] <- c("EMS")
fakeLinRelations[["Ea"]][["grandparent"]] <- c("EMS")
fakeLinRelations[["Ex"]][["grandparent"]] <- c("EMS")
fakeLinRelations[["Cx1"]][["grandparent"]] <- c("P2")
fakeLinRelations[["Cx2"]][["grandparent"]] <- c("P2")
fakeLinRelations[["Cx"]][["grandparent"]] <- c("P2")
fakeLinRelations[["D"]][["grandparent"]] <- c("P2")
fakeLinRelations[["P4"]][["grandparent"]] <- c("P2")

# GREAT GRAND DAUGHTERS
fakeLinRelations[["P0"]][["greatgranddaughters"]] <- c("ABal", "ABar", "ABpl", "ABpr", "E", "MS", "C", "P3")

fakeLinRelations[["P1"]][["greatgranddaughters"]] <- c("Ep", "Ea", "Ex", "MSx1", "MSx2", "MSx", "Cx1", "Cx2", "Cx", "D", "P4")
fakeLinRelations[["AB"]][["greatgranddaughters"]] <- c("ABalx", "ABarx", "ABplx", "ABprx")

# GREAT GRAND PARENTS
fakeLinRelations[["ABal"]][["greatgrandparent"]] <- c("P0")
fakeLinRelations[["ABar"]][["greatgrandparent"]] <- c("P0")
fakeLinRelations[["ABpl"]][["greatgrandparent"]] <- c("P0")
fakeLinRelations[["ABpr"]][["greatgrandparent"]] <- c("P0")
fakeLinRelations[["MS"]][["greatgrandparent"]] <- c("P0")
fakeLinRelations[["E"]][["greatgrandparent"]] <- c("P0")
fakeLinRelations[["C"]][["greatgrandparent"]] <- c("P0")
fakeLinRelations[["P3"]][["greatgrandparent"]] <- c("P0")

fakeLinRelations[["ABalx"]][["greatgrandparent"]] <- c("AB")
fakeLinRelations[["ABarx"]][["greatgrandparent"]] <- c("AB")
fakeLinRelations[["ABplx"]][["greatgrandparent"]] <- c("AB")
fakeLinRelations[["ABprx"]][["greatgrandparent"]] <- c("AB")
fakeLinRelations[["MSx1"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["MSx2"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["MSx"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["Ep"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["Ea"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["Ex"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["Cx1"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["Cx2"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["Cx"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["D"]][["greatgrandparent"]] <- c("P1")
fakeLinRelations[["P4"]][["greatgrandparent"]] <- c("P1")

# GREAT GREAT GRAND DAUGHTERS
fakeLinRelations[["P0"]][["greatgreatgranddaughters"]] <- c("ABalx", "ABarx", "ABplx", "ABprx", "Ep", "Ea", "Ex", "MSx1", "MSx2", "MSx", "Cx1", "Cx2", "Cx", "D", "P4")

# GREAT GREAT GRAND PARENTS
fakeLinRelations[["ABalx"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["ABarx"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["ABplx"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["ABprx"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["MSx1"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["MSx2"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["MSx"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["Ep"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["Ea"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["Ex"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["Cx1"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["Cx2"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["Cx"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["D"]][["greatgreatgrandparent"]] <- c("P0")
fakeLinRelations[["P4"]][["greatgreatgrandparent"]] <- c("P0")


# SISTERS
fakeLinRelations[["P1"]][["sister"]] <- c("AB")
fakeLinRelations[["AB"]][["sister"]] <- c("P1")

fakeLinRelations[["ABa"]][["sister"]] <- c("ABp")
fakeLinRelations[["ABp"]][["sister"]] <- c("ABa")
fakeLinRelations[["EMS"]][["sister"]] <- c("P2")
fakeLinRelations[["P2"]][["sister"]] <- c("EMS")

fakeLinRelations[["ABal"]][["sister"]] <- c("ABar")
fakeLinRelations[["ABar"]][["sister"]] <- c("ABal")
fakeLinRelations[["ABpl"]][["sister"]] <- c("ABpr")
fakeLinRelations[["ABpr"]][["sister"]] <- c("ABpl")
fakeLinRelations[["MS"]][["sister"]] <- c("E")
fakeLinRelations[["E"]][["sister"]] <- c("MS")
fakeLinRelations[["C"]][["sister"]] <- c("P3")
fakeLinRelations[["P3"]][["sister"]] <- c("C")

fakeLinRelations[["ABalx"]][["sister"]] <- c("ABarx")
fakeLinRelations[["ABarx"]][["sister"]] <- c("ABalx")
fakeLinRelations[["ABplx"]][["sister"]] <- c("ABprx")
fakeLinRelations[["ABprx"]][["sister"]] <- c("ABplx")
fakeLinRelations[["MSx1"]][["sister"]] <- c("MSx2")
fakeLinRelations[["MSx2"]][["sister"]] <- c("MSx1")
fakeLinRelations[["MSx"]][["sister"]] <- c("MSx")
fakeLinRelations[["Ep"]][["sister"]] <- c("Ea")
fakeLinRelations[["Ea"]][["sister"]] <- c("Ep")
fakeLinRelations[["Ex"]][["sister"]] <- c("Ex")
fakeLinRelations[["Cx1"]][["sister"]] <- c("Cx2")
fakeLinRelations[["Cx2"]][["sister"]] <- c("Cx1")
fakeLinRelations[["Cx"]][["sister"]] <- c("Cx")
fakeLinRelations[["D"]][["sister"]] <- c("P4")
fakeLinRelations[["P4"]][["sister"]] <- c("D")

# COUSINS #

fakeLinRelations[["ABa"]][["cousins"]] <- c("EMS", "P2")
fakeLinRelations[["ABp"]][["cousins"]] <- c("EMS", "P2")
fakeLinRelations[["EMS"]][["cousins"]] <- c("ABa", "ABp")
fakeLinRelations[["P2"]][["cousins"]] <- c("ABa", "ABp")

fakeLinRelations[["ABal"]][["cousins"]] <- c("ABpl", "ABpr")
fakeLinRelations[["ABar"]][["cousins"]] <- c("ABpl", "ABpr")
fakeLinRelations[["ABpl"]][["cousins"]] <- c("ABal", "ABar")
fakeLinRelations[["ABpr"]][["cousins"]] <- c("ABal", "ABar")
fakeLinRelations[["MS"]][["cousins"]] <- c("C", "P3")
fakeLinRelations[["E"]][["cousins"]] <- c("C", "P3")
fakeLinRelations[["C"]][["cousins"]] <- c("E", "MS")
fakeLinRelations[["P3"]][["cousins"]] <- c("E", "MS")

fakeLinRelations[["ABalx"]][["cousins"]] <- c("ABplx", "ABprx")
fakeLinRelations[["ABarx"]][["cousins"]] <- c("ABplx", "ABprx")
fakeLinRelations[["ABplx"]][["cousins"]] <- c("ABalx", "ABarx")
fakeLinRelations[["ABprx"]][["cousins"]] <- c("ABalx", "ABarx")
fakeLinRelations[["MSx1"]][["cousins"]] <- c("Ep", "Ea")
fakeLinRelations[["MSx2"]][["cousins"]] <- c("Ep", "Ea")
fakeLinRelations[["MSx"]][["cousins"]] <- c("Ex")
fakeLinRelations[["Ep"]][["cousins"]] <- c("MSx1", "MSx2")
fakeLinRelations[["Ea"]][["cousins"]] <- c("MSx1", "MSx2")
fakeLinRelations[["Ex"]][["cousins"]] <- c("MSx")
fakeLinRelations[["Cx1"]][["cousins"]] <- c("D", "P4")
fakeLinRelations[["Cx2"]][["cousins"]] <- c("D", "P4")
fakeLinRelations[["Cx"]][["cousins"]] <- c("D", "P4")
fakeLinRelations[["D"]][["cousins"]] <- c("Cx1", "Cx2", "Cx")
fakeLinRelations[["P4"]][["cousins"]] <- c("Cx1", "Cx2", "Cx")


