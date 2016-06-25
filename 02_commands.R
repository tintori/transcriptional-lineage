#~#~#~#~#
#   
#   All commands used to generate analysis figures in Tintori et al. 2016
#   Code by Sophia Tintori sophia.tintori@gmail.com
#


# Load scripts
source("01_functions.R")


# Load libraries
require('reshape2')
install.packages("ggplot2")
require('ggplot2')
require('gplots')
require('scales')
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
require('ggmap')


# Load RNA-seq data
tlRPKMs <- as.matrix(read.table("03_dataIn/0_seqData/tlRPKMs.txt", stringsAsFactors = F))
tlCounts <- as.matrix(read.table("03_dataIn/0_seqData/tlCounts.txt", stringsAsFactors = F))
tlDesign <- read.table("03_dataIn/0_seqData/tlDesign.txt", stringsAsFactors = F)
tlDesign$stage <- factor(tlDesign$stage, levels = c("1-cell", "2-cell", "4-cell", "8-cell", "16-cell"))
tlDesign$safeID <- factor(tlDesign$safeID, levels = c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3", "ABalx", "ABarx", "ABplx", "ABprx", "MSx", "Ex", "Cx", "D", "P4"))
repColors=c("dark blue", "red", "orange", "purple", "brown", "darkcyan", "darksalmon", "darkgreen", "violetred")

# Generate data tables by cell type rather than replicate
synthData <- function(){
    safeCellTypes <- as.character(unique(tlDesign$safeID)[!unique(tlDesign$safeID) %in% c("tossed", "ABxxx", NA)])
    tlSafeCellAvg <- matrix(0, nrow=nrow(tlRPKMs), ncol=length(safeCellTypes), 
                            dimnames = list(rownames(tlRPKMs), safeCellTypes))
    for(cT in safeCellTypes){
        tlSafeCellAvg[,cT] <- apply(tlRPKMs[,typeCells(cT)], 1, mean) 
    }
    
    tlSafeCellMed <- matrix(0, nrow=nrow(tlRPKMs), ncol=length(safeCellTypes), 
                            dimnames = list(rownames(tlRPKMs), safeCellTypes))
    
    for(cT in safeCellTypes){
        tlSafeCellMed[,cT] <- apply(tlRPKMs[,typeCells(cT)], 1, median) 
    }
    rm(cT)
}
synthData()


##
# Figure 1A
#   Synthesis of Sulston et al 1983 data

if(!file.exists("04_analysesOut/Figure1/")){
    dir.create("04_analysesOut/Figure1/")
}

makeSulstonPlot <- function(){
    
    founders <- read.table(file = "03_dataIn/1A_sulston/founder_cell_fates.txt", 
                           row.names = 1, header = T, sep="\t")
    # (remove total line)
    founders <- founders[1:dim(founders)[1]-1,]
    # subscript p4
    colnames(founders)[which(colnames(founders)=="P4")] <- expression(P[4])
    founders <- rbind(founders, "neurons" = apply(founders, 2, function(x){sum(x["misc neurons"], x["ventral"], x["ring ganglia"])}))
    founders["head meso",] <- apply(founders, 2, function(x){sum(x["head meso"], x["arcade"])})
    founders <- founders[which(!rownames(founders)%in%c("ring ganglia", "ventral", "misc neurons", "arcade", "valve")),]
    
    for (i in 1:16) {founders[,i] = founders[,i]/sum(founders[,i])}
    
    founders.m <- melt(founders)
    founders.m <- cbind(founders.m, rep(row.names(founders), 16))
    founders.m <- data.frame(founders.m[founders.m$value != 0,], row.names = NULL)
    colnames(founders.m) = c("Founder_Cell", "Proportion", "Cell_Fate")
    founders.m = transform (founders.m,
                            Cell_Fate = factor(Cell_Fate,
                                               levels = c("apoptosis", "neurons", "germ line", "intestine", "hypodermis", "pharynx", "muscle", "distal tip", "excretory", "ceolomocyte", "head meso"),
                                               labels = c("Apoptosis", "Neurons", "Germ Line", "Intestine", "Epidermis", "Pharynx", "Muscle", "Gonad", "Excretory", "Ceolomocyte", "Head Mesoderm")
                            )
    )
    
    qual_11 <- c("#a6cee3", "#1f78b4", "#b2bf8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
    
    ggplot(founders.m, aes (x=Founder_Cell, y=Proportion, color= Cell_Fate, fill=Cell_Fate))+
        geom_bar(stat="identity", colour = "#FFFFFF00")+ scale_fill_manual(values = qual_11)+ 
        scale_y_continuous(breaks=NULL)+
        guides(fill=guide_legend(title="Descendant\nFates"))+
        theme_minimal()+
            scale_x_discrete(labels=c(colnames(founders)[1:15],expression(P[4])))+
        theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              axis.text.y = element_blank(), text=element_text(size=20),
              legend.text = element_text(size=10),
              legend.title = element_text(size=12))+
        labs(#title="Descendant fates of the 16 cell stage", 
            x="Ancestor Cell", y="Proportion\nof descendants")
    
    ggsave(file = "04_analysesOut/Figure1/1A_sulstonData.pdf", 
           width=10, height = 4)
    
    # clear out workspace
    rm(founders, founders.m, qual_11, i)
}

makeSulstonPlot()


##
# Figure 1C-E
#   QC of single-cell RNA-seq datasets

makeQCplots <- function(){

    # Figure 1C
    # Picograms mRNA in each whole embryo
    
    # Load in RNA-seq data for ERCC reads only
    erccCounts <- read.table("03_dataIn/1C_ERCC/erccCounts.txt", stringsAsFactors = F)
    
    # Load in concentration data (from ERCC documentation)
    ERCCconc <- read.table("03_dataIn/1C_ERCC/ERCC_Controls_Analysis.txt", 
                               header=T, sep="\t")[,c(1,3)]
    colnames(ERCCconc) <- c("ERCC.ID", "attomoles.per.ul")
    temp <- read.table("03_dataIn/1C_ERCC/ERCC92.gff2")
    colnames(temp)[c(1,5)] <- c("ERCC.ID", "length")
    ERCCconc <- merge(ERCCconc, temp[,c(1,5)], by="ERCC.ID")
    
    # attomoles are 10^-18 of a mole
    # a mole is 6.02214179 * 10^23 transcripts
    ERCCconc <- cbind(ERCCconc, ERCCconc[,2]*602214.179)
    ERCCconc <- cbind(ERCCconc, ERCCconc[,4]/2000000)
    colnames(ERCCconc)[c(4,5)] <- c("molecules.per.ul", "molecules.2millDilution")
    rownames(ERCCconc) <- ERCCconc[,1]
    
    # Calculate the pg expected from each cell
    erccBasesExp <- 0
    for(gene in rownames(erccCounts)){
        erccBasesExp <- erccBasesExp + (ERCCconc[gene, "molecules.2millDilution"]*ERCCconc[gene, "length"])
    }
    
    tlDesign$pgExpected <- 0
    baseToPgConstant <- (499.5*10^12)/(6.02214179*10^23)
    for(cell in rownames(tlDesign)){
        erccBasesObs <- sum(erccCounts[,cell])*50
        ExpObsRatio <- erccBasesExp/erccBasesObs
        ce10BasesObs <- sum(tlCounts[,cell])*50
        ce10BasesExp <- ce10BasesObs * ExpObsRatio
        tlDesign[cell, "pgExpected"] <- ce10BasesExp*baseToPgConstant
    }
    # The first four cell replicate was different though, 
    #   because it only had a 1:500,000 dilution instead of 2 million
    for(cell in rownames(tlDesign[which(tlDesign$stage=="4cell"&tlDesign$replicate==1),])){
        erccBasesObs <- sum(erccCounts[,cell])*50
        ExpObsRatio <- erccBasesExp/erccBasesObs
        ce10BasesObs <- sum(tlCounts[,cell])*50
        ce10BasesExp <- ce10BasesObs * ExpObsRatio
        tlDesign[cell, "pgExpected"] <- ce10BasesExp*baseToPgConstant*4
    }
    
    # Calculate the pg expected from each embryo
    
    tlDesign$embryoPgExpected <- 0
    for(embryo in unique(paste(tlDesign$stage, tlDesign$replicate, sep="."))){
        embSplit = strsplit(embryo, "\\.")[[1]]
        stage=embSplit[1]
        repl = embSplit[2]
        outNumb = 0
        celset <- intersect(colnames(tlCounts[,which(tlDesign$stage==stage)]), 
                            colnames(tlCounts[,which(tlDesign$replicate==repl)]))
        for(cell in celset){
            outNumb = outNumb + tlDesign[cell, "pgExpected"]
        }
        tlDesign[which(tlDesign$stage==stage & tlDesign$replicate==repl), "embryoPgExpected"] = outNumb
    }
    
    #   Plot pg material detected in whole embryos
    
    pdf(file = "04_analysesOut/Figure1/1C_wholeEmbryoPg.pdf", 
        width=8.5, height=4.5)
    par(mar=c(5.1, 6.1, 4.1, 6.1))
    embMean = mean(unique(log10(tlDesign[, "embryoPgExpected"]+1)))
    embSD = sd(unique(log10(tlDesign[, "embryoPgExpected"]+1)))
    plot(-1,-1, xlim=c(1,31), ylim=c(-.05,1.8),
         xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
    rect(-3,embMean-embSD,34,embMean+embSD, col="light gray", border = NA)
    abline(h=embMean, col=alpha("black", .2))
    points(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop=TRUE)), log10(tlDesign$embryoPgExpected+1), 
           pch=18, col=repColors[tlDesign$replicate], cex=2
    )
    abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
    axis(1, at=c(3.5,9,14,19.5,27), tick=FALSE, cex.axis=1.5,
         labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
    axis(2, at=c(0,1.5), cex.axis=1.5)
    title(ylab = expression(paste("log"[10], "(pg RNA)", sep="")), cex.lab = 1.5, line = 2.5)
    title(xlab = "Embryos (aggregated from single-cell data)", cex.lab=1.5)
    dev.off()
    
    #   Change quality scores for outlying embryos, based on this plot
    
    tlDesign[which(tlDesign$stage=="1-cell" & tlDesign$replicate==3),"quality"] <- "Toss"
    tlDesign[which(tlDesign$stage=="8-cell" & tlDesign$replicate==5),"quality"] <- "Toss"
    tlDesign[which(tlDesign$stage=="16-cell" & tlDesign$replicate%in%c(4,5,9)),"quality"] <- "Toss"
    
    #   Figure 1D
    # Whole embryo vs number of genes
    
    thresh = 25
    
    # calculate genesDetected
    tlDesign$genesDetected <- 0
    for(cell in rownames(tlDesign)){
        binLog <- (tlRPKMs[,cell]>thresh)
        tlDesign[cell,"genesDetected"] = length(which(binLog))
    }
    
    # calculate whole embryo genes detected
    tlDesign$embryoGenesDetected <- 0
    for(embryo in unique(paste(tlDesign$stage, tlDesign$replicate, sep="."))){
        embSplit = strsplit(embryo, "\\.")[[1]]
        stage=embSplit[1]
        repl = embSplit[2]
        tempMatr = as.matrix(tlRPKMs[,which(tlDesign$stage==stage & tlDesign$replicate==repl)])
        tempLog <- rep(FALSE, nrow(tempMatr))
        for(tempRow in 1:nrow(tempMatr)){
            for(tempCol in 1:ncol(tempMatr)){
                if(tempMatr[tempRow,tempCol] > thresh){
                    tempLog[tempRow] = TRUE
                }
            }
        }
        outNum = length(which(tempLog))
        tlDesign[which(tlDesign$stage==stage & tlDesign$replicate==repl), "embryoGenesDetected"] = outNum
    }
    
    
    pdf(file = "04_analysesOut/Figure1/1D_wholeEmbryoNumGenes_thresh.pdf", 
        width=8.5, height=4.5)
    par(mar=c(5.1, 6.1, 4.1, 6.1))
    plot(-1,-1, xlim=c(1,31), ylim=c(3.6,3.95),
         xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
    tempColors <- repColors[tlDesign$replicate]
    tempColors[which(tlDesign$quality=="Toss")] <- "gray"
    points(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop=TRUE)), 
           log10(tlDesign$embryoGenesDetected+1), 
           pch=18, col=tempColors, cex=2
    )
    axis(1, at=c(3.5,9,14,19.5,27),tick=FALSE, cex.axis=1.5,
         labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
    axis(2, cex.axis=1.5, at=c(3.6,3.7, 3.8, 3.9)
    )
    abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
    title(ylab = expression(paste("  log"[10], "(genes detected)", sep="")), cex.lab = 1.5, line = 2.5)
    title(xlab = "Embryos (aggregated from single-cell data)", cex.lab=1.5)
    dev.off()
    
    #   Figure 1E
    # Single cells vs RPKMs
    
    pdf(file = "04_analysesOut/Figure1/1E_singleCellNumGenes_thresh.pdf", 
        width=8.5, height=4.5)
    par(mar=c(5.1, 6.1, 4.1, 6.1))
    library('scales')
    # set up color scheme for cells that got dropped
    cellCol <- repColors[tlDesign$replicate]
    cellCol[c(which(tlDesign$stage=="1-cell" & tlDesign$replicate==3),
              which(tlDesign$stage=="8-cell" & tlDesign$replicate==5),
              which(tlDesign$stage=="16-cell" & tlDesign$replicate%in%c(4,5,9)))] <- "light gray"
    plot(jitter(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop = TRUE)), factor = .7), (log10(tlDesign$genesDetected+1)), 
         cex=((25*tlDesign$pgExpected)^.1),
         pch=16, col=alpha(cellCol, .6), xlim=c(1,31), ylim=c(2.6,3.8),
         xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
    #title(main = "Number of Genes Detected Above 25 RPKM in each cell")
    abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
    axis(1, at=c(3.5,9,14,19.5,27), tick=FALSE, cex.axis=1.5,
         labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
    axis(2, cex.axis=1.5)
    title(ylab = expression(paste("  log"[10], "(genes detected)", sep="")), cex.lab = 1.5, line = 2.5)
    title(xlab = "Individual Cells", cex.lab=1.5)
    legend(legend = c(".2", "2", "20"), title = "pg mRNA", #horiz = T,
           #x = par('usr')[2]+1, y = par('usr')[4]-.15,
           x = .5, y=3.4, box.col = "dark gray", bg = "white",
           #col = "dark gray", border="dark gray", text.col="dark gray",
           pch=16, pt.cex=c(5, 50, 500)^(.1), cex=1.2, 
           y.intersp=1.2, x.intersp=.8,  xpd=TRUE)
    dev.off()
    
    # Remove all these new values
    rm(baseToPgConstant, binLog, ce10BasesExp, ce10BasesObs, cell, cellCol, 
       celset, embMean, embryo, embSD, embSplit, erccBasesExp, erccBasesObs, 
       ExpObsRatio, gene, outNum, outNumb, repl, stage, tempCol, tempColors, 
       tempLog, tempRow, thresh, ERCCconc, erccCounts, temp, tempMatr)
}

makeQCplots()

##
# Figure 2
#   PCA plots to match transcriptomes up with others of similar cell type

if(!file.exists("04_analysesOut/Figure2/")){
    dir.create("04_analysesOut/Figure2/")
}

makePCAs <- function(){
        
    # TWO CELL
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="2-cell" & tlDesign$quality=="Use")], falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2B_2cell.pdf")
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="2-cell" & tlDesign$quality=="Use")], falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2B_2cell_lab.pdf", labelSamples = T)
    
    # FOUR CELL
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="4-cell" & tlDesign$quality=="Use")], falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2E_4cell1.pdf")
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="4-cell" & tlDesign$quality=="Use")], falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2E_4cell1_lab.pdf", labelSamples = T)
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="4-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P2", "EMS")))],  falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2F_4cell2.pdf")
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="4-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P2", "EMS")))],  falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2F_4cell2_lab.pdf", labelSamples = T)
    
    # EIGHT CELL
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use")], falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2I_8cell1.pdf")
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use")], falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2I_8cell1_lab.pdf", labelSamples = T)
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P3", "C")))],  falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2J_8cell2.pdf")
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P3", "C")))],  falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2J_8cell2_lab.pdf", labelSamples = T)
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P3", "C", "E", "MS")))],  falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2K_8cell3.pdf")
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P3", "C", "E", "MS")))],  falseNeg = 0,
                    outPDF="04_analysesOut/Figure2/2K_8cell3_lab.pdf", labelSamples = T)
    
    
    
    # SIXTEEN CELL
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use")], topN=50, falseNeg=0,
                    outPDF="04_analysesOut/Figure2/2N_16cell1.pdf")
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use")], topN=50, falseNeg=0,
                    outPDF="04_analysesOut/Figure2/2N_16cell1_lab.pdf", labelSamples = T)
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" 
                                   & tlDesign$quality=="Use" 
                                   &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea")) 
                                   #   &(!rownames(tlDesign)%in%c("st315", "st353"))
                                   )], topN=50, falseNeg=0,
                                   outPDF="04_analysesOut/Figure2/2O_16cell2.pdf")
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" 
                                   & tlDesign$quality=="Use" 
                                   &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea")) 
                                   #   &(!rownames(tlDesign)%in%c("st315", "st353"))
                                   )], topN=50, falseNeg=0,
                    outPDF="04_analysesOut/Figure2/2O_16cell2_lab.pdf", labelSamples = T)
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" 
                                   & tlDesign$quality=="Use" 
                                   &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2")) 
                                   #   &(!rownames(tlDesign)%in%c("st315", "st353"))
                                   )], topN=50,    
                    outPDF="04_analysesOut/Figure2/2P_16cell3_1.pdf")
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" 
                                   & tlDesign$quality=="Use" 
                                   &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2")) 
                                   #   &(!rownames(tlDesign)%in%c("st315", "st353"))
                                   )], topN=50,
                    outPDF="04_analysesOut/Figure2/2P_16cell3_1_lab.pdf", labelSamples = T)
    
    # three of the replicates seem to be segregating Cx's very strongly. 
    # use those to ratchet out the Cs in the other replicates
    tempCs <- c("st322", "st321", "st357", "st358", "st530", "st531")
    tempNonCs <- rownames(tlDesign[which(tlDesign$stage=="16-cell" 
                                         & tlDesign$replicate%in%c(2,3,6)
                                         & !tlDesign$ID%in%c("P4", "D", "Ep", "Ea", "MSx1", "MSx2")
                                         & !rownames(tlDesign)%in%tempCs),])
    CxVsNonCGenes <- ERdiffExV1(tlRPKMs[,c(tempCs, tempNonCs)], 
                                sig = .05,
                                group = factor(c(rep("Cs", length(tempCs)), rep("nonCs", length(tempNonCs)))),
                                to="list")
    
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use" 
                                   &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2"))
                                   #   &(!rownames(tlDesign)%in%c("st315", "st353"))
    )], 
    forceGenes = union(rownames(CxVsNonCGenes[[1]]), rownames(CxVsNonCGenes[[2]])), falseNeg = 0,
    outPDF="04_analysesOut/Figure2/2P_16cell3_2.pdf")
    
    unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use" 
                                   &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2"))
                                   #   &(!rownames(tlDesign)%in%c("st315", "st353"))
    )], 
    forceGenes = union(rownames(CxVsNonCGenes[[1]]), rownames(CxVsNonCGenes[[2]])), falseNeg = 0,
    outPDF="04_analysesOut/Figure2/2P_16cell3_2_lab.pdf", labelSamples = T)
    
    # Remove all the functions and objects
    rm(CxVsNonCGenes, tbxSeedTemp, tempCs, tempNonCs)
}

makePCAs()


##
# Figure 4
# Overview of scRNA-seq dataset

if(!file.exists("04_analysesOut/Figure4/")){
    dir.create("04_analysesOut/Figure4/")
}

# Load in files to make cell maps
cellMapShapes <- list(
    "nameRef" = read.table("03_dataIn/4A_pictogramming/cellMapNameRef.txt", stringsAsFactors = F),
    "mapAllCells" = read.table("03_dataIn/4A_pictogramming/mapAllCells.txt", stringsAsFactors = F),
    "mapAllSafeCells" = read.table("03_dataIn/4A_pictogramming/mapAllSafeCells.txt", stringsAsFactors = F),
    "mapLineage" = read.table("03_dataIn/4A_pictogramming/mapLineage.txt", stringsAsFactors = F),
    "mapCellsHorz" = read.table("03_dataIn/4A_pictogramming/mapCellsHorz.txt", stringsAsFactors = F),
    "mapSafeHorz" = read.table("03_dataIn/4A_pictogramming/mapSafeHorz.txt", stringsAsFactors = F)
)

# Figure 4ACD
# Pictograms of gene expression

expectedGenes <- function(){
    mapExpression("sdz-38", mapType = "horz",
                  to = "04_analysesOut/Figure4/4A1_sdz-38.pdf")
    mapExpression("tbx-37",  mapType = "horz",
                  to = "04_analysesOut/Figure4/4A2_tbx-37.pdf")
    mapExpression("ceh-51",  mapType = "horz",
                  to = "04_analysesOut/Figure4/4A3_ceh-51.pdf")
    mapExpression("elt-7", mapType = "horz",
                  to = "04_analysesOut/Figure4/4A4_elt-7.pdf")
    mapExpression("cwn-1", mapType = "horz",
                  to = "04_analysesOut/Figure4/4A5_cwn-1.pdf")
    mapExpression("cey-2", mapType = "horz",
                  to = "04_analysesOut/Figure4/4A6_cey-2.pdf")
    
    mapExpression("daz-1", mapType = "horz",
                  to = "04_analysesOut/Figure4/4C_daz-1.pdf")
    
    mapExpression("skr-10", mapType = "horz",
                  to = "04_analysesOut/Figure4/4D_skr10.pdf")
}

expectedGenes()

# Figure 4B
# Export table for big blue heatmap

bigBlueHeatmapData <- function(){

    # isolate detectable genes
    subMatr <- tlRPKMs[unique(rownames(which(tlSafeCellAvg>25, arr.ind = T))),
                       rownames(tlDesign[which(tlDesign$quality=="Use"),])]

    # replace stNNN names with descriptive names
    newColNames = NULL
    for(colmn in 1:dim(subMatr)[2]){
        sample = colnames(subMatr)[colmn]
        newName = paste(tlDesign[sample, "stage"], tlDesign[sample, "replicate"], tlDesign[sample, "safeID"],sep="_")
        newColNames[colmn] = newName
    }
    colnames(subMatr) = newColNames
    
    # reorder columns by cell type
    goodCellTypes = tlDesign[which(tlDesign$quality=="Use"), "safeID"]

    colOrder = NULL
    for(i in 1:length(levels(tlDesign$safeID))){
        colOrder = c(colOrder, which(goodCellTypes==levels(tlDesign$safeID)[i]))
    }

    subMatr <- subMatr[,colOrder]
    
    # log matrix
    subMatr <- log10(subMatr+1)
    
    # export this matrix
    write.table(subMatr, "04_analysesOut/Figure4/4B1_allGenes.txt", sep="\t", col.names = NA, quote=F)

    # Put this dataset into cluster, and generate a cdt.
    # Open that in java tree view to flip branches around a bit for a nice layout
    
    # import cdt of this full .cdt dataset
    tempAllCells <-  read.csv("03_dataIn/4B_RPKMheatmap/allGenes.cdt", sep = "\t")
    # import the gene order list
    tempTest <- read.csv("03_dataIn/4B_RPKMheatmap/geneOrder.txt")
    tempReorder <- tempTest[,2]
    names(tempReorder) <- tempTest[,1]
    tempOut <- tempAllCells[c(1,tempReorder),]
    write.table(tempOut, "04_analysesOut/Figure4/4B2_allGenesOrdered.cdt", 
                sep="\t", quote=F, row.names = F)

    # Open this up in java tree view to generate a PNG

    # Clear out all these temporary files
    rm(subMatr, newColNames, sample, newName, goodCellTypes, colOrder, 
       tempAllCells, tempTest, tempReorder, tempOut)
}

bigBlueHeatmapData()

# Figure 4E-K
# Calculate and map statistics on each cell type

calcPlotStats <- function(){
    tlCellStats <- data.frame(matrix(nrow = length(safeCellTypes), ncol = 0))
    rownames(tlCellStats) <- safeCellTypes
    
    calcStats <- function(){}
    
    # Load in lineage relationship information
    source("03_dataIn/4E_linRelatinships/cellRelations.r")
    
    
    # Calculate stats
    tlCellStats$upReg <- 0
    for(cell in safeCellTypes[-which(safeCellTypes=="P0")]){
        if(cell %in% tlDesign$safeID){
            print(cell)
            # identify parent
            parent <- linRelations[[cell]][["parent"]]
            # find genes that are expressed in cell
            genesInCell <- names(which(tlSafeCellAvg[,cell]>25))
            # just keep genes that are twice the expression as in the parent
            newGenes <- NULL
            for(gene in genesInCell){
                if(tlSafeCellAvg[gene, parent]<(.5*tlSafeCellAvg[gene, cell])){newGenes <- c(newGenes, gene)}
            }
            tempNum <- length(newGenes)
            tlCellStats[cell,"upReg"] <- tempNum
        }
    }
    
    tlCellStats$downReg <-0
    for(cell in safeCellTypes[-which(safeCellTypes=="P0")]){
        if(cell %in% tlDesign$safeID){
            print(cell)
            # identify parent
            parent <- linRelations[[cell]][["parent"]]
            # find genes that are expressed in the parent cell
            genesInParent <- names(which(tlSafeCellAvg[,cell]>25))
            # just keep genes that are half the expression of the parent
            newGenes <- NULL
            for(gene in genesInParent){
                if(tlSafeCellAvg[gene, parent]>(2*tlSafeCellAvg[gene, cell])){newGenes <- c(newGenes, gene)}
            }
            tempNum <- length(newGenes)
            tlCellStats[cell,"downReg"] <- tempNum
        }
    }
    
    tlCellStats$uniqueGenes <- 0
    for(cell in safeCellTypes){
        if(cell %in% tlDesign$safeID){
            print(cell)
            genesInCell <- names(which(tlSafeCellAvg[,cell]>25))
            # exclude genes that are expressed anywhere else    
            otherMax <- apply(tlSafeCellAvg[genesInCell,-which(colnames(tlSafeCellAvg)==cell)], 1, max)
            newGenes <- NULL
            for(gene in genesInCell){
                if(tlSafeCellAvg[gene,cell]>(10*otherMax[gene])){newGenes <- c(newGenes, gene)}
            }
            tempNum <- length(newGenes)
            tlCellStats[cell,"uniqueGenes"] <- tempNum
        }
    }
    
    tlCellStats$uniqueGenesPercTF <- 0
    TFs <- read.table("03_dataIn/4H_transcriptionFactors/transcriptionFactors.txt", stringsAsFactors = F)[,1]
    for(cell in safeCellTypes){
        if(cell %in% tlDesign$safeID){
            print(cell)
            genesInCell <- names(which(tlSafeCellAvg[,cell]>25))
            # exclude genes that are expressed anywhere else    
            otherMax <- apply(tlSafeCellAvg[genesInCell,-which(colnames(tlSafeCellAvg)==cell)], 1, max)
            newGenes <- NULL
            for(gene in genesInCell){
                if(tlSafeCellAvg[gene,cell]>(10*otherMax[gene])){newGenes <- c(newGenes, gene)}
            }
            tempNum <- length(newGenes)
            tempTFnum <- length(intersect(newGenes, TFs))
            if(tempNum==0){tlCellStats[cell,"uniqueGenesPercTF"] <- 0}
            else{
                tempPerc <- (tempTFnum/tempNum)*100
                tlCellStats[cell,"uniqueGenesPercTF"] <- tempPerc
            }
        }
    }
    
    tlCellStats$pgRNA <- 0
    for(cell in safeCellTypes){
        if(cell %in% tlDesign$safeID){
            tempNum <- mean(tlDesign[which(tlDesign$safeID==cell),"input"])
            tlCellStats[cell,"pgRNA"] <- tempNum
        }
    }
    
    tlCellStats$presentOver25 <- 0
    for(cell in safeCellTypes){
        if(cell %in% tlDesign$safeID){
            tempNum <- length(which(tlSafeCellAvg[,cell]>25))
            tlCellStats[cell,"presentOver25"] <- tempNum
        }
    }
    
    tlCellStats$cycleLen <- 0
    tlCellStats[c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", 
                  "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                  "ABalx", "ABarx", "ABplx", "ABprx", "MSx", "Ex", "Cx", "D", "P4"), "cycleLen"] <- c(
                      40, 17, 18, 17, 17, 22, 25,
                      18, 18, 18, 18, 20, 23, 24, 33,
                      22, 22, 22, 22, 27, 45, 30, 40, 67)
    
    # Export statistic heatmap pictograms
    
    mapStats("upReg", legendTitle = "# genes\nupregulated\n",
             to = "04_analysesOut/Figure4/4E_upreg_legend.pdf")
    mapStats("downReg", legendTitle = "# genes\ndownregulated\n",
             to = "04_analysesOut/Figure4/4F_downreg_legend.pdf")
    mapStats("uniqueGenes", legendTitle = "# of unique\ngenes expressed\n",
             to = "04_analysesOut/Figure4/4G_uniqueGenes_legend.pdf")
    mapStats("uniqueGenesPercTF", legendTitle = "percentage of\nunique genes\nthat are\ntxn factors\n",
             to = "04_analysesOut/Figure4/4H_uniqueGenesPercTF_legend.pdf")
    mapStats("pgRNA", legendTitle = "pg mRNA\ndetected\n",
             to = "04_analysesOut/Figure4/4I_pgRNA_legend.pdf")
    mapStats("presentOver25", legendTitle = "# genes\ndetected\n", fillRange = c(0,max(tlCellStats$presentOver25)),
             to = "04_analysesOut/Figure4/4J_detectedOver25_legend.pdf")
    mapStats("cycleLen", legendTitle = "length of\ncell cycle\n(minutes)\n", fillRange = c(0,70),
             to = "04_analysesOut/Figure4/4K_cellCycle_legend.pdf")
    
    rm(cell, colmn, gene, genesInCell, genesInParent, i, linRelations, newGenes, 
       otherMax, parent, tempNum, tempPerc, tempTFnum, TFs)
    
}

calcPlotStats()

# Figure 4L
# Correlations of cell statistics

cellStatCors <- function(){
    nonP4 <- safeCellTypes[which(safeCellTypes!=c("P0", "P1", "P2", "P3", "P4"))]
    tempMatr <- tlCellStats[nonP4, c("upReg", "downReg", "uniqueGenes", "uniqueGenesPercTF", "pgRNA", "presentOver25", "cycleLen")]

    pdf("04_analysesOut/Figure4/4L_corrHeatmap.pdf", 5,5)
    heatmap.2(cor(tempMatr), Rowv = F, Colv = F, trace = 'none', 
              col=colorRampPalette(c('black', 'gray', 'white')), 
              labRow = c("E", "F", "G", "H", "I", "J", "K"),
              labCol = c("E", "F", "G", "H", "I", "J", "K"),
              density.info = 'none', colsep=4, rowsep=4, sepcolor = "red")
    dev.off()
    
    rm(nonP4)
}

cellStatCors()

# Figure 4M
# Correlations between cell type transcriptomes

transcriptomeCors <- function(){
    
    pdf("04_analysesOut/Figure4/4M1_transcriptomeCorrs.pdf", 5, 5)
    heatmap.2(cor(tlSafeCellAvg), trace = "none", density.info = "none", 
              col=colorRampPalette(c('#ffffff', '#fcbaaa', '#eb765a', '#d01e0a')))
    dev.off()
    
    # make a blank cell map to color code
    mapExpression("2L52.2", fillRange = c(0,100), legend = F, 
                  to="04_analysesOut/Figure4/4M2_blankCellMap.pdf")
    
}

transcriptomeCors()


##
# Figure 5
# Spatially dynamic gene expression

if(!file.exists("04_analysesOut/Figure5/")){
    dir.create("04_analysesOut/Figure5/")
}

tbx32AndCo <- function(){
    mapExpression("tbx-32", to="04_analysesOut/Figure5/5A1_tbx-32_vert.pdf")
    mapExpression("tbx-32", mapType = "lineage", to="04_analysesOut/Figure5/5A2_tbx-32_lin.pdf")
    mapExpression("tbx-31", to="04_analysesOut/Figure5/5C1_tbx-31.pdf")
    mapExpression("tbx-40", to="04_analysesOut/Figure5/5C2_tbx-40.pdf")
    mapExpression("Y43D4A.6", to="04_analysesOut/Figure5/5C3_Y43D4A.6.pdf")
    mapExpression("Y116A8C.20", to="04_analysesOut/Figure5/5C4_Y116A8C.20.pdf")
    mapExpression("ZK666.1", to="04_analysesOut/Figure5/5C5_ZK666.1.pdf")
}

tbx32AndCo()


##
# Figure 6
# RNAi experiments

if(!file.exists("04_analysesOut/Figure6/")){
    dir.create("04_analysesOut/Figure6/")
}

# Figure 6A-C
# Synexpressed Paralogs

# Load homology data
ceHoms <- read.table("03_dataIn/6A_homology/homologyGroups.txt", 
                     header = F, stringsAsFactors = F)
colnames(ceHoms) <- c("gene", "CHOM")

# Figure 6A
# Correlations of homology groups from true data

trueHomCors <- function(){
    # Reduce the list to only those that are expressed
    expressed6266 <- rownames(tlSafeCellMed[which(apply(tlSafeCellMed, 1, max)>25),])
    subCeHoms <- ceHoms[which(ceHoms[,"gene"] %in% expressed6266),]
    
    # Turn that into a list organized by homology group (CHOM)
    ceHomList <- list()
    n=0
    for(tempCHOM in unique(subCeHoms$CHOM)){
        n<- n+1
        print(n)
        ceHomList[[tempCHOM]] <- subCeHoms[which(subCeHoms$CHOM==tempCHOM),"gene"]
    }
    
    # Get rid of the homology groups that have less than two genes or more than 5 genes
    for(tempCHOM in names(ceHomList)){
        if(length(ceHomList[[tempCHOM]])<2){
            ceHomList[[tempCHOM]] <- NULL
        }
        if(length(ceHomList[[tempCHOM]])>5){
            ceHomList[[tempCHOM]] <- NULL
        }
    }
    
    # Remove all the sets that are equivalent but have different names
    for(tempCHOM in names(ceHomList)){
        aVec <- ceHomList[[tempCHOM]]
        i <- which(names(ceHomList)==tempCHOM)
        print(paste(i, "of", length(ceHomList)))
        for(j in (i+1):length(ceHomList)){
            bVec <- ceHomList[[j]]
            if(all(aVec %in% bVec)){
                ceHomList[[tempCHOM]] <- NULL
                break
            }
        }
    }
    
    # Make a table that evaluates the number of genes and average correlation
    chomCors <- matrix(0, nrow=length(ceHomList), ncol = 2, 
                       dimnames = list(names(ceHomList), c("numGenes", "avgCor")))
    for(tempCHOM in names(ceHomList)){
        # make a table of correlations
        tempCors <- as.vector(cor(t(tlSafeCellAvg[ceHomList[[tempCHOM]],])))
        # calculate the average of all non-1 correlations
        tempAvg <- mean(tempCors[which(tempCors!=1)])
        # push this information to the table
        chomCors[tempCHOM, c("numGenes", "avgCor")] <- c(length(ceHomList[[tempCHOM]]), tempAvg)
    }
    
    # How does it look?
    pdf("04_analysesOut/Figure6/6A_corrsOfSets.pdf", 3,4)
    plot(jitter(chomCors[,"numGenes"], factor = .75), 
         chomCors[,"avgCor"], pch=16, col="#00000030", xlim=c(1.5,5.5), ylim=c(-1,1),
         xlab = "number of genes in\nparalogous set", ylab="correlation of genes' expression patterns", main = "true data")
    abline(h = 0)
    dev.off()
    
    rm(chomCors, subCeHoms, aVec, bVec, ceHomList, 
       expressed6266, i, j, n, tempAvg, tempCHOM, tempCors)
}

trueHomCors()

# Figure 6B
# Correlations of homology groups from scrambled data

scramHomCors <- function(){
    # Make a scrambled dataset
    set.seed(1)
    scramGeneNames <- sample(row.names(tlSafeCellMed))
    scramSafeCellMed <- tlSafeCellMed
    rownames(scramSafeCellMed) <- scramGeneNames
    scramSafeCellAvg <- tlSafeCellAvg
    rownames(scramSafeCellAvg) <- scramGeneNames
    
    # Reduce the list to only those that are expressed
    expressed6266 <- rownames(scramSafeCellMed[which(apply(scramSafeCellMed, 1, max)>25),])
    subCeHoms <- ceHoms[which(ceHoms[,"gene"] %in% expressed6266),]
    
    # Turn that into a list organized by homology group (CHOM)
    ceHomList <- list()
    n=0
    for(tempCHOM in unique(subCeHoms$CHOM)){
        ceHomList[[tempCHOM]] <- subCeHoms[which(subCeHoms$CHOM==tempCHOM),"gene"]
    }
    
    # Get rid of the homology groups that have less than two genes or more than 5 genes
    for(tempCHOM in names(ceHomList)){
        if(length(ceHomList[[tempCHOM]])<2){
            ceHomList[[tempCHOM]] <- NULL
        }
        if(length(ceHomList[[tempCHOM]])>5){
            ceHomList[[tempCHOM]] <- NULL
        }
    }
    
    for(tempCHOM in names(ceHomList)){
        aVec <- ceHomList[[tempCHOM]]
        i <- which(names(ceHomList)==tempCHOM)
        print(paste(i, "of", length(ceHomList)))
        for(j in (i+1):length(ceHomList)){
            bVec <- ceHomList[[j]]
            if(all(aVec %in% bVec)){
                ceHomList[[tempCHOM]] <- NULL
                break
            }
        }
    }
    
    # Make a table that evaluates the number of genes and average correlation
    scramChomCors <- matrix(0, nrow=length(ceHomList), ncol = 2, 
                            dimnames = list(names(ceHomList), c("numGenes", "avgCor")))
    for(tempCHOM in names(ceHomList)){
        # make a table of correlations
        tempCors <- as.vector(cor(t(scramSafeCellAvg[ceHomList[[tempCHOM]],])))
        # calculate the average of all non-1 correlations
        tempAvg <- mean(tempCors[which(tempCors!=1)])
        # push this information to the table
        scramChomCors[tempCHOM, c("numGenes", "avgCor")] <- c(length(ceHomList[[tempCHOM]]), tempAvg)
    }
    
    # Plot this scrambled version
    pdf("04_analysesOut/Figure6/6B_scrambledCorrsOfSets.pdf", 3,4)
    plot(jitter(scramChomCors[,"numGenes"], factor = .75), 
         scramChomCors[,"avgCor"], pch=16, col="#00000030", xlim=c(1.5,5.5), ylim=c(-1,1), 
         xlab = "number of genes in\nparalogous set", ylab="correlation of genes' expression patterns", main = "scrambled data")
    abline(h = 0)
    dev.off()
    
    rm(scramChomCors, scramSafeCellAvg, scramSafeCellMed, subCeHoms,
       aVec, bVec, ceHomList, expressed6266, i, j, n, tempAvg, tempCHOM, tempCors)
}

scramHomCors()

# Figure 6C
# Make histogram of 100 scrambled datasets

histo100scrambs <- function(){
    
    randoRedundTrials <-c()
    
    for(seedn in 1:100){
        set.seed(seedn)
        scramGeneNames <- sample(row.names(tlSafeCellMed))
        scramSafeCellMed <- tlSafeCellMed
        rownames(scramSafeCellMed) <- scramGeneNames
        scramSafeCellAvg <- tlSafeCellAvg
        rownames(scramSafeCellAvg) <- scramGeneNames
        
        # Reduce the list to only those that are expressed
        expressed6266 <- rownames(scramSafeCellMed[which(apply(scramSafeCellMed, 1, max)>25),])
        subCeHoms <- ceHoms[which(ceHoms[,"gene"] %in% expressed6266),]
    
        # Turn that into a list organized by homology group (CHOM)
        ceHomList <- list()
        for(tempCHOM in unique(subCeHoms$CHOM)){
            ceHomList[[tempCHOM]] <- subCeHoms[which(subCeHoms$CHOM==tempCHOM),"gene"]
        }
        
        # Get rid of the homology groups that have less than two genes or more than 5 genes
        for(tempCHOM in names(ceHomList)){
            if(length(ceHomList[[tempCHOM]])<2){
                ceHomList[[tempCHOM]] <- NULL
            }
            if(length(ceHomList[[tempCHOM]])>5){
                ceHomList[[tempCHOM]] <- NULL
            }
        }
    
        for(tempCHOM in names(ceHomList)){
            aVec <- ceHomList[[tempCHOM]]
            i <- which(names(ceHomList)==tempCHOM)
            len <- length(ceHomList)
            print(paste(i, "of", length(ceHomList)))
            if(i<len){
                for(j in (i+1):len){
                    bVec <- ceHomList[[j]]
                    if(all(aVec %in% bVec)){
                        ceHomList[[tempCHOM]] <- NULL
                        break
                    }
                }
            }
        }
    
        # Make a table that evaluates the number of genes and average correlation
        scramChomCors <- matrix(0, nrow=length(ceHomList), ncol = 2, 
                                dimnames = list(names(ceHomList), c("numGenes", "avgCor")))
        for(tempCHOM in names(ceHomList)){
            # make a table of correlations
            tempCors <- as.vector(cor(t(scramSafeCellAvg[ceHomList[[tempCHOM]],])))
            # calculate the average of all non-1 correlations
            tempAvg <- mean(tempCors[which(tempCors!=1)])
            # push this information to the table
            scramChomCors[tempCHOM, c("numGenes", "avgCor")] <- c(length(ceHomList[[tempCHOM]]), tempAvg)
        }
        
        # determine number of sets of genes
        randoRedundTrials <- c(randoRedundTrials, dim(scramChomCors[which(scramChomCors[,"avgCor"]>.25),])[1])  #295 groups of genes in real data
        
        # make a plot to show progress
        hist(c(randoRedundTrials, 295), breaks=25,
             xlab="number of paralogous sets", ylab="frequency",
             main = paste("paralogous gene set analysis\ndone with randomized data : ", seedn, sep=""))
        
    }
    
    # plot this final set of trails. include the result 
    #   of the real data, which is 295 genes
    
    pdf("04_analysesOut/Figure6/6C_histo100scrambles.pdf", 3,3)
    hist(c(randoRedundTrials,295), breaks=25, xlim=c(0,300), ylim=c(0,30), 
         xlab="number of paralogous sets found in a trial", 
         ylab="frequency", axes = F,
         main = "paralogous gene set analysis\ndone with randomized data", 
         col = c(rep("gray", 15), rep("red", 15)))
    axis(1, at=c(0,100, 200, 300))
    axis(2, at=c(0,10,20,30))
    dev.off()
    
    rm(scramChomCors, scramSafeCellAvg, scramSafeCellMed, subCeHoms,
       aVec, bVec, ceHomList, expressed6266, i, j, n, tempAvg, tempCHOM, 
       tempCors, len, randoRedundTrials, scramGeneNames, seedn)
    
}

histo100scrambs()

# Figure 6D
# RNAi targeting T12E24.1 and T12E24.13 in N2

RNAi_N2 <- function(){

    RNAi.T24 <- data.frame(hours = rep(c("12-24", "24-36", "36-48"),3), 
                           lethality=c(.3642,.7692,.6455, .2354, .5353, .699, .827, .8349, .7034),
                           total = c(1219,1109,347, 1028,934,289, 1578, 1169, 263),
                           dead = rep(0,9), max95 = rep(0,9), min95=rep(0,9),
                           condition = c(rep("both", 3), rep("T24E12.1", 3), rep("T24E12.13", 3)))
    
    # Calculate 95% confidence intervals
    RNAi.T24$dead <- round(RNAi.T24$lethality * RNAi.T24$total)
    for(rown in 1:dim(RNAi.T24)[1]){
        temp <- RNAi.T24[rown,"lethality"] + c(-qnorm(0.975),qnorm(0.975))*sqrt((rown/RNAi.T24[1,"total"])*RNAi.T24[rown,"lethality"]*(1- RNAi.T24[rown,"lethality"]))
        RNAi.T24[rown,c("min95", "max95")] <- temp
    }
    
    c <- ggplot(RNAi.T24[which(RNAi.T24$condition=="both"),], 
                aes(x=factor(hours), y=lethality))+
        geom_bar(stat="identity", width=.75, color="black", fill="dark gray")+
        xlab("Hours after dsRNA injection\nof T24E12.1 and T24E12.13")+
        ylab("Embryonic lethality")+
        scale_x_discrete(breaks=c(1,2,3),  labels=c("12", "24", "36"))+
        geom_errorbar(aes(ymax=RNAi.T24[which(RNAi.T24$condition=="both"),"max95"], 
                          ymin=RNAi.T24[which(RNAi.T24$condition=="both"), "min95"]), 
                      position=position_dodge(width=0.9), width=0.15)+
        theme_classic()+
        theme(axis.line.x = element_line(color = "black"), 
              axis.line.y = element_line(color = "black"))
    
    c
    
    ggsave("04_analysesOut/Figure6/6D_RNAiInN2.pdf", width = 2, height = 3)
    rm(RNAi.T24, c)
}

RNAi_N2()

# Figure 6E
# Pictograms of T24.. genes

mapExpression("T24E12.1", to="04_analysesOut/Figure6/6E1_T24E12.1.pdf")
mapExpression("T24E12.13", to="04_analysesOut/Figure6/6E2_T24E12.12.pdf")


##
# Supplementary Figure 2
# Extra pieces for IDing cell types of transcriptomes

if(!file.exists("04_analysesOut/SuppFig2/")){
    dir.create("04_analysesOut/SuppFig2/")
}

makeIDextras <- function(){
    
    # SF1A
    # Parsing D from P4 transcriptomes
    
    unsupervisedPCA(tlRPKMs[,c("st263", "st264", "st327", "st328", "st359", "st360", "st532", "st533", "st478", "st479", "st519", "st520")], falseNeg = 0,
                    outPDF="04_analysesOut/SuppFig2/2A_D_P4.pdf")
    tlDesign[c("st263", "st328", "st360", "st532", "st478", "st519"), "ID"] <- "P4"
    tlDesign[c("st264", "st327", "st359", "st533", "st479", "st520"), "ID"] <- "D"
    
    # SF2B
    # Parsing ABal from ABar transcriptomes
    
    tbxSeedTemp <- ERdiffExV1(matr=tlCounts[,c("st226", "st227", "st370", "st372")], 
                              to = "list",  sig=.05,
                              group = factor(rep(c("tbx+", "tbx-"), 2), levels = c("tbx+", "tbx-")))
    
    tbxSeedTemp <- c(rownames(tbxSeedTemp[[1]]), rownames(tbxSeedTemp[[2]]))
    
    unsupervisedPCA(tlRPKMs[,c("st226", "st227", "st337", "st338", "st370", "st372", "st393", "st396", "st482", "st484")],  
                    falseNeg = 6, forceGenes = tbxSeedTemp, labelSamples = T,
                    outPDF="04_analysesOut/SuppFig2/2B_ABal_ABar_lab.pdf")
    
    unsupervisedPCA(tlRPKMs[,c("st226", "st227", "st337", "st338", "st370", "st372", "st393", "st396", "st482", "st484")],  
                    falseNeg = 6, forceGenes = tbxSeedTemp, 
                    outPDF="04_analysesOut/SuppFig2/2B_ABal_ABar.pdf")
    
    # Figure S2C
    # Ratchet back from ABal+ABar vs ABpl+ABpr, to ABa vs ABp
    
    ABaxVsABpxGenes <- ERdiffExV1(matr=tlCounts[,typeCells(c("ABal", "ABar", "ABpl", "ABpr"))], 
                                  group = factor(c(rep("ABax", length(typeCells(c("ABal", "ABar")))), 
                                                   rep("ABpx", length(typeCells(c("ABpl", "ABpr"))))), levels = c("ABax", "ABpx")), 
                                  postThresh = 25, to = "list", sig=.5)
    
    ABaVsABpGenes <- ERdiffExV1(matr=tlCounts[,typeCells(c("ABa", "ABp"))], 
                                postThresh = 25, to = "list", sig=.5)
    
    # make a heatmap of these overlaps:
    pdf("04_analysesOut/SuppFig2/2C_ratchetHeatmap.pdf", 5, 5)
    heatmap.2(matrix(c(6.996,6.049,9.021,9.462), nrow = 2), trace = "none", Rowv=F, Colv=F, 
              labRow = c("ABa", "ABp"), labCol = c("ABpx", "ABax"), density.info = "none")
    dev.off()
    
    # Figure S2D
    # PCAs of ABxxx cells using notch target transcripts
    
    ABxxxPCA(geneList = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", "tbx-38", "tbx-39", "ref-1"), 
             outPDF = "04_analysesOut/SuppFig2/2D_all_ABxxx.pdf", splitReps = T)
    
    ABxxxPCA(geneList = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", "tbx-38", "tbx-39", "ref-1"), 
             outPDF = "04_analysesOut/SuppFig2/2D_all_ABxxx_lab.pdf", splitReps = T, labelSamps = T)
    
}

makeIDextras()


##
# Supplementary Figure 3
# Comparisons of present data to previous studies

if(!file.exists("04_analysesOut/SuppFig3/")){
    dir.create("04_analysesOut/SuppFig3/")
}

# Supplementary Figure 3A
# AB/P1 enrichment index comparisons (present, Osborne-Nishimura, Hashimshony)

makeAllABP1compPlots <- function(){

    # Import fold change and avg expression data (calculated using edgeR)
    all2cellFCs <- read.table("03_dataIn/S3_prevStudies/ABP1_FC.txt")
    all2cellCPMs <- read.table("03_dataIn/S3_prevStudies/ABP1_CPMs.txt")
    all2cellFCgenes <- rownames(all2cellFCs)
    
    # Establish the genes each studied ID'ed as significantly enriched
    erinIDgenes <- c("erm-1", "Y69H2.3", "C50E3.13", "F32D1.6", "afd-1", "F40G12.11", "lem-3", "C50E3.12", "cdc-25.3", "W02F12.3", "F52B5.2", "his-24", "C44B9.3", "C55C3.5", "F14H3.4", "R09A8.2", "ZK507.6", "T09B4.1", "F23F12.9", "pkn-1", "F13H10.4", "E02H4.6", "age-1", "F33D11.12", "F14D7.2", "ape-1", "M116.5", "F14H3.6", "C17H11.4", "C50D2.9", "T05H4.6a", "T21B10.1", "F08F3.6", "Y39G10AR.9", "C50B6.3", "cyk-1", "Y5F2A.4", "F02H6.2", "hlh-2", "F16A11.1", "K02C4.3", "mex-3", "syn-16", "F22G12.5", "T02G5.12", "epg-2", "sdz-30", "D1069.3", "aman-2", "apa-2", "pcaf-1", "T24E12.1", "Y32H12A.8", "hmp-2", "frm-7", "unc-73", "F23H11.4", "sma-1", "cyy-1", "dnc-1", "Y48G10A.1", "dpf-3", "Y73F8A.24", "C17G10.1", "cht-1", "akt-1", "R53.2", "exoc-7", "K01G12.3", "C37C3.9", "C14C11.2", "lsd-1", "nob-1", "let-2", "tag-224", "gon-2", "C52B9.8", "Y54F10AM.11", "Y53C10A.6", "aph-2", "chs-1", "cyp-31A3", "pgal-1", "clu-1", "egg-2", "F17C11.2", "cyp-31A2", "sqv-4", "rme-2", "T01H3.4", "T22F3.3", "cpar-1", "T27C4.1", "sac-1", "fbf-2", "R10H10.7", "K08F4.3", "T04A8.7", "F12A10.8", "lbp-9", "daf-22", "zig-7", "tdc-1", "Y113G7B.16", "F57B10.3", "NA", "rab-5", "T24H10.1", "T25G3.4", "cpg-2", "H05L03.3", "pgl-3", "C05C8.6", "T23B12.6", "spdl-1", "sbt-1", "R05F9.6", "C30F12.3", "C53C7.5", "T24D1.3", "F54E12.2", "smgl-1", "unc-60", "rab-1", "pyc-1", "flp-12", "nlp-12", "cpg-3", "cogc-2", "K07D4.9", "W08F4.3", "rme-1", "Y39E4B.5", "cab-1", "F36F2.2", "egl-21", "F23A7.8", "sqv-1", "F54D5.2", "egg-1", "T02E9.5", "F40B5.2", "oac-39", "F22B3.4", "F16D3.4", "dhs-28", "H36L18.2", "nlp-20", "B0272.4", "acl-2", "W05H9.4", "rskn-1", "Y75B8A.3", "T14G10.5", "puf-3", "T12G3.4", "nep-2", "K08D12.3", "Y119D3B.21", "clk-2", "acs-7", "C28H8.11", "W06F12.2", "C25D7.10", "K07H8.3", "C40H1.6", "F09E5.3", "hip-1", "Y24D9A.8", "dhs-20", "yop-1", "F52F12.7", "ppw-2", "acdh-13", "T21H3.1", "vars-1", "T14B4.1", "F59B2.2", "F57F4.4", "mdh-1", "dhs-18", "tin-13", "nlp-21", "Y39A3CL.1", "btb-19", "W01A11.2", "rps-4", "acs-4", "C07D8.6", "B0334.5", "sar-1", "pmp-4", "F32A5.4", "rpl-31", "R05F9.9", "M01E11.1", "F26A10.2", "pod-2", "sqv-2", "nmt-1", "rps-14", "W05F2.3", "F07A11.2", "lips-15", "M60.4", "oxi-1", "Y82E9BR.14", "rpt-6", "atf-7", "C46F11.3", "rpl-14", "rps-16", "rpl-28", "mtx-1", "mppa-1", "F25H2.5", "Y54F10AL.1", "him-14", "nhr-114", "F22F7.1", "rab-8", "prx-13", "F57B10.14", "Y53C12B.1", "spsb-2", "NA", "cpn-1", "Y39B6A.42", "vha-7", "sox-4", "sek-1", "T12B3.4", "rps-6", "rps-23", "rps-13", "B0001.2", "Y69A2AR.3", "rpl-10", "F16A11.2", "C48B4.6", "C37A2.7", "snt-4", "haf-2", "T20B12.3", "Y110A7A.19", "rla-1", "F49E2.5", "Y51F10.2", "bpl-1", "CC8.2", "F42C5.9", "D2023.6", "F25G6.8", "jtr-1", "rps-9", "plk-3", "Y38F2AR.9", "ubc-12", "mlc-2", "rmo-1", "rpl-32", "rps-8", "hsp-60", "Y37E3.8", "rps-30", "rps-22", "arx-6", "K05C4.2", "gsk-3", "Y71H2AM.7", "nol-6", "flp-22", "C34G6.1", "ZC410.5", "cutl-3", "rps-12", "C03H5.2", "pfn-1", "vps-33.2", "pph-4.1")
    # 280 genes
    hashimIDgenes <- c("mex-3", "par-3", "sca-1", "F32D1.6", "mig-5", "T21B10.4", "Y59A8B.12", "E02H4.6", "erm-1", "hsp-4", "ZK1127.6", "hmp-2", "Y73F8A.24", "T23B12.6", "pyc-1", "F07A11.2", "Y18H1A.7")
    # 17 genes
    my2cellIDgenes <- ERdiffExV1(tlRPKMs[,typeCells(c("AB", "P1"))], postThresh = 25, to="list", sig = .05)
    my2cellIDgenes <- rownames(rbind(my2cellIDgenes[[1]], my2cellIDgenes[[2]]))
    # 236 genes
    
    # make plots
    enrichMatr <- makeEnrichTable(constant=4)
    
    # Tintori vs Hashimshony
    pdf("04_analysesOut/SuppFig3/A1_tintoriHashimshony_c4.pdf",4,4)
    par(mar=c(5.1, 5.1, 2.1, 2.1))
    ABcompPlot("Tintori", "Hashimshony", axisLim = 20000)
    dev.off()
    
    # Hashimshony v Nishimura
    pdf("04_analysesOut/SuppFig3/A2_hashimshonyNishimura_c4.pdf",4,4)
    par(mar=c(5.1, 5.1, 2.1, 2.1))
    ABcompPlot("Hashimshony", "Nishimura", axisLim = 20000)
    dev.off()
    
    # Nishimura v Tintori
    pdf("04_analysesOut/SuppFig3/A3_nishimuraTintori_c4.pdf",4,4)
    par(mar=c(5.1, 5.1, 2.1, 2.1))
    ABcompPlot("Nishimura", "Tintori", axisLim = 20000)
    dev.off()
    
    rm(all2cellFCgenes, all2cellFCs, all2cellCPMs, enrichMatr, 
       erinIDgenes, hashimIDgenes, my2cellIDgenes)
}

makeAllABP1compPlots()

# Supplementary Figure 3B
# Upregulated and Downregulated transcripts (present, Baugh study)

makeBaughHeatmaps <- function(){
    
    # Import Ryan Baugh's 2005 microarray data timecourse for the relevant time points
    
    PAbaugh2005 <- read.table("03_dataIn/S3_prevStudies/baugh_avgScores.txt", sep="\t", header = T, stringsAsFactors = F)
    GNbaugh2005 <- read.table("03_dataIn/S3_prevStudies/baugh_probeNameLookup.txt", sep="\t", header=T, stringsAsFactors = F)
    GN2baugh2005 <- read.table("03_dataIn/S3_prevStudies/baugh_geneNames.txt", sep="\t", header=T, stringsAsFactors = F)
    GNbaugh2005[GNbaugh2005=="---"] <- NA
    GN2baugh2005[GN2baugh2005=="---"] <- NA
    GNbaugh2005 <- merge(GNbaugh2005, GN2baugh2005, by.x = "Probe.Set.ID", by.y="PROBEID")
    GNbaugh2005 <- cbind(GNbaugh2005, tlMatchID=NA)
    
    # Scour alternative gene names for potentially out-of-date names
    # Set up alt gene name database
    tempGeneNameMatr <- read.table("03_dataIn/S3_prevStudies/altGeneNames.txt", sep = "\t", stringsAsFactors = F)
    geneNameMatr <- matrix(NA, nrow = length(unique(tempGeneNameMatr$V1)), ncol = length(unique(tempGeneNameMatr$V2)))
    rownames(geneNameMatr) <- unique(tempGeneNameMatr$V1)
    for(rown in 1:dim(tempGeneNameMatr)[1]){
        print(paste(rown, "in", dim(tempGeneNameMatr)[1]))
        tempGene <- tempGeneNameMatr[rown,1]
        tempCol <- tempGeneNameMatr[rown,2]
        tempAlt <- tempGeneNameMatr[rown,3]
        geneNameMatr[tempGene, tempCol] <- tempAlt
    }
    # Change baugh's gene names, where appropriate
    geneNameMatr <- rearrangeNameMatrix(geneNameMatr, rownames(tlRPKMs))
    for(rowNum in 1:length(GNbaugh2005[,1])){
        print(rowNum) # should go up to 26,433
        for(attempt in GNbaugh2005[rowNum, unique(c("display", "ALIAS", "Gene.Symbol", "Ensembl", "WORMBASE"))]){
            if(attempt %in% rownames(tlRPKMs)){
                GNbaugh2005[rowNum,"tlMatchID"] <- attempt
                break
            } else {
                transformedOption <- renameGeneVect(attempt, geneNameMatr)
                if(transformedOption %in% rownames(tlRPKMs)){
                    GNbaugh2005[rowNum, "tlMatchID"] <- transformedOption
                    break
                }
            }
        }
    }
    # 88 of which are dead/weirdly named and won't match up. Not that bad, out of dozens of thousands.
    
    baugh2005 <- merge(GNbaugh2005[,c("Probe.Set.ID", "tlMatchID")], 
                       PAbaugh2005[,c("QUALIFIER", "PC.6", "X3C", "PC32", "X0", "X23", "X41", "X53", "X66")], 
                       by.x = "Probe.Set.ID", by.y = "QUALIFIER")
    
    # Remove all the duplicates
    for(gene in unique(baugh2005$tlMatchID)){
        geneValueSet <- which(baugh2005$tlMatchID == gene)
        if(length(geneValueSet) > 1){
            newLine <- baugh2005[geneValueSet[1],]
            for(samp in c("PC.6", "X3C", "PC32", "X0", "X23", "X41", "X53", "X66")){
                newLine[,samp] <- mean(baugh2005[geneValueSet, samp])
            }
            baugh2005[geneValueSet[1],] <- newLine
            baugh2005 <- baugh2005[-geneValueSet[2:length(geneValueSet)],]
        }
    }
    
    baugh2005 <- baugh2005[which(baugh2005$tlMatchID != "NA"),]
    temp <- baugh2005[,c("PC.6", "X3C", "PC32", "X0", "X23", "X41", "X53", "X66")]    
    rownames(temp) <- baugh2005[,"tlMatchID"]
    baugh2005 <- temp
    
    rm(PAbaugh2005, GNbaugh2005, GN2baugh2005, mergeBaugh2005, temp, 
       tempGeneNameMatr, tempGene, tempCol, tempAlt, rown, rowNum, attempt, 
       transformedOption, gene, geneValueSet, newLine, samp)
    
    # Curate lists of transcribed and degraded genes, from both my and Ryan's data:
    thresh = 25
    
    # what genes have increasing expression in baugh's data
    baughTransGenes <- NULL
    for(gene in rownames(baugh2005[which(apply(baugh2005, 1, max)>thresh),])){
        # find the high point in expression
        tempMax <- max(baugh2005[gene,])
        tempWhichMax <- which(baugh2005[gene,]==tempMax)
        tempWhichMax <- tempWhichMax[length(tempWhichMax)]
        # for all the timepoints before that, is the minimum considerably lower?
        if(tempWhichMax==1){next}
        tempMin <- min(baugh2005[gene, 1:tempWhichMax-1])
        if(tempMin*2 < tempMax){
            baughTransGenes <- c(baughTransGenes, gene)
        }
    }
    
    # what genes have decreasing expression in baugh's dataset
    baughDegGenes <- NULL
    for(gene in rownames(baugh2005[which(apply(baugh2005, 1, max)>thresh),])){
        # find the low point in expression
        tempMin <- min(baugh2005[gene,])
        tempWhichMin <- which(baugh2005[gene,]==tempMin)
        tempWhichMin <- tempWhichMin[length(tempWhichMin)]
        # for all the timepoints before that, is the maximum considerably higher?
        if(tempWhichMin == 1){ next }
        tempMax <- max(baugh2005[gene,1:tempWhichMin-1])
        if(tempMin*2 < tempMax){
            baughDegGenes <- c(baughDegGenes, gene)
        }
    }
    
    source("03_dataIn/4E_linRelatinships/cellRelations.r")
    
    # what genes are transcribed in my dataset
    myTransGenes <- NULL
    for(gene in rownames(tlSafeCellAvg[which(apply(tlSafeCellAvg, 1, max)>thresh),])){
        # for each cell that has expression above threshold..
        for(cellWith in names(sort(which(tlSafeCellAvg[gene,]>thresh), decreasing = T))){
            # check the minimum of ancestor cells, to see if any of them are significantly lower than the cellWith
            cellAncestors <- as.vector(unlist(linRelations[[cellWith]][c("parent", "grandparent", "greatgrandparent", "greatgreatgrandparent")]))
            if(length(cellAncestors)==0){next}
            tempMin <- min(tlSafeCellAvg[gene, cellAncestors])
            if(tempMin*2 < tlSafeCellAvg[gene, cellWith]){
                myTransGenes <- c(myTransGenes, gene)
                break
            }
        }
    }
    
    # what genes are degraded in my dataset?
    myDegGenes <- NULL
    for(gene in rownames(tlSafeCellAvg[which(apply(tlSafeCellAvg, 1, max)>thresh),])){
        # for each cell that has expression above threshold..
        for(cellWith in names(which(tlSafeCellAvg[gene,]>thresh))){
            # check the minimum of descendant cells, to see if any of them are significantly lower than the cellWith
            cellDescendants <- as.vector(unlist(linRelations[[cellWith]][c("daughters", "granddaughters", "greatgranddaughters", "greatgreatgranddaughters")]))
            cellDescendants <- intersect(cellDescendants, safeCellTypes)
            if(length(cellDescendants)==0){next}
            tempMin <- min(tlSafeCellAvg[gene, cellDescendants])
            if(tempMin*2 < tlSafeCellAvg[gene, cellWith]){
                myDegGenes <- c(myDegGenes, gene)
                break
            }
        }
    }
    
    # what are the genes that both studies had access to?
    tintoBaughCommonAll <- intersect(rownames(tlRPKMs), rownames(baugh2005))
    
    # export these groups of genes to make matrices
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B1_transBaughOnly.txt", 
                        geneSet=setdiff(baughTransGenes, myTransGenes))
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B2_transBoth.txt", 
                        geneSet=intersect(baughTransGenes, myTransGenes))
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B3_transMine_baughAccess.txt", 
                        geneSet=intersect(setdiff(myTransGenes, baughTransGenes), tintoBaughCommonAll))
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B4_transMine_baughNoAccess.txt", 
                        geneSet=setdiff(setdiff(myTransGenes, baughTransGenes), tintoBaughCommonAll))
    
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B5_degBaughOnly.txt", 
                        geneSet=setdiff(baughDegGenes, myDegGenes))
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B6_degBoth.txt", 
                        geneSet=intersect(baughDegGenes, myDegGenes))
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B7_degMine_baughAccess.txt", 
                        geneSet=intersect(setdiff(myDegGenes, baughDegGenes), tintoBaughCommonAll))
    exportGeneExpMatrix(outFile="04_analysesOut/SuppFig3/B8_degMine_baughNoAccess.txt", 
                        geneSet=setdiff(setdiff(myDegGenes, baughDegGenes), tintoBaughCommonAll))
    
    rm(baugh2005, thresh, baughTransGenes, baughDegGenes, 
       myTransGenes, myDegGenes, tintoBaughCommonAll)

}

makeBaughHeatmaps()

##
# Supplementary Figure 4
# RNAi on rrf-3 worms

if(!file.exists("04_analysesOut/SuppFig4/")){
    dir.create("04_analysesOut/SuppFig4/")
}

RNAi_rrf3 <- function(){

    RNAiResults <- matrix(0, nrow = 11, ncol=3,
                          dimnames = list(c("no-inj", "par-6", "set1", "set2", "set3", 
                                            "set4", "set5", "set6", "set7", "set8", "set9"),
                                          c("12-24 hrs", "24-36 hrs", "36-48 hrs")))
    
    RNAiResults[1:6,"12-24 hrs"] <- c(7.92,98.95,9.24,11.85,12.84,12.79)
    RNAiResults[1:6,"24-36 hrs"] <- c(14.29,100,6.51,18.75,14.85,12.07)
    RNAiResults[1:6,"36-48 hrs"] <- c(28.57,100,12.5,13.56,17.03,5.41)
    
    RNAiResults[7:11,"12-24 hrs"] <- c(16.4,9.26,10.8,50.12,17.48)
    RNAiResults[7:11,"24-36 hrs"] <- c(15.4,14.51,13.12,94.23,10.62)
    RNAiResults[7:11,"36-48 hrs"] <- c(14.29,25.29,6.25,91.3,12.86)
    
    library(reshape2)
    tempRNAiPlot <- melt(RNAiResults)
    colnames(tempRNAiPlot) <- c("condition", "timepoint", "lethality")
    timepoints <- as.character(unique(tempRNAiPlot[,"timepoint"]))
    
    ggplot(tempRNAiPlot, aes(x=timepoint, y=lethality, fill=factor(condition)))+
        geom_bar(stat="identity", color="black", position="dodge")+
        scale_fill_discrete(name="condition",#breaks=c(1,2),
                            labels=c("no inj", "par-6", "set1", "set2", "set3", 
                                     "set4", "set5", "set6", "set7", "set8", "set9"))+
        ylab("percent lethality")+
        theme_classic()
    
    ggsave("04_analysesOut/SuppFig4/4_RNAi_rrf3.pdf",
           width = 6, height = 3)

    rm(RNAiResults)
    
}

RNAi_rrf3()

