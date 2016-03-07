######
# 2015.12.31
#   All commands for transcriptional lineage paper

# import custom functions
source('~/Desktop/whoopsgraduateschool archive/20151231_transLinPape_functions.R')

# import data tables
source('~/Desktop/whoopsgraduateschool archive/20150715_makeCellMaps.R')


# downloaded all RPKM files to Desktop/whoopsgraduateschool archive/20150904_finalRPKMs/
# downloaded all alignent statistics to Desktop/whoopsgraduateschool archive/20150904_finalAlignmentStats/

#~#~#~#
#
#   Import and organize data

#
#  Load RPKMs and count files
#

orderOfAllSamples <- c("st411", "st413", "st415", "st441", "st449", "st451", "st265", "st266", "st301", "st302", "st311", "st312", "st361", "st362", "st409", "st410", "st193", "st194", "st195", "st196", "st305", "st306", "st307", "st308", "st333", "st334", "st335", "st336", "st363", "st364", "st365", "st366", "st405", "st406", "st407", "st408", "st225", "st226", "st227", "st228", "st229", "st230", "st231", "st232", "st337", "st338", "st339", "st340", "st341", "st342", "st343", "st344", "st369", "st370", "st371", "st372", "st373", "st374", "st375", "st376", "st393", "st394", "st395", "st396", "st397", "st398", "st399", "st400", "st417", "st418", "st419", "st420", "st421", "st422", "st423", "st424", "st481", "st482", "st483", "st484", "st485", "st486", "st487", "st488", "st249", "st250", "st251", "st252", "st253", "st254", "st255", "st256", "st257", "st258", "st259", "st260", "st261", "st262", "st263", "st264", "st313", "st314", "st315", "st316", "st317", "st318", "st319", "st320", "st321", "st322", "st323", "st324", "st325", "st327", "st328", "st345", "st346", "st347", "st348", "st349", "st350", "st351", "st352", "st353", "st354", "st355", "st356", "st357", "st358", "st359", "st360", "st377", "st378", "st379", "st380", "st381", "st382", "st383", "st385", "st386", "st387", "st388", "st389", "st390", "st391", "st392", "st425", "st426", "st427", "st428", "st429", "st430", "st431", "st432", "st433", "st434", "st435", "st436", "st437", "st438", "st439", "st440", "st527", "st528", "st529", "st530", "st531", "st532", "st533", "st535", "st536", "st537", "st538", "st539", "st540", "st541", "st542", "st465", "st466", "st467", "st468", "st469", "st470", "st471", "st472", "st473", "st474", "st475", "st476", "st477", "st478", "st479", "st505", "st506", "st507", "st508", "st509", "st510", "st513", "st514", "st515", "st516", "st517", "st518", "st519", "st520", "st554", "st555", "st556", "st557", "st558", "st561", "st562", "st563", "st564", "st565", "st566", "st567", "st568")

tlRPKMs <- feedNewRPKMs("Desktop/whoopsgraduateschool archive/20150904_finalRPKMs/", by="rpkm")
tlRPKMs <- tlRPKMs[,orderOfAllSamples]
tlCounts <- feedNewRPKMs("Desktop/whoopsgraduateschool archive/20150904_finalRPKMs/", by="count")
tlCounts <- tlCounts[,orderOfAllSamples]

#
#  Set up Design meta data table
#

# The beginning of tlDesign table is all input by hand
tlDesign <- v3Design[,1:4]
tlDesign[,"quality"] <- "Use"
tlDesign[,"ID"] <- as.character(tlDesign[,"ID"])

# import all the data from the alignment files
alignStats = data.frame()
for(file in list.files("Desktop/whoopsgraduateschool archive/20150904_finalAlignmentStats/")){
    tempIn <- read.table(paste("Desktop/whoopsgraduateschool archive/20150904_finalAlignmentStats/", file, sep=""), sep="\t", row.names="V1")
    alignStats = rbind(alignStats, tempIn)
}
alignStats

# push those values into the design file
for(cell in row.names(tlDesign)){
    tlDesign[cell, "uniqueCounts"] = round(alignStats[cell, "V4"])
    tlDesign[cell, "input"] = alignStats[cell, "V6"]
    tlDesign[cell, "ERCCcounts"] = round(alignStats[cell, "V5"])
}

thresh = 25

# calculate genesDetected
for(cell in rownames(tlDesign)){
    binLog <- (tlRPKMs[,cell]>thresh)
    tlDesign[cell,"genesDetected"] = length(which(binLog))
}

# calculate whole embryo genes detected
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

# calculate embryo input
for(embryo in unique(paste(tlDesign$stage, tlDesign$replicate, sep="."))){
    embSplit = strsplit(embryo, "\\.")[[1]]
    stage=embSplit[1]
    repl = embSplit[2]
    outNumb = 0
    celset <- intersect(colnames(tlCounts[,which(tlDesign$stage==stage)]), 
                        colnames(tlCounts[,which(tlDesign$replicate==repl)]))
    for(cell in celset){
        outNumb = outNumb + tlDesign[cell, "input"]
    }
    tlDesign[which(tlDesign$stage==stage & tlDesign$replicate==repl), "embryoInput"] = outNumb
}


#~#~#~#
#
#   Make Sulston founder cell plot

founders <- read.table(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/A_sulston/founder_cell_fates.txt", 
                    row.names = 1, header = T, sep="\t")
# (remove total line)
founders2 <- founders[1:dim(founders)[1]-1,]
# subscript p4
colnames(founders2)[which(colnames(founders2)=="P4")] <- expression(P[4])
founders2 <- rbind(founders2, "neurons" = apply(founders2, 2, function(x){sum(x["misc neurons"], x["ventral"], x["ring ganglia"])}))
founders2["head meso",] <- apply(founders2, 2, function(x){sum(x["head meso"], x["arcade"])})
founders2 <- founders2[which(!rownames(founders2)%in%c("ring ganglia", "ventral", "misc neurons", "arcade", "valve")),]

for (i in 1:16) {founders2[,i] = founders2[,i]/sum(founders2[,i])}

founders2.m <- melt(founders2)
cellFates <- row.names(founders2)
founders2.m <- cbind(founders2.m, rep(cellFates, 16))
founders2.m <- data.frame(founders2.m[founders2.m$value != 0,], row.names = NULL)
colnames(founders2.m) = c("Founder_Cell", "Proportion", "Cell_Fate")
founders2.m = transform (founders2.m,
                         Cell_Fate = factor(Cell_Fate,
                                            levels = c("apoptosis", "neurons", "germ line", "intestine", "hypodermis", "pharynx", "muscle", "distal tip", "excretory", "ceolomocyte", "head meso"),
                                            labels = c("Apoptosis", "Neurons", "Germ Line", "Intestine", "Epidermis", "Pharynx", "Muscle", "Gonad", "Excretory", "Ceolomocyte", "Head Mesoderm")
                         )
)

qual_11 <- c("#a6cee3", "#1f78b4", "#b2bf8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")

founderPlot <- ggplot(founders2.m, aes (x=Founder_Cell, y=Proportion, color= Cell_Fate, fill=Cell_Fate))+
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

ggsave(founderPlot, file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/A_sulston/20150831_founderCellSulston.pdf", 
       width=10, height = 4)

# Use ERCC counts to calculate the number of transcripts for each gene, and the pg detected in each cell
source("~/Desktop/whoopsgraduateschool archive/20150918_ERCCanalysis.R")


#~#~#~#
#
#   Make bubble party plots 

#
#   Whole embryo vs pg RNA
#

pdf(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/B_bubbleParty/01_summaryStatsAllReps_wholeEmbryo_pgRNA_v2.pdf", 
    width=8.5, height=4.5)
par(mar=c(5.1, 6.1, 4.1, 6.1))
embMean = mean(unique(log10(tlDesign[, "embryoPgExpected"]+1)))
embSD = sd(unique(log10(tlDesign[, "embryoPgExpected"]+1)))
#abline(h=c(embMean, embMean-embSD, embMean+embSD), lty=c(1,2,2), col=alpha("black", .4))
plot(-1,-1, xlim=c(1,31), ylim=c(-.05,1.8),
     xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
rect(-3,embMean-embSD,34,embMean+embSD, col="light gray", border = NA)
abline(h=embMean, col=alpha("black", .2))
points(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop=TRUE)), log10(tlDesign$embryoPgExpected+1), 
       pch=18, col=repColors[tlDesign$replicate], #cex=((tlDesign$embryoGenesDetected)^.6)/100
       cex=2
)
#title(main = "Picograms mRNA recovered from sampled embryos")
abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
axis(1, at=c(3.5,9,14,19.5,27), tick=FALSE, cex.axis=1.5,
     labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
axis(2, at=c(0,1.5), cex.axis=1.5)
title(ylab = expression(paste("log"[10], "(pg RNA)", sep="")), cex.lab = 1.5, line = 2.5)
title(xlab = "Embryos (aggregated from single-cell data)", cex.lab=1.5)
#legend(legend = c("5,000 genes", "10,000 genes"), title = "Number of Genes \nDetected above 25 RPKM\n in whole embryos",
#       x = par('usr')[2]+1, y = par('usr')[4]-1.25,
#       pch=c(18, 18), pt.cex=c(5000,10000)^.6/100, cex=.6, 
#       y.intersp = 2, x.intersp=1.5, bty="n", xpd=TRUE)
dev.off()


#
#   Change quality scores for outlying embryos, based on this plot
#

tlDesign[which(tlDesign$stage=="1-cell" & tlDesign$replicate==3),"quality"] <- "Toss"
tlDesign[which(tlDesign$stage=="8-cell" & tlDesign$replicate==5),"quality"] <- "Toss"
tlDesign[which(tlDesign$stage=="16-cell" & tlDesign$replicate%in%c(4,5,9)),"quality"] <- "Toss"

#
#   Whole embryo vs number of genes
#

pdf(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/B_bubbleParty/02_summaryStatsAllReps_wholeEmbryo_genesDetected_v2_25RPKM.pdf", 
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
#title(main = "Number of Genes Detected Above 10 RPKM in each embryo")
axis(1, at=c(3.5,9,14,19.5,27),tick=FALSE, cex.axis=1.5,
     labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
axis(2, cex.axis=1.5, at=c(3.6,3.7, 3.8, 3.9)
     )
abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
title(ylab = expression(paste("  log"[10], "(genes detected)", sep="")), cex.lab = 1.5, line = 2.5)
title(xlab = "Embryos (aggregated from single-cell data)", cex.lab=1.5)
dev.off()


#
#   Single cells vs RPKMs
#

pdf(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/B_bubbleParty/03_summaryStatsAllReps_singleCell_RPKMs_v2_25RPKM.pdf", 
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

# Is total RPKM increasing in later embryos because there are larger numbers of 
#   lower complexity libraries being combined?
# ID the number of genes found in all the replicates for each stage, 
#   and then ID what percentage of that total is found in each embryo
#   If 
tlDesign$percRepUnique <- -20
thresh <- 25
for(stage in unique(tlDesign$stage)){
    print(stage)
    qualReps <- unique(tlDesign[which(tlDesign$stage==stage 
                                      & tlDesign$quality=="Use"),"replicate"])
    for(rep in 1:length(qualReps)){
        # select 5 replicates
        useReps <- (((rep-1):(rep+3))%%length(qualReps))+1
        useReps <- qualReps[useReps]
        subMatr <- tlRPKMs[,rownames(tlDesign[which(tlDesign$stage==stage 
                                                     & tlDesign$quality=="Use"
                                                     & tlDesign$replicate %in% useReps[2:length(useReps)]),])]
        # come up with a number of how many genes are detected over threshhold
        stageGenesNames <- unique(rownames(which(subMatr > thresh, arr.ind = T)))
        # for this particular rep, how many genes are detected?
        subSubMatr <- tlRPKMs[,rownames(tlDesign[which(tlDesign$stage==stage 
                                                    & tlDesign$quality=="Use"
                                                    & tlDesign$replicate==useReps[1]),])]
        repGenesNum <- 0
        if(stage=="1-cell"){
            repGenesNames <- names(which(subSubMatr > thresh))
        } else {
            repGenesNames <- unique(rownames(which(subSubMatr > thresh, arr.ind = T)))
        }
        # input that percentage for each cell of that replicate
        tempPerc <- (length(setdiff(repGenesNames, stageGenesNames))/length(repGenesNames))*100
        tlDesign[which(tlDesign$stage==stage 
                       & tlDesign$replicate==useReps[1] 
                       & tlDesign$quality=="Use"),"percRepUnique"] <- tempPerc
    }
}


#
#   Whole embryo reproducibility
#

pdf(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/B_bubbleParty/04_summaryStatsAllReps_wholeEmbryo_reproducibility_byRep.pdf", 
    width=8.5, height=4.5)
par(mar=c(5.1, 6.1, 4.1, 6.1))
plot(-1,-1, xlim=c(1,31), ylim=c(-5,105),
     xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
tempColors <- repColors[tlDesign$replicate]
tempColors[which(tlDesign$quality=="Toss")] <- "gray"
points(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop=TRUE)), 
       tlDesign$percRepUnique, 
       pch=18, col=tempColors, cex=2
)
#title(main = "Number of Genes Detected Above 10 RPKM in each embryo")
axis(1, at=c(3.5,9,14,19.5,27),tick=FALSE, cex.axis=1.5,
     labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
axis(2, cex.axis=1.5, at=c(0,20,40,60,80,100))
abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
title(ylab = "Percentage of genes detected that are\nreplicate-specific compared to stage", cex.lab = 1.5, line = 2.5)
title(xlab = "Embryos (aggregated from single-cell data)", cex.lab=1.5)
dev.off()



tlDesign$percCellUnique <- -20
thresh <- 25
for(stage in c("2-cell", "4-cell", "8-cell", "16-cell")){
    print(stage)
    goodReps <- unique(tlDesign[which(tlDesign$stage==stage 
                                      & tlDesign$quality=="Use"),"replicate"])
    for(repl in goodReps){
        allCells <- rownames(tlDesign[which(tlDesign$stage==stage &
                                                tlDesign$replicate==repl),])
        for(tempCell in allCells){
            # find genes enriched in this cell
            cellGenes <- names(which(tlRPKMs[,tempCell] > thresh))
            # find genes enriched in other cells
            subMatr <- tlRPKMs[,allCells[-which(allCells==tempCell)]]
            embryoGenes <- NULL
            if(stage=="2-cell"){
                embryoGenes <- names(which(subMatr > thresh))
            } else {
                embryoGenes <- unique(rownames(which(subMatr > thresh, arr.ind = T)))
            }
            # come up with percentage unique
            tempPerc <- (length(setdiff(cellGenes, embryoGenes))/length(cellGenes))*100
            #put that in the design dataframe
            tlDesign[tempCell, "percCellUnique"] <- tempPerc
        }
    }
}


#
#   Single cells reproducibility
#

pdf(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/B_bubbleParty/05_summaryStatsAllReps_singleCell_reproducibility_25RPKM.pdf", 
    width=8.5, height=4.5)
par(mar=c(5.1, 6.1, 4.1, 6.1))
library('scales')
# set up color scheme for cells that got dropped
cellCol <- repColors[tlDesign$replicate]
cellCol[c(which(tlDesign$stage=="1-cell" & tlDesign$replicate==3),
          which(tlDesign$stage=="8-cell" & tlDesign$replicate==5),
          which(tlDesign$stage=="16-cell" & tlDesign$replicate%in%c(4,5,9)))] <- "light gray"
plot(jitter(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop = TRUE)), factor = .7), 
     tlDesign$percCellUnique, 
     pch=16, col=alpha(cellCol, .6), xlim=c(1,31), ylim=c(0,100),
     xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
axis(1, at=c(3.5,9,14,19.5,27), tick=FALSE, cex.axis=1.5,
     labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
axis(2, cex.axis=1.5)
title(ylab = "Percentage of genes detected that are\ncell-specific compared to embryo", cex.lab = 1.5, line = 2.5)
title(xlab = "Individual Cells", cex.lab=1.5)
dev.off()

    

#
# A second attempt at this reproducibility business:
#   Each score consists of an average of every one-to-one comparison possible
#

tlDesign$AvgdRepUnique <- -20
thresh <- 25
for(stage in unique(tlDesign$stage)){
    print(stage)
    qualReps <- unique(tlDesign[which(tlDesign$stage==stage 
                                      & tlDesign$quality=="Use"),"replicate"])
    for(rep in 1:length(qualReps)){
        # select 5 replicates
        thisRep <- qualReps[rep]
        otherReps <- qualReps[-which(qualReps==thisRep)]
        scores <- c()
        thisGeneNames <- NULL
        thisMatr <- tlRPKMs[,rownames(tlDesign[which(tlDesign$stage==stage
                                                     & tlDesign$replicate==thisRep
                                                     & tlDesign$quality=="Use"),])]
        if(stage=="1-cell"){
            thisGeneNames <- names(which(thisMatr > thresh))
        }else{
            thisGeneNames <- unique(rownames(which(thisMatr > thresh, arr.ind = T)))
        }
        for(othRep in otherReps){
            #make a subMatrix
            othMatr <- tlRPKMs[,rownames(tlDesign[which(tlDesign$stage==stage
                                                        & tlDesign$replicate==othRep
                                                        & tlDesign$quality=="Use"),])]
            # find the genes that are detected in that subMatrs
            othGeneNames <- NULL
            if(stage=="1-cell"){
                # 1 cell
                othGeneNames <- names(which(othMatr > thresh))
            }else{
                # others
                othGeneNames <- unique(rownames(which(othMatr > thresh, arr.ind = T)))
            }
            #compare thisMatr to othMatr
            tempPerc <- (length(setdiff(thisGeneNames, othGeneNames))/length(thisGeneNames))*100
            #push value to scores vector
            scores <- c(scores, tempPerc)
        }
        # average scores vector
        avgScores <- mean(scores)
        #assign the average to all the cells of that replicate embryo
        tlDesign[which(tlDesign$stage==stage & tlDesign$replicate==thisRep),"AvgdRepUnique"] <- avgScores
    }
}



#
#   Whole embryo reproducibility
#

pdf(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/B_bubbleParty/06_summaryStatsAllReps_wholeEmbryo_avgdReproducibility_byRep.pdf", 
    width=8.5, height=4.5)
par(mar=c(5.1, 6.1, 4.1, 6.1))
plot(-1,-1, xlim=c(1,31), ylim=c(-5,105),
     xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
tempColors <- repColors[tlDesign$replicate]
tempColors[which(tlDesign$quality=="Toss")] <- "gray"
points(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop=TRUE)), 
       tlDesign$AvgdRepUnique, 
       pch=18, col=tempColors, cex=2
)
#title(main = "Number of Genes Detected Above 10 RPKM in each embryo")
axis(1, at=c(3.5,9,14,19.5,27),tick=FALSE, cex.axis=1.5,
     labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
axis(2, cex.axis=1.5, at=c(0,20,40,60,80,100))
abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
title(ylab = "Percentage of genes detected that are\nreplicate-specific compared to stage", cex.lab = 1.5, line = 2.5)
title(xlab = "Embryos (aggregated from single-cell data)", cex.lab=1.5)
dev.off()


tlDesign$avgdCellUnique <- -20
thresh <- 25
for(stage in c("2-cell", "4-cell", "8-cell", "16-cell")){
    print(stage)
    goodReps <- unique(tlDesign[which(tlDesign$stage==stage 
                                      & tlDesign$quality=="Use"),"replicate"])
    for(repl in goodReps){
        allCells <- rownames(tlDesign[which(tlDesign$stage==stage &
                                                tlDesign$replicate==repl),])
        for(tempCell in allCells){
            # find genes enriched in this cell
            cellGenes <- names(which(tlRPKMs[,tempCell] > thresh))
            scores <- c()
            # find genes enriched in other cells
            for(othCell in allCells[-which(allCells==tempCell)]){
                # find genes enriched in the other cell
                othGenes <- names(which(tlRPKMs[,othCell] > thresh))
                # calculate a percentage of unique genes
                tempPerc <- (length(setdiff(cellGenes, othGenes))/length(cellGenes))*100
                # push that score to the list
                scores <- c(scores, tempPerc)
            }
            # make an average score
            avgScores <- mean(scores)
            # push that average score to the design file
            tlDesign[tempCell, "avgdCellUnique"] <- avgScores
        }
    }
}


#
#   Single cells reproducibility
#

pdf(file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig1_sulstonAndQC/B_bubbleParty/07_summaryStatsAllReps_ingleCell_avgdReproducibility_25RPKM.pdf", 
    width=8.5, height=4.5)
par(mar=c(5.1, 6.1, 4.1, 6.1))
library('scales')
# set up color scheme for cells that got dropped
cellCol <- repColors[tlDesign$replicate]
cellCol[c(which(tlDesign$stage=="1-cell" & tlDesign$replicate==3),
          which(tlDesign$stage=="8-cell" & tlDesign$replicate==5),
          which(tlDesign$stage=="16-cell" & tlDesign$replicate%in%c(4,5,9)))] <- "light gray"
plot(jitter(as.numeric(interaction(tlDesign$replicate, tlDesign$stage, drop = TRUE)), factor = .7), 
     tlDesign$avgdCellUnique, 
     pch=16, col=alpha(cellCol, .6), xlim=c(1,31), ylim=c(0,100),
     xaxt = "n", yaxt = "n", ann=FALSE, oma = c(5.1, 4.1, 4.1, 8.1))
abline(v = c(6.5, 11.5, 16.5, 22.5), lty=3)
axis(1, at=c(3.5,9,14,19.5,27), tick=FALSE, cex.axis=1.5,
     labels=c("1 cell", "2 cell", "4 cell", "8 cell", "16 cell"))
axis(2, cex.axis=1.5)
title(ylab = "Percentage of genes detected that are\ncell-specific compared to embryo", cex.lab = 1.5, line = 2.5)
title(xlab = "Individual Cells", cex.lab=1.5)
dev.off()

# this looks pretty interesting actually ... think about it tomorrow



#~#~#~#
#
#  PCAs to indentify cell types in replicates

# TWO CELL
unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="2-cell" & tlDesign$quality=="Use")], falseNeg = 0,
                   outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/01_2cell_43genes.pdf")
tlDesign[c("st265", "st301", "st311", "st361", "st409"),"ID"] <- "AB"
tlDesign[c("st266", "st302", "st312", "st362", "st410"),"ID"] <- "P1"

# FOUR CELL
unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="4-cell" & tlDesign$quality=="Use")], falseNeg = 0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/02_4cell_150genes.pdf")
tlDesign[c("st196", "st308", "st336", "st366", "st408"),"ID"] <- "P2"
tlDesign[c("st193", "st305", "st335", "st364", "st406"),"ID"] <- "EMS"

unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="4-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P2", "EMS")))],  falseNeg = 0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/03_4cell_40genes.pdf")
tlDesign[c("st194", "st306", "st334", "st363", "st407"),"ID"] <- "ABa"
tlDesign[c("st195", "st307", "st333", "st365", "st405"),"ID"] <- "ABp"

# EIGHT CELL
unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use")], falseNeg = 0,
           outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/04_8cell_453genes.pdf")
tlDesign[c("st232", "st344", "st376", "st400", "st488"),"ID"] <- "P3"
tlDesign[c("st230", "st341", "st373", "st397", "st487"),"ID"] <- "C"

unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P3", "C")))],  falseNeg = 0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/05_8cell_295genes.pdf")
tlDesign[c("st231", "st343", "st375", "st398", "st486"),"ID"] <- "E"
tlDesign[c("st229", "st342", "st374", "st399", "st485"),"ID"] <- "MS"

unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="8-cell" & tlDesign$quality=="Use" &(!tlDesign$ID %in% c("P3", "C", "E", "MS")))],  falseNeg = 0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/06_8cell_169genes.pdf")

# sorted ABxx's by hand. hlh-27 and ref-1 took care of the two ABpx's. 
# of the rest, c("st226", "st227", "st337", "st338", "st370", "st372", "st393", "st396", "st482", "st484")
    # some have tbx-38 (c("st226", "st370", "st337")) and their matches dont (c("st227", "st372", "st338"))

# generate list of top informative genes between these sets
ERdiffExV1(matr=tlCounts[,c("st226", "st227", "st370", "st372")], plotWidth = 5, plotHeight = 5,
           to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/05_tbx38MAplot.pdf",  
           sig=.05, group = factor(rep(c("tbx+", "tbx-"), 2), levels = c("tbx+", "tbx-")))

temp <- ERdiffExV1(matr=tlCounts[,c("st226", "st227", "st370", "st372")], 
                   to = "list",  sig=.05,
           group = factor(rep(c("tbx+", "tbx-"), 2), levels = c("tbx+", "tbx-")))

temp <- c(rownames(temp[[1]]), rownames(temp[[2]]))

unsupervisedPCA(tlRPKMs[,c("st226", "st227", "st337", "st338", "st370", "st372", "st393", "st396", "st482", "st484")],  
                    falseNeg = 6, forceGenes = temp, labelSamples = T,
                    outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/06_tbx38BasedPCA_119genes_lab.pdf")

unsupervisedPCA(tlRPKMs[,c("st226", "st227", "st337", "st338", "st370", "st372", "st393", "st396", "st482", "st484")],  
                falseNeg = 6, forceGenes = temp, 
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/06_tbx38BasedPCA_119genes.pdf")

    

# According to smFISH on notch pathway genes, ABal and ABar have no hlh-27, whereas ABpl and ABpr do.
#   ref-1 should be expressed in only one of those; ABpl.
tlDesign[c("st225", "st340", "st371", "st394", "st483"),"ID"] <- "ABpl"
tlDesign[c("st228", "st339", "st369", "st395", "st481"),"ID"]  <- "ABpr"
#   and tbx-38 should be expressed in only one from the other set; ABal
tlDesign[c("st226", "st337", "st370", "st396", "st484"),"ID"] <- "ABal"
tlDesign[c("st227", "st338", "st372", "st393", "st482"),"ID"]  <- "ABar"


# perhaps it's possible to ratchet back from ABal+ABar vs ABpl+ABpr here, to ABa vs ABp in the four cell stage...
ABaxVsABpxGenes <- ERdiffExV1(matr=tlCounts[,typeCells(c("ABal", "ABar", "ABpl", "ABpr"))], 
           group = factor(c(rep("ABax", length(typeCells(c("ABal", "ABar")))), 
                            rep("ABpx", length(typeCells(c("ABpl", "ABpr"))))), levels = c("ABax", "ABpx")), 
           postThresh = 25, to = "list", sig=.5)
str(ABaxVsABpxGenes)
# List of 2
# $ ABax:'data.frame':	1054 obs. of  3 variables:
# $ ABpx:'data.frame':	1665 obs. of  3 variables:
ABaVsABpGenes <- ERdiffExV1(matr=tlCounts[,typeCells(c("ABa", "ABp"))], 
                            postThresh = 25, to = "list", sig=.5)
str(ABaVsABpGenes)
# List of 2
# $ ABa:'data.frame':	873 obs. of  3 variables:
# $ ABp:'data.frame':	857 obs. of  3 variables:
    
length(intersect(rownames(ABaxVsABpxGenes[["ABax"]]), rownames(ABaVsABpGenes[["ABa"]]))) # 126
length(union(rownames(ABaxVsABpxGenes[["ABax"]]), rownames(ABaVsABpGenes[["ABa"]]))) # 126
(126/1801)*100 # 6.996

length(intersect(rownames(ABaxVsABpxGenes[["ABax"]]), rownames(ABaVsABpGenes[["ABp"]]))) # 109
length(union(rownames(ABaxVsABpxGenes[["ABax"]]), rownames(ABaVsABpGenes[["ABp"]]))) # 109
(109/1802)*100 # 6.049

length(intersect(rownames(ABaxVsABpxGenes[["ABpx"]]), rownames(ABaVsABpGenes[["ABa"]]))) # 210
length(union(rownames(ABaxVsABpxGenes[["ABpx"]]), rownames(ABaVsABpGenes[["ABa"]]))) # 210
(210/2328)*100 # 9.021

length(intersect(rownames(ABaxVsABpxGenes[["ABpx"]]), rownames(ABaVsABpGenes[["ABp"]]))) # 218
length(union(rownames(ABaxVsABpxGenes[["ABpx"]]), rownames(ABaVsABpGenes[["ABp"]]))) # 218
(218/2304)*100 # 9.462

# yea, by the tiniest bit, that looks like it supports my hypothesis.
# this result remains:

#tlDesign[c("st194", "st306", "st334", "st363", "st407"),"ID"] <- "ABa"
#tlDesign[c("st195", "st307", "st333", "st365", "st405"),"ID"] <- "ABp"

# make a heatmap of these overlaps:
heatmap.2(matrix(c(126,109,210,218), nrow = 2), )

pdf("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/07_ABapABaxpx_heatmap.pdf", 5, 5)
heatmap.2(matrix(c(6.996,6.049,9.021,9.462), nrow = 2), trace = "none", Rowv=F, Colv=F, 
          labRow = c("ABa", "ABp"), labCol = c("ABpx", "ABax"), density.info = "none")
dev.off()


# SIXTEEN CELL
unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use")], topN=50, falseNeg=0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/07_16cell_50topN_518genes.pdf")
tlDesign[c("st263", "st264", "st327", "st328", "st359", "st360", "st532", "st533", "st478", "st479", "st519", "st520"),"ID"] #<- "MS"

unsupervisedPCA(tlRPKMs[,c("st263", "st264", "st327", "st328", "st359", "st360", "st532", "st533", "st478", "st479", "st519", "st520")], falseNeg = 0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/08_16cell_73genes.pdf")
tlDesign[c("st263", "st328", "st360", "st532", "st478", "st519"), "ID"] <- "P4"
tlDesign[c("st264", "st327", "st359", "st533", "st479", "st520"), "ID"] <- "D"
# putting this in the supplement
unsupervisedPCA(tlRPKMs[,c("st263", "st264", "st327", "st328", "st359", "st360", "st532", "st533", "st478", "st479", "st519", "st520")], falseNeg = 0,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/01_DP4_73genes.pdf")
unsupervisedPCA(tlRPKMs[,c("st263", "st264", "st327", "st328", "st359", "st360", "st532", "st533", "st478", "st479", "st519", "st520")], falseNeg = 0, labelSamples = T,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/01_DP4_73genes_lab.pdf")


# Ex's with singlets
unsupervisedPCA(tlRPKMs[,c("st325", "st260", "st261", "st355", "st356", "st527", "st473", "st475", "st514", "st516")], 
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/09_16cell_98genes.pdf")
# the pairing is ambiguous, because of the singlets
# Ex's, only doublets
unsupervisedPCA(tlRPKMs[,c("st260", "st261", "st355", "st356", "st473", "st475", "st514", "st516")], falseNeg = 0, 
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/10_16cell_34genes.pdf")
# gonna stick this in the supplement
unsupervisedPCA(tlRPKMs[,c("st260", "st261", "st355", "st356", "st473", "st475", "st514", "st516")], falseNeg = 0, labelSamples = T,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/03_Ex_34genes_lab.pdf")
unsupervisedPCA(tlRPKMs[,c("st260", "st261", "st355", "st356", "st473", "st475", "st514", "st516")], falseNeg = 0, 
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/03_Ex_34genes.pdf")
# difference between those sets:
Ex1Ex2Genes <- ERdiffExV1(tlRPKMs[,c("st261", "st355", "st473", "st516", "st260", "st356", "st475", "st514")], 
           group = factor(c(rep("Ep", 4), rep("Ea", 4))), label = rownames(tlRPKMs), sig=.05, to="list")
# Ex's with singlets, using genes from Ex doublets
unsupervisedPCA(tlRPKMs[,c("st325", "st260", "st261", "st355", "st356", "st527", "st473", "st475", "st514", "st516")], 
            forceGenes = union(rownames(Ex1Ex2Genes[[1]]), rownames(Ex1Ex2Genes[[2]])), falseNeg = 0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/11_16cell_69genes.pdf")
tlDesign[c("st261", "st325", "st355", "st527", "st473", "st516"), "ID"] <- "Ep"
tlDesign[c("st260", "st356", "st475", "st514"), "ID"] <- "Ea"

unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" 
                           & tlDesign$quality=="Use" 
                           &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea")) 
                        #   &(!rownames(tlDesign)%in%c("st315", "st353"))
                        )], topN=50, falseNeg=0,
            outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/12_16cell_50topN_370genes.pdf")

unsupervisedPCA(tlRPKMs[,c("st257", "st258", "st323", "st324", "st354", "st528", "st529", "st474", "st476", "st513", "st517")], 
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/13_16cell_72genes.pdf")
tlDesign[c("st257", "st323", "st529", "st476", "st513"), "ID"] <- "MSx1"                        
tlDesign[c("st258", "st324", "st354", "st528", "st474", "st517"), "ID"] <- "MSx2"                        

unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" 
                               & tlDesign$quality=="Use" 
                               &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2")) 
                               #   &(!rownames(tlDesign)%in%c("st315", "st353"))
                               )], topN=50,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/14_16cell_321genes.pdf")
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
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/15_16cell_499genes.pdf")

# ID Cx1 and Cx2 for cells in doublets
unsupervisedPCA(tlRPKMs[,c("st259", "st262", "st321", "st322", "st357", "st358", "st530", "st531",  "st515", "st518")], 
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/16_16cell_28genes.pdf")
c("st262", "st322", "st357", "st531", "st515")
c("st259", "st321", "st358", "st530", "st518")
# ID genes separating Cx1 and Cx2 doublets
Cx1Cx2Genes <- ERdiffExV1(tlRPKMs[,c("st262", "st322", "st357", "st531", "st515", "st259", "st321", "st358", "st530", "st518")], 
                            group = factor(c(rep("Cx1", 5), rep("Cx2", 5))), sig = .05,
                            to="list")
unsupervisedPCA(tlRPKMs[,c("st259", "st262", "st321", "st322", "st357", "st358", "st530", "st531",  "st477", "st515", "st518")], 
                forceGenes = union(rownames(Cx1Cx2Genes[[1]]), rownames(Cx1Cx2Genes[[2]])), falseNeg = 2,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/17_16cell_79genes.pdf")
tlDesign[c("st262", "st322", "st357", "st531", "st515"), "ID"] <- "Cx1"
tlDesign[c("st259", "st321", "st358", "st530", "st477", "st518"), "ID"] <- "Cx2"
# Looking over all replicates to make sure all P1 descendants are accounted for
tlDesign["st353", "ID"] <- "MSx1"
unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use" 
                               &(!tlDesign$ID %in% c("P4", "D", "Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2")))], 
outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/18_16cell_557genes.pdf")


# Use candidate genes to cluster AB descendants

# can tell by hlh rpkms that ABxx1 and ABxx2 are a group
# can tell by ref-1 that ABxx1 and ABxx2 are separate
tlDesign[c("st225", "st340", "st371", "st394", "st483"), "ID"] <- "ABxx1"
tlDesign[c("st228", "st339", "st369", "st395", "st481"), "ID"] <- "ABxx2"
# tbx-38 sets three (st226, st337, st370) of the other ABxx's apart. Use those to ratchet the others to one side or the other
ABxx3ABxx4Genes <- ERdiffExV1(tlRPKMs[,c("st226", "st337", "st370", "st227", "st338", "st372")], 
                          group = factor(c(rep("ABxx3", 3), rep("ABxx4", 3))), sig = .05,
                          to="list")
unsupervisedPCA(tlRPKMs[,c("st226", "st337", "st370", "st393", "st482", "st227", "st338", "st372", "st396", "st484")], 
                forceGenes = union(rownames(ABxx3ABxx4Genes[[1]]), rownames(ABxx3ABxx4Genes[[2]])), falseNeg = 3,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/19_8cell_genes.pdf")
tlDesign[c("st226", "st337", "st370", "st396", "st484"), "ID"] <- "ABxx3"
tlDesign[c("st227", "st338", "st372", "st393", "st482"), "ID"] <- "ABxx4"

unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use"
                               & !tlDesign$ID%in%c("P4", "D", "Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2"))], 
                forceGenes = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", 
                               "tbx-38", "tbx-39", "ref-1"), falseNeg = 3,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/20_16cell_notchGenes.pdf")
# that makes 2 clusters:
# not having hlh's -
# c("st253", "st254", "st255", "st256", "st316", "st320", "st314", "st319", 
#   "st346", "st348", "st345", "st351", "st535", "st538", "st536", "st539", 
#   "st466", "st470", "st467", "st469", "st505", "st509")
# having hlh's - 
# c("st249", "st252", "st250", "st251", "st317", "st313", "st315", "st318", 
#   "st347", "st349", "st350", "st352", "st541", "st542", "st537", "st540", 
#   "st468", "st471", "st465", "st472", "st506", "st507", "st508", "st510")
unsupervisedPCA(tlRPKMs[,c("st253", "st254", "st255", "st256", "st316", "st320", "st314", "st319", 
                           "st346", "st348", "st345", "st351", "st535", "st538", "st536", "st539", 
                           "st466", "st470", "st467", "st469" #, "st505", "st509"
                           )],  falseNeg = 2,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/21_16cell_noHlh_134genes.pdf")
unsupervisedPCA(tlRPKMs[,c("st249", "st252", "st250", "st251", "st317", "st313", "st315", "st318", 
                           "st347", "st349", "st350", "st352", "st541", "st542", "st537", "st540", 
                           "st468", "st471", "st465", "st472" #, "st506", "st507", "st508", "st510"
                           )], #falseNeg = 2,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/PCA_images/22_16cell_wHlh_139genes.pdf")

tlDesign[c("st255", "st256", "st314", "st319", "st345", 
           "st351", "st536", "st539", "st467", "st469"), c("ID", "safeID")] <- "ABalx"

tlDesign[c("st250", "st251", "st313", "st318", "st350", "st352", 
           "st537", "st540", "st465", "st472", "st505", "st509"), c("ID", "safeID")] <- "ABprx"

tlDesign[c("st253", "st254", "st316", "st320", "st346", "st348", 
           "st535", "st538", "st466", "st470", "st506", "st508"), c("ID", "safeID")] <- "ABarx"

tlDesign[c("st249", "st252", "st317", "st347", "st349", 
           "st541", "st542", "st468", "st471", "st507", "st510"), c("ID", "safeID")] <- "ABplx"

# putting those in the supps:
unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use" & 
                                   !tlDesign$ID%in%c("P4", "D", "Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2"))], 
                forceGenes = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", "tbx-38", "tbx-39", "ref-1"), falseNeg = 6, labelSamples = T,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/08_16cell_notchPCA1_lab.pdf")
unsupervisedPCA(tlRPKMs[,which(tlDesign$stage=="16-cell" & tlDesign$quality=="Use" & 
                                   !tlDesign$ID%in%c("P4", "D", "Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2"))], 
                forceGenes = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", "tbx-38", "tbx-39", "ref-1"), falseNeg = 6,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/08_16cell_notchPCA1.pdf")

# then split that bolus into two halves 
   # # # on the right:
unsupervisedPCA(tlRPKMs[,c("st253", "st254", "st255", "st256", "st316", "st320", "st314", "st319", 
                           "st346", "st348", "st345", "st351", "st535", "st538", "st536", "st539", 
                           "st466", "st470", "st467", "st469" , "st505", "st509")],  falseNeg = 2, labelSamples = T,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/09_rightHalfNotch_labs.pdf")

unsupervisedPCA(tlRPKMs[,c("st253", "st254", "st255", "st256", "st316", "st320", "st314", "st319", 
                           "st346", "st348", "st345", "st351", "st535", "st538", "st536", "st539", 
                           "st466", "st470", "st467", "st469" , "st505", "st509")],  falseNeg = 2,
                outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/09_rightHalfNotch.pdf")

    # # # on the left:
unsupervisedPCA(tlRPKMs[,c("st249", "st252", "st250", "st251", "st317", "st313", #"st315", 
                           "st318", "st347", "st349", "st350", "st352", "st541", "st542", "st537", "st540", 
                           "st468", "st471", "st465", "st472" , "st506", "st507", "st508", "st510"
)], forceGenes = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", "tbx-38", "tbx-39", "ref-1"), falseNeg = 6, 
labelSamples = T, #falseNeg = 2,
outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/10_leftHalfNotch_labs.pdf")

unsupervisedPCA(tlRPKMs[,c("st249", "st252", "st250", "st251", "st317", "st313", "st315", 
                           "st318", "st347", "st349", "st350", "st352", "st541", "st542", "st537", "st540", 
                           "st468", "st471", "st465", "st472" , "st506", "st507", "st508", "st510"
)], forceGenes = c("hlh-27","ref-1"), falseNeg = 6, 
labelSamples = T, #falseNeg = 2,
outPDF="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/10_leftHalfNotch_labs.pdf")

ABxxxPCA(geneList = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", "tbx-38", "tbx-39", "ref-1"), 
         outPDF = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/11_16cell_notch_minPCA_lab.pdf",
         labelSamps = T)
ABxxxPCA(geneList = c("hlh-25", "hlh-26", "hlh-27", "hlh-29", "tbx-38", "tbx-39", "ref-1"), 
         outPDF = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS1_IDingCells/11_16cell_notch_minPCA.pdf",
         width = 5, height = 5, splitReps = T)



#
#   Make MA plots that show where the handful of candidate genes are showing up
#

# ABa vs ABp
ERdiffExV1(matr=tlCounts[,typeCells(c("ABa", "ABp"))], 
           overlay = c("ref-1", "daf-4", "vab-2"), 
           label=c("ref-1", "daf-4", "vab-2"), 
           plotWidth = 5, plotHeight = 5,
           to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/newCandidateMAplots/ABaABp.pdf")

# MSx1 vs MSx2
ERdiffExV1(matr=tlCounts[,typeCells(c("MSx1", "MSx2"))], 
           overlay = "C43H6.4", 
           label="C43H6.4", 
           plotWidth = 5, plotHeight = 5,
           to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/newCandidateMAplots/MSx1MSx2.pdf")

# Cx1 vs Cx2
ERdiffExV1(matr=tlCounts[,typeCells(c("Cx1", "Cx2"))], 
           overlay = "par-5", 
           label="par-5", 
           plotWidth = 5, plotHeight = 5,
           to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig2_PCAs/newCandidateMAplots/Cx1Cx2.pdf")

# 
#   Add and extra column to the design matrix which offers more conservative ID assignments
#

safeID = tlDesign[,"ID"]
temp <- cbind(tlDesign[,1:4], safeID, tlDesign[,5:ncol(tlDesign)])
temp$safeID <- as.character(temp$safeID)
temp[which(temp$safeID%in%c("Ea", "Ep")), "safeID"] <- "Ex"
temp[which(temp$safeID%in%c("Cx1", "Cx2")), "safeID"] <- "Cx"
temp[which(temp$safeID%in%c("MSx1", "MSx2")), "safeID"] <- "MSx"
tlDesign <- temp




#
#  Set up average and median RPKMs matrices
#


cellTypes <- sort(factor(unique(tlDesign$ID), 
                    levels=c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2",
                             "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                             "ABxxx1", "ABalx", "ABala", "ABalp", 
                             "ABxxx2", "ABarx", "ABara", "ABarp", 
                             "ABxxx3", "ABplx", "ABpla", "ABplp", 
                             "ABxxx4", "ABprx", "ABpra", "ABprp", 
                             "MSx", "MSx1", "MSx2", "MSa", "MSp",
                             "Ex", "Ex1", "Ex2", "Ea", "Ep",
                             "Cx", "Cx1", "Cx2", "Ca", "Cp",
                             "D", "P4")))
cellTypes <- as.character(cellTypes[!cellTypes %in% c("tossed", "ABxxx", NA)])

tlCellAvg <- matrix(0, nrow=nrow(tlRPKMs), ncol=length(cellTypes), dimnames = list(rownames(tlRPKMs), cellTypes))
for(i in 1:dim(tlCellAvg)[1]){
    print(i)
    for(cellType in cellTypes){
        tlCellAvg[i, cellType] <- mean(tlRPKMs[i, typeCells(cellType)])
    }
}

tlCellMed <- matrix(0, nrow=nrow(tlRPKMs), ncol=length(cellTypes), dimnames = list(rownames(tlRPKMs), cellTypes))
for(i in 1:dim(tlCellMed)[1]){
    print(i)
    for(cellType in cellTypes){
        tlCellMed[i, cellType] <- median(tlRPKMs[i, typeCells(cellType)])
    }
}

safeCellTypes <- sort(factor(unique(tlDesign$safeID), 
                             levels=c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2",
                                      "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                                      "ABxxx1", "ABalx", "ABala", "ABalp", 
                                      "ABxxx2", "ABarx", "ABara", "ABarp", 
                                      "ABxxx3", "ABplx", "ABpla", "ABplp", 
                                      "ABxxx4", "ABprx", "ABpra", "ABprp", 
                                      "MSx", "MSx1", "MSx2", "MSa", "MSp",
                                      "Ex", "Ex1", "Ex2", "Ea", "Ep",
                                      "Cx", "Cx1", "Cx2", "Ca", "Cp",
                                      "D", "P4")))
safeCellTypes <- as.character(safeCellTypes[!safeCellTypes %in% c("tossed", "ABxxx", NA)])

tlSafeCellAvg <- matrix(0, nrow=nrow(tlRPKMs), ncol=length(safeCellTypes), 
                        dimnames = list(rownames(tlRPKMs), safeCellTypes))
for(i in 1:dim(tlSafeCellAvg)[1]){
    print(i)
    for(cellType in safeCellTypes){
        tlSafeCellAvg[i, cellType] <- mean(tlRPKMs[i, typeCells(cellType)])
    }
}

tlSafeCellMed <- matrix(0, nrow=nrow(tlRPKMs), ncol=length(safeCellTypes), 
                        dimnames = list(rownames(tlRPKMs), safeCellTypes))
for(i in 1:dim(tlSafeCellMed)[1]){
    print(i)
    for(cellType in safeCellTypes){
        tlSafeCellMed[i, cellType] <- median(tlRPKMs[i, typeCells(cellType)])
    }
}

# How many genes do we have expression patterns for?
length(unique(rownames(which(tlSafeCellAvg>25, arr.ind = T))))
# 8575


#~#~#~#
#
#   Redundant genes proof of concept

#
# EMS plot with three homologous genes labelled, and each gene's pattern
#

singleCellMAplots(geneList = c(), 
                  outDir = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/", 
                  geneTypeName = "3emsHoms_noLabLighter", cellTypeList = "EMS", 
                  plotHeight = 4, plotWidth = 4, sig=0, opac = .1)

# alt version with enrichment significance: 
singleCellMAplots(geneList = c("Y75B12A.2", "Y50E8A.11", "AH10.2"), 
                  outDir = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/", 
                  geneTypeName = "3emsHoms_05", cellTypeList = "EMS", 
                  plotHeight = 5, plotWidth = 5, sig=.05)

# pictograms of 3 ems homologs
mapExpression("Y75B12A.2", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/EMS3_pictograms/Y75B12A.2.pdf")
mapExpression("Y50E8A.11", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/EMS3_pictograms/Y50E8A.11.pdf")
mapExpression("AH10.2", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/EMS3_pictograms/AH10.2.pdf")

#
#   RNAi plot
#

# Input RNAi data
RNAi.Y75.Y50 <- data.frame(hours = c("12-24", "24-36", "36-48"), 
                           unhatched = c(140, 403, 105),
                           total = c(704, 572, 361),
                           lethality=c(19.9, 70.4, 29.1))


# Input 95% confidence intervals (ward method), generated with graph pad stats calculator
limits <- aes(ymax=c(23,74,34), ymin=c(17,67,25))
dodge <- position_dodge(width=0.9)

c <- ggplot(RNAi.Y75.Y50, aes(x=factor(hours), y=lethality))+
    geom_bar(stat="identity", width=.75, color="black", fill="dark gray")+
    xlab("Hours after dsRNA injection\nof Y75B12A.2 and Y50E8A.11")+
    ylab("Percent embryonic lethality")+
    scale_x_discrete(breaks=c(1,2,3),  labels=c("12", "24", "36"))+
    geom_errorbar(limits, position=dodge, width=0.15)+
    theme_classic()

c

ggsave("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/RNAiBarplot_smaller.pdf",
       width = 2, height = 3)

#
#   Identifying groups of redundant genes
#

# Load in table
ceHoms <- read.table("Desktop/whoopsgraduateschool archive/20150904_homology_000000000000001.txt", header = F)
colnames(ceHoms) <- c("gene", "CHOM")

# Make gene names fixable
geneNameMatr <- read.csv("Desktop/whoopsgraduateschool archive/20150904_finalGeneNameMatr.txt", header=F, sep="\t")
geneNameMatr <- as.matrix(geneNameMatr)
# Fix gene names
geneNameMatr <- rearrangeNameMatrix(geneNameMatr, rownames(tlRPKMs))
ceHoms[,1] <- renameGeneVect(ceHoms[,1], geneNameMatr)

# Make all the sets of plots
# legacyCHOMsetsP05 <- CHOMsetsP05

CHOMsetsP05 <- list()
CHOMsetsP05 <- c(CHOMsetsP05, batchHomplots(cellList=list(c("Ep", "Ea"), c("MSx1", "MSx2"), c("Cx1", "Cx2"), c("ABa", "ABp"),
                                                    c("ABal", "ABar", "ABpl", "ABpr"), c("ABalx", "ABarx", "ABplx", "ABprx")), 
                            cellNameVect = c("Ex", "MSx", "Cx", "ABx", "ABxx", "ABxxx"), 
                            outDir = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/MAplotsOfCHOMs_P05/", 
                            returnSets = T, sig = .05, elim = "overlap"))
CHOMsetsP05 <- c(CHOMsetsP05, 
              batchHomplots(cellList = "all", 
                            outDir = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/MAplotsOfCHOMs_P05/", 
                            returnSets=T, sig=.05, elim="overlap"))
length(unlist(CHOMsetsP05))
CHOMsetsP05 <- clearRedundCHOMs(CHOMsetsP05, elim = "overlap")
length(unlist(CHOMsetsP05))

# make a table that summarizes all these plots

CHOMtable <- data.frame(cellType=character(),
                        CHOM = character(),
                        geneSet = character())

CHOMtable$cellType <- as.character(CHOMtable$cellType)
CHOMtable$CHOM <- as.character(CHOMtable$CHOM)
CHOMtable$geneSet <- as.character(CHOMtable$geneSet)

n=1
for(cellT in c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABx", 
               "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3", "ABxx",
               "ABalx", "ABarx", "ABplx", "ABprx", "ABxxx",
               "MSx", "MSx1", "MSx2", "Ex", "Ea", "Ep", "Cx", "Cx1", "Cx2", "D", "P4")){
    print(cellT)
    for(tempCHOM in names(CHOMsetsP05[[cellT]])){
        # make a character vector of all the genes
        tempGeneSet <- paste(CHOMsetsP05[[cellT]][[tempCHOM]], collapse = ", ")
        # print all the datas to the table
        CHOMtable[n,c("cellType", "CHOM", "geneSet")] <- c(cellT, tempCHOM, tempGeneSet)
        # one step forward
        n <- n+1
    }
}

write.table(x = CHOMtable, 
            file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/CHOMtable.txt", 
            sep = "\t", quote = F, col.names = NA)

# fill in the cellStats column for paralogous sets
#  ~~~ go to line ~1261

# Negative control: 
#   Make a fake scramble tlRPKM set and run the algorithm again to see how many sets I get
randoRedundTrials <-c()

for(i in 1:100){
    print(paste("ITERATION NUMBER", i, sep=" "))
    # Make fake-out matrix:
    randoRPKMs <- tlRPKMs
    set.seed(i)
    rownames(randoRPKMs) <- sample(rownames(randoRPKMs))
    
    randoCounts <- tlCounts
    set.seed(i)
    rownames(randoCounts) <- sample(rownames(randoCounts))
    
    # Stow real matrix elsewhere, swap out tlRPKMs
    #backupTlRPKMs <- tlRPKMs
    #backupTlCounts <- tlCounts
    tlRPKMs <- randoRPKMs
    tlCounts <- randoCounts
    
    # Run the same script, but store everything elsewhere
    randoCHOMsetsP05 <- list()
    randoCHOMsetsP05 <- c(randoCHOMsetsP05, batchHomplots(cellList=list(c("Ep", "Ea"), c("MSx1", "MSx2"), c("Cx1", "Cx2"), c("ABa", "ABp"),
                                                              c("ABal", "ABar", "ABpl", "ABpr"), c("ABalx", "ABarx", "ABplx", "ABprx")), 
                                                cellNameVect = c("Ex", "MSx", "Cx", "ABx", "ABxx", "ABxxx"), 
                                                outDir = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/MAplotsOfCHOMs_P05_rando/", 
                                                returnSets = T, sig = .05, elim = "overlap"))
    randoCHOMsetsP05 <- c(randoCHOMsetsP05, 
                     batchHomplots(cellList = c("AB", "P1", "ABa", "ABp", "EMS", "P2", 
                                                "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3",
                                                "ABalx", "ABarx", "ABplx", "ABprx", "D", "P4"), 
                                   outDir = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/MAplotsOfCHOMs_P05_rando/", 
                                   returnSets=T, sig=.05, elim="overlap"))
    length(unlist(randoCHOMsetsP05))
    randoCHOMsetsP05 <- clearRedundCHOMs(randoCHOMsetsP05, elim = "overlap")
    length(unlist(randoCHOMsetsP05))
    
    # make a table that summarizes all these plots
    
    randoCHOMtable <- data.frame(cellType=character(),
                            CHOM = character(),
                            geneSet = character())
    
    randoCHOMtable$cellType <- as.character(randoCHOMtable$cellType)
    randoCHOMtable$CHOM <- as.character(randoCHOMtable$CHOM)
    randoCHOMtable$geneSet <- as.character(randoCHOMtable$geneSet)
    
    n=1
    for(cellT in c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABx", 
                   "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3", "ABxx",
                   "ABalx", "ABarx", "ABplx", "ABprx", "ABxxx",
                   "MSx", "MSx1", "MSx2", "Ex", "Ea", "Ep", "Cx", "Cx1", "Cx2", "D", "P4")){
        print(cellT)
        for(tempCHOM in names(randoCHOMsetsP05[[cellT]])){
            # make a character vector of all the genes
            tempGeneSet <- paste(randoCHOMsetsP05[[cellT]][[tempCHOM]], collapse = ", ")
            # print all the datas to the table
            randoCHOMtable[n,c("cellType", "CHOM", "geneSet")] <- c(cellT, tempCHOM, tempGeneSet)
            # one step forward
            n <- n+1
        }
    }
    
    write.table(x = randoCHOMtable, 
                file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/CHOMtable_rando.txt", 
                sep = "\t", quote = F, col.names = NA)
    
    # figure out number of sets:
    tempNum <- 0
    for(cellT in c("AB", "P1", "ABa", "ABp", "ABx", "EMS", "P2", 
                   "ABal", "ABar", "ABpl", "ABpr", "ABxx", "MS", "E", "C", "P3",
                   "ABalx", "ABarx", "ABplx", "ABprx", "ABxxx", "MSx", "Ex", "Cx", "D", "P4")){
        tempNum <- tempNum + length(randoCHOMsetsP05[[cellT]])
    }
    # put number of sets to vector
    randoRedundTrials <-c(randoRedundTrials, tempNum)
    # print vector, as is
    print(paste(randoRedundTrials, collapse = ", "))
}

# Restore the original tlRPKM file
tlRPKMs <- backupTlRPKMs
tlCounts <- backupTlCounts

pdf("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/paralogousRandos.pdf", 
    3,3)
hist(c(randoRedundTrials, 99), breaks=25, xlim=c(0,100), ylim=c(0,40),
     xlab="number of paralogous sets\nfound enriched in any cell type", 
     ylab="frequency", axes = F,
     main = "paralogous gene set analysis\ndone with randomized data", 
     col = c(rep("gray", 15), rep("red", 10)))
axis(1, at=c(0,25, 50, 75, 100))
axis(2, at=c(0,20,40))
dev.off()


#~#~#~#
#
#   Overview of data

#
#   Export a big heatmap of all genes x all timepoints
#

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
for(i in 1:length(cellTypes)){
    colOrder = c(colOrder, which(goodCellTypes==cellTypes[i]))
}

subMatr <- subMatr[,colOrder]

# log matrix
subMatr <- log10(subMatr+1)

# export this matrix
write.table(subMatr, "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/01_overview_bigBlue/heatmap_pieces/allGenes.txt", sep="\t", col.names = NA, quote=F)

# import cdt of this full dataset
tempAllCells <-  read.csv("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/01_overview_bigBlue/heatmap_pieces/allGenes.cdt", sep = "\t")
# import the gene order list
tempTest <- read.csv("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/01_overview_bigBlue/heatmap_pieces/10_geneOrder.txt")
tempReorder <- tempTest[,2]
names(tempReorder) <- tempTest[,1]
tempOut <- tempAllCells[c(1,tempReorder),]
write.table(tempOut, "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/01_overview_bigBlue/heatmap_pieces/09_allGenes.cdt", 
            sep="\t", quote=F, row.names = F)
# then open this up in java tree view to generate a PNG


#
#   My/Tamar's/Erin's data
#

# make differential expression tables for each of our datasets

erinFC <- ERdiffExV1(erinsRPKMs, erinDesign, mastMatr = erinsRPKMs, to = "list", sig = 2)
erinFC <- rbind(erinFC[[1]], erinFC[[2]])
#13119 genes, call everything else no fold change

hashimFC <- ERdiffExV1(hashimRPKMs[,1:6], hashimDesign, mastMatr = hashimRPKMs, to = "list", sig = 2)
hashimFC <- rbind(hashimFC[[1]], hashimFC[[2]])
# 8421 genes, call everything else no fold change

my2cellFC <- ERdiffExV1(tlRPKMs[,typeCells(c("AB", "P1"))], to = "list", sig = 2)
my2cellFC <- rbind(my2cellFC[[1]], my2cellFC[[2]])

# concatenate list of genes
all2cellFCgenes <- union(union(rownames(erinFC), rownames(hashimFC)), rownames(my2cellFC))
all2cellFCgenes <- intersect(all2cellFCgenes, rownames(tlRPKMs))
all2cellFCgenes <- intersect(all2cellFCgenes, rownames(erinsRPKMs))

# set up tables that will be used for plotting

# legacyRef <- ref
# legacyAll2cellFCs <- all2cellFCs
# legacyAll2cellCPMs <- all2cellCPMs
# legacyEnrichMatr <- enrichMatr

ref = list("Nishimura" = erinFC, "Hashimshony" = hashimFC, "Tintori" = my2cellFC)

all2cellFCs <- matrix(0, ncol=3, nrow=length(all2cellFCgenes), 
                      dimnames = list(all2cellFCgenes, c("Nishimura", "Hashimshony", "Tintori")))
all2cellCPMs <- matrix(0, ncol=3, nrow=length(all2cellFCgenes), 
                       dimnames = list(all2cellFCgenes, c("Nishimura", "Hashimshony", "Tintori")))

for(gene in all2cellFCgenes){
    for(study in c("Nishimura", "Hashimshony", "Tintori")){
        all2cellFCs[gene, study] <- ref[[study]][gene,"logFC"]
        all2cellCPMs[gene, study] <- ref[[study]][gene, "logCPM"]
    }
}
all2cellFCs[is.na(all2cellFCs)] <- 0
all2cellCPMs[is.na(all2cellCPMs)] <- 0

# make plots

enrichMatr <- makeEnrichTable(constant=4)

# Tintori vs Hashimshony
pdf("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS2_overviewAndComparisons/01_hashimNishimTinto_comparisons/tintoriHashimshony_c4.pdf",4,4)
par(mar=c(5.1, 5.1, 2.1, 2.1))
ABcompPlot("Tintori", "Hashimshony", axisLim = 20000)
dev.off()

# Hashimshony v Nishimura
pdf("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS2_overviewAndComparisons/01_hashimNishimTinto_comparisons/hashimshonyNishimura_c4.pdf",4,4)
par(mar=c(5.1, 5.1, 2.1, 2.1))
ABcompPlot("Hashimshony", "Nishimura", axisLim = 20000)
dev.off()

# Nishimura v Tintori
pdf("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS2_overviewAndComparisons/01_hashimNishimTinto_comparisons/nishimuraTintori_c4.pdf",4,4)
par(mar=c(5.1, 5.1, 2.1, 2.1))
ABcompPlot("Nishimura", "Tintori", axisLim = 20000)
dev.off()

#
#   Big blue heatmaps
#

cellTypes <- c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3", 
               "ABxxx1", "ABxxx2", "ABxxx3", "ABxxx4", "MSx1", "MSx2", "Ep", "Ea", "Cx1", "Cx2", "D", "P4")

# select just genes with a reasonable amount of variation
tempAvgLogged <- log10(tlCellAvg+1)
tempSDgenes <- rownames(tlCellAvg[which(apply(tempAvgLogged, 1, sd)>.5),])

exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_overviewAndComparisons/geneMatrices/testScript.txt", 
                                geneSet=tempSDgenes)


# print out some example genes

# expected:
mapExpression("ceh-51", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/ceh-51.pdf")
mapExpression("cey-2", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/cey-2.pdf")
mapExpression("end-1", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/end-1.pdf")
mapExpression("tbx-37", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/tbx-37.pdf")
mapExpression("sdz-38", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/sdz-38.pdf")
mapExpression("elt-7", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/elt-7.pdf")
mapExpression("cwn-1", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/cwn-1.pdf")

# with legends
mapExpression("ceh-51", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/ceh-51_legend.pdf")
mapExpression("cey-2",
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/cey-2_legend.pdf")
mapExpression("end-1",
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/end-1_legend.pdf")
mapExpression("tbx-37", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/tbx-37_legend.pdf")
mapExpression("sdz-38",
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/sdz-38_legend.pdf")
mapExpression("elt-7",
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/elt-7_legend.pdf")
mapExpression("cwn-1",
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/cwn-1_legend.pdf")

# ranges, for figuring out how to group legends together:
t(apply(tlSafeCellAvg[c("sdz-38", "tbx-37", "ceh-51", "elt-7", "end-1", "cwn-1", "daz-1", "cey-2", "skr-10"), safeCellTypes], 1, range))

#mapExpression("M106.3",  mapType = "horz", width = 6, height = 1.5, fillRange = c(0,1204),
#              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/M106.3.pdf")
# for legends:
mapExpression("end-1", mapType = "horz", width = 6, height = 1.5, fillRange = c(0,1204),
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/end-1_legend.pdf")
mapExpression("tbx-37", mapType = "horz", width = 6, height = 1.5, fillRange = c(0,300),
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/tbx-37_legend.pdf")


# somatic vs germ dichotomy:
mapExpression("daz-1", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/daz-1.pdf")
mapExpression("skr-10", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/skr-10.pdf")
# for legends
mapExpression("daz-1", mapType = "horz", width = 6, height = 1.5,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/daz-1_legend.pdf")
mapExpression("skr-10", mapType = "horz", width = 6, height = 1.5,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/exampleCellMaps/skr-10_legend.pdf")


#
#   Baugh comparison
#

# Import Ryan Baugh's 2005 microarray data timecourse for the relevant time points

PAbaugh2005 <- read.table("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS2_overviewAndComparisons/02_baughData/avgScores.txt", sep="\t", header = T)
GNbaugh2005 <- read.table("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS2_overviewAndComparisons/02_baughData/probeNameLookup.txt", sep="\t", header=T)
GN2baugh2005 <- read.table("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS2_overviewAndComparisons/02_baughData/geneNames.txt", sep="\t", header=T)
GNbaugh2005[GNbaugh2005=="---"] <- NA
GN2baugh2005[GN2baugh2005=="---"] <- NA
GNbaugh2005 <- merge(GNbaugh2005, GN2baugh2005, by.x = "Probe.Set.ID", by.y="PROBEID")
GNbaugh2005 <- cbind(GNbaugh2005, tlMatchID=NA)

for(colmn in colnames(GNbaugh2005)){
    GNbaugh2005[,colmn] <- as.character(GNbaugh2005[,colmn])
}

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

# have to remove all the duplicates
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
rm(PAbaugh2005, GNbaugh2005, GN2baugh2005, mergeBaugh2005, temp)

###
# Curate lists of transcribed and degraded genes, from both my and Ryan's data:
#

thresh = 25

# come up with a list of genes that have increasing expression in baugh's data
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

# come up with a list of genes that have decreasing expression in baugh's dataset
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


# what genes are transcribed in my dataset
myTransGenes <- NULL
for(gene in rownames(tlCellAvg[which(apply(tlCellAvg, 1, max)>thresh),])){
    # for each cell that has expression above threshold..
    for(cellWith in names(sort(which(tlCellAvg[gene,]>thresh), decreasing = T))){
        # check the minimum of ancestor cells, to see if any of them are significantly lower than the cellWith
        cellAncestors <- as.vector(unlist(fakeLinRelations[[cellWith]][c("parent", "grandparent", "greatgrandparent", "greatgreatgrandparent")]))
        if(length(cellAncestors)==0){next}
        tempMin <- min(tlCellAvg[gene, cellAncestors])
        if(tempMin*2 < tlCellAvg[gene, cellWith]){
            myTransGenes <- c(myTransGenes, gene)
            break
        }
    }
}


# what genes are degraded in my dataset?
myDegGenes <- NULL
for(gene in rownames(tlCellAvg[which(apply(tlCellAvg, 1, max)>thresh),])){
    # for each cell that has expression above threshold..
    for(cellWith in names(which(tlCellAvg[gene,]>thresh))){
        # check the minimum of descendant cells, to see if any of them are significantly lower than the cellWith
        cellDescendants <- as.vector(unlist(fakeLinRelations[[cellWith]][c("daughters", "granddaughters", "greatgranddaughters", "greatgreatgranddaughters")]))
        if(length(cellDescendants)==0){next}
        tempMin <- min(tlCellAvg[gene, cellDescendants])
        if(tempMin*2 < tlCellAvg[gene, cellWith]){
            myDegGenes <- c(myDegGenes, gene)
            break
        }
    }
}

# what are the genes that both studies had access to?
tintoBaughCommonAll <- intersect(rownames(tlRPKMs), rownames(baugh2005))

# FOR THRESH = 25
length(baughTransGenes)     # 1935
length(myTransGenes)        # 8239
length(intersect(baughTransGenes, myTransGenes))    # 1760
length(intersect(myTransGenes, tintoBaughCommonAll))  # 7218
setdiff(myTransGenes, baughTransGenes)  # 6479
length(setdiff(myTransGenes, tintoBaughCommonAll))  # 1021

length(baughDegGenes)   # 2164
length(myDegGenes)      # 6944
length(intersect(baughDegGenes, myDegGenes))        # 2099
length(intersect(myDegGenes, tintoBaughCommonAll))    # 6324
setdiff(myDegGenes, baughDegGenes))     # 4845
length(setdiff(myDegGenes, tintoBaughCommonAll))    # 620

# export these groups of genes to make matrices
exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/transBaughOnly.txt", 
                    geneSet=setdiff(baughTransGenes, myTransGenes))
exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/transBoth.txt", 
                    geneSet=intersect(baughTransGenes, myTransGenes))
exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/transMine_baughAccess.txt", 
                    geneSet=intersect(setdiff(myTransGenes, baughTransGenes), tintoBaughCommonAll))
exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/transMine_baughNoAccess.txt", 
                    geneSet=setdiff(setdiff(myTransGenes, baughTransGenes), tintoBaughCommonAll))

exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_TranscrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/degBaughOnly.txt", 
                    geneSet=setdiff(baughDegGenes, myDegGenes))
exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_TranscrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/degBoth.txt", 
                    geneSet=intersect(baughDegGenes, myDegGenes))
exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_TranscrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/degMine_baughAccess.txt", 
                    geneSet=intersect(setdiff(myDegGenes, baughDegGenes), tintoBaughCommonAll))
exportGeneExpMatrix(outFile="Desktop/whoopsgraduateschool archive/20151231_TranscrLineagePaper/Figures/v2/fig4_overviewAndComparisons/03_baughTintoComparisonMatrices/degMine_baughNoAccess.txt", 
                    geneSet=setdiff(setdiff(myDegGenes, baughDegGenes), tintoBaughCommonAll))




#~#~#~#
#
#   MapStats

# tlCellStatsLegacy <- tlCellStats

cellTypes <- unique(c(tlDesign$ID, tlDesign$safeID, tlDesign$ABsafe))
cellTypes <- cellTypes[-which(cellTypes %in% c(NA, "tossed"))]

tlCellStats <- matrix(0, ncol=1, nrow=length(cellTypes))
rownames(tlCellStats) <- cellTypes

# first column, how many genes are present over 25 rpkm
colnames(tlCellStats) <- "presentOver25"

for(cell in cellTypes){
    if(cell %in% tlDesign$safeID){
        tempNum <- length(which(tlSafeCellAvg[,cell]>25))
        tlCellStats[cell,"presentOver25"] <- tempNum
    }
}

# second column, avg pg rna found in each cell
tlCellStats <- cbind(tlCellStats, pgRNA=0)
for(cell in cellTypes){
    if(cell %in% tlDesign$safeID){
        tempNum <- mean(tlDesign[which(tlDesign$safeID==cell),"input"])
        tlCellStats[cell,"pgRNA"] <- tempNum
    }
}

# third column, other metric for pg, I forget how this was done...
tlCellStats <- cbind(tlCellStats, pgExpected=0)
for(cell in cellTypes){
    if(cell %in% tlDesign$safeID){
        tempNum <- mean(tlDesign[which(tlDesign$safeID==cell),"pgExpected"])
        tlCellStats[cell,"pgExpected"] <- tempNum
    }
}

# fourth column, newly expressed genes
tlCellStats <- cbind(tlCellStats, newlyExpr=0)

for(cell in cellTypes[-which(cellTypes=="P0")]){
    if(cell %in% tlDesign$safeID){
        print(cell)
        # identify all previous cells
        relatsAvail <- intersect(names(fakeLinRelations[[cell]]), 
                                 c("greatgreatgrandparent", "greatgrandparent", "grandparent", "parent"))
        relatNames <- NULL
        for(relat in relatsAvail){relatNames <- c(relatNames, fakeLinRelations[[cell]][[relat]])}
        genesInCell <- names(which(tlSafeCellAvg[,cell]>25))
        newGenes <- NULL
        for(gene in genesInCell){
            if(max(tlSafeCellAvg[gene, relatNames])>10){next}
            else{newGenes <- c(newGenes, gene)}
        }
        tempNum <- length(newGenes)
        tlCellStats[cell,"newlyExpr"] <- tempNum
    }
}

# fifth column, upregulated genes
tlCellStats <- cbind(tlCellStats, upReg=0)

for(cell in cellTypes[-which(cellTypes=="P0")]){
    if(cell %in% tlDesign$safeID){
        print(cell)
        # identify parent
        parent <- fakeLinRelations[[cell]][["parent"]]
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

# sixth column, downregulated genes
tlCellStats <- cbind(tlCellStats, downReg=0)

for(cell in cellTypes[-which(cellTypes=="P0")]){
    if(cell %in% tlDesign$safeID){
        print(cell)
        # identify parent
        parent <- fakeLinRelations[[cell]][["parent"]]
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

# seventh column, average expression level of expressed genes
tlCellStats <- cbind(tlCellStats, avgRPKMs=0)

for(cell in cellTypes){
    if(cell %in% tlDesign$safeID){
        print(cell)
        genesInCell <- names(which(tlSafeCellAvg[,cell]>25))
        tempNum <- mean(tlSafeCellAvg[genesInCell, cell])
        tlCellStats[cell,"avgRPKMs"] <- tempNum
    }
}

# eight column, how many genes are unique to that cell?
tlCellStats <- cbind(tlCellStats, uniqueGenes=0)

for(cell in cellTypes){
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

# add a column for sets of homologous genes
tlCellStats <- cbind(tlCellStats, CHOMgroups=0)
for(cell in names(CHOMsetsP05)){
    if(cell %in% rownames(tlCellStats)){
        tlCellStats[cell,"CHOMgroups"] <- length(CHOMsetsP05[[cell]])
    }
}
for(cell in c("ABal", "ABar", "ABpl", "ABpr")){
    tlCellStats[cell, "CHOMgroups"] <- tlCellStats[cell, "CHOMgroups"] + length(CHOMsetsP05[["ABxx"]])
}
for(cell in c("ABalx", "ABarx", "ABplx", "ABprx")){
    tlCellStats[cell, "CHOMgroups"] <- tlCellStats[cell, "CHOMgroups"] + length(CHOMsetsP05[["ABxxx"]])
}

# map cell cycle length

tlCellStats <- cbind(tlCellStats, cycleLen = 0)
tlCellStats[c("P0", "AB", "P1", "ABa", "ABp", "EMS", "P2", "ABx", 
              "ABal", "ABar", "ABpl", "ABpr", "MS", "E", "C", "P3", "ABxx",
              "ABalx", "ABarx", "ABplx", "ABprx", "ABxxx",
              "MSx", "MSx1", "MSx2", "Ex", "Ea", "Ep", "Cx", "Cx1", "Cx2", "D", "P4"), "cycleLen"] <- c(
                  40, 17, 18, 17, 17, 22, 25, 17, 
                  18, 18, 18, 18, 20, 23, 24, 33, 18,
                  22, 22, 22, 22, 22,
                  27, 27, 27, 45, 44, 46, 30, 30, 30, 40, 67)

# what percentage of genes expressed are transcription factors
tlCellStats <- cbind(tlCellStats, presentPercTF=0)
for(cell in cellTypes){
    if(cell %in% tlDesign$safeID){
        tempGenes <- names(which(tlSafeCellAvg[,cell]>25))
        tempNum <- length(tempGenes)
        tempTFnum <- length(intersect(tempGenes, TFs_all))
        tempPerc <- (tempTFnum/tempNum)*100
        tlCellStats[cell,"presentPercTF"] <- tempPerc
    }
}

mapStats("presentPercTF", fillRange = c(7,9))
# this looks... like nothing. everything is all right around 8 percent for some reason

# mapStats("presentOver25", legendTitle = "# genes\ndetected",
#          to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/detectedOver25.pdf")


# what percentage of new genes are transcription factors
tlCellStats <- cbind(tlCellStats, upRegPercTF=0)

for(cell in cellTypes[-which(cellTypes=="P0")]){
    if(cell %in% tlDesign$safeID){
        print(cell)
        # identify parent
        parent <- fakeLinRelations[[cell]][["parent"]]
        # find genes that are expressed in cell
        genesInCell <- names(which(tlSafeCellAvg[,cell]>25))
        # just keep genes that are twice the expression as in the parent
        newGenes <- NULL
        for(gene in genesInCell){
            if(tlSafeCellAvg[gene, parent]<(.5*tlSafeCellAvg[gene, cell])){newGenes <- c(newGenes, gene)}
        }
        tempGeneNum <- length(newGenes)
        tempTFnum <- length(intersect(newGenes, TFs_all))
        tempPerc <- (tempTFnum/tempGeneNum)*100
        tlCellStats[cell,"upRegPercTF"] <- tempPerc
    }
}

mapStats("upRegPercTF")
# more static

tlCellStats <- cbind(tlCellStats, newlyExprPercTF=0)

for(cell in cellTypes[-which(cellTypes=="P0")]){
    if(cell %in% tlDesign$safeID){
        print(cell)
        # identify all previous cells
        relatsAvail <- intersect(names(fakeLinRelations[[cell]]), 
                                 c("greatgreatgrandparent", "greatgrandparent", "grandparent", "parent"))
        relatNames <- NULL
        for(relat in relatsAvail){relatNames <- c(relatNames, fakeLinRelations[[cell]][[relat]])}
        genesInCell <- names(which(tlSafeCellAvg[,cell]>25))
        newGenes <- NULL
        for(gene in genesInCell){
            if(max(tlSafeCellAvg[gene, relatNames])>10){next}
            else{newGenes <- c(newGenes, gene)}
        }
        tempNum <- length(newGenes)
        tempTFnum <- length(intersect(newGenes, TFs_all))
        tempPerc <- (tempTFnum/tempNum)*100
        tlCellStats[cell,"newlyExprPercTF"] <- tempPerc
    }
}

mapStats("newlyExprPercTF")

# mapStats("newlyExpr", legendTitle = "# genes\nnewly expressed",
#          to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/newlyExpressed.pdf")


# what percentage of unique genes are transcription factors
# eight column, how many genes are unique to that cell?
tlCellStats <- cbind(tlCellStats, uniqueGenesPercTF=0)

for(cell in cellTypes){
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
        tempTFnum <- length(intersect(newGenes, TFs_all))
        if(tempNum==0){tlCellStats[cell,"uniqueGenesPercTF"] <- 0}
        else{
            tempPerc <- (tempTFnum/tempNum)*100
            tlCellStats[cell,"uniqueGenesPercTF"] <- tempPerc
        }
    }
}



mapStats("uniqueGenesPercTF")

# hm. this one is cool. 

# Make all map stats
mapStats("upReg", legend = F,
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/upreg.pdf")
mapStats("downReg", legendTitle = "# genes\ndownregulated\n", legend = F,
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/downreg.pdf")
mapStats("uniqueGenes", legendTitle = "# of unique\ngenes expressed\n", legend = F,
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/uniqueGenes.pdf")
mapStats("uniqueGenesPercTF", legendTitle = "percentage of\nunique genes\nthat are\ntxn factors\n", legend = F,
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/uniqueGenesPercTF.pdf")
mapStats("presentOver25", legendTitle = "# genes\ndetected\n", legend = F,
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/detectedOver25.pdf")
mapStats("pgRNA", legendTitle = "pg mRNA\ndetected\n", legend = F,
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/pgRNA.pdf")
mapStats("cycleLen", legendTitle = "length of\ncell cycle\n(minutes)\n", fillRange = c(0,70), legend = F,
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/cellCycle.pdf")

# get legends
mapStats("upReg", legendTitle = "# genes\nupregulated\n",
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/upreg_legend.pdf")
mapStats("downReg", legendTitle = "# genes\ndownregulated\n",
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/downreg_legend.pdf")
mapStats("uniqueGenes", legendTitle = "# of unique\ngenes expressed\n",
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/uniqueGenes_legend.pdf")
mapStats("uniqueGenesPercTF", legendTitle = "percentage of\nunique genes\nthat are\ntxn factors\n",
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/uniqueGenesPercTF_legend.pdf")
mapStats("presentOver25", legendTitle = "# genes\ndetected\n",
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/detectedOver25_legend.pdf")
mapStats("pgRNA", legendTitle = "pg mRNA\ndetected\n",
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/pgRNA_legend.pdf")
mapStats("cycleLen", legendTitle = "length of\ncell cycle\n(minutes)\n", fillRange = c(0,70),
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/cellCycle_legend.pdf")


# mapStats("newlyExpr", legendTitle = "# genes\nnewly expressed\n",
#          to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/newlyExpressed_legend.pdf")
# mapStats("avgRPKMs", legendTitle = "avg RPKM of\ngenes expressed\n",
#          to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/avgRPKM_legend.pdf")
# mapStats("pgExpected", legendTitle = "pg mRNA detected",
#          to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/mapStats/pgRNAexpected_legend.pdf")


# this one doesn't have to be horizontal. for a different figure
mapStats("CHOMgroups", legendTitle = "# of sets of\nparalogous\ngenes enriched\nin each cell",
         to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig3_redundantGenes/paralogs.pdf")

# Make some plots of how amount of how pgRNA, # of genes, and cell cycle length
#   corresponds to upreg, downreg, unique, %TF unique

nonP4 <- c("ABalx", "ABarx", "ABplx", "ABprx", "MSx", "Ex", "Cx", "D")
tempMatr <- tlCellStats[nonP4, c("upReg", "downReg", "uniqueGenes", "uniqueGenesPercTF", "pgRNA", "presentOver25", "cycleLen")]

pdf("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig5_newlySomatic/statMap_corrs/corrHeatmap.pdf", 5,5)
heatmap.2(cor(tempMatr), Rowv = F, Colv = F, trace = 'none', 
          col=colorRampPalette(c('black', 'gray', 'white')), 
          labRow = c("E", "F", "G", "H", "I", "J", "K"),
          labCol = c("E", "F", "G", "H", "I", "J", "K"),
          #margins=c(12,12), 
          density.info = 'none', colsep=4, rowsep=4, sepcolor = "red")
dev.off()



#~#~#~#
#
#   Spatially dynamic gene expression

#
#   tbx-32 lineage map
#

mapExpression("tbx-32", 
              to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/tbx-32_plots/tbx-32_pictogram.pdf")
mapExpression("tbx-32", mapType = "lineage", safeCell=F,width = 20, height=12,
              to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/tbx-32_plots/tbx-32_lineage.pdf")
# for the outline
mapExpression("2L52.2", mapType = "lineage", safeCell=F,width = 20, height=12,
              to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/tbx-32_plots/outline_lineage.pdf")
mapExpression("tbx-32", mapType = "horz",
              to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/tbx-32_plots/tbx-32_horizontal.pdf")

# Pictograms of other genes with tbx-32 expression patterns
mapExpression("tbx-32", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/tbx-32.pdf")
mapExpression("tbx-31", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/tbx-31.pdf")
mapExpression("tbx-40", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/tbx-40.pdf")
mapExpression("Y43D4A.6", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/Y43D4A.6.pdf")
mapExpression("Y116A8C.20", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/Y116A8C.20.pdf")
mapExpression("ZK666.1", mapType = "horz", width = 6, height = 1.5, legend = F,
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/ZK666.1.pdf")

# with legends
mapExpression("tbx-32", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/tbx-32_legend.pdf")
mapExpression("tbx-31", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/tbx-31_legend.pdf")
mapExpression("tbx-40", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/tbx-40_legend.pdf")
mapExpression("Y43D4A.6", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/Y43D4A.6_legend.pdf")
mapExpression("Y116A8C.20", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/Y116A8C.20_legend.pdf")
mapExpression("ZK666.1", 
              to = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/figS4_tbx32_smFISH_stepBy/similarPictos/ZK666.1_legend.pdf")


# #
# #   ID dynamic genes
# #
# 
# # First, establish a list of lists, where each gene points to the cells that express it
# expressionHits <- list()
# for(gene in rownames(tlSafeCellAvg)){
#     tmMax <- max(tlSafeCellAvg[gene,])
#     if(tmMax < 25) {next}
#     # establish basic stats on gene
#     tmMean <- mean(tlSafeCellAvg[gene,])
#     # if it doesn't pass the differential expression threshold, skip it
#     if((2*tmMean) >= tmMax) {next}
#     # for genes that pass, push a list of expressing cells to the list
#     expressionHits[[gene]] <- names(which(tlSafeCellAvg[gene,]>(2*tmMean)))
# }
# length(expressionHits) # 7546
# 
# # For each gene, if it's pattern spans multiple islands of cells, add it to a new pruned list
# spatDynList <- list()
# for(gene in names(expressionHits)){
#     # ID gene with highest expression
#     topCell <- names(which(tlSafeCellAvg[gene,]==max(tlSafeCellAvg[gene,])))
#     tmThresh = (2*mean(tlSafeCellAvg[gene,]))
#     startGroup = whichNeighborsHaveGene(topCell, gene, tmThresh, safeCells = T)
#     newGroup = NULL
#     for(groupCell in startGroup){
#         newGroup <- c(newGroup, whichNeighborsHaveGene(groupCell, gene, tmThresh, safeCells = T))
#     }
#     newGroup <- unique(newGroup)
#     while(suppressWarnings(all.equal(sort(startGroup), sort(newGroup))) != TRUE){
#         startGroup <- newGroup
#         tempGroup <- NULL
#         for(groupCell in startGroup){
#             tempGroup <- c(tempGroup, whichNeighborsHaveGene(groupCell, gene, tmThresh, safeCells = T))
#         }
#         newGroup = unique(tempGroup)
#     }
#     # if there are more expressing cells left that are not encompassed in this list, 
#     #   the gene is spatially dynamic
#     if(suppressWarnings(all.equal(sort(startGroup), sort(expressionHits[[gene]]))) != TRUE){
#         spatDynList[[gene]] <- expressionHits[[gene]]
#     }
# }
# length(spatDynList) # 4185
# 
# # Finally, figure out how to group these by similar or exactly expression patterns
# # For now, starting with exact only
# spatDynPattList <- list()
# for(gene in names(spatDynList)){
#     # make a compound name for the expression pattern
#     pattName <- patternToPattName(gene)
#     # add gene to existing or not existing pattern
#     spatDynPattList[[pattName]] <- c(spatDynPattList[[pattName]], gene)
# } 
# length(spatDynPattList) # 2369
# 
# # Get a sense of how many patterns are represented by more than one gene
# spatDynPattScores <- sapply(spatDynPattList, length)
# spatDynPattScores <- sort(spatDynPattScores, decreasing = T)
# hist(spatDynPattScores, ylim=c(0,10), xlim=c(10,40), breaks=100)
# 
# # Exploring the most common expression patterns
# 
# spatDynPattScores[1]
# mapExpression(spatDynPattList[[names(spatDynPattScores[1])]],  safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/01.pdf")
# # ABxxx1_ABxxx2_ABxxx3_ABxxx4_MSx 
# # 40 
# 
# spatDynPattScores[2]
# mapExpression(spatDynPattList[[names(spatDynPattScores[2])]],  safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/02.pdf")
# # ABxxx1_ABxxx2_ABxxx3_ABxxx4_MSx_Ex 
# # 39 
# 
# spatDynPattScores[3]
# mapExpression(spatDynPattList[[names(spatDynPattScores[3])]],  safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/03.pdf")
# # Ex_Cx 
# # 29 
# 
# spatDynPattScores[4]
# mapExpression(spatDynPattList[[names(spatDynPattScores[4])]],  safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/04.pdf")
# # MS_D 
# # 28 
# 
# spatDynPattScores[5]
# mapExpression(spatDynPattList[[names(spatDynPattScores[5])]],  safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/05.pdf")
# # ABxxx2_Ex 
# # 22 
# 
# spatDynPattScores[6]
# mapExpression(spatDynPattList[[names(spatDynPattScores[6])]], safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/06.pdf")
# # C_D 
# # 19 
# 
# spatDynPattScores[7]
# mapExpression(spatDynPattList[[names(spatDynPattScores[7])]], safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/07.pdf")
# # ABxxx1_ABxxx2_ABxxx3_ABxxx4_Ex 
# # 17 
# 
# spatDynPattScores[8]
# mapExpression(spatDynPattList[[names(spatDynPattScores[8])]], safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/08.pdf")
# # MS_Cx 
# # 17 
# 
# spatDynPattScores[9]
# mapExpression(spatDynPattList[[names(spatDynPattScores[9])]], safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/09.pdf")
# # Ex_D 
# # 17 
# 
# spatDynPattScores[10]
# mapExpression(spatDynPattList[[names(spatDynPattScores[10])]], safeCell = T,
#               to="Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/fig4_spatiallyDynamicGenes/newSpatDynGroups/10.pdf")
# # ABxxx2_MSx 
# # 17 
# 
# 
# 
# #
# #   Dynamic vs null hypothesis histogram
# #
# 
# #
# #   Hottest cells? Hottest patterns? Hottest genes?
# #
# 

### for design table with details about all samples
tempDes <- tlDesign[,c("stage", "replicate", "quality", "safeID", "input", "genesDetected", "uniqueCounts", "ERCCcounts")]
#tempDes <- cbind(tempDes, "V2" = rownames(tempDes))
fromNCBI <- read.table("Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/sampleTable/fromNCBI.txt", 
                       sep="\t")
tempDes <- merge(x=fromNCBI[,c(2,3,1)], 
                 y=tempDes[,c("quality", "input", "genesDetected", "uniqueCounts", "ERCCcounts")], 
                 by.x="V2", by.y="row.names", sort = F)
tempDes$input <- round(tempDes$input, digits = 2)
colnames(tempDes) <- c("short name", "sample", "NCBI GEO accession", "Usable Quality?", "estimated pg mRNA", "Number of Genes Detected", "Number of Unique Worm Reads", "Number of ERCC reads")
tempDes[which(tempDes$`Usable Quality?`=="Use"), "Usable Quality?"] <- "Yes"
tempDes[which(tempDes$`Usable Quality?`=="Toss"), "Usable Quality?"] <- "No"

write.table(x = tempDes[,2:length(colnames(tempDes))], 
            file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/sampleTable/finalTable.txt", 
            sep = "\t", quote = F, row.names = F)

### for RPKM tables
write.table(x = tlRPKMs, 
            file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/RPKM tables/wormRPKMs.txt", 
            sep = "\t", quote = F, row.names = F)
write.table(x = erccRPKMs, 
            file = "Desktop/whoopsgraduateschool archive/20151231_transcrLineagePaper/Figures/v2/supps/RPKM tables/erccRPKMs.txt", 
            sep = "\t", quote = F, row.names = F)


### for interactive data viz
cellTypes <- unique(tlDesign$ID)
cellTypes <- cellTypes[!cellTypes %in% c("tossed", "ABxxx", NA)]
singleCellMAplots(outDir = "Desktop/whoopsgraduateschool archive/20150915_singleCellMAplotsForGenomeBrowser/", cellTypeList = "all")

ERdiffExV1(tlCounts[,c(typeCells(c("ABxx1", "ABxx2", "ABxx3", "ABxx4")), typeCells(c("E", "MS", "C", "P3")))],
           group=factor(c(rep("ABxx's", length(typeCells(c("ABxx1", "ABxx2", "ABxx3", "ABxx4")))), 
                          rep("P1desc's", length(typeCells(c("E", "MS", "C", "P3")))))),
           label=rownames(tlCounts), to="Desktop/whoopsgraduateschool archive/20150915_singleCellMAplotsForGenomeBrowser/8cell_ABvP1.pdf" )

ERdiffExV1(tlCounts[,c(typeCells(c("ABxxx1", "ABxxx2", "ABxxx3", "ABxxx4")), typeCells(c("Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2", "D", "P4")))],
           group=factor(c(rep("ABxxx's", length(typeCells(c("ABxxx1", "ABxxx2", "ABxxx3", "ABxxx4")))), 
                          rep("P1desc's", length(typeCells(c("Ep", "Ea", "MSx1", "MSx2", "Cx1", "Cx2", "D", "P4")))))),
           label=rownames(tlCounts), to="Desktop/whoopsgraduateschool archive/20150915_singleCellMAplotsForGenomeBrowser/16cell_ABvP1.pdf" )


