#------------------------------------------------------------------------------#
#							raw per gene exp data							   #	
#------------------------------------------------------------------------------#

setwd("~/TCE/GSE89307_Scheitinger2017")

geneExp <- get(load("Scheitinger2017_exp.RData"))

M <- rowMeans(geneExp)													# 1570 all zero genes
geneExp <- geneExp[which(M != 0), ]										# 18584 x 38

colnames(geneExp) <- gsub("^L", "TD", colnames(geneExp))


#------------------------------------------------------------------------------#
#							4K most var genes								   #	
#------------------------------------------------------------------------------#

SDs <- apply(log2(geneExp), 1, sd)
topSDs <- names(sort(SDs, decreasing=TRUE)[1:4000])
top4kExp <- geneExp[match(topSDs, rownames(geneExp)), ]


#------------------------------------------------------------------------------#
#								GSE41867 geneSet							   #	
#------------------------------------------------------------------------------#

setwd("~/TCE/_Data_infection/GSE41867")
markersGSE41867 <- read.delim("module1.3.7.8.11.12.MarkerGenesFromPaper.txt",
					  sep="\t", header=FALSE, as.is=TRUE)$V1
markersGSE41867 <- toupper(markersGSE41867)
shared <- intersect(markersGSE41867, rownames(geneExp))
markersGSE41867Exp <- geneExp[shared, ]									# 141 x 26


#------------------------------------------------------------------------------#
#							Schietinger geneSets							   #	
#------------------------------------------------------------------------------#

peaks.L5.v.E5 <- c("FOS", "SMARCC1", "NFATC1", "NFAT5", "FOSB", "JUND", "NFATC2", 
				   "NFATC3", "BATF", "NRF4A2", "JUNB", "BACH1", "PPARGC1", "ESRRA", 
				   "NFE2", "RREB1", "E2F4", "KLF16", "ETS1", "GABPA", "E2F3", 
				   "MAZ", "KLF7", "MBD2", "ZFX", "ZBTB7A", "KLF6", "ELF2", "ELF4", 
				   "NR2C1", "NFE2L1", "ZFP148", "SP1", "KLF", "ELK", "E2F1", 
				   "FLI1", "ETV3", "SP3")

exp.L5.v.E5 <- c("S1PR1", "ARL4C", "SPICE1", "ADRB2", "RASA3", "CARD11", "ITGAE", 
				 "CTLA4", "SPRY1", "MTX2", "ENDOD1", "GAS2", "CHST2", "HIVEP3", 
				 "EYA2", "SESN1","NFKB1")

upDnL14 <- c("IRF4", "ST6GAL1", "DUSP4", "CCR7", "MAK11", "LEF1", "CD200", "CCR4", "BTLA", 
		   "POU2AF1", "LIF", "TCF7", "TNFSF11", "SELL", "CD40LG", "ID3", "SAMF6", "NR4A3",
		   "EGR2", "ARMCX2", "SATB1", "TNFS4", "CNN3", "TNFSF13B", "MBL2", "IL1A", "CD86",
		   "IDO1", "CD86", "CCRL2", "ENTPD1", "PRDM1", "ADAM8", "MCAM", "CD244", "SOCS2",
		   "CAR2", "GZMA", "CD200R1", "CCL5", "CD7", "KLRB1", "ADRB2", "ITGAE", "BCAS3", 
		   "EXT1", "PXYLP1")

MarkersSchiet <- unique(c(peaks.L5.v.E5, exp.L5.v.E5, upDnL14)) 		# 100 genes
shared <- intersect(MarkersSchiet, rownames(geneExp))
markersSchietExp <- geneExp[shared, ]									# 84 x 26


#------------------------------------------------------------------------------#
#							TCE network geneSet								   #	
#------------------------------------------------------------------------------#

setwd("~/TCE/CyJS/July2017")
edges <- get(load("cancerEdgeAnnot_17July2017.RData"))					# 505
Q <- "Adhesion|proliferation|glycolysis|ZN|CA"
i1 <- !grepl(Q, edges[ , 1])
i2 <- !grepl(Q, edges[ , 2])
edges <- edges[(i1 & i2), ]												# 478

nodes <- unique(c(edges[ , 1], edges[ , 2]))
nodesInGeneExp <- intersect(nodes, rownames(geneExp))					# 108 AOK

tceExp <- geneExp[nodesInGeneExp, ]


#------------------------------------------------------------------------------#
#							all markers geneSet								   #	
#------------------------------------------------------------------------------#

allMarkers <- unique(c(nodesInGeneExp, MarkersSchiet, markersGSE41867))
shared <- intersect(allMarkers, rownames(geneExp))						# 323
allMarkerExp <- geneExp[shared, ]

#------------------------------------------------------------------------------#
#							get fold changes								   #	
#------------------------------------------------------------------------------#

setwd("~/TCE/_Data_cancer")
T.D5.2N.89307 <- get(load("T.D5.2N.89307.RData"))
T.D7.2N.89307 <- get(load("T.D7.2N.89307.RData"))
T.D14.2N.89307 <- get(load("T.D14.2N.89307.RData"))
T.D21.2N.89307 <- get(load("T.D21.2N.89307.RData"))
T.D28.2N.89307 <- get(load("T.D28.2N.89307.RData"))
T.D35.2N.89307 <- get(load("T.D35.2N.89307.RData"))
T.D60.2N.89307 <- get(load("T.D60.2N.89307.RData"))

allFC <- cbind(T.D5.2N.89307, T.D7.2N.89307, T.D14.2N.89307,
			   T.D21.2N.89307, T.D28.2N.89307, T.D35.2N.89307,
			   T.D60.2N.89307)
colnames(allFC) <- gsub(".2N.89307", "", colnames(allFC))
colnames(allFC) <- gsub("T.D", "TD", colnames(allFC))


#------------------------------------------------------------------------------#
#							calculate edge scores							   #	
#------------------------------------------------------------------------------#

mkTbl <- function(expRatiosTbl, edgesTbl) {

	N <- dim(expRatiosTbl)[2]
	edgeScores = NULL
	
	SS <- edgesTbl[ , "Source"]
	TT <- edgesTbl[ , "Target"]
		
	fcTbl = fcRownames = NULL
	for (i in 1:length(SS)) {
		iS <- match(rownames(expRatiosTbl), SS[i])
		iS <- which(!is.na(iS))
		
		iT <- match(rownames(expRatiosTbl), TT[i])
		iT <- which(!is.na(iT))
		
		posNeg <- ifelse(edgesTbl[i , "edgeType"] == "promotes", 1, -1)
		
		if (length(iS) > 0 & length(iT) > 0) {
			fcS <- expRatiosTbl[iS, ]
			fcT <- expRatiosTbl[iT, ]
			fc <- c(fcS, fcT, posNeg)
			fcTbl <- rbind(fcTbl, fc)
			fcRownames <- c(fcRownames, 
							paste(SS[i], TT[i], posNeg, sep=":"))
		}
	}
	rownames(fcTbl) <- fcRownames
	
	edgeScores <- apply(fcTbl, 1, function(x) {
						log(x[1:N]) * log(x[(N+1):(2*N)]) * x[2*N + 1] })
	colnames(edgeScores) <- fcRownames
	
	return(t(edgeScores))
}

rawScores <- mkTbl(allFC, edges)
asinhScore <- asinh(rawScores)

# boxplot(rawScores, las=2, cex.axis=0.4, cex=0.5, col=rgb(0,0.5,1,0.5), 
		# outcol=rgb(0,0,0,0.2), border="blue", lwd =0.5)
# dev.off()


#------------------------------------------------------------------------------#
#						Cosine Edge Concordance similarity					   #
#------------------------------------------------------------------------------#

setwd("~/TCE")

colCosDist <- function(TBL) {

	cosSim <- function(vector1, vector2) {
		sim <- sum(vector1 * vector2) # / 
			   ( (sum(vector1^2)^0.5) * (sum(vector2^2)^0.5) )
		return(sim)
	}

	getSims <- function(y, mtrx) {
		simTbl <- apply(mtrx, 2, function(x) cosSim(x, y))
		return(simTbl)
	}

	cosSimsTbl <- apply(TBL, 2, function(y) getSims(y, TBL))
	normSims <- cosSimsTbl
	# if (min(normSims) < 0) {
		normSims <- normSims - min(normSims)
		normSims <- normSims / max(normSims)
	# }

	D <- max(normSims) - normSims
	diag(D) <- 0
	return(D)
}

cosDistRaw <- colCosDist(rawScores)


#------------------------------------------------------------------------------#
#						Boolean Edge Concordance similarity					   #
#------------------------------------------------------------------------------#

boolScores <- apply(rawScores, 2, function(y) sign(y))

grps <- sapply(colnames(boolScores), function(x) strsplit(x, "[.]"))
grps <- unlist(lapply(grps, function(x) x[[1]]))

grps <- split(colnames(boolScores), grps)
grps <- grps[paste0("TD", sort(as.numeric(gsub("TD", "", names(grps)))))]

avBoolTbl = NULL
for (G in grps) {
	avBool <- apply(boolScores[ , G], 1, function(x) sum(x)/length(x))
	avBoolTbl <- cbind(avBoolTbl, avBool)
}
colnames(avBoolTbl) <- names(grps)

cosSim <- function(vector1, vector2) {
	sim <- sum(vector1 * vector2) /
		   ( (sum(vector1^2)^0.5) * (sum(vector2^2)^0.5) )
	return(sim)
}

getSims <- function(y, mtrx) {
	simTbl <- apply(mtrx, 2, function(x) cosSim(x, y))
	return(simTbl)
}

indx <- which(abs(rowSums(avBoolTbl)) != ncol(avBoolTbl))
avBoolTbl <- avBoolTbl[indx, ]

cosSimsBool <- apply(avBoolTbl, 2, function(y) getSims(y, avBoolTbl))
cosDistBool <- round(max(cosSimsBool) - cosSimsBool, 5) # stop dist(self)==1E-16

# ZERO edges that are perfectly con/dis-cordant at all times.
# names(which(rowSums(avBoolTbl)== 7))
# names(which(rowSums(avBoolTbl)== -7 )) 


#------------------------------------------------------------------------------#
#									plot									   #
#------------------------------------------------------------------------------#

setwd("~/TCE")

expID <- "GSE89307"
if (expID == "GSE89307") {
	cellTypes <- sapply(colnames(allFC), function(x) strsplit(x, "[.]"))
	cellTypes <- c(unlist(unique(lapply(cellTypes, 
				   function(x) x[1]))), "N")

	A <- 0.75
	library(dichromat)
	myBlues <- adjustcolor(colorschemes$LightBluetoDarkBlue.10, alpha.f=A)[4:10]
	myColors <- c(myBlues, adjustcolor("black", alpha.f=A))	
	names(myColors) <- cellTypes
}


mkPlot <- function(dataType="geneExp") {
	plotData <- get(dataType)

	if (grepl("Dist", dataType, ignore.case=TRUE)) {
		cellTypes <- cellTypes[-(length(cellTypes))]
		
		XY <- cmdscale(plotData, k=2) 
		colnames(XY) <- c("MDS1", "MDS2")
		fileTxt <- "MDS"
		propVar <- "NA"

	} else {
		PCA <- prcomp(t(plotData), center=TRUE, scale=TRUE)

		propVar <- summary(PCA)$importance[2, 1:2]
		XY <- PCA$x
		fileTxt <- "PCA"

		if (!grepl("exp", dataType, ignore.case=TRUE)) {
			cellTypes <- cellTypes[-(length(cellTypes))]	
			PCA <- prcomp(t(plotData), center=TRUE, scale=TRUE)
			XY <- PCA$x
		} 
	}
	
	XY <- XY[grep(paste(cellTypes, collapse="|"), rownames(XY)), 1:2]

	centroids = NULL
	for (Q in cellTypes) {
		indx <- grep(Q, rownames(XY))
		centroidX <- median(XY[indx, 1])
		centroidY <- median(XY[indx, 2])
		centroids <- rbind(centroids, c(x=centroidX, y=centroidY))
	}
	rownames(centroids) <- cellTypes

	myColors <- myColors[cellTypes]

	setwd("~/TCE")
	pdf(paste0(expID, "_", dataType, "_", fileTxt, ".pdf"))
		S = 0.75

		if (length(propVar) == 2) {
			propVar <- sapply(propVar, function(x) round(x, 2))
			txt <- paste0("Proportions of Variance=", 
						  round(propVar[1], 2), ", ", 
						  round(propVar[2], 2))
		} else txt=""
		
		plot(XY, xaxt='n', yaxt='n', cex=S, col="gray60",
			 main=txt, cex.main=0.75)

		for (i in 1:length(cellTypes)) {
			indx <- grep(cellTypes[i], rownames(XY))
			points(XY[indx, 1], XY[indx, 2], pch=19, col=myColors[i], cex=S)
		}

		for (i in cellTypes) {
			points(centroids[i, 1], centroids[i, 2], 
			col=myColors[i], cex=2.5, lwd=4)
		}

		legend("topleft", legend=cellTypes, 
				pch=19, col=myColors, bty="n", cex=0.7)
	dev.off()
}	

# "geneExp" 			# "top4kExp"		# tceExp		
# "markersGSE41867Exp" 	# markersSchietExp	# allMarkerExp
# "rawScores" 			# asinhScore
# "cosDistRaw"			# cosDistBool	# avBoolTbl

mkPlot("rawScores")
mkPlot("asinhScore")
mkPlot("avBoolTbl")
mkPlot("cosDistBool")
mkPlot("cosDistRaw")
mkPlot("markersGSE41867Exp")
mkPlot("markersSchietExp")
mkPlot("tceExp")
mkPlot("allMarkerExp")
mkPlot("top4kExp")
mkPlot("geneExp")


#------------------------------------------------------------------------------#
#						heatmap allmarkers -> select marker					   #
#------------------------------------------------------------------------------#

tbl <- log2(geneExp[intersect(allMarkers, rownames(geneExp)), 
			   grep("^N|^T", colnames(geneExp))] + 1)

conditions <- gsub("N1|N2|N3", "N", colnames(tbl))
cellTypes <- unique(conditions)

library(dichromat)
aColors <- c("black", colorschemes$LightBluetoDarkBlue.10[-(1:3)])
names(aColors) <- cellTypes

DR <- dist(tbl, method="euclidean")
hr <- hclust(DR, method="ward.D2")

DC <- dist(t(tbl), method="euclidean")
hc <- hclust(DC, method="ward.D2")

library(gplots)
library("dichromat")
myColors <- colorschemes$DarkRedtoBlue.18

colAnnot <- conditions
for (i in 1:length(cellTypes)) colAnnot[colAnnot == cellTypes[i]] <- aColors[i]

setwd("/home/hbolouri/TCE")
pdf('GSE89307_allMarkersHeatmap2.pdf')
par(oma=c(2,0,0,0))
par(family="mono")
H = heatmap.2(tbl, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc),
	scale="row",keysize=1,symkey=FALSE,
	key.title="",key.xlab="z(log2(expr))",key.ylab="",na.rm=TRUE,
	ylab="",xlab="",trace="none",col=myColors,
	cexRow=0.1, cexCol=0.6, ColSideColors=colAnnot) # 
dev.off()











































