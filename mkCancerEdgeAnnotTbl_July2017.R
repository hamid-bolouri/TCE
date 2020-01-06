#------------------------ make table ------------------------------------------#

setwd("~/TCE/CyJS/July2017")
# setwd("/Volumes/hbolouri/TCE/CyJS/July2017")
edges <- read.delim("cancerEdgeAnnot_17July2017.txt", 
					sep="\t", header=TRUE, as.is=TRUE)
nodes <- sort(unique(c(edges[ ,1], edges[ ,2])))


calcScores <- function(expRatiosTbl) {
	indx <- match(nodes, rownames(expRatiosTbl))

	selExp <- expRatiosTbl[indx[!is.na(indx)], ]

	i1 <- match(edges[ ,1], rownames(selExp))
	i2 <- match(edges[ ,2], rownames(selExp))
	selEdges <- edges[intersect(which(!is.na(i1)),which(!is.na(i2))), ]

	edgeTbl <- selEdges[ ,c("Source", "Target", "edgeType")]

	N <- dim(expRatiosTbl)[2]
	edgeScores = NULL
	for (i in 1:dim(edgeTbl)[1]) {
		SS <- edgeTbl[i, "Source"]
		TT <- edgeTbl[i, "Target"]
		
		fcS <- expRatiosTbl[match(SS, rownames(expRatiosTbl)), ]
		fcT <- expRatiosTbl[match(TT, rownames(expRatiosTbl)), ]
		edgeDir <- edgeTbl[i, "edgeType"]
		
		if (edgeDir == "promotes") {
			dirMatch <- length(which(round(log(fcS), 2)*round(log(fcT), 2) >= 0)) # both high or both low
			misMatch <- length(which(round(log(fcS), 2)*round(log(fcT), 2) < 0))
		} else {
			dirMatch <- length(which(round(log(fcS), 2)*round(log(fcT), 2) <= 0))
			misMatch <- length(which(round(log(fcS), 2)*round(log(fcT), 2) > 0))
		}
		
		score <- dirMatch / (dirMatch + misMatch)
		
		res <- c(Source=SS, Target=TT, 
				 score=round(score, 3))
		edgeScores <- rbind(edgeScores, res)
	}
	return(edgeScores)
}

setwd("~/TCE/_Data_cancer")
# setwd("/Volumes/hbolouri/TCE/_Data_cancer")
aPDL1.v.TIL.GSE93007 <- get(load("aPDL1.v.TIL.GSE93007.RData"))
PDL1plus.v.PDL1minus.GSE39205 <- get(load("PDL1plus.v.PDL1minus.GSE39205.RData"))
Tim3Plus2Tim3Minus.GSE84072 <- get(load("Tim3Plus2Tim3Minus.GSE84072.RData"))
TILN2CMV.GSE24536 <- get(load("TILN2CMV.GSE24536.RData"))
TILN2naive.GSE24536 <- get(load("TILN2naive.GSE24536.RData"))
TILN2EBV.GSE24536 <- get(load("TILN2EBV.GSE24536.RData"))
TILN2tumorPBMC.GSE24536 <- get(load("TILN2tumorPBMC.GSE24536.RData"))

# Scheitinger data:
acuteD5.2N.89307 <- get(load("acuteD5.2N.89307.RData"))
acuteD7.2N.89307 <- get(load("acuteD7.2N.89307.RData"))
T.D5.2N.89307 <- get(load("T.D5.2N.89307.RData"))
T.D7.2N.89307 <- get(load("T.D7.2N.89307.RData"))
T.D14.2N.89307 <- get(load("T.D14.2N.89307.RData"))
T.D21.2N.89307 <- get(load("T.D21.2N.89307.RData"))

aPDL1.v.TIL.GSE93007.scores <- calcScores(aPDL1.v.TIL.GSE93007)
PDL1plus.v.PDL1minus.GSE39205.scores <- calcScores(PDL1plus.v.PDL1minus.GSE39205)
Tim3Plus2Tim3Minus.GSE84072.scores <- calcScores(Tim3Plus2Tim3Minus.GSE84072)
TILN2CMV.GSE24536.scores <- calcScores(TILN2CMV.GSE24536)
TILN2naive.GSE24536.scores <- calcScores(TILN2naive.GSE24536)
TILN2EBV.GSE24536.scores <- calcScores(TILN2EBV.GSE24536)
TILN2tumorPBMC.GSE24536.scores <- calcScores(TILN2tumorPBMC.GSE24536)

acuteD5.2N.89307.scores <- calcScores(acuteD5.2N.89307)
acuteD7.2N.89307.scores <- calcScores(acuteD7.2N.89307)
T.D5.2N.89307.scores <- calcScores(T.D5.2N.89307)
T.D7.2N.89307.scores <- calcScores(T.D7.2N.89307)
T.D14.2N.89307.scores <- calcScores(T.D14.2N.89307)
T.D21.2N.89307.scores <- calcScores(T.D21.2N.89307)


litAnnotTbl <- edges[ ,1:4]

i <- grep("edgeScore", colnames(litAnnotTbl))
if (length(i) > 0) litAnnotTbl <- litAnnotTbl[ ,-i]

mkScoreTbl <- function(edgeScores) {
	scoreOut <- litAnnotTbl[ ,c("Source", "Target")]
	scoreOut$score <- 0; IST = NULL
	for (i in 1:dim(edgeScores)[1]) {
		iS <- match(scoreOut$Source, edgeScores[i , "Source"])
		iT <- match(scoreOut$Target, edgeScores[i , "Target"])
		iS <- which(!is.na(iS))
		iT <- which(!is.na(iT))
		iST <- intersect(iS, iT)
		if (length(iST) > 0) {
			scoreOut$score[iST] <- edgeScores[i, "score"]
			# print(paste0(i, "   ", iST, "   ", edgeScores[i, "score"]))
		} else print(paste0("Error   ", i))
	}
	return(scoreOut)
}

aPDL1TILScores <- mkScoreTbl(aPDL1.v.TIL.GSE93007.scores)
PDL1plusPDL1minusScores <- mkScoreTbl(PDL1plus.v.PDL1minus.GSE39205.scores)
Tim3Plus2Tim3MinusScores <- mkScoreTbl(Tim3Plus2Tim3Minus.GSE84072.scores)
TILN2CMVScores <- mkScoreTbl(TILN2CMV.GSE24536.scores)
TILN2naiveScores <- mkScoreTbl(TILN2naive.GSE24536.scores)
TILN2EBVScores <- mkScoreTbl(TILN2EBV.GSE24536.scores)
TILN2tumorPBMCScores <- mkScoreTbl(TILN2tumorPBMC.GSE24536.scores)

acuteD5.2NScores <- mkScoreTbl(acuteD5.2N.89307.scores)
acuteD7.2NScores <- mkScoreTbl(acuteD7.2N.89307.scores)
T.D5.2NScores <- mkScoreTbl(T.D5.2N.89307.scores)
T.D7.2NScores <- mkScoreTbl(T.D7.2N.89307.scores)
T.D14.2NScores <- mkScoreTbl(T.D14.2N.89307.scores)
T.D21.2NScores <- mkScoreTbl(T.D21.2N.89307.scores)


litAnnotTbl$aPDL1TILScores <- aPDL1TILScores$score
litAnnotTbl$PDL1plusPDL1minusScores <- PDL1plusPDL1minusScores$score
litAnnotTbl$Tim3Plus2Tim3MinusScores <- Tim3Plus2Tim3MinusScores$score
litAnnotTbl$TILN2CMVScores <- TILN2CMVScores$score
litAnnotTbl$TILN2naiveScores <- TILN2naiveScores$score
litAnnotTbl$TILN2EBVScores <- TILN2EBVScores$score
litAnnotTbl$TILN2tumorPBMCScores <- TILN2tumorPBMCScores$score

litAnnotTbl$acuteD5.2NScores <- acuteD5.2NScores$score
litAnnotTbl$acuteD7.2NScores <- acuteD7.2NScores$score
litAnnotTbl$T.D5.2NScores <- T.D5.2NScores$score
litAnnotTbl$T.D7.2NScores <- T.D7.2NScores$score
litAnnotTbl$T.D14.2NScores <- T.D14.2NScores$score
litAnnotTbl$T.D21.2NScores <- T.D21.2NScores$score

setwd("~/TCE/CyJS/July2017")
# setwd("/Volumes/hbolouri/TCE/CyJS/July2017")

save(litAnnotTbl, file="cancerEdgeAnnot_17July2017.RData")
write.table(litAnnotTbl, file="cancerEdgeAnnot_17July2017.txt", sep="\t",
			col.names=TRUE, row.names=FALSE, quote=FALSE)














############################# NOT USE ############################################

# #------------------------- scale fold changes to {-1, +1} ---------------------#

# scaleFC <- function(fcFile) {
	# MAX <- max(unlist(fcFile))
	# MIN <- min(unlist(fcFile))
	# scaled <- apply(fcFile, 2, function(y) {
					# y[y > 1] <- (y[y > 1] - 1)/(MAX - 1)
					# # y[y <= 1] <- y[y <= 1] - 1 
					# } )
	# # scaled <- scaled - 1
	# return(scaled)
# }

# setwd("~/TCE/_Data_cancer")

# aPDL1.v.TIL.GSE93007 <- get(load("aPDL1.v.TIL.GSE93007.RData"))
# aPDL1.v.TIL.GSE93007 <- scaleFC(aPDL1.v.TIL.GSE93007)


# PDL1plus.v.PDL1minus.GSE39205 <- get(load("PDL1plus.v.PDL1minus.GSE39205.RData"))
# PDL1plus.v.PDL1minus.GSE39205 <- scaleFC(PDL1plus.v.PDL1minus.GSE39205)


# Tim3Plus2Tim3Minus.GSE84072 <- get(load("Tim3Plus2Tim3Minus.GSE84072.RData"))
# Tim3Plus2Tim3Minus.GSE84072 <- scaleFC(Tim3Plus2Tim3Minus.GSE84072)


# TILN2CMV.GSE24536 <- get(load("TILN2CMV.GSE24536.RData"))
# TILN2CMV.GSE24536 <- scaleFC(TILN2CMV.GSE24536)


# TILN2naive.GSE24536 <- get(load("TILN2naive.GSE24536.RData"))
# TILN2naive.GSE24536 <- scaleFC(TILN2naive.GSE24536)


# TILN2EBV.GSE24536 <- get(load("TILN2EBV.GSE24536.RData"))
# TILN2EBV.GSE24536 <- scaleFC(TILN2EBV.GSE24536)


# TILN2tumorPBMC.GSE24536 <- get(load("TILN2tumorPBMC.GSE24536.RData"))
# TILN2tumorPBMC.GSE24536 <- scaleFC(TILN2tumorPBMC.GSE24536)


