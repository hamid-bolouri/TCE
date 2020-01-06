setwd("~/TCE")
pos2Naive <- get(load("CXCR5pos2naive.RData"))
tex2Naive <- get(load("tex2Naive.RData"))
acute2Naive <- get(load("acute2Naive.RData"))
tex2Acute <- get(load("tex2Acute.RData"))
pos2neg <- get(load("CXCR5pos2neg.RData"))
neg2Naive <- get(load("CXCR5neg2naive.RData"))

#### NB T-bet (TBX21) is NOT in the above tables ####

setwd("~/TCE/CyJS/July2017")
edges <- read.delim("edges.txt", sep="\t", header=TRUE, as.is=TRUE)
nodes <- sort(unique(c(edges[ ,1], edges[ ,2])))
nodes <- setdiff(nodes, c("Adhesion", "CA", "Exhaustion", 
				"interferons", "proliferation", "ZN", "glycolysis"))


mkTbl <- function(expRatiosTbl) {
	indx <- match(nodes, rownames(expRatiosTbl))

	# nodes[which(is.na(indx))]
	 # [1] "CD160"      "CD3D"       "CD3E"       "CD8A"       "glycolysis"
	 # [6] "IL12RB2"    "LEF1"       "MAP3"       "SLC2A2"     "SLC2A4"
	# [11] "SLC2A5"     "TBX21"      "TIGIT"      "TNFRSF9"    "ZAP70"	
	# # Confirmed that these genes are NOT in the expression data.

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
		posNeg <- edgeTbl[i, "edgeType"]
		
		if (posNeg == "promotes") {
			score <- abs(log(fcS * fcT))
			dirMatch <- length(which(log(fcS)*log(fcT) > 0))
		} else {
			if (posNeg == "inhibits") score <- abs(log(fcS / fcT))
			dirMatch <- length(which(log(fcS)*log(fcT) < 0))
		}
		
		fracGd <- dirMatch / N
		
		res <- c(Source=SS, Target=TT, 
				 score=paste(round(score, 3), sep=","), 
				 fracGd=round(fracGd, 3))
		edgeScores <- rbind(edgeScores, res)
	}
	return(edgeScores)
}

pos2NaiveEdgeScores <- mkTbl(pos2Naive)
neg2NaiveEdgeScores <- mkTbl(neg2Naive)
tex2NaiveEdgeScores <- mkTbl(tex2Naive)
acute2NaiveEdgeScores <- mkTbl(acute2Naive)
tex2AcuteEdgeScores <- mkTbl(tex2Acute)
pos2NegEdgeScores <- mkTbl(pos2neg)

litAnnotTbl <- edges

i <- grep("edgeScore", colnames(litAnnotTbl))
if (length(i) > 0) litAnnotTbl <- litAnnotTbl[ ,-i]

getScores <- function(edgeScores) {
	scoreOut <- litAnnotTbl[ ,c("Source", "Target")]
	scoreOut$score <- 0; IST = NULL
	for (i in 1:dim(edgeScores)[1]) {
		iS <- match(scoreOut$Source, edgeScores[i , "Source"])
		iT <- match(scoreOut$Target, edgeScores[i , "Target"])
		iS <- which(!is.na(iS))
		iT <- which(!is.na(iT))
		iST <- intersect(iS, iT)
		if (length(iST) > 0) {
			scoreOut$score[iST] <- edgeScores[i, "fracGd"]
			# print(paste0(i, "   ", iST, "   ", edgeScores[i, "fracGd"]))
		} else print(paste0("Error   ", i))
	}
	return(scoreOut)
}

pos2NaiveScores <- getScores(pos2NaiveEdgeScores)
neg2NaiveScores <- getScores(neg2NaiveEdgeScores)
tex2NaiveScores <- getScores(tex2NaiveEdgeScores)
acute2NaiveScores <- getScores(acute2NaiveEdgeScores)
tex2AcuteScores <- getScores(tex2AcuteEdgeScores)
pos2NegScores <- getScores(pos2NegEdgeScores)

litAnnotTbl$pos2NaiveScores <- pos2NaiveScores$score
litAnnotTbl$neg2NaiveScores <- neg2NaiveScores$score
litAnnotTbl$tex2NaiveScores <- tex2NaiveScores$score
litAnnotTbl$acute2NaiveScores <- acute2NaiveScores$score
litAnnotTbl$tex2AcuteScores <- tex2AcuteScores$score
litAnnotTbl$pos2NegScores <- pos2NegScores$score

setwd("~/TCE/CyJS/July2017")

write.table(litAnnotTbl, file="cytoscapeEdgeAnnot_11July2017.txt", sep="\t",
			col.names=TRUE, row.names=FALSE, quote=FALSE)






















