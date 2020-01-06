setwd("~/TCE/_Data_cancer")
acuteD5.2N.89307 <- get(load("acuteD5.2N.89307.RData"))
acuteD7.2N.89307 <- get(load("acuteD7.2N.89307.RData"))
T.D5.2N.89307 <- get(load("T.D5.2N.89307.RData"))
T.D7.2N.89307 <- get(load("T.D7.2N.89307.RData"))
T.D14.2N.89307 <- get(load("T.D14.2N.89307.RData"))
T.D21.2N.89307 <- get(load("T.D21.2N.89307.RData"))


setwd("~/TCE/CyJS/July2017")
edges <- get(load("cancerEdgeAnnot_17July2017.RData"))					# 505
i1 <- !grepl("proliferation|glycolysis|ZN|CA", edges[ , 1])
i2 <- !grepl("Adhesion|proliferation|glycolysis|ZN|CA", edges[ , 2])
edges <- edges[(i1 & i2), ]												# 478

nodes <- unique(c(edges[ , 1], edges[ , 2]))
nodes <- intersect(nodes, rownames(acuteD5.2N.89307))					# 108


mkTbl <- function(expRatiosTbl) {

	N <- dim(expRatiosTbl)[2]
	edgeScores = NULL
	
	SS <- edges[ , "Source"]
	TT <- edges[ , "Target"]
		
	fcTbl = fcRownames = NULL
	for (i in 1:length(SS)) {
		iS <- match(rownames(expRatiosTbl), SS[i])
		iS <- which(!is.na(iS))
		
		iT <- match(rownames(expRatiosTbl), TT[i])
		iT <- which(!is.na(iT))
		
		posNeg <- ifelse(edges[i , "edgeType"] == "promotes", 1, -1)
		
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


acuteD5.2N.89307.edgeScores <- mkTbl(acuteD5.2N.89307)
acuteD7.2N.89307.edgeScores <- mkTbl(acuteD7.2N.89307)
T.D5.2N.89307.edgeScores <- mkTbl(T.D5.2N.89307)
T.D7.2N.89307.edgeScores <- mkTbl(T.D7.2N.89307)
T.D14.2N.89307.edgeScores <- mkTbl(T.D14.2N.89307)
T.D21.2N.89307.edgeScores <- mkTbl(T.D21.2N.89307)

allScores  <- cbind(T.D5.2N.89307.edgeScores, T.D7.2N.89307.edgeScores,
					T.D14.2N.89307.edgeScores, T.D21.2N.89307.edgeScores,
					acuteD5.2N.89307.edgeScores, acuteD7.2N.89307.edgeScores)

scoreTbl <- allScores
scoreTbl[scoreTbl > 0] <- 1
scoreTbl[scoreTbl < 0] <- 0

library(pheatmap)
setwd("/home/hbolouri/TCE")
pheatmap(scoreTbl, border_color=NA, scale="none", 
		 fontsize_row=1, fontsize_col=5, 
		 clustering_method="ward.D2", legend=TRUE,
		 filename="EdgeConcordanceOverTime_25May2018.pdf")
dev.off()

length(which(apply(scoreTbl, 1, function(x) any(x == 1)))) # 461 of 478

setdiff(rownames(scoreTbl), 
		rownames(scoreTbl)[which(apply(scoreTbl, 1, function(x) any(x == 1)))])

 [1] AKT1:MTOR:1      BATF:IRF4:-1     CD160:IFNG:-1    CD160:IL2:-1
 [5] CTLA4:AKT1:-1    DNMT3A:DNMT3A:-1 FOXO1:PDCD1:1    IFNG:SLC2A4:-1
 [9] LEF1:EOMES:1     TBX21:CD160:-1   TBX21:LAG3:-1    TBX21:PDCD1:-1
[13] TCF7:EOMES:1     TCF7:IFNGR2:-1   TET1:PDCD1:1     TET3:PDCD1:1
[17] ZEB2:IL2:-1


#------------------------------------------------------------------------------#
#						edge scores of ranodomized network					   #
#------------------------------------------------------------------------------#

#### randEdgesTbl:
# Source = Target = NULL
# for (i in rep(nodes, 1000)) {
	# Source <- c(Source, i)
	# Target <- c(Target, sample(setdiff(nodes, i), 1))
# }

# edgeType <- sample(edges$edgeType, size=length(Source), replace=TRUE)
# randEdgesTbl <- cbind(Source, Target, edgeType)

setwd("/home/hbolouri/TCE")
# save(randEdgesTbl, file="randEges100K.RData")
randEdges <- get(load("randEges100K.RData"))


setwd("~/TCE/CyJS/July2017")
edges <- get(load("cancerEdgeAnnot_17July2017.RData"))					# 505
i1 <- !grepl("proliferation|glycolysis|ZN|CA", edges[ , 1])
i2 <- !grepl("Adhesion|proliferation|glycolysis|ZN|CA", edges[ , 2])
edgesTCE <- edges[(i1 & i2), ]											# 478

nodes <- unique(c(edges[ , 1], edges[ , 2]))
nodes <- intersect(nodes, rownames(acuteD5.2L60.89307))					# 113


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

acuteD5.2N.89307.scoresRand <- mkTbl(acuteD5.2N.89307, randEdges)
acuteD7.2N.89307.scoresRand <- mkTbl(acuteD7.2N.89307, randEdges)
L5.2N.89307.scoresRand <- mkTbl(T.D5.2N.89307, randEdges)
L7.2N.89307.scoresRand <- mkTbl(T.D7.2N.89307, randEdges)
L14.2N.89307.scoresRand <- mkTbl(T.D14.2N.89307, randEdges)
L21.2N.89307.scoresRand <- mkTbl(T.D21.2N.89307, randEdges)


scoresRand.89307  <- cbind(L5.2N.89307.scoresRand, L7.2N.89307.scoresRand,
					L14.2N.89307.scoresRand, L21.2N.89307.scoresRand,
					# L35.2N.89307.scoresRand, 
					acuteD5.2N.89307.scoresRand, 
					acuteD7.2N.89307.scoresRand)

setwd("~/TCE")
save(scoresRand.89307, file="scoresRand_wrt2N.89307.RData")


randScoreTbl <- scoresRand.89307
randScoreTbl[randScoreTbl > 0] <- 1
randScoreTbl[randScoreTbl < 0] <- 0

length(which(apply(randScoreTbl, 1, function(x) any(x == 1)))) # 102369 of 108000


colnames(scoreTbl) <- strtrim(gsub("T.|acuteD", "", colnames(scoreTbl)), 3)
res = NULL
for (j in unique(colnames(scoreTbl))) {
	tmpTbl <- scoreTbl[ ,grep(j, colnames(scoreTbl))]
	res <- c(res, unlist(apply(tmpTbl, 1, function(x) sum(x)/length(x))))
}
length(which(res == 1)) # == 741/2868 ~ 26%

# table(colnames(randScoreTbl))
# 5.2 7.2 D14 D21 D5. D7.
  # 9   9   9   9   9   9
(0.5^9) # 0.195% ~ 0.2%
(0.5^3) # 0.125


colnames(randScoreTbl) <- strtrim(gsub("T.|acuteD", "", colnames(randScoreTbl)), 3)
resRand = NULL
for (j in unique(colnames(randScoreTbl))) {
	tmpTbl <- randScoreTbl[ ,grep(j, colnames(randScoreTbl))]
	resRand <- c(resRand, unlist(apply(tmpTbl, 1, function(x) sum(x)/length(x))))
}

H.TCE = hist(res, breaks=50, col="blue", border="white", freq=FALSE)
H.rand = hist(resRand, breaks=50, col=rgb(1,0,0,0.5), 
			  border="white", add=TRUE, freq=FALSE)
dev.off()


library(pheatmap)
setwd("/home/hbolouri/TCE")
pheatmap(randScoreTbl, border_color=NA, scale="none", 
		 fontsize_row=1, fontsize_col=5, 
		 clustering_method="ward.D2", legend=TRUE,
		 cluster_rows=FALSE, cluster_cols=FALSE,
		 filename="EdgeConcordanceRan100K_25May2018.pdf")
dev.off()

length(which(resRand > 0.5)) # 







