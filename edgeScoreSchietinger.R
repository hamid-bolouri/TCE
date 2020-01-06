setwd("~/TCE/_Data_cancer") 
acuteD5.2N.89307 <- get(load("acuteD5.2N.89307.RData"))
acuteD7.2N.89307 <- get(load("acuteD7.2N.89307.RData"))
T.D5.2N.89307 <- get(load("T.D5.2N.89307.RData"))
T.D7.2N.89307 <- get(load("T.D7.2N.89307.RData"))
T.D14.2N.89307 <- get(load("T.D14.2N.89307.RData"))
T.D21.2N.89307 <- get(load("T.D21.2N.89307.RData"))

setwd("~/TCE/CyJS/July2017")
edges <- get(load("cancerEdgeAnnot_17July2017.RData"))

nodes <- unique(c(edges[ , 1], edges[ , 2]))
nodes <- intersect(nodes, rownames(acuteD5.2N.89307))					# 113 AOK

#### Using only genes in our network, even after shuffling, ~ half edges get score == 1
# Source = Target = NULL
# for (i in rep(nodes, 1000)) {
	# Source <- c(Source, i)
	# Target <- c(Target, sample(setdiff(nodes, i), 1))
# }

# edgeVal <- runif(length(Source))
# edgeType <- sapply(edgeVal, function(x) {
				   # if (x >= 0.5) "promotes" else "inhibits"})
# edgesTbl <- cbind(Source, Target, edgeType)

# save(edgesTbl, file="SchietingerRandEges100K.RData")

genes <- setdiff(rownames(acuteD5.2N.89307), nodes)

Source <- sample(rownames(acuteD5.2N.89307), 1E5, replace=TRUE)
Target <- sample(rownames(acuteD5.2N.89307), 1E5, replace=TRUE)
posNeg <- sample(c(+1, -1), 1E5, replace=TRUE)
edgesTbl <- cbind(Source=Source, Target=Target, edgeType=posNeg)

library(caret)
mkTbl <- function(expRatiosTbl) {

	getScore <- function(x) {
				score <- log(x[1:N] * x[(N+1):(2*N)])
				Match <- length(which(score * x[2*N + 1] > 0))
				fracGd <- Match / N
				return(fracGd)
				} 

	nearZero <- apply(expRatiosTbl, 1, nearZeroVar)
	zeroI <- which(nearZero == 1)
	expRatiosTbl <- expRatiosTbl[-zeroI, ]
	
	N <- dim(expRatiosTbl)[2]
	edgeScores = NULL

	SS <- edgesTbl[ , "Source"]
	TT <- edgesTbl[ , "Target"]
	
	iS <- match(SS, rownames(expRatiosTbl))
	iT <- match(TT, rownames(expRatiosTbl))
	
	fcS <- expRatiosTbl[iS, ]
	fcT <- expRatiosTbl[iT, ]
	posNeg <- as.numeric(edgesTbl[ , "edgeType"])
	
	fcTbl <- cbind(fcS, fcT, posNeg)

	edgeScores <- apply(fcTbl, 1, function(x) getScore(x))
	return(edgeScores)
}


acuteD5.2N.89307.edgeScores <- mkTbl(acuteD5.2N.89307)
acuteD7.2N.89307.edgeScores <- mkTbl(acuteD7.2N.89307)
T.D5.2N.89307.edgeScores <- mkTbl(T.D5.2N.89307)
T.D7.2N.89307.edgeScores <- mkTbl(T.D7.2N.89307)
T.D14.2N.89307.edgeScores <- mkTbl(T.D14.2N.89307)
T.D21.2N.89307.edgeScores <- mkTbl(T.D21.2N.89307)

scores500K <- c(T.D5.2N.89307.edgeScores, T.D7.2N.89307.edgeScores,
				T.D14.2N.89307.edgeScores, T.D21.2N.89307.edgeScores,
				acuteD5.2N.89307.edgeScores, acuteD7.2N.89307.edgeScores)

H = hist(scores500K, breaks=50, border="white", 
		 col="skyblue", xlab="Edge Concordance Score")
dev.off()

# H$counts
 # [1] 277110      0      0      0      0  22177      0      0      0      0
# [11]      0  30849      0      0      0      0  40880      0      0      0
# [21]      0      0  26785      0      0      0      0  26201      0      0
# [31]      0      0      0  39327      0      0      0      0  27282      0
# [41]      0      0      0      0  20691      0      0      0      0  88698
# >
# H$breaks
 # [1] 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18 0.20 0.22 0.24 0.26 0.28
# [16] 0.30 0.32 0.34 0.36 0.38 0.40 0.42 0.44 0.46 0.48 0.50 0.52 0.54 0.56 0.58
# [31] 0.60 0.62 0.64 0.66 0.68 0.70 0.72 0.74 0.76 0.78 0.80 0.82 0.84 0.86 0.88
# [46] 0.90 0.92 0.94 0.96 0.98 1.00

# FDR:
# 88698/sum(H$counts)	== 0.14783

# 25 May 2018: 
# H$counts[length(H$counts)]/sum(H$counts) = 0.1465183


tex2NaiveEdgeScores <- mkTbl(tex2Naive)
acute2NaiveEdgeScores <- mkTbl(acute2Naive)
tex2AcuteEdgeScores <- mkTbl(tex2Acute)
pos2NegEdgeScores <- mkTbl(pos2neg)

save(pos2NaiveEdgeScores, file="pos2NaiveEdgeScoresRand100K.RData")
save(tex2NaiveEdgeScores, file="pos2NaiveEdgeScoresRand100K.RData")
save(acute2NaiveEdgeScores, file="pos2NaiveEdgeScoresRand100K.RData")
save(tex2AcuteEdgeScores, file="pos2NaiveEdgeScoresRand100K.RData")
save(pos2NegEdgeScores, file="pos2NaiveEdgeScoresRand100K.RData")

allScores <- c(pos2NaiveEdgeScores[ ,"fracGd"],
			   tex2NaiveEdgeScores[ ,"fracGd"],
			   acute2NaiveEdgeScores[ ,"fracGd"],
			   tex2AcuteEdgeScores[ ,"fracGd"],
			   pos2NegEdgeScores[ ,"fracGd"])							# 554,100
			   
allScores <- sapply(allScores, as.numeric)
pdf("scoreDistributionHistRand100K.pdf")
H = hist(allScores, col=rgb(0,0.3,1,0.5), breaks=100,
		 border="white", xlab="edge score", main="") # breaks=20, 
dev.off()

# # # Scores > 0.94 and < 0.06 have <5% chance of being erroneous
# # # sum(H$counts[95:100])/sum(H$counts) # 0.04320339
# # # sum(H$counts[94:100])/sum(H$counts) # 0.04720628

# # # sum(H$counts[1:6])/sum(H$counts)	# 0.04821693

