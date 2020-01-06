setwd("~/TCE/CyJS/July2017")

# edges <- read.delim("cytoscapeEdgeAnnot_11July2017.txt", 
					# sep="\t", header=TRUE, as.is=TRUE)
# allGenes <- unique(c(edges[ ,1], edges[ ,2]))

setwd("~/TCE/CyJS/July2017")
edges <- get(load("cancerEdgeAnnot_17July2017.RData"))					# 505
i1 <- !grepl("proliferation|glycolysis|ZN|CA", edges[ , 1])
i2 <- !grepl("Adhesion|proliferation|glycolysis|ZN|CA", edges[ , 2])
edges <- edges[(i1 & i2), ]												# 478

nodes <- unique(c(edges[ , 1], edges[ , 2]))
allGenes <- intersect(nodes, rownames(acuteD5.2N.89307))				# 108


setwd("~/TCE/GSE89307_Scheitinger2017")
TBL <- get(load("Scheitinger2017_exp.RData"))

iA <- grep("^E5", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iA) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("acuteD5.2N.89307.", 1:dim(allFC)[2])
acuteD5.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(acuteD5.2N.89307, file="acuteD5.2N.89307.RData")					# 9 comparisons


iA <- grep("^E7", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iA) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("acuteD7.2N.89307.", 1:dim(allFC)[2])
acuteD7.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(acuteD7.2N.89307, file="acuteD7.2N.89307.RData")					# 9 comparisons


#### Liver cancer:

iT <- grep("^L5", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iT) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("T.D5.2N.89307.", 1:dim(allFC)[2])
T.D5.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(T.D5.2N.89307, file="T.D5.2N.89307.RData")							# 9 comparisons


iT <- grep("^L7", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iT) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("T.D7.2N.89307.", 1:dim(allFC)[2])
T.D7.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(T.D7.2N.89307, file="T.D7.2N.89307.RData")							# 9 comparisons


iT <- grep("^L14", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iT) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("T.D14.2N.89307.", 1:dim(allFC)[2])
T.D14.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(T.D14.2N.89307, file="T.D14.2N.89307.RData")						# 9 comparisons


iT <- grep("^L21", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iT) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("T.D21.2N.89307.", 1:dim(allFC)[2])
T.D21.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(T.D21.2N.89307, file="T.D21.2N.89307.RData")						# 9 comparisons


iT <- grep("^L28", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iT) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("T.D28.2N.89307.", 1:dim(allFC)[2])
T.D28.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(T.D28.2N.89307, file="T.D28.2N.89307.RData")						# 9 comparisons


iT <- grep("^L35", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iT) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("T.D35.2N.89307.", 1:dim(allFC)[2])
T.D35.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(T.D35.2N.89307, file="T.D35.2N.89307.RData")						# 9 comparisons


iT <- grep("^L60", colnames(TBL))
iN <- grep("^N", colnames(TBL))

allFC = NULL
for (T in iT) {
	for (C in iN) {
		FC <- (TBL[ , T] + 0.1) / (TBL[ , C]	+ 0.1)
		allFC <- cbind(allFC, FC)
	}
}
colnames(allFC) <- paste0("T.D60.2N.89307.", 1:dim(allFC)[2])
T.D60.2N.89307 <- allFC

setwd("~/TCE/_Data_cancer")
save(T.D60.2N.89307, file="T.D60.2N.89307.RData")						# 9 comparisons






































