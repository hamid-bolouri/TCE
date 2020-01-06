# setwd("/Volumes/hbolouri/TCE/CyJS/July2017")
setwd("~/TCE/CyJS/July2017")

nodesTbl <- read.delim("cytoscapeEdgeAnnot_11July2017.txt", sep="\t",
						header=TRUE, quote="", as.is=TRUE)

nodesList <- unique(c(nodesTbl[ ,1], nodesTbl[ ,2]))

nodesTbl <- data.frame(Source=nodesList)


# setwd("/Volumes/hbolouri/TCE/_Data_cancer")
setwd("~/TCE/_Data_cancer")

aPDL1.v.TIL <- get(load("aPDL1.v.TIL.GSE93007.RData"))
PDL1plusPDL1minus <- get(load("PDL1plus.v.PDL1minus.GSE39205.RData"))
Tim3PlusTim3Minus <- get(load("Tim3Plus2Tim3Minus.GSE84072.RData"))
TILN2CMV <- get(load("TILN2CMV.GSE24536.RData"))
TILN2naive <- get(load("TILN2naive.GSE24536.RData"))
TILN2EBV <- get(load("TILN2EBV.GSE24536.RData"))
TILN2tumorPBMC <- get(load("TILN2tumorPBMC.GSE24536.RData"))

T.D21.2N.89307 <- get(load("T.D21.2N.89307.RData"))
T.D14.2N.89307 <- get(load("T.D14.2N.89307.RData"))
T.D7.2N.89307 <- get(load("T.D7.2N.89307.RData"))
T.D5.2N.89307 <- get(load("T.D5.2N.89307.RData"))
acuteD7.2N.89307 <- get(load("acuteD7.2N.89307.RData"))
acuteD5.2N.89307 <- get(load("acuteD5.2N.89307.RData"))

# setwd("/Volumes/hbolouri/TCE/CyJS/July2017")
setwd("~/TCE/CyJS/July2017")

indx <- match(nodesTbl$Source, rownames(aPDL1.v.TIL))
avExp.aPDL1.v.TIL <- apply(aPDL1.v.TIL[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.aPDL1.v.TIL <- 0
A <- log2(avExp.aPDL1.v.TIL)
nodesTbl$avExp.aPDL1.v.TIL[which(!is.na(indx))] <- 
	2*A/(max(A) + abs(min(A)))

indx <- match(nodesTbl$Source, rownames(PDL1plusPDL1minus))
avExp.PDL1plusPDL1minus <- apply(PDL1plusPDL1minus[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.PDL1plusPDL1minus <- 0
A <- log2(avExp.PDL1plusPDL1minus)
nodesTbl$avExp.PDL1plusPDL1minus[which(!is.na(indx))] <- 
	2*A/(max(A) + abs(min(A)))


indx <- match(nodesTbl$Source, rownames(Tim3PlusTim3Minus))
avExp.Tim3PlusTim3Minus <- apply(Tim3PlusTim3Minus[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.Tim3PlusTim3Minus <- 0
A <- log2(avExp.Tim3PlusTim3Minus)
nodesTbl$avExp.Tim3PlusTim3Minus[which(!is.na(indx))] <- 
	2*A/(max(A) + abs(min(A)))


indx <- match(nodesTbl$Source, rownames(TILN2CMV))
avExp.TILN2CMV <- apply(TILN2CMV[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.TILN2CMV <- 0
A <- log2(avExp.TILN2CMV)
nodesTbl$avExp.TILN2CMV[which(!is.na(indx))] <- 
	2*A/(max(A) + abs(min(A)))


indx <- match(nodesTbl$Source, rownames(TILN2naive))
avExp.TILN2naive <- apply(TILN2naive[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.TILN2naive <- 0
A <- log2(avExp.TILN2naive)
nodesTbl$avExp.TILN2naive[which(!is.na(indx))] <- 
	2*A/(max(A) + abs(min(A)))


indx <- match(nodesTbl$Source, rownames(TILN2EBV))
avExp.TILN2EBV <- apply(TILN2EBV[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.TILN2EBV <- 0
A <- log2(avExp.TILN2EBV)
nodesTbl$avExp.TILN2EBV[which(!is.na(indx))] <- 
	2*A/(max(A) + abs(min(A)))


indx <- match(nodesTbl$Source, rownames(TILN2tumorPBMC))
avExp.TILN2tumorPBMC <- apply(TILN2tumorPBMC[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.TILN2tumorPBMC <- 0
A <- log2(avExp.TILN2tumorPBMC)
nodesTbl$avExp.TILN2tumorPBMC[which(!is.na(indx))] <- 
	2*A/(max(A) + abs(min(A)))

#### Sheitinger:

indx <- match(nodesTbl$Source, rownames(acuteD5.2N.89307))
avExp.acuteD5.2N.89307 <- apply(acuteD5.2N.89307[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.acuteD5.2N.89307 <- 0
A <- log2(avExp.acuteD5.2N.89307)
nodesTbl$avExp.acuteD5.2N.89307[which(!is.na(indx))] <- 
	2*A/(max(A) - min(A))

indx <- match(nodesTbl$Source, rownames(acuteD7.2N.89307))
avExp.acuteD7.2N.89307 <- apply(acuteD7.2N.89307[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.acuteD7.2N.89307 <- 0
A <- log2(avExp.acuteD7.2N.89307)
nodesTbl$avExp.acuteD7.2N.89307[which(!is.na(indx))] <- 
	2*A/(max(A) - min(A))

indx <- match(nodesTbl$Source, rownames(T.D5.2N.89307))
avExp.T.D5.2N.89307 <- apply(T.D5.2N.89307[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.T.D5.2N.89307 <- 0
A <- log2(avExp.T.D5.2N.89307)
nodesTbl$avExp.T.D5.2N.89307[which(!is.na(indx))] <- 
	2*A/(max(A) - min(A))

indx <- match(nodesTbl$Source, rownames(T.D7.2N.89307))
avExp.T.D7.2N.89307 <- apply(T.D7.2N.89307[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.T.D7.2N.89307 <- 0
A <- log2(avExp.T.D7.2N.89307)
nodesTbl$avExp.T.D7.2N.89307[which(!is.na(indx))] <- 
	2*A/(max(A) - min(A))

indx <- match(nodesTbl$Source, rownames(T.D14.2N.89307))
avExp.T.D14.2N.89307 <- apply(T.D14.2N.89307[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.T.D14.2N.89307 <- 0
A <- log2(avExp.T.D14.2N.89307)
nodesTbl$avExp.T.D14.2N.89307[which(!is.na(indx))] <- 
	2*A/(max(A) - min(A))

indx <- match(nodesTbl$Source, rownames(T.D21.2N.89307))
avExp.T.D21.2N.89307 <- apply(T.D21.2N.89307[indx[!is.na(indx)], ], 1, median)
nodesTbl$avExp.T.D21.2N.89307 <- 0
A <- log2(avExp.T.D21.2N.89307)
nodesTbl$avExp.T.D21.2N.89307[which(!is.na(indx))] <- 
	2*A/(max(A) - min(A))


# setwd("/Volumes/hbolouri/TCE/CyJS/July2017")
setwd("~/TCE/CyJS/July2017")

save(nodesTbl, file="cancerNodeAnnot_17July2017.RData")
write.table(nodesTbl, file="cancerNodeAnnot_17July2017.txt", sep="\t",
			col.names=TRUE, row.names=FALSE, quote=FALSE)











