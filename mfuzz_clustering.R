setwd("~/TCE/CyJS/July2017")
# setwd("/Volumes/hbolouri/TCE/CyJS/July2017")
edges <- read.delim("cancerEdgeAnnot_17July2017.txt", 
					sep="\t", header=TRUE, as.is=TRUE)
nodes <- sort(unique(c(edges[ ,1], edges[ ,2])))



setwd("~/TCE/GSE89307_Scheitinger2017")
TBL <- get(load("Scheitinger2017_exp.RData"))

SDs <- apply(TBL, 1, sd)
iBad <- which(SDs == 0)

TBL <- TBL[-iBad, ]
TBL <- TBL[!is.na(rownames(TBL)), ]

TBL <- log10(TBL + 1)

colnames(TBL)[grep("^N", colnames(TBL))] <- "N"
colnames(TBL)[grep("^M1|^M2|^M3", colnames(TBL))] <- "M"
colnames(TBL)[grep("^ML7", colnames(TBL))] <- "ML7"
colnames(TBL)[grep("^ML35", colnames(TBL))] <- "ML35"

library(limma)
TBL <- t(avereps(t(TBL), colnames(TBL)))

# library(caret)
# nearZero <- apply(TBL, 1, nearZeroVar)							# nothing
# zeroI <- which(nearZero == 1)
# TBL <- TBL[-zeroI, ]

zeroCount <- apply(TBL, 1, function(x) length(which(x < 0.5)))
nearZero <- which(zeroCount > ncol(TBL)*0.75)

TBL <- TBL[-nearZero, ]												# 13740 x 13

corTbl <- cor(t(TBL))

triTbl <- corTbl
triTbl[lower.tri(triTbl, diag=TRUE)] <- NA

indx <- which(abs(triTbl) > 0.9, arr.ind=TRUE)						# 614,177 x 2

tmp <- rev(sort(table(c(rownames(corTbl)[indx[ ,1]], colnames(corTbl)[indx[ ,2]]))))

selCors <- corTbl
selCors[selCors < 0.9] <- NA

countsNAs <- apply(selCors, 1, function(x) length(is.na)

library(pheatmap)
pheatmap(selCors)
dev.off()


#### Fuzzy time course clustering:
library(Mfuzz)

tTbl <- TBL[ ,grep("^L", colnames(TBL))]

eset <- ExpressionSet(tTbl)
# M <- mestimate(eset)												# 1.053

eset <- standardise(eset)

cl <- mfuzz(eset, c=49, m=1.5)

pdf('mfuzzClusters49M1.5.Z.pdf', height=21, width=21)
mfuzz.plot(eset, cl=cl, mfrow=c(7, 7), new.window=FALSE,
		   time.labels=c("5", "7", "14", "21", "28", "35", "60"),
		   min.mem=0.5) 
# mfuzzColorBar(horizontal=TRUE)
dev.off()

acore.list <- acore(eset, cl=cl, min.acore=0.5)

sumHits = 0
for (i in 1:length(acore.list)) {
	clustI <- rownames(acore.list[[i]])
	genes <- intersect(clustI, nodes)
	X <- length(genes)
	print(paste0(i, "    ", X, "   ", paste(genes, collapse=", ")))
	sumHits <- sumHits + X
}

5     3   CD28, HIF1A, TCF7
9     1   MTOR
13    1   RELB
14    3   ID2, IL2RB, PDCD1
19    6   TIGIT, PRDM1, GZMA, IL12RB2, CD244, TNFRSF1A
24    4   ID3, LEF1, BTLA, ICOS
25    1   IFNGR1
27    1   HAVCR2
35    1   PRKCQ
37    1   RICTOR
38    1   NFKB1
39    3   EGR3, IRF4, MAPK8
43    1   REL
46    1   SLC39A14
49    1   FOXO1

save(acore.list, file="mfuzzClusters.RData")

















