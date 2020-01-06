setwd("~/TCE/CyJS/July2017")

litAnnotTbl <- read.delim("cytoscapeEdgeAnnot_11July2017.txt", sep="\t",
						  header=TRUE, quote="", as.is=TRUE)

nodesFromEdges <- unique(c(litAnnotTbl[ ,1], litAnnotTbl[ ,2]))

nodesTbl <- data.frame(Symbol = unique(c(litAnnotTbl[ ,1], litAnnotTbl[ ,2])))

setwd("~/TCE")
pos2Naive <- get(load("CXCR5pos2naive.RData")) 
tex2Naive <- get(load("tex2Naive.RData"))
acute2Naive <- get(load("acute2Naive.RData"))
tex2Acute <- get(load("tex2Acute.RData"))
pos2neg <- get(load("CXCR5pos2neg.RData"))
neg2Naive <- get(load("CXCR5neg2naive.RData"))
setwd("~/TCE/CyJS/July2017")
#### NB T-bet (TBX21) is NOT in the above tables ####

indx <- match(nodesTbl$Symbol, rownames(pos2Naive))
avExpPos2Naive <- apply(pos2Naive[indx[!is.na(indx)], ], 1, mean)
nodesTbl$avExpPos2Naive <- 0
nodesTbl$avExpPos2Naive[which(!is.na(indx))] <- log2(avExpPos2Naive + 1E-4)


indx <- match(nodesTbl$Symbol, rownames(neg2Naive))
avExpNeg2Naive <- apply(neg2Naive[indx[!is.na(indx)], ], 1, mean)
nodesTbl$avExpNeg2Naive <- 0
nodesTbl$avExpNeg2Naive[which(!is.na(indx))] <- log2(avExpNeg2Naive + 1E-4)


indx <- match(nodesTbl$Symbol, rownames(tex2Naive))
avExpTex2Naive <- apply(tex2Naive[indx[!is.na(indx)], ], 1, mean)
nodesTbl$avExpTex2Naive <- 0
nodesTbl$avExpTex2Naive[which(!is.na(indx))] <- log2(avExpTex2Naive + 1E-4)


indx <- match(nodesTbl$Symbol, rownames(acute2Naive))
avExpAcute2Naive <- apply(acute2Naive[indx[!is.na(indx)], ], 1, mean)
nodesTbl$avExpAcute2Naive <- 0
nodesTbl$avExpAcute2Naive[which(!is.na(indx))] <- log2(avExpAcute2Naive + 1E-4)


indx <- match(nodesTbl$Symbol, rownames(tex2Acute))
avExpTex2Acute <- apply(tex2Acute[indx[!is.na(indx)], ], 1, mean)
nodesTbl$avExpTex2Acute <- 0
nodesTbl$avExpTex2Acute[which(!is.na(indx))] <- log2(avExpTex2Acute + 1E-4)


indx <- match(nodesTbl$Symbol, rownames(pos2neg))
avExpPos2Neg <- apply(pos2neg[indx[!is.na(indx)], ], 1, mean)
nodesTbl$avExpPos2Neg <- 0
nodesTbl$avExpPos2Neg[which(!is.na(indx))] <- log2(avExpPos2Neg + 1E-4)

write.table(nodesTbl, file="CytoscapeNodeAnnot_11July2017.txt", sep="\t",
			col.names=TRUE, row.names=FALSE, quote=FALSE)











