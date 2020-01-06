setwd("~/TCE/CyJS")

cytoEdges <- get(load("cancerEdgeAnnot_14March2017.RData"))
cytoEdges <- cytoEdges[ ,1:2]

simNodes <- c("AKT","JUN","FOS", "BATF","BCL6","PRDM1","FOXO1","IFNg","IL12","IL2","IL21","IRF4","mTOR","NFATC1","NFATC2","NFkB","NR4A1","PDCD1","STAT3","TCF1","TCR")
simNodes <- toupper(simNodes)

selNodes = NULL; j = 0
allRes = list()
for (i in 1:length(simNodes)) {
	i1 <- grep(simNodes[i], cytoEdges[ ,1])
	i2 <- grep(simNodes[i], cytoEdges[ ,2])
	
	if (length(i1)>0 | length(i2)>0) {		
		j <- j + 1
		iCytoS <- i1[!is.na(i1)]
		iCytoT <- i2[!is.na(i2)]
	
		selNodes <- c(selNodes, simNodes[i])
		# allRes <- rbind(allRes, cytoEdges[c(iCytoS, iCytoT), ])
		allRes[[j]] <- cytoEdges[c(iCytoS, iCytoT), ]
	}
}

names(allRes) <- selNodes













