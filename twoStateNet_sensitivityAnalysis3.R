MPEC <- c("FOXO1","AKT","EZH2/PRC2","PIK3K","MTOR","NFATC2","NFkB","HIF1A","PGC1A","MAPK","RAS","JUN","FOS","EGR2/3","MYC","proliferation","glycolysis","MT biogenesis","NR4A1","TCF-1","BCL6","BACH2","ID3","E2A","CXCR5")
SLEC <- c("IL2-R","PPARa","NFATC1","IRF4","TBET","ZEB2","EOMES","GZMA/B","PRF1","IFNG","RUNX3","KLRG1","BLIMP-1","IL15-R","IL2-R","ID2","BATF","IL12-R","IL21-R")
other <- c("CD247","CD274","CD244","CD244/2B4","CD160","TIGIT","FAS/FASL","Tregs","TIM3","LAG3","BTLA","PD-1","CTLA-4","CTLA4","degranulation","CD3, TCR, CD8, CD28","4-1BB","CD137L")

setwd("~/TCE")
net <- read.delim("twoStateNetworkInteractionsForR.txt", 
				  sep="\t", as.is=TRUE, header=TRUE)

# Exclude feedback paths via receptors:
iBad <- grep("CD3|TCR|CD8|CD28", net$Target)
net <- net[-iBad, ]
iBad <- grep("CD3|TCR|CD8|CD28", net$Source)
net <- net[-iBad, ]

iBad <- which(net$Source == "NFATC2" & net$Target == "IL12-R")
net <- net[-iBad, ]

iBad <- which(net$Source == "ID3" & net$Target == "IL2-R")
net <- net[-iBad, ]

allS <- unique(net$Source)
allT <- unique(net$Target)
allNodes <- unique(c(allS, allT))

library(igraph)

G <- graph_from_data_frame(net)

# plot.igraph(tmpG, vertex.shape= "circle", edge.arrow.size=.3, 
	 # edge.color=ifelse(E(tmpG)$Interaction =="promotes", "forestgreen", "red"), 
	 # vertex.color="gold", vertex.size=5, vertex.frame.color="gray", 
	 # vertex.label.color="black", vertex.label.cex=0.5, 
	 # vertex.label.dist=0.5, edge.curved=0.5)
# dev.off()


allPaths = list(); iHits = 0
for (root in allS) {
	iHits <- iHits + 1
	
	paths <- get.all.shortest.paths(G, V(G)[root], allNodes)$res
	paths <- lapply(paths, function(v) V(G)[v]$name)
	allDirs = NULL
	pathDirs <-lapply(paths, function(x) {
				   for (i in 1:(length(x) - 1)) {
						if (length(x) > 1) {
							iS <- grep(x[i], net$Source)
							iT <- grep(x[i + 1], net$Target)
							edgeDir <- net[intersect(iS, iT), "Interaction"]
							names(edgeDir) <- x[i+1]
							allDirs <- c(allDirs, edgeDir)
						}
				   
				   }
				   return(allDirs)
			   })
	pathDirs <- pathDirs[!(unlist(lapply(pathDirs, is.null)))]
	names(pathDirs) <- rep(root, length(pathDirs)) 
	allPaths <- c(allPaths, pathDirs)
}

targets <- lapply(allPaths, function(x) names(x)[max(which(!is.na(names(x))))])
signT <- lapply(allPaths, function(x) {
				ifelse(any(x == "inhibits"), 
					ifelse(table(x)['inhibits'] %% 2 != 0, "-", "#"), 
					"+")
				})
names(targets) <- names(allPaths)
names(signT) <- names(allPaths)

signsTbl <- t(mapply(c, targets, signT))
colnames(signsTbl) <- c("target", "sign")


target.state <- rep("-", nrow(signsTbl))
iM <- which(!is.na(match(signsTbl[, 'target'], MPEC)))
target.state[iM] <- "MPEC"
iS <- which(!is.na(match(signsTbl[, 'target'], SLEC)))
target.state[iS] <- "SLEC"

signsTbl <- cbind(signsTbl, targetState=target.state)
signsTbl <- signsTbl[grepl("MPEC|SLEC", signsTbl[ ,"targetState"]), ]

write.csv(signsTbl, file="signsTbl.CSV")

tuples <- paste0(rownames(signsTbl), ":", signsTbl[ ,1])
i1 <- which(duplicated(tuples, fromLast=FALSE))
i2 <- which(duplicated(tuples, fromLast=TRUE))
iBad <- c(i1,i2)
signsTbl <- signsTbl[-iBad, ]

signsList <- apply(signsTbl, 1, function(x) paste(x[3], x[2], sep=""))
signsList <- split(signsList, names(signsList))
signsList <- lapply(signsList, table)

proSLEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("MPEC\\-|SLEC\\+", names(x))])))
proMPEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("MPEC\\+|SLEC\\-", names(x))])))

pdf("targetSelectionfrom2StateNetBarPlot_geneEffects.pdf", width=6, height=6)
par(mar=c(10,5,2,2))
M <- max(c(proMPEC, proSLEC))
barplot(hitCounts[ ,'proMPEC'], names.arg=rownames(hitCounts), 
		col=rgb(0.8,0.1, 0,0.35), 
		border="white", ylim=c(0, M), 
		las=2, cex.names=0.5, cex.lab=0.65, cex.axis=0.65,
		ylab="Number of downstream genes contributing to effect")
barplot(hitCounts[ ,'proSLEC'], names.arg="", col=rgb(0,0.5,0.25,0.5), 
		border="white", add=TRUE, yaxt='n')
legend("topright", legend=c("mem/prolif", "eff/exh"), 
	   col=c(rgb(0.8,0.1, 0,0.35), rgb(0,0.5,0.25,0.5)), 
	   pch=15, bty='n', cex=0.7, pt.cex=1.5)
dev.off()

#### 9Aug2018: Switch to show effects of KO, rather than effects of gene:
proSLEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("MPEC\\+|SLEC\\-", names(x))])))
proMPEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("MPEC\\-|SLEC\\+", names(x))])))

ratios <- proMPEC/(proSLEC + 0.1)
hitCounts <- cbind(proMPEC, proSLEC, ratios)
hitCounts <- hitCounts[rev(order(ratios, proMPEC, -1*proSLEC)), ]

# write.csv(hitCounts, file="hitCounts_iGraphV.CSV")

pdf("targetSelectionfrom2StateNetBarPlot_KO.effects.pdf", width=6, height=6)
par(mar=c(10,5,2,2))
M <- max(c(proMPEC, proSLEC))
barplot(hitCounts[ ,'proMPEC'], names.arg=rownames(hitCounts), 
		col=rgb(0.8,0.1, 0,0.35), 
		border="white", ylim=c(0, M), 
		las=2, cex.names=0.5, cex.lab=0.65, cex.axis=0.65,
		ylab="Number of downstream genes contributing to effect")
barplot(hitCounts[ ,'proSLEC'], names.arg="", col=rgb(0,0.5,0.25,0.5), 
		border="white", add=TRUE, yaxt='n')
legend("topright", legend=c("mem/prolif up", "eff/exh up"), 
	   col=c(rgb(0.8,0.1, 0,0.35), rgb(0,0.5,0.25,0.5)), 
	   pch=15, bty='n', cex=0.7, pt.cex=1.5)
dev.off()

upMPEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("MPEC\\-", names(x))])))
dnSLEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("SLEC\\+", names(x))])))

hitCounts <- cbind(upMPEC, dnSLEC)
hitCounts <- hitCounts[rev(order(upMPEC)), ]

pdf("targetSelectionfrom2StateNetBarPlot_upMMPECdnSLEC.pdf",
	height=7, width=7)
par(mar=c(10,5,2,2), cex=1.1)
M <- max(c(upMPEC, dnSLEC))
barplot(hitCounts[ ,'upMPEC'], names.arg=rownames(hitCounts), 
		col=rgb(0.8,0.1, 0,0.35), 
		border="white", ylim=c(0, M), 
		las=2, cex.names=0.5, cex.lab=0.65, cex.axis=0.65,
		ylab="Number of downstream genes affected")
barplot(hitCounts[ ,'dnSLEC'], names.arg="", col=rgb(0,0.5,0.25,0.5), 
		border="white", add=TRUE, yaxt='n')
legend("topright", legend=c("mem/prolif activated", "eff/exh inhibited"), 
	   col=c(rgb(0.8,0.1, 0,0.35), rgb(0,0.5,0.25,0.5)), 
	   pch=15, bty='n', cex=0.7, pt.cex=1.5)
dev.off()


upSLEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("SLEC\\-", names(x))])))
dnMPEC <- unlist(lapply(signsList, 
				  function(x) sum(x[grep("MPEC\\+", names(x))])))

hitCounts <- cbind(dnMPEC, upSLEC)
hitCounts <- hitCounts[order(upSLEC), ]

pdf("targetSelectionfrom2StateNetBarPlot_upSLEC.dnMPEC.pdf",
	height=1.25*7, width=7)
par(mar=c(10,5,2,2), cex=1.25)
M <- max(c(dnMPEC, upSLEC))
barplot(hitCounts[ ,'dnMPEC'], names.arg=rownames(hitCounts), 
		col=rgb(0.8,0.1, 0,0.35), 
		border="white", ylim=c(0, M), 
		las=2, cex.names=0.5, cex.lab=0.65, cex.axis=0.65,
		ylab="Number of downstream genes affected")
barplot(hitCounts[ ,'upSLEC'], names.arg="", col=rgb(0,0.5,0.25,0.5), 
		border="white", add=TRUE, yaxt='n')
legend("topleft", legend=c("mem/prolif inhibited", "eff/exh activated"), 
	   col=c(rgb(0.8,0.1, 0,0.35), rgb(0,0.5,0.25,0.5)), 
	   pch=15, bty='n', cex=0.7, pt.cex=1.5)
dev.off()


#### Source of multiple drug-target interaction lists:
# http://www.guidetopharmacology.org/lists.jsp

#### For now, using the integrated data from "dgiDB.org":
library(readxl)
setwd("~/TCE/drug_gene_interactions_guidetopharmacology.orglists.jsp")
drug2Gene <- read_excel("dgiDB.xlsx")
drug2Gene <- drug2Gene[drug2Gene$gene_name != "-" & 
					   !is.na(drug2Gene$gene_name), ]
drug2Gene <- drug2Gene[drug2Gene$interaction_types != "-" & 
					   !is.na(drug2Gene$interaction_types), ]			# 22755 x 8

setwd("~/TCE")

grep("EGR2|EGR3", drug2Gene$gene_name) 									# 0
# drug2Gene[grep("RAS", drug2Gene$gene_name), ]
# 1 HRAS      3265      GuideToPharmacologyI… inhibitor        LONAFARNIB
# 2 KRAS      3845      MyCancerGenomeClinic… inhibitor        REGORAFENIB
# 3 KRAS      3845      GuideToPharmacologyI… inhibitor        LONAFARNIB
# 4 KRAS      3845      MyCancerGenomeClinic… inhibitor        SELUMETINIB
# 5 KRAS      3845      MyCancerGenomeClinic… inhibitor        SIMVASTIN
# 6 KRAS      3845      CancerCommons         inhibitor        REOLYSIN
# 7 NRAS      4893      GuideToPharmacologyI… inhibitor        LONAFARNIB
# 8 KRAS      3845      TALC                  vaccine          RAS PEPTIDE CANCE…
drug2Gene[grep("EZH2", drug2Gene$gene_name), ]
# 1 EZH2      2146      GuideToPharmacologyI… inhibitor        EI1
# 2 EZH2      2146      GuideToPharmacologyI… inhibitor        UNC1999
# 3 EZH2      2146      GuideToPharmacologyI… inhibitor        JQEZ5
# 4 EZH2      2146      GuideToPharmacologyI… inhibitor        EPZ005687
# 5 EZH2      2146      GuideToPharmacologyI… inhibitor        GSK126
# 6 EZH2      2146      GuideToPharmacologyI… inhibitor        CPI-1205
# 7 EZH2      2146      GuideToPharmacologyI… inhibitor        TAZEMETOSTAT
# 8 EZH2      2146      GuideToPharmacologyI… inhibitor        EPZ011989
# 9 EZH2      2146      GuideToPharmacologyI… inhibitor        GSK343


#------------------------------------------------------------------------------#
################################ Cytoscape network #############################
#------------------------------------------------------------------------------#\

# Get Schietinger's data:
setwd("~/TCE/GSE89307_Scheitinger2017")

geneExp <- get(load("Scheitinger2017_exp.RData"))
geneExp <- geneExp[ which(!is.na(rownames(geneExp))), ]

M <- rowMeans(geneExp)													# 1570 all zero genes
geneExp <- geneExp[which(M != 0), ]										# 18584 x 38

grp1 = grep("^L5|^L7", colnames(geneExp))
grp2 = grep("^L14|^L28", colnames(geneExp))
types <- c(rep(0, length(grp1)), rep(1, length(grp2)))

library(limma)
library(edgeR)

dge <- DGEList(counts=geneExp[ ,c(grp1, grp2)])
dge <- calcNormFactors(dge, method="TMM")

limmaD <- model.matrix(~ types)
rownames(limmaD) <- colnames(dge)

v <- voom(counts=dge, design=limmaD, plot=TRUE)
fit <- lmFit(v, design=limmaD)
fit <- eBayes(fit)

DE <- topTable(fit, coef=ncol(limmaD), number=1E4, 
			   sort.by="logFC", p.value=1E-1)
# logFC > 0 means higher in grp2 (D14, 28)
earlyGenes <- rownames(DE)[DE$logFC < 0]								# 3426
lateGenes <- rownames(DE)[DE$logFC > 0]									# 3707


MPEC <- c("FOXO1","AKT1","AKT2","AKT3","EZH2","DNMT3A","DNMT3B","DNMT3L","KDM6B","PIK3CA","PIK3CB","PIK3CG","MTOR","NFATC2","NFKB1","NFKB2","REL","RELB","RELA","HIF1A","PGC1A","MAPK1","MAPK3","MAPK8","NRAS","KRAS","HRAS","JUN","JUNB","JUND","FOS","EGR1","EGR2","EGR3","MYC","PROLIFERATION","GLYCOLYSIS","MT BIOGENESIS","NR4A1","TCF7","BCL6","BACH2","ID3","E2A","CXCR5")
SLEC <- c("IL2","IL2RA","IL2RB","PPARA","NFATC1","IRF4","TBX21","ZEB2","EOMES","GZMA","GZMB","PRF1","IFNG","IFNGR1","IFNGR2","RUNX3","KLRG1","PRDM1","IL15RA","ID2","BATF","IL12RB1","IL12RB2","IL21R")

receptors <- c("TNFRSF1","TNFRSF2","TNFRSF9","CD247","CD3D","CD3E","CD3G","CD8A","CD28","TCR","IFNAR1","IFNAR2","IFNGR1","IFNGR2","CD244","CD160","TIGIT","FASL","HAVCR2","LAG3","BTLA","PDCD1","CTLA4","CD137","LCK","ZAP70")
ligands <- c("CD274","CD137L","FAS$","IL2$","IL12$","IL15$","IL21$","TNF$","ZN","^CA$") # "IFNG$",

setwd("~/TCE/CyJS/July2017")
net <- read.delim("cancerEdgeAnnot_17July2017.txt", 
				   sep="\t", header=TRUE, as.is=TRUE)

net <- net[ ,c("Source", "Target", "edgeType")]
colnames(net)[3] <- "Interaction"

# Cut off cell-level feedback paths
badQ <- paste(ligands, collapse="|")
iBad <- grep(badQ, net$Target)
net <- net[-iBad, ]

net[ , 1:2] <- apply(net[ ,1:2], 2, function(x) toupper(x))


#### Do 2 passes: (1) Early genes only, (2) late genes only ####

keepM <- unique(c(earlyGenes, MPEC))
keepS <- unique(c(lateGenes, SLEC))

# intersect(keepM, keepS)
 # [1] "IRF4"   "IL2RA"  "EOMES"  "IFNG"   "AKT3"   "DNMT3B" "PIK3CA" "PIK3CB"
 # [9] "PIK3CG" "MTOR"   "RELB"   "MAPK3"
# See Youngblood ~ DNMT3B

misclassified <- intersect(keepM, keepS)
keepM <- setdiff(keepM, misclassified)
keepS <- setdiff(keepS, misclassified)

setwd("~/TCE")
pdf("targetSelectionfromReduced505CytoscapeNet.pdf")
par(mar=c(10,5,2,2))

for (selKeep in list(keepM, keepS)) {
	iS <- match(net$Source, selKeep)
	iT <- match(net$Target, selKeep)
	iKeep <- unique(intersect(which(!is.na(iS)), which(!is.na(iT))))
	reducedNet <- net[iKeep, ]
	
	allS <- unique(reducedNet$Source)
	allT <- unique(reducedNet$Target)
	allNodes <- unique(c(allS, allT))

	library(igraph)
	G <- graph_from_data_frame(reducedNet)

	allPaths = list(); iHits = 0
	for (root in allS) {
		iHits <- iHits + 1
		
		paths <- get.all.shortest.paths(G, V(G)[root], allNodes)$res
		paths <- lapply(paths, function(v) V(G)[v]$name)
		
		pathLengths <- lapply(paths, length)
		paths <- paths[pathLengths > 1]
		
		allDirs = NULL
		pathDirs <-lapply(paths, function(x) {
						if (length(x) > 1) {
							for (i in 1:(length(x) - 1)) {
								iS <- grep(x[i], reducedNet$Source)
								iT <- grep(x[i + 1], reducedNet$Target)
								edgeDir <- reducedNet[intersect(iS, iT), "Interaction"]
								names(edgeDir) <- x[i+1]
								allDirs <- c(allDirs, edgeDir)
							}
					   }
					   return(allDirs)
				   })
		pathDirs <- pathDirs[!(unlist(lapply(pathDirs, is.null)))]
		names(pathDirs) <- rep(root, length(pathDirs)) 
		allPaths <- c(allPaths, pathDirs)
	}

	Q <- "ADHESION|PROLIFERATION|EXHAUSTION|INTERFERONS|GLYCOLYSIS|ZN|CA|BIOGENESIS|DEGRANULATION|TREGS"
	iBad <- grepl(Q, names(allPaths))
	allPaths <- allPaths[!iBad]

	nonReceptors <- which(is.na(match(names(allPaths), receptors)))
	allPaths <- allPaths[nonReceptors]

	# Should even-length repression chains be considered activation?
	targets <- lapply(allPaths, function(x) names(x)[max(which(!is.na(names(x))))])
	signT <- lapply(allPaths, function(x) {
					ifelse(any(x == "inhibits"), 
						ifelse(table(x)['inhibits'] %% 2 != 0, "-", "#"), 
						"+")
					})
	names(targets) <- names(allPaths)
	names(signT) <- names(allPaths)

	signsTbl <- t(mapply(c, targets, signT))
	colnames(signsTbl) <- c("target", "sign")


	target.state <- rep("-", nrow(signsTbl))
	iM <- which(!is.na(match(signsTbl[, 'target'], MPEC)))
	target.state[iM] <- "MPEC"
	iS <- which(!is.na(match(signsTbl[, 'target'], SLEC)))
	target.state[iS] <- "SLEC"

	signsTbl <- cbind(signsTbl, targetState=target.state)
	signsTbl <- signsTbl[grepl("MPEC|SLEC", signsTbl[ ,"targetState"]), ]

	tuples <- paste0(rownames(signsTbl), ":", signsTbl[ ,1])
	i1 <- which(duplicated(tuples, fromLast=FALSE))
	i2 <- which(duplicated(tuples, fromLast=TRUE))
	iBad <- c(i1,i2)
	signsTbl <- signsTbl[-iBad, ]

	signsList <- apply(signsTbl, 1, function(x) paste(x[3], x[2], sep=""))
	signsList <- split(signsList, names(signsList))
	signsList <- lapply(signsList, table)


	#### Effect of KO:
	proSLEC <- unlist(lapply(signsList, 
					  function(x) sum(x[grep("MPEC\\+|SLEC\\-", names(x))])))
	proMPEC <- unlist(lapply(signsList, 
					  function(x) sum(x[grep("MPEC\\-|SLEC\\+", names(x))])))
	ratios <- proMPEC/(proSLEC + 0.1)
	hitCounts <- cbind(proMPEC, proSLEC, ratios)
	hitCounts <- hitCounts[rev(order(ratios, proMPEC, -1*proSLEC)), ]

	M <- max(c(proMPEC, proSLEC))
	barplot(hitCounts[ ,'proMPEC'], names.arg=rownames(hitCounts), 
			col=rgb(0.8,0.1, 0,0.35), 
			border="white", ylim=c(0, M), 
			las=2, cex.names=0.5, cex.lab=0.65, cex.axis=0.65,
			ylab="Number of downstream genes affected")
	barplot(hitCounts[ ,'proSLEC'], names.arg="", col=rgb(0,0.5,0.25,0.5), 
			border="white", add=TRUE, yaxt='n')

	print("####")
}
dev.off()

plot(1:10, type='n')
legend("top", legend=c("mem/prolif up", "eff/exh up"), 
	   col=c(rgb(0.8,0.1, 0,0.35), rgb(0,0.5,0.25,0.5)), 
	   pch=15, bty='n', cex=0.7, pt.cex=1.5)
dev.off()


	# setwd("~/TCE")
	# save(allPaths, file="allPaths505NodeCytoscapeNet.RData")

	# setwd("~/TCE")
	# allPaths <- get(load("allPaths505NodeCytoscapeNet.RData"))

	# write.csv(signsTbl, file="signsTbl505NodeCytoscapeNet.CSV")
	
	# write.csv(hitCounts, file="hitCounts505NodeCytoscapeNet.CSV")

# save(hitCounts, file="reduced505net_MPEC_KO.pdf")


	# proSLEC <- unlist(lapply(signsList, 
					  # function(x) sum(x[grep("MPEC\\-|SLEC\\+", names(x))])))
	# proMPEC <- unlist(lapply(signsList, 
					  # function(x) sum(x[grep("MPEC\\+|SLEC\\-", names(x))])))
	# ratios <- proMPEC/(proSLEC + 0.1)
	# hitCounts <- cbind(proMPEC, proSLEC, ratios)
	# hitCounts <- hitCounts[order(ratios, proMPEC, -1*proSLEC), ]

	# M <- max(c(proMPEC, proSLEC))
	# barplot(hitCounts[ ,'proMPEC'], names.arg=rownames(hitCounts), 
			# col=rgb(0.8,0.1, 0,0.35), 
			# border="white", ylim=c(0, M), 
			# las=2, cex.names=0.35, cex.lab=0.65, cex.axis=0.65,
			# ylab="Number of downstream genes affected")
	# barplot(hitCounts[ ,'proSLEC'], names.arg="", col=rgb(0,0.5,0.25,0.5), 
			# border="white", add=TRUE, yaxt='n')
	# legend("topright", legend=c("pro mem/prolif", "pro eff/exh"), 
		   # col=c(rgb(0.8,0.1, 0,0.35), rgb(0,0.5,0.25,0.5)), 
		   # pch=15, bty='n', cex=0.7, pt.cex=1.5)








    Source Target Interaction
402  PRDM1   BCL6    inhibits
403  PRDM1  CXCR5    inhibits
404  PRDM1    ID3    inhibits
406  PRDM1    MYC    inhibits
407  PRDM1 NFATC1    inhibits
408  PRDM1  PDCD1    inhibits
409  PRDM1   TCF7    inhibits



