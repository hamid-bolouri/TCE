# From 'mkPosNegExpTbl.R'
setwd("/home/hbolouri/TCE/_Data_infection")

iqrGSE76279 <- log2(get(load("iqrGSE76279.RData")) + 1)
iqrGSE74148 <- log2(get(load("iqrGSE74148.RData")) + 1)
iqrGSE84105 <- log2(get(load("iqrGSE84105.RData")) + 1)


#### For iqrGSE76279:
colnames(iqrGSE76279) <- c(rep("CXCR5+", 3), rep("CXCR5-", 3), rep("Naive",2))

library(limma)
avGSE76279 <- t(avereps(t(iqrGSE76279), colnames(iqrGSE76279)))
avGSE76279 <- avGSE76279[ ,c(3,1,2)]

#### For iqrGSE74148:
avGSE74148 <- iqrGSE74148
colnames(avGSE74148) <- c("Naive", "Effector", "CXCR5+", "CXCR5-")

#### For iqrGSE84105:
tmp <-sapply(colnames(iqrGSE84105), function(x) strsplit(x, "_"))
colnames(iqrGSE84105) <- unlist(lapply(tmp, function(x) x[1]))

library(limma)
avGSE84105 <- t(avereps(t(iqrGSE84105), colnames(iqrGSE84105)))


#------------------------------------------------------------------------------#
#					Immunological Genome Project MOUSE data 				   #	
#------------------------------------------------------------------------------#

setwd("~/TCE")
immProj <- read.delim("ImmunologicalGenomeProject_CD8ExpTimeCourse_2013.txt",
					  sep="\t", header=TRUE, as.is=TRUE)
iExp <- immProj[ , 3:ncol(immProj)]
library(limma)
iExp <- avereps(iExp, immProj[ , 2])
iExp <- iExp[ , grep("Sp.OT1", colnames(iExp))]							# 50 samples
rownames(iExp) <- toupper(rownames(iExp))
iExp <- log2(iExp + 1)

library(limma)
colnames(iExp) <- gsub("^T.8|.1$|.2$|.3$|.4$|.5$", "", 
							colnames(iExp))
avImmProj <- t(avereps(t(iExp)))
iRm <- grep("Mem", colnames(avImmProj))
avImmProj <- avImmProj[ ,-iRm]


#------------------------------------------------------------------------------#
#						Shietinger raw per gene exp 						   #	
#------------------------------------------------------------------------------#

setwd("~/TCE/GSE89307_Scheitinger2017")

Scheitinger <- log2(get(load("Scheitinger2017_exp.RData")) +1)

M <- rowMeans(Scheitinger)												# 1570 all zero genes
Scheitinger <- Scheitinger[which(M != 0), ]								# 18584 x 38

colnames(Scheitinger) <- gsub("^L", "TD", colnames(Scheitinger))

colnames(Scheitinger) <- gsub("N1|N2|N3", "N", colnames(Scheitinger))
colnames(Scheitinger) <- gsub("M1|M2|M3", "M", colnames(Scheitinger))

library(limma)
avScheitinger <- t(avereps(t(Scheitinger), colnames(Scheitinger)))


#------------------------------------------------------------------------------#
#						Wherry raw per gene exp 							   #	
#------------------------------------------------------------------------------#

setwd("~/TCE/_Data_infection/GSE41867")
Wherry <- log2(get(load("perGeneExp.GSE41867.RData")) +1)
setwd("~/TCE")															# 36 samples

rownames(Wherry) <- toupper(rownames(Wherry))

outliers <- c("Naive_4","AcuteD6_5","AcuteD8_9","AcuteD15_14","AcuteD30_18",
			  "ChronicD6_21","ChronicD8_25","ChronicD15_29","ChronicD15_30",
			  "ChronicD30_35")
keep <- setdiff(colnames(Wherry), outliers)
Wherry <- Wherry[ , keep]												# 26 samples

tmp <- sapply(colnames(Wherry), function(x) strsplit(x, "_"))
colnames(Wherry) <- unlist(lapply(tmp, function(x) x[1]))

library(limma)
avWherry <- t(avereps(t(Wherry), colnames(Wherry)))


#------------------------------------------------------------------------------#
#								metabolic pathways							   #
#------------------------------------------------------------------------------#

setwd("/home/hbolouri/TCE/Becky_28July2017")
beckyCD8 <- get(load("expCD8_corrected955.RData"))

library(limma)
avBecky <- t(avereps(t(beckyCD8), colnames(beckyCD8)))

j <- c(1,6,2:3,5,4,7,17,18,8,grep(":gluL", colnames(avBecky)),
	   grep(":L:", colnames(avBecky)), 
	   grep(":1:", colnames(avBecky)), 
	   grep(":6:", colnames(avBecky)))
avBecky <- avBecky[ , j]

#------------------------------------------------------------------------------#
#								metabolic pathways							   #
#------------------------------------------------------------------------------#

hemeBioGenesis <- c("ALAD","ALAS1","ALAS2","CPOX",
					"FECH","HMBS","PPOX","UROD","UROS") # from wikipathways

lyposis <- C("ATGL")

FAO <- c("AMPK", "PPARG", "PPARGC1A", "CPT1A")

aaMetabolism <- c("SLC7A5", "SLC3A2", "SLC38A2", "SLC1A5", "SLC38A1")

glycolysis <- c("PFKP", "PGAM1", "HK1", "HK2", "PKM2", "PKM",
				"LDHA", "SLC2A1", "SLC16A1", "SLC16A3")

gluconeogenUP <- c("FOXO1", "G6PC", "G6PC3", "PCK1", "PCK2", 
				   "PPARGC1A", "CREB1")
gluconeogenDN <- "PPARG"

lipidBioGenUP <- c("SREBP1", "SREBF1", "LDLR", "HMGCR")
lipidBioGenDN <- "ABCA1"

oxphos <- c("NDUFAF2", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4",
			"NDUFA5", "NDUFA6", "NDUFA7", "NDUFA8", "NDUFA9", 
			"NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "COX6C")

TCA <- c("ACSS1", "CPT1A", "ACO1","ACO2","CS","DLAT","DLD","DLST",
		 "FH1","IDH2","IDH3A","IDH3B","IDH3G","MDH2","OGDH","OGDHL",
		 "PCK1","PCK2","PDHA1","PDHA2","PDHB","SDHA","SDHB","SDHC",
		 "SDHD","SUCLA2","SUCLG1","SUCLG2")

glutaminoLysis <- c("SLC1A5", "SLC7A5", "GLS", "GLUD1")


IFNG <- c("IFNG", "GZMB", "PRF1", "PDCD1")


allGenes <- c(glycolysis, gluconeogenUP, gluconeogenDN,
			  lipidBioGenUP, lipidBioGenDN, oxphos,
			  TCA, glutaminoLysis, IFNG)

allMetGenes <- c(glycolysis, gluconeogenUP, gluconeogenDN,
				 lipidBioGenUP, lipidBioGenDN, oxphos,
				 TCA, glutaminoLysis)

metPathways <- list(glycolysis, gluconeogenUP, gluconeogenDN,
					lipidBioGenUP, lipidBioGenDN, oxphos,
					TCA, glutaminoLysis)
names(metPathways) <- c("glycolysis","gluconeogenUP","gluconeogenDN",
						"lipidBioGenUP","lipidBioGenDN","oxphos",
						"TCA","glutaminoLysis")

#------------------------------------------------------------------------------#
#						calculate metabolic signatures						   #
#------------------------------------------------------------------------------#

#### NB gene sets have both mouse and human symbols !

calcSigs <- function(expTbl) {
	iGly <- match(glycolysis, rownames(expTbl))
	iGly <- iGly[!is.na(iGly)]

	glySIG <- apply(expTbl[iGly, ], 2, mean)
	glySIG <- glySIG - min(glySIG)
	glySIG <- glySIG / max(abs(glySIG))


	iGluUP <- match(gluconeogenUP, rownames(expTbl))
	iGluUP <- iGluUP[!is.na(iGluUP)]
	iGluDN <- match("PPARG", rownames(expTbl))

	gluSIG.1 <- apply(expTbl[iGluUP, ], 2, mean)
	# gluSIG.2 <- apply(expTbl[iGluDN, ], 2, mean)
	gluSIG.2 <- mean(unlist(expTbl[iGluDN, ]))
	glucSIG <- gluSIG.1 / gluSIG.2
	glucSIG <- glucSIG - min(glucSIG)
	glucSIG <- glucSIG / max(abs(glucSIG))


	iLipidUP <- match(lipidBioGenUP, rownames(expTbl))
	iLipidUP <- iLipidUP[!is.na(iLipidUP)]
	iLipidDN <- match(lipidBioGenDN, rownames(expTbl))

	lipSIG.1 <- apply(expTbl[iLipidUP, ], 2, mean)
	# lipSIG.2 <- apply(expTbl[iLipidDN, ], 2, mean)
	lipSIG.2 <- mean(unlist(expTbl[iLipidDN, ]))
	lipSIG <- lipSIG.1 / lipSIG.2
	lipSIG <- lipSIG - min(lipSIG)
	lipSIG <- lipSIG / max(abs(lipSIG))


	iOxphos <- match(oxphos, rownames(expTbl))
	iOxphos <- iOxphos[!is.na(iOxphos)]

	oxSIG <- apply(expTbl[iOxphos, ], 2, mean)
	oxSIG <- oxSIG - min(oxSIG)
	oxSIG <- oxSIG / max(abs(oxSIG))


	iTCA <- match(TCA, rownames(expTbl))
	iTCA <- iTCA[!is.na(iTCA)]

	tcaSIG <- apply(expTbl[iTCA, ], 2, mean)
	tcaSIG <- tcaSIG - min(tcaSIG)
	tcaSIG <- tcaSIG / max(abs(tcaSIG))


	iGlutamin <- match(glutaminoLysis, rownames(expTbl))
	iGlutamin <- iGlutamin[!is.na(iGlutamin)]

	glutSIG <- apply(expTbl[iGlutamin, ], 2, mean)
	glutSIG <- glutSIG - min(glutSIG)
	glutSIG <- glutSIG / max(abs(glutSIG))

	ifngSIG <- apply(expTbl[intersect(rownames(expTbl), IFNG), ], 2, mean)
	ifngSIG <- ifngSIG - min(ifngSIG)
	ifngSIG <- ifngSIG / max(abs(ifngSIG))

	sigMtrx <- rbind(glycolysis=glySIG, gluconeogenesis=glucSIG,
					 lipidSynthesis=lipSIG, oxphos=oxSIG,
					 TCA=tcaSIG, glutaminoLysis=glutSIG, IFNG=ifngSIG)

	return(sigMtrx)
}


calcRawSigs <- function(expTbl) {
	iGly <- match(glycolysis, rownames(expTbl))
	iGly <- iGly[!is.na(iGly)]

	glySIG <- apply(expTbl[iGly, ], 2, mean)


	iGluUP <- match(gluconeogenUP, rownames(expTbl))
	iGluUP <- iGluUP[!is.na(iGluUP)]
	iGluDN <- match("PPARG", rownames(expTbl))

	gluSIG.1 <- apply(expTbl[iGluUP, ], 2, mean)
	# gluSIG.2 <- apply(expTbl[iGluDN, ], 2, mean)
	gluSIG.2 <- mean(unlist(expTbl[iGluDN, ]))
	glucSIG <- gluSIG.1 # / gluSIG.2


	iLipidUP <- match(lipidBioGenUP, rownames(expTbl))
	iLipidUP <- iLipidUP[!is.na(iLipidUP)]
	iLipidDN <- match(lipidBioGenDN, rownames(expTbl))

	lipSIG.1 <- apply(expTbl[iLipidUP, ], 2, mean)
	# lipSIG.2 <- apply(expTbl[iLipidDN, ], 2, mean)
	lipSIG.2 <- mean(unlist(expTbl[iLipidDN, ]))
	lipSIG <- lipSIG.1 # / lipSIG.2


	iOxphos <- match(oxphos, rownames(expTbl))
	iOxphos <- iOxphos[!is.na(iOxphos)]

	oxSIG <- apply(expTbl[iOxphos, ], 2, mean)


	iTCA <- match(TCA, rownames(expTbl))
	iTCA <- iTCA[!is.na(iTCA)]

	tcaSIG <- apply(expTbl[iTCA, ], 2, mean)


	iGlutamin <- match(glutaminoLysis, rownames(expTbl))
	iGlutamin <- iGlutamin[!is.na(iGlutamin)]

	glutSIG <- apply(expTbl[iGlutamin, ], 2, mean)
	
	ifngSIG <- apply(expTbl[intersect(rownames(expTbl), IFNG), ], 2, mean)

	sigMtrx <- rbind(glycolysis=glySIG, gluconeogenesis=glucSIG,
					 lipidSynthesis=lipSIG, oxphos=oxSIG,
					 TCA=tcaSIG, glutaminoLysis=glutSIG, IFNG=ifngSIG)

	return(sigMtrx)
}


#------------------------------------------------------------------------------#
#						Construct SIGNATURE table and makePlot				   #
#------------------------------------------------------------------------------#

setwd("/home/hbolouri/TCE/")

#### Current data sets:
# Scheitinger Wherry iqrGSE76279 iqrGSE74148 iqrGSE84105

setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolicAvScore_heatmap_GSE84105.pdf", width=3.5, height=8)
pheatmap(calcSigs(avGSE84105), legend_breaks=seq(0, 1, by=0.2),
		 fontsize_col=10, fontsize_row=10, scale="none",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()


setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolicAvScore_heatmap_GSE74148.pdf", width=4, height=8)
pheatmap(calcSigs(avGSE74148), legend_breaks=seq(0, 1, by=0.2),
		 fontsize_col=10, fontsize_row=10, scale="none",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()


setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolicAvScore_heatmap_GSE76279.pdf", width=3.5, height=8)
pheatmap(calcSigs(avGSE76279), legend_breaks=seq(0, 1, by=0.2),
		 fontsize_col=10, fontsize_row=10, scale="none",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()


#### ImmuneGenomeProject
setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolicAvScore_heatmap_avImmProj.pdf", width=6, height=10)
colnames(avImmProj) <- gsub("Eff.Sp.OT1.", "", colnames(avImmProj))
pheatmap(calcSigs(avImmProj), legend_breaks=seq(0, 1, by=0.2),
		 fontsize_col=15, fontsize_row=15, scale="none",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off() #idth=4, height=8


#### For Wherry:
setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolicAvScore_heatmap_Wherry.pdf", width=15, height=10)
pheatmap(calcSigs(avWherry), 
		 fontsize_col=15, fontsize_row=15, scale="none",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()


setwd("/home/hbolouri/TCE/")
# mtrx <- calcRawSigs(avWherry)
# mtrx <- mtrx/mtrx[ ,1]
mtrx <- avWherry[intersect(allMetGenes, rownames(avWherry)), c(1:4, 6:8)]
PCA <- prcomp(t(mtrx), center=FALSE, scale=FALSE)
pdf("PCA_6MetabolicPathwayGenes_Wherry.pdf")
plot(PCA$x, pch=19, col=rgb(0,0.3,1,0.7))
text(PCA$x[ ,1], PCA$x[ ,2], rownames(PCA$x), pos=1, cex=0.7, co="blue")
dev.off()


#### For Scheitinger:
setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolicAvScore_heatmap_ScheitingerIFNG_v2.pdf", width=6, height=10) 
pheatmap(calcSigs(avScheitinger[ , 1:11]), 								
		 fontsize_col=15, fontsize_row=15, scale="none",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off() # was 'width=15' # added '[1:11]'

setwd("/home/hbolouri/TCE/")
# mtrx <- calcRawSigs(avScheitinger)
# mtrx <- mtrx/mtrx[ ,1]
mtrx <- avScheitinger[intersect(allMetGenes, rownames(avScheitinger)), 1:8]
PCA <- prcomp(t(mtrx), center=FALSE, scale=FALSE)
pdf("PCA_6MetabolicPathwayGenes_Shietinger.pdf")
plot(PCA$x, pch=19, col=rgb(0,0.3,1,0.7))
text(PCA$x[ ,1], PCA$x[ ,2], rownames(PCA$x), pos=1, cex=0.7, co="blue")
dev.off()


#### For Becky:
setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolicAvScore_heatmap_Becky.pdf", width=15, height=8)
pheatmap(calcSigs(avBecky), 
		 fontsize_col=15, fontsize_row=15, scale="none",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()


setwd("/home/hbolouri/TCE/")
# mtrx <- calcRawSigs(avBecky)
# mtrx <- mtrx/mtrx[ ,1]
mtrx <- avBecky[intersect(allGenes, rownames(avBecky)), ]
PCA <- prcomp(t(mtrx), center=FALSE, scale=FALSE)
pdf("PCA_6MetabolicPathwayGenes_Becky.pdf")
plot(PCA$x, pch=19, col=rgb(0,0.3,1,0.7), ylim=c(-12,5), xlim=c(-89, -74))
points(PCA$x[6, 1], PCA$x[6, 2], pch=19, col="red")
text(PCA$x[ ,1], PCA$x[ ,2], rownames(PCA$x), pos=1, cex=0.7, co=rgb(0,0.3,0.7))
text(PCA$x[6 ,1], PCA$x[6 ,2], rownames(PCA$x)[6], pos=1, cex=0.7, co="red")
dev.off()


#------------------------------------------------------------------------------#
#					metabolic profile of immune-active tumors				   #
#------------------------------------------------------------------------------#

setwd("~/TCGA.TARGET.GTEX")
bigTbl <- readRDS("expTbl.RDS", refhook = NULL)


setwd("/home/hbolouri/TCE/")
library(pheatmap)
# png("heatmap_metGenes_19KBigTbl_rowScaled.png", width=4800, height=1600)
tmpTbl <- bigTbl[intersect(allMetGenes, rownames(bigTbl)), ]
tmpTbl <- t(apply(tmpTbl, 1, scale))
tmpTbl[tmpTbl > 2] <- 2
tmpTbl[tmpTbl <  -2] <- (-2)

SDs <- apply(tmpTbl, 1, sd)
tmpTbl <- tmpTbl[SDs > 0, ]

colnames(tmpTbl) <- colnames(bigTbl)

annotCol = data.frame(sampleType=colnames(bigTbl), stringsAsFactors=FALSE)
rownames(annotCol) <- colnames(bigTbl)
annotCol$sampleType[grep("GTEX", annotCol$sampleType)] <- "GTEX"
annotCol$sampleType[grep("K.562", annotCol$sampleType)] <- "K562"
annotCol$sampleType[grep("TCGA", annotCol$sampleType)] <- "TCGA"
annotCol$sampleType[grep("TARGET", annotCol$sampleType)] <- "TARGET"

annotRow = data.frame(pathway=rownames(tmpTbl), stringsAsFactors=FALSE)
rownames(annotRow) <- rownames(tmpTbl)
for (i in 1:length(metPathways)) {
	indx <- match(metPathways[[i]], rownames(annotRow))
	indx <- indx[!is.na(indx)]
	annotRow$pathway[indx] <- 
			  names(metPathways)[i]
}

colors22 <- rev(c("#fffac8","#000000","#0082c8","#aa6e28","#ffd8b1","#46f0f0","#3cb44b","#808080","#e6beff","#d2f53c","#f032e6","#800000","#aaffc3","#000080","#808000","#f58231","#fabebe","#911eb4","#e6194b","#008080","#FFFFFF","#ffe119"))
colColors <- colors22[1:length(unique(annotCol$sampleType))]
names(colColors) <- unique(annotCol$sampleType)
rowColors <- colors22[1:length(unique(annotRow$pathway))]
names(rowColors) <- unique(annotRow$pathway)

myColors <- list(sampleType=colColors, 
				 pathway=rowColors)

pheatmap(tmpTbl, 
		 fontsize_col=0.1, fontsize_row=5, scale="none",
		 cluster_cols=TRUE, cluster_rows=FALSE,
		 annotation_col=annotCol, annotation_row=annotRow,
		 annotation_colors=rev(myColors),
		 clustering_method="ward.D2",
		 filename="heatmap_metGenes_19KBigTbl_rowScaled.png")
dev.off()


setwd("/home/hbolouri/TCE/")
mtrx <- bigTbl[intersect(allGenes, rownames(bigTbl)), ]
PCA <- prcomp(t(mtrx), center=FALSE, scale=FALSE)
png("PCA_6MetabolicPathwayGenes_19KBigTbl.png", width=1000, height=1000)
plot(PCA$x, pch=19, col=rgb(0,0.3,1,0.5))
# i <- 
# points(PCA$x[6, 1], PCA$x[6, 2], pch=19, col="red")
dev.off()










#------------------------------------------------------------------------------#
#					Construct EXPRESSION table and makePlot					   #
#------------------------------------------------------------------------------#

tmpMtrx <- expTbl[c(iGly, iGluUP, iGluDN, iLipidUP, iLipidDN, iOxphos, iTCA, iGlutamin), ]

# CXCR5.pos.Im <- apply(tmpMtrx[ ,3:5], 1, mean)-apply(tmpMtrx[ ,1:2], 1, mean)
# CXCR5.neg.Im <- apply(tmpMtrx[ ,6:8], 1, mean)-apply(tmpMtrx[ ,1:2], 1, mean)

# CXCR5.pos.He <- tmpMtrx[ , 11]-(tmpMtrx[ ,9] + 0.1)
# CXCR5.neg.He <- tmpMtrx[ , 12]-(tmpMtrx[ ,9] + 0.1)
# Teff.He <- tmpMtrx[ , 10]-(tmpMtrx[ ,9] + 0.1)

# CXCR5.pos.Leong <- apply(tmpMtrx[ ,13:15], 1, mean)-apply(tmpMtrx[ ,19:20], 1, mean)
# CXCR5.neg.Leong <- apply(tmpMtrx[ ,16:18], 1, mean)-apply(tmpMtrx[ ,19:20], 1, mean)

# pltMtrx <- cbind(CXCR5.pos.Im, CXCR5.neg.Im, 
				 # CXCR5.pos.He, CXCR5.neg.He, Teff.He,
				 # CXCR5.pos.Leong, CXCR5.neg.Leong))


setwd("/home/hbolouri/TCE/")
library(pheatmap)
pdf("metabolic_heatmap_Leong.pdf", width=10, height=15)
j <- grep("Leong", colnames(tmpMtrx))[c(7:8, 1:3, 4:6)]
pheatmap(tmpMtrx[ , j], 
		 fontsize_col=15, fontsize_row=15, scale="row",
		 cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()













