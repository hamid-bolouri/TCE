# Set this before run:
initState = "acute"	# TCF1hi Tex naive

TCR = NFkB = NFATC2 = AP1 = IRF4 = IFNg = AKT = mTOR = glycolysis = NFATC1.med = NFATC1.lo = NFATC1 = NR4A1 = FOXO1 = AP1.DNA = BLIMP1 = BCL6 = TCF1 = BATF = IL2 = IL12 = IL21 = PD1 = BATF.IRF4 = STAT3 = 0.0
startState <- sort(c(TCR,NFkB,NFATC2,AP1,IRF4,IFNg,AKT,mTOR,glycolysis,NFATC1.med,NFATC1.lo,NFATC1,NR4A1,FOXO1,AP1.DNA,BLIMP1,BCL6,TCF1,BATF,IL2,IL12,IL21,PD1,BATF.IRF4,STAT3))
names(startState) <- sort(c("TCR","NFkB","NFATC2","AP1","IRF4","IFNg","AKT","mTOR","glycolysis","NFATC1.med","NFATC1.lo","NFATC1","NR4A1","FOXO1","AP1.DNA","BLIMP1","BCL6","TCF1","BATF","IL2","IL12","IL21","PD1","BATF.IRF4","STAT3"))

prevS2 = prevS1 = NULL
prevS2 = prevS1 = S = startState

if (initState == "acute") {
	startState[c("TCR", "IL2", "IL12", "IL21")] <- 1
}

texS <- rep(0, length(startState))
names(texS) <- names(startState)
texS[c("PD1","NFATC1","NR4A1","BLIMP1")] <- 1


setwd("~/TCE"); allSS = goodIndx = NULL
for (N in 1:10000) {
	TIME = 0
	S <- startState
	while(!identical(prevS2, S) & TIME < 100) {

		TIME <- TIME + 1
		if (TIME == 2) prevS1 = S else if (TIME > 2) {
			prevS2 <- prevS1
			prevS1 <- S
		}
		
		indx <- sample(1:22, 22)
		for (i in indx) {
			nextCMD <- rev(readLines("TexSimCmds.R", n=indx[i]))[1]
			eval(parse(text=nextCMD))
		}
	}
	goodSS <- 0
	if (identical(texS, S)) {
		goodSS <- 1
		goodIndx <- rbind(goodIndx, indx)
	}
	
	allSS <- c(allSS, goodSS)
	print(paste0(N, "   ", TIME, "   ", Sys.time()))	# (paste(S, collapse=" : "))
}
sum(allSS) # 28
setwd("~/TCE")
save(goodIndx, file="28updateOrdersGivingTexSS_in10Kruns.RData")

similarity = matrix(0, nrow=dim(goodIndx)[1], ncol=dim(goodIndx)[1])
for (i in 1:dim(goodIndx)[1]) {
	for (j in 1:dim(goodIndx)[2]) {
		if (i != j) {
			similarity[i,j] <- length(which(goodIndx[i, ] == 
											goodIndx[j, ]))
		}
	}
}
max(similarity)
table(unlist(similarity))
#   0   1   2   3   4   5
# 396 203 136  45   2   2


setwd("~/TCE")
goodIndx <- get(load("28updateOrdersGivingTexSS_in10Kruns.RData"))
rownames(goodIndx) <- paste0("sim", 1:dim(goodIndx)[1])
colnames(goodIndx) <- paste0("Rule", 1:dim(goodIndx)[2])

write.table(goodIndx, file="updateOrdersOf28RunsWithTexSS.txt",
			sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)

myColors <- rep(c("red", "blue", "green", "black", "cyan", 
			  "pink", "forestgreen", "pink", "gray",
			  "purple", "tan1"), 5)

library(circlize)
set.seed(123); n=22
data = data.frame(
    factor = sample(paste0("Rule", 1:n), 1000, replace=TRUE),
    x=runif(1000), y=runif(1000))
 
# Initialize the plot.
par(mar = c(1, 1, 1, 1) ) 
circos.initialize(factors=data$factor, x=data$x)
 
# Build the regions of track
circos.trackPlotRegion(factors=data$factor, y=data$y, 
					   bg.col=rep(rgb(0,0.25,1,0.4), n), 
					   bg.border=NA,
					   panel.fun = function(x, y) {
						  xlim = get.cell.meta.data("xlim")
						  ylim = get.cell.meta.data("ylim")
						  sector.name = get.cell.meta.data("sector.index")
						  circos.text(mean(xlim), ylim[1] + 0.25, cex=0.75,
						  sector.name, facing = "clockwise", 
						  niceFacing = TRUE, adj=c(0, 0.5), col="blue")
					   } )					   

# Add links
for (i in 1:(dim(goodIndx)[1])) {
	S <- colnames(goodIndx)[grep(paste0("e", i, "$"), colnames(goodIndx))]
	T <- colnames(goodIndx)[grep(paste0("e", i+1, "$"), colnames(goodIndx))]
	for (j in 1:(dim(goodIndx)[2] - 1)) {
		S <- colnames(goodIndx)[match(j, goodIndx[i, ])]
		T <- colnames(goodIndx)[match(j+1, goodIndx[i, ])]
		circos.link(S, 0.5, T, 0.5, h = (i/dim(goodIndx)[1])*0.5, col=myColors[i])
	}
}
dev.off()





