# Set this before run:
initState = "acute"	# TCF1hi Tex naive

TCR = NFkB = NFATC2 = AP1 = IRF4 = IFNg = AKT = mTOR = glycolysis = NFATC1.med = NFATC1.lo = NFATC1 = NR4A1 = FOXO1 = AP1.DNA = BLIMP1 = BCL6 = TCF1 = BATF = IL2 = IL12 = IL21 = PD1 = BATF.IRF4 = STAT3 = 0.0
startState <- sort(c(TCR,NFkB,NFATC2,AP1,IRF4,IFNg,AKT,mTOR,glycolysis,NFATC1.med,NFATC1.lo,NFATC1,NR4A1,FOXO1,AP1.DNA,BLIMP1,BCL6,TCF1,BATF,IL2,IL12,IL21,PD1,BATF.IRF4,STAT3))
names(startState) <- sort(c("TCR","NFkB","NFATC2","AP1","IRF4","IFNg","AKT","mTOR","glycolysis","NFATC1.med","NFATC1.lo","NFATC1","NR4A1","FOXO1","AP1.DNA","BLIMP1","BCL6","TCF1","BATF","IL2","IL12","IL21","PD1","BATF.IRF4","STAT3"))

prevS2 = prevS1 = NULL
prevS2 = prevS1 = S = startState

if (initState == "acute") {
	S[c("TCR", "IL2", "IL12", "IL21")] <- 1
}

allStates = S; TIME = 0
while(!identical(prevS2, S)) {

	TIME <- TIME + 1
	if (TIME == 2) prevS1 = S else if (TIME > 2) {
		prevS2 <- prevS1
		prevS1 <- S
	}
	
	# NB evaluation order matters!
	
	S['TCR'] <- ( S['TCR']  & !S['PD1'] )
	S['IL2'] = S['IL21'] = S['IL12'] <- S['TCR']
	
	S['FOXO1'] <- ( !S['IL2'] & !S['AKT'] ) & 
				  ( prevS2['TCR'] | prevS2['STAT3'] )
	S['AKT'] <- ( S['TCR'] & !S['PD1'] )

	S['NR4A1'] <- S['NFATC1'] | S['PD1']
	S['PD1'] <- ( ( (prevS2['AP1'] & prevS1['AP1'] & S['AP1']) | 
				   S['FOXO1'] ) | S['NFATC1'] )
	
	S['AP1.DNA'] <- ( prevS1['AP1'] | prevS2['AP1'] | S['AP1'] ) & 
					  !S['BCL6'] 
	S['BLIMP1'] <- ( S['BATF.IRF4'] | S['IL2'] | S['AP1.DNA'] ) 

	S['IRF4'] <- ( S['NFkB'] & !S['NR4A1'] )
	S['BCL6'] <- ( ( S['FOXO1'] | S['BATF']) ) & 
				   ! S['IRF4']
	S['BCL6'] -> S['TCF1']

	S['BATF.IRF4'] <- ( S['BATF'] & S['TCR'] )
	S['BATF'] <- ( S['IL12'] | S['IL21'] )

	S['IL2'] <- ( S['IL2'] & !S['FOXO1'] )

	S['NFATC1'] <- ( prevS2['NFATC1.med'] | prevS1['NFATC1'] )
	S['NFATC1.med'] <- prevS2['NFATC1.lo']
	S['NFATC1.lo'] <- prevS2['NFATC2'] 

	S['IFNg'] <- ( S['NFATC2'] & S['AP1'] & !S['TCF1'] )
	
	S['TCR'] -> S[c( 'NFkB', 'NFATC2', 'AP1' )]
	S['AP1'] <- S['TCR']  | S['AP1'] # AP-1 stays on once turned on
	S['AKT'] -> S['mTOR'] -> S['glycolysis']
	S['STAT3'] <- prevS1['IL21'] 

	allStates <- rbind(allStates, S)
}

setwd("~/TCE")
write.table(allStates, file="allStates.txt", sep="\t", quote=FALSE, row.names=FALSE)
