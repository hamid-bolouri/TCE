S['TCR'] <- ( S['TCR']  & !S['PD1'] )
S['IL2'] = S['IL21'] = S['IL12'] <- S['TCR']
S['FOXO1'] <- ( !S['IL2'] & !S['AKT'] ) & ( prevS2['TCR'] | prevS2['STAT3'] )
S['AKT'] <- ( S['TCR'] & !S['PD1'] )
S['NR4A1'] <- S['NFATC1'] | S['PD1']
S['PD1'] <- ( ( (prevS2['AP1'] & prevS1['AP1'] & S['AP1']) | S['FOXO1'] ) | S['NFATC1'] )
S['AP1.DNA'] <- ( prevS1['AP1'] | prevS2['AP1'] | S['AP1'] ) & !S['BCL6'] 
S['BLIMP1'] <- ( S['BATF.IRF4'] | S['IL2'] | S['AP1.DNA'] ) 
S['IRF4'] <- ( S['NFkB'] & !S['NR4A1'] )
S['BCL6'] <- ( ( S['FOXO1'] | S['BATF']) ) & ! S['IRF4']
S['BCL6'] -> S['TCF1']
S['BATF.IRF4'] <- ( S['BATF'] & S['TCR'] )
S['BATF'] <- ( S['IL12'] | S['IL21'] )
S['IL2'] <- ( S['IL2'] & !S['FOXO1'] )
S['NFATC1'] <- ( prevS2['NFATC1.med'] | prevS1['NFATC1'] )
S['NFATC1.med'] <- prevS2['NFATC1.lo']
S['NFATC1.lo'] <- prevS2['NFATC2'] 
S['IFNg'] <- ( S['NFATC2'] & S['AP1'] & !S['TCF1'] )
S['TCR'] -> S[c( 'NFkB', 'NFATC2', 'AP1' )]
S['AP1'] <- 1 # AP-1 stays on once turned on
S['AKT'] -> S['mTOR'] -> S['glycolysis']
S['STAT3'] <- prevS1['IL21'] 
