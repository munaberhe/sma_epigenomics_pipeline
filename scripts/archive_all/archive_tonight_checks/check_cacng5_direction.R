df1 <- read.csv('results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv')
hit1 <- df1[df1$SYMBOL=='CACNG5' & !is.na(df1$SYMBOL),]
hit1 <- hit1[order(hit1$pValue),][1,]
cat('ASO_CTRL_vs_Scramble_CTRL (cond_a=ASO_CTRL, cond_b=Scramble_CTRL):\n')
cat('  proportion1 (ASO_CTRL)      =', hit1$proportion1, '\n')
cat('  proportion2 (Scramble_CTRL) =', hit1$proportion2, '\n')
cat('  regionType =', hit1$regionType, '(gain=hypo in cond_a, loss=hyper in cond_a)\n')
cat('  Interpretation: ASO_CTRL is',
    ifelse(hit1$proportion1 < hit1$proportion2, 'LOWER (hypomethylated)', 'HIGHER (hypermethylated)'),
    'than Scramble_CTRL\n\n')

df2 <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
hit2 <- df2[df2$SYMBOL=='CACNG5' & !is.na(df2$SYMBOL),]
hit2 <- hit2[order(hit2$pValue),][1,]
cat('ASO_VPA_vs_Scramble_CTRL (cond_a=ASO_VPA, cond_b=Scramble_CTRL):\n')
cat('  proportion1 (ASO_VPA)       =', hit2$proportion1, '\n')
cat('  proportion2 (Scramble_CTRL) =', hit2$proportion2, '\n')
cat('  regionType =', hit2$regionType, '\n')
cat('  Interpretation: ASO_VPA (combination) is',
    ifelse(hit2$proportion1 < hit2$proportion2, 'LOWER (hypomethylated)', 'HIGHER (hypermethylated)'),
    'than Scramble_CTRL\n')
