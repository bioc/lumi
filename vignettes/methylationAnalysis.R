### R code from vignette source 'methylationAnalysis.Rnw'

###################################################
### code chunk number 1: load library
###################################################
library(lumi)


###################################################
### code chunk number 2: load example dataset
###################################################
## load example data (a methyLumiM object)
data(example.lumiMethy)
## summary of the example data
example.lumiMethy
## print sample Names
sampleNames(example.lumiMethy)


###################################################
### code chunk number 3: sampleRelation
###################################################
plotSampleRelation(example.lumiMethy, method='mds', cv.Th=0) 


###################################################
### code chunk number 4: sampleRelationTree
###################################################
plotSampleRelation(example.lumiMethy, method='cluster', cv.Th=0) 


###################################################
### code chunk number 5: load example titration dataset
###################################################
## load the tritration data (a methyLumiM object)
data(example.methyTitration)
## summary of the example data
example.methyTitration
## print sample Names
sampleNames(example.methyTitration)


###################################################
### code chunk number 6: sampleRelationTitration
###################################################
plotSampleRelation(example.methyTitration, method='mds', cv.Th=0) 


###################################################
### code chunk number 7: densityMTitration
###################################################
## plot the density
density(example.methyTitration, xlab="M-value")	


###################################################
### code chunk number 8: densityM
###################################################
## specify the colors of control and treatment samples
sampleColor <- rep(1, ncol(example.lumiMethy))
sampleColor[grep("Treat", sampleNames(example.lumiMethy))] <- 2 

density(example.lumiMethy, col=sampleColor, xlab="M-value")	## plot the density


###################################################
### code chunk number 9: boxplotM
###################################################
## Because the distribution of M-value has two modes, we use a boxplot different from regular ones
boxplot(example.lumiMethy)


###################################################
### code chunk number 10: densityColorBiasBoth
###################################################
plotColorBias1D(example.lumiMethy)


###################################################
### code chunk number 11: densityColorBiasMethy
###################################################
plotColorBias1D(example.lumiMethy, channel='methy')


###################################################
### code chunk number 12: densityColorBiasUnmethy
###################################################
plotColorBias1D(example.lumiMethy, channel='unmethy')


###################################################
### code chunk number 13: boxplotColorBiasMethy
###################################################
boxplotColorBias(example.lumiMethy, channel='methy')


###################################################
### code chunk number 14: boxplotColorBiasUnmethy
###################################################
boxplotColorBias(example.lumiMethy, channel='unmethy')


###################################################
### code chunk number 15: color balance summary
###################################################
## summary of color balance information of individual samples
colorBiasSummary(example.lumiMethy[,1:8], channel='methy')


###################################################
### code chunk number 16: scatterColorBias1
###################################################
plotColorBias2D(example.lumiMethy, selSample=1, cex=2)


###################################################
### code chunk number 17: boxplotColorBiasSum
###################################################
boxplotColorBias(example.lumiMethy, channel='sum')


###################################################
### code chunk number 18: densityColorBiasSum
###################################################
plotColorBias1D(example.lumiMethy, channel='sum')


###################################################
### code chunk number 19: densityIntensity
###################################################
density(estimateIntensity(example.lumiMethy), xlab="log2(CpG-site Intensity)")


###################################################
### code chunk number 20: boxplotIntensity
###################################################
boxplot(estimateIntensity(example.lumiMethy))


###################################################
### code chunk number 21: pairsColor
###################################################
## get the color channel information
colorChannel <- as.character(pData(featureData(example.lumiMethy))[, "COLOR_CHANNEL"])
## replace the "Red" and "Grn" as color names defined in R
colorChannel[colorChannel == 'Red'] <- 'red'
colorChannel[colorChannel == 'Grn'] <- 'green'
## select a subet of sample for pair plot
selSample <- c( "Ctrl1", "Ctrl1.rep", "Treat1", "Treat1.rep")
## plot pair plot with the dots in scatter plot colored based on the color channels
pairs(estimateIntensity(example.lumiMethy[, selSample]), dotColor= colorChannel, main="Pair plot of CpG-site Intensity")


###################################################
### code chunk number 22: color balance adjustment
###################################################
## summary of color balance information of individual samples
lumiMethy.c.adj <- lumiMethyC(example.lumiMethy)


###################################################
### code chunk number 23: densityColorBiasSumAdj
###################################################
plotColorBias1D(lumiMethy.c.adj, channel='sum')


###################################################
### code chunk number 24: boxplotColorBiasSumAdj
###################################################
boxplotColorBias(lumiMethy.c.adj, channel='sum')


###################################################
### code chunk number 25: scatterColorBias1Adj
###################################################
## plot the color balance adjusted scatter plot of two color channels
plotColorBias2D(lumiMethy.c.adj, selSample=1, cex=2)


###################################################
### code chunk number 26: pairsColorAdj
###################################################
## plot pairwise plot after color balance adjustment 
pairs(estimateIntensity(lumiMethy.c.adj[, selSample]), dotColor= colorChannel, main="Pair plot of CpG-site Intensity after color balance adjustment")


###################################################
### code chunk number 27: background adjustment
###################################################
##separately adjust backgrounds of two color channels
lumiMethy.b.adj <- lumiMethyB(example.lumiMethy, method="bgAdjust2C", separateColor=TRUE)

##background adjustment of individual samples
lumiMethy.bc.adj <- lumiMethyB(lumiMethy.c.adj, method="bgAdjust2C")


###################################################
### code chunk number 28: bgDensityMethy
###################################################
## plot the background mode of methylated probe data of first five example samples
plotColorBias1D(example.lumiMethy[,1:5], channel='methy', xlim=c(-1000,5000), logMode=FALSE)


###################################################
### code chunk number 29: bgAdjDensityMethy
###################################################
## plot the background mode of methylated probe data of first five example samples
plotColorBias1D(lumiMethy.b.adj [,1:5], channel='methy', xlim=c(-1000,5000), logMode=FALSE)


###################################################
### code chunk number 30: bcAdjDensityMethy
###################################################
## plot the background mode of methylated probe data of first five example samples
plotColorBias1D(lumiMethy.bc.adj [,1:5], channel='methy', xlim=c(-1000,5000), logMode=FALSE)


###################################################
### code chunk number 31: Normalization
###################################################
## Perform quantile normalization based on color balance adjusted data
lumiMethy.c.quantile <- lumiMethyN(lumiMethy.c.adj, method='quantile')
## Perform SSN normalization based on color balance adjusted data
lumiMethy.c.ssn <- lumiMethyN(lumiMethy.c.adj, method='ssn')


###################################################
### code chunk number 32: sampleRelationTreeNormalized
###################################################
plotSampleRelation(lumiMethy.c.quantile, method='cluster', cv.Th=0) 


###################################################
### code chunk number 33: densityQuantile
###################################################
## plot the density of M-values after quantile normalization
density(lumiMethy.c.quantile, col= sampleColor, main="Density plot after quantile normalization")


###################################################
### code chunk number 34: densityIntensityNormalizedQ
###################################################
density(estimateIntensity(lumiMethy.c.quantile), col= sampleColor,  xlab="log2(CpG-site Intensity)")


###################################################
### code chunk number 35: boxplotColorBiasNormalized
###################################################
boxplotColorBias(lumiMethy.c.quantile, channel='sum')


###################################################
### code chunk number 36: pairMNormalize
###################################################
## select a subet of sample for pair plot
selSample <- c( "Ctrl1", "Ctrl1.rep", "Treat1", "Treat1.rep")
## plot pair plot with the dots in scatter plot colored based on the color channels
pairs(lumiMethy.c.quantile[, selSample], dotColor= colorChannel, main="Pair plot of M-value after normalization")


###################################################
### code chunk number 37: user defined preprocessing functions (eval = FALSE)
###################################################
## ## suppose "userB" is a user defined background adjustment method
##  lumiMethy.b.adj <- lumiMethyB(example.lumiMethy, method=userB, separateColor=TRUE) # not Run
## 
## ## suppose "userC" is a user defined color balance adjustment method
##  lumiMethy.c.adj <- lumiMethyC(example.lumiMethy, method=userC, separateColor=TRUE) # not Run
## 
## ## suppose "userN" is a user defined probe level normalization method
##  lumiMethy.c.n <- lumiMethyN(lumiMethy.c.adj, method= userN, separateColor=TRUE) # not Run


###################################################
### code chunk number 38: Color channel information
###################################################
## retrieve the featureData information
ff <- pData(featureData(example.lumiMethy))

## show the color channel information
head(ff)

## add user provided color channel information if it is not existed in the featureData
# example.lumiMethy <- addAnnotationInfo(example.lumiMethy, lib="IlluminaHumanMethylation27k.db")


###################################################
### code chunk number 39: options to separately process each color channel (eval = FALSE)
###################################################
## ## suppose "userB" is a user defined background adjustment method
##  lumiMethy.b.adj <- lumiMethyB(example.lumiMethy,  separateColor=TRUE) 
## 
## ## suppose "userC" is a user defined color balance adjustment method
##  lumiMethy.c.adj <- lumiMethyC(example.lumiMethy,  separateColor=TRUE) 
## 
## ## suppose "userN" is a user defined probe level normalization method
##  lumiMethy.c.n <- lumiMethyN(lumiMethy.c.adj,  separateColor=TRUE) 


###################################################
### code chunk number 40: Estimate detection call of a CpG site (eval = FALSE)
###################################################
## ## Estimate the detection call of a CpG site
## presentCount <- detectionCall(example.lumiMethy)


###################################################
### code chunk number 41: Preprocessing methylation microarray (eval = FALSE)
###################################################
## 
## library(lumi)
## ## specify the file name
## # fileName <- 'Example_Illumina_Methylation_profile.txt' 
## ## load the data
## # example.lumiMethy <- lumiMethyR(fileName, lib="IlluminaHumanMethylation27k.db")		# Not Run
## 
## ## Quality and color balance assessment
## data(example.lumiMethy)
## ## summary of the example data
## example.lumiMethy
## 
## ## preprocessing and quality control after normalization
## plotColorBias1D(example.lumiMethy, channel='sum')
## boxplotColorBias(example.lumiMethy, channel='sum')
## ## select interested sample to further check color balance in 2D scatter plot
## plotColorBias2D(example.lumiMethy, selSample=1)
## 
## ## Color balance adjustment between two color channels
## lumiMethy.c.adj <- lumiMethyC(example.lumiMethy)
## ## Check color balance after color balance adjustment
## boxplotColorBias(lumiMethy.c.adj, channel='sum')
## 
## ## Background adjustment is skipped because the SSN normalization includes background adjustment
## 
## ## Normalization
## ## Perform SSN normalization based on color balance adjusted data
## lumiMethy.c.ssn <- lumiMethyN(lumiMethy.c.adj, method='ssn')
## 
## ## Or we can perform quantile normalization based on color balance adjusted data
## # lumiMethy.c.q <- lumiMethyN(lumiMethy.c.adj, method='quantile')
## 
## ## plot the density of M-values after SSN normalization
## density(lumiMethy.c.ssn, main="Density plot of M-value after SSN normalization")
## ## comparing with the density of M-values before normalization
## density(example.lumiMethy, main="Density plot of M-value of the raw data")
## 
## ## output the normlized M-value as a Tab-separated text file
## write.exprs(lumiMethy.c.ssn, file='processedMethylationExampleData.txt')


###################################################
### code chunk number 42: sessionInfo
###################################################
toLatex(sessionInfo())


