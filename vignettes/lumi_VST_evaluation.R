### R code from vignette source 'lumi_VST_evaluation.Rnw'

###################################################
### code chunk number 1: libraries
###################################################
library("lumi")
library("vsn")
library("genefilter")
library("RColorBrewer")
library("limma")
library("lumiBarnes")
set.seed(0xbadbeef)
## Load the Barnes data set
data("lumiBarnes")


###################################################
### code chunk number 2: figFittingLinear
###################################################
temp <- lumiT(lumiBarnes[,1], fitMethod='linear', ifPlot=TRUE)


###################################################
### code chunk number 3: figFittingQuad
###################################################
temp <- lumiT(lumiBarnes[,1], fitMethod='quadratic', ifPlot=TRUE)


###################################################
### code chunk number 4: load
###################################################
## Select the blood and placenta samples
selChip = !is.na(lumiBarnes$pctBlood)
x.lumi <- lumiBarnes[, selChip]
presentCount <- detectionCall(x.lumi)

## Since the Barnes data was not background removed, we will do background adjustment first.
## The background estimation will be based on the control probe information. 
##  As the old version lumiBarnes library does not include controlData slot, we will check it first.
if (nrow(controlData(x.lumi)) == 0) {
	## We will use the control probe information in the example.lumi in the updated lumi package
	data(example.lumi)
	controlData(x.lumi) <- controlData(example.lumi)
}

x.lumi <- lumiB(x.lumi, method='bgAdjust')

repl1 <- which(x.lumi$replicate=="A")
repl2 <- which(x.lumi$replicate=="B")
stopifnot(sum(selChip)==12L, length(repl1)==6L, length(repl2)==6L)


###################################################
### code chunk number 5: preprocess
###################################################
## VST transform and Quantile normalization
x.lumi.vst <- lumiT(x.lumi)
x.lumi.vst.quantile <- lumiN(x.lumi.vst, method='quantile')


## log2 transform and Quantile normalization
x.lumi.log <- lumiT(x.lumi, method='log2')
x.lumi.log.quantile <- lumiN(x.lumi.log, method='quantile')


## VSN normalization: use lts.quantile=0.5 since in the blood/placenta
##  comparison more genes are differentially expressed than what is
##   expected by the default of 0.9.
x.lumi.vsn <- lumiN(x.lumi, method='vsn', lts.quantile=0.5)


## Add the vsn based on technical replicates
vsn.pair <- exprs(x.lumi)
cor.i <- NULL
for(i in 1:length(repl1)) {
	vsn.pair[, c(i, i+length(repl1))] <- exprs(vsn2(vsn.pair[, c(repl1[i], repl2[i])], verbose=FALSE))
}
# vsn.quantile <- normalize.quantiles(vsn.pair)
# rownames(vsn.quantile) <- rownames(vsn.pair)
# colnames(vsn.quantile) <- colnames(vsn.pair)


normDataList <- list('VST-Quantile'=exprs(x.lumi.vst.quantile), 
                    'Log2-Quantile'=exprs(x.lumi.log.quantile),
	            'VSN'=exprs(x.lumi.vsn))		# , 'VSN-Quantile'=vsn.quantile)

## scatter plots: 
## pairs(exprs(x.lumi.vsn), panel=function(...){par(new=TRUE);smoothScatter(..., nrpoints=0)})


###################################################
### code chunk number 6: chipCorList
###################################################
## Check the correlation between technique replicates
tempDataList <- c(normDataList, list(vsn.pair))
names(tempDataList) <- c(names(normDataList), 'VSN-techReplicate')
chipCorList <- matrix(as.numeric(NA), nrow=length(repl1), ncol=length(tempDataList))
colnames(chipCorList) <- names(tempDataList)
for (i in seq(along= tempDataList))
  for (j in seq(along=repl1))
    chipCorList[j,i] = cor(tempDataList[[i]][, c(repl1[j], repl2[j])])[1,2]


###################################################
### code chunk number 7: figboxplot
###################################################
labels <- colnames(chipCorList)
## set the margin of the plot
mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
oldpar = par(xaxt='n', mar=mar)
boxplot(chipCorList ~ col(chipCorList),  xlab='', 
        ylab='Correlation between technique replicate chips',
        col='skyblue')
par(xaxt='s')
axis(1, at=1:ncol(chipCorList), labels=labels, tick=TRUE, las=2)
par(oldpar)


###################################################
### code chunk number 8: figmeanSdPlot
###################################################
## select the technique replicates
selChip <-  c(repl1[1],repl2[1]) 
oldpar <- par(mfrow=c(length(normDataList) + 1,1))
for (i in 1:length(normDataList)) {
  meanSdPlot(normDataList[[i]][, selChip], ylab='Standard deviation')
}
meanSdPlot(vsn.pair[, selChip], ylab='Standard deviation')
par(oldpar)


###################################################
### code chunk number 9: fstatistic
###################################################
fac <- factor(paste(x.lumi$pctBlood, x.lumi$pctPlacenta, sep=":"))
rf <- lapply(normDataList, function(x) {
filtered.x = x[presentCount > 0,]
ftest.x = rowFtests(filtered.x, fac=fac)
ftest.x$IDs <- rownames(filtered.x)
return(ftest.x)
})
ef <- sapply( rf, function(x) ecdf(x$p.value))


###################################################
### code chunk number 10: figfstat
###################################################
pcol <- seq(along= normDataList) 
plty <- (1:(1 + length(normDataList))) [-3] 
plwd <- 1.5
x <- seq(0, 0.05, by=0.0001); x = x[-1]
plot(x, ef[[1]](x), type='l', lwd=plwd, lty=plty[1], col=pcol[1],
	main="Cumulative distribution of F-test p-value", xlab="F-test p-value", ylab="Empirical probability", log='x')
for (i in 2:length(ef)) {
	lines(x, ef[[i]](x), lwd=plwd, lty=plty[i], col=pcol[i])
}
legend(0.01, 0.25, names(normDataList), lwd=plwd, lty=plty, col=pcol)


###################################################
### code chunk number 11: Expression and dilution profile correlation
###################################################
modelProfile1 <- c(100, 95, 75, 50, 25, 0, 100, 95, 75, 50, 25, 0)
corrList <- lapply(normDataList, function(x) {
	x <- x[presentCount > 0, ]
	corr1 <- apply(x, 1, cor, y=modelProfile1)
	return(corr1)
	} )


###################################################
### code chunk number 12: histCorrelation
###################################################
freqMatrix <- NULL
breaks <- NULL
for (i in 1:length(corrList)) {
	hist.i <- hist(abs(corrList[[i]]), 30, plot=FALSE)
	breaks <- cbind(breaks, hist.i$breaks)
	freqMatrix <- cbind(freqMatrix, hist.i$counts)
}
freqMatrix <- rbind(freqMatrix, freqMatrix[nrow(freqMatrix),])
matplot(breaks, freqMatrix, type='s', lty=plty, col=pcol, lwd=plwd, ylab='Frequency', xlab='Absolute values of correlation coefficients')
legend(x=0.1, y=2800, legend=names(normDataList), lty=plty, col=pcol, lwd=plwd)


###################################################
### code chunk number 13: fstatistic concordance percentage
###################################################
topNumList <- seq(50, 3000, by=100)
corTh <- 0.8
highCorrNumMatrix <- NULL
for (i in 1:length(rf)) {
	probeList <- rf[[i]]$IDs
	ordProbe.i <- probeList[order(abs(rf[[i]]$p.value), decreasing=FALSE)]
	corr1 <- corrList[[i]]
	matchNum.j <- NULL
	for (topNum.j in topNumList) {
		topProbe.j <- ordProbe.i[1:topNum.j]
		matchNum.j <- c(matchNum.j, length(which(abs(corr1[topProbe.j]) > corTh)))
	}
	highCorrNumMatrix <- cbind(highCorrNumMatrix, matchNum.j)
}
rownames(highCorrNumMatrix) <- topNumList
colnames(highCorrNumMatrix) <- names(rf)


###################################################
### code chunk number 14: figfstatCor
###################################################
matplot(topNumList, (100 * highCorrNumMatrix/(topNumList %*% t(rep(1,ncol(highCorrNumMatrix))))),
 	type='l', xlab='Number of most significant probes by ranking their p-values (F-test)', 
	ylab='Percentage of concordant probes (%)', 
        lty=plty, col=pcol, lwd=plwd, ylim=c(50,100))
legend(x=2000, y=70, legend=colnames(highCorrNumMatrix), lty=plty, col=pcol, lwd=plwd)


###################################################
### code chunk number 15: fitList.limma
###################################################
## Select the comparing chip index
sampleInfo <- pData(phenoData(x.lumi))
sampleType <- paste(sampleInfo[,'pctBlood'], sampleInfo[,'pctPlacenta'], sep=':')
sampleType <- paste('c', sampleType, sep='')
## Comparing index
## used in the paper (the most challenging comparison):
compareInd <- c(repl1[1:2], repl2[1:2])  	
compareType <- sampleType[compareInd]
fitList.limma <- NULL
for (i in 1:length(normDataList)) {
	selDataMatrix <- normDataList[[i]]
	selDataMatrix <- selDataMatrix[presentCount > 0, ]
	selProbe <- rownames(selDataMatrix)
	compareMatrix <- selDataMatrix[, compareInd]
	
	design <- model.matrix(~ 0 + as.factor(compareType))
	colnames(design) <- c('A', 'B')
	fit1 <- lmFit(compareMatrix, design)
	contMatrix <- makeContrasts('A-B'=A - B, levels=design)
	fit2 <- contrasts.fit(fit1, contMatrix)
	fit <- eBayes(fit2)
	fitList.limma <- c(fitList.limma, list(fit))
}
names(fitList.limma) <- names(normDataList)


###################################################
### code chunk number 16: highCorrNumMatrix
###################################################
## Check the correlation of the top differentiated probes based on the limma results 
## rank the probes based on the p-values of limma result
fitList <- fitList.limma
topNumList <- c(30, seq(35, 1000, by=30))
corTh <- 0.8
highCorrNumMatrix <- NULL
for (i in 1:length(fitList)) {
	probeList <- rownames(fitList[[i]]$p.value)
	ordProbe.i <- probeList[order(abs(fitList[[i]]$p.value[,1]), decreasing=FALSE)]
	profileMatrix <- normDataList[[i]][ordProbe.i, ]

	modelProfile1 <- c(100, 95, 75, 50, 25, 0, 100, 95, 75, 50, 25, 0)
	corr1 <- apply(profileMatrix, 1, cor, y=modelProfile1)
	names(corr1) <- ordProbe.i
	matchNum.j <- NULL
	for (topNum.j in topNumList) {
		topProbe.j <- ordProbe.i[1:topNum.j]
		matchNum.j <- c(matchNum.j, length(which(abs(corr1[topProbe.j]) > corTh)))
	}
	highCorrNumMatrix <- cbind(highCorrNumMatrix, matchNum.j)
}
rownames(highCorrNumMatrix) <- topNumList
colnames(highCorrNumMatrix) <- names(fitList)


###################################################
### code chunk number 17: figLimmaConcordance
###################################################
matplot(topNumList, (100 * highCorrNumMatrix/(topNumList %*% t(rep(1,ncol(highCorrNumMatrix))))),
 	type='l', xlab='Number of most significant probes by ranking their p-values', 
	ylab='Percentage of concordant probes (%)', 
        lty=plty, col=pcol, lwd=plwd, ylim=c(0,100))
legend(x=700, y=50, legend=colnames(highCorrNumMatrix), lty=plty, col=pcol, lwd=plwd)


###################################################
### code chunk number 18: sessionInfo
###################################################
toLatex(sessionInfo())


