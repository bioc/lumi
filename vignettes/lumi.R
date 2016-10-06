### R code from vignette source 'lumi.Rnw'

###################################################
### code chunk number 1: load library
###################################################
library(lumi)


###################################################
### code chunk number 2: load example data
###################################################
## load example data (a LumiBatch object)
data(example.lumi)
## summary of the example data
example.lumi


###################################################
### code chunk number 3: summary of example data
###################################################
## summary of the quality control
summary(example.lumi, 'QC')


###################################################
### code chunk number 4: densityPlot
###################################################
plot(example.lumi, what='density')	## plot the density


###################################################
### code chunk number 5: densityPlot
###################################################
plot(example.lumi, what='density')	## plot the density


###################################################
### code chunk number 6: CDF
###################################################
 plotCDF(example.lumi)


###################################################
### code chunk number 7: pairPlot
###################################################
plot(example.lumi, what='pair')		## pairwise plots


###################################################
### code chunk number 8: lumi.Rnw:233-234
###################################################
plot(example.lumi, what='pair', smoothScatter=T)


###################################################
### code chunk number 9: MAplot
###################################################
## pairwise MAplot
plot(example.lumi, what='MAplot')	


###################################################
### code chunk number 10: lumi.Rnw:265-267
###################################################
## pairwise MAplot
plot(example.lumi, what='MAplot', smoothScatter=T)	


###################################################
### code chunk number 11: lumi.Rnw:279-281
###################################################
## density plot of coefficient of varience
plot(example.lumi, what='cv')	


###################################################
### code chunk number 12: lumi.Rnw:292-293
###################################################
plot(example.lumi, what='sampleRelation')


###################################################
### code chunk number 13: lumi.Rnw:310-311
###################################################
plot(example.lumi, what='sampleRelation', method='mds', color=c("01", "02", "01", "02"))


###################################################
### code chunk number 14: VST transform
###################################################
## Do default VST variance stabilizing transform
lumi.T <- lumiT(example.lumi)


###################################################
### code chunk number 15: plotVST
###################################################
trans <- plotVST(lumi.T)


###################################################
### code chunk number 16: lumi.Rnw:356-358
###################################################
matplot(log2(trans$untransformed), trans$transformed, type='l', main='compare VST and log2 transform', xlab='log2 transformed', ylab='vST transformed')
abline(a=0, b=1, col=2)


###################################################
### code chunk number 17: Default normalization
###################################################
## Do quantile between microarray normaliazation
lumi.N <- lumiN(lumi.T)


###################################################
### code chunk number 18: Customized normalization (eval = FALSE)
###################################################
## ## Do RSN between microarray normaliazation
## lumi.N <- lumiN(lumi.T, method='rsn')


###################################################
### code chunk number 19: QC after normalization
###################################################
## Do quality control estimation after normalization
lumi.N.Q <- lumiQ(lumi.N)

## summary of the quality control
summary(lumi.N.Q, 'QC')		## summary of QC


###################################################
### code chunk number 20: lumi.Rnw:405-406
###################################################
plot(lumi.N.Q, what='density')		## plot the density


###################################################
### code chunk number 21: lumi.Rnw:415-417
###################################################
plot(lumi.N.Q, what='boxplot')		## box plot
# boxplot(lumi.N.Q)


###################################################
### code chunk number 22: lumi.Rnw:426-427
###################################################
plot(lumi.N.Q, what='pair')			## pairwise plots


###################################################
### code chunk number 23: lumi.Rnw:436-437
###################################################
plot(lumi.N.Q, what='MAplot')		## plot the pairwise MAplot


###################################################
### code chunk number 24: lumi.Rnw:446-448
###################################################
## plot the sampleRelation using hierarchical clustering
plot(lumi.N.Q, what='sampleRelation')


###################################################
### code chunk number 25: lumi.Rnw:457-459
###################################################
## plot the sampleRelation using MDS
plot(lumi.N.Q, what='sampleRelation', method='mds', color=c("01", "02", "01", "02"))


###################################################
### code chunk number 26: Encapsulate the processing steps
###################################################
## Do all the default preprocessing in one step
lumi.N.Q <- lumiExpresso(example.lumi)


###################################################
### code chunk number 27: lumi.Rnw:476-478
###################################################
## Do all the preprocessing with customized settings
# lumi.N.Q <- lumiExpresso(example.lumi, normalize.param=list(method='rsn'))


###################################################
### code chunk number 28: inverse VST
###################################################
## Parameters of VST transformed LumiBatch object
names(attributes(lumi.T))
## VST parameters: "vstParameter"  and  "transformFun"
attr(lumi.T, 'vstParameter')
attr(lumi.T, 'transformFun')
## Parameters of VST transformed and RSN normalized LumiBatch object
names(attributes(lumi.N.Q))
## VSN "targetArray" , VST parameters: "vstParameter"  and  "transformFun"
attr(lumi.N.Q, 'vstParameter')
attr(lumi.N.Q, 'transformFun')
## After doing statistical analysis of the data, users can recover to the raw scale for the fold-change estimation.
## Inverse VST to the raw scale
lumi.N.raw <- inverseVST(lumi.N.Q)


###################################################
### code chunk number 29: retrieve normalized data
###################################################
dataMatrix <- exprs(lumi.N.Q)


###################################################
### code chunk number 30: filtering
###################################################
presentCount <- detectionCall(example.lumi)
selDataMatrix <- dataMatrix[presentCount > 0,]
if (require(lumiHumanAll.db) & require(annotate)) {
  selDataMatrix <- selDataMatrix[!is.na(getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')),]
}
probeList <- rownames(selDataMatrix)


###################################################
### code chunk number 31: Identify differentially expressed genes
###################################################
## Specify the sample type
sampleType <- c('100:0', '95:5', '100:0', '95:5')
if (require(limma)) {
	##  compare '95:5' and '100:0'
	design <- model.matrix(~ factor(sampleType))
	colnames(design) <- c('100:0', '95:5-100:0')
	fit <- lmFit(selDataMatrix, design)
	fit <- eBayes(fit)
	## Add gene symbols to gene properties
	if (require(lumiHumanAll.db) & require(annotate)) {
               geneSymbol <- getSYMBOL(probeList, 'lumiHumanAll.db')
               geneName <- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
               fit$genes <- data.frame(ID= probeList, geneSymbol=geneSymbol, geneName=geneName, stringsAsFactors=FALSE)
          }
	## print the top 10 genes
	print(topTable(fit, coef='95:5-100:0', adjust='fdr', number=10))
	
	## get significant gene list with FDR adjusted p.values less than 0.01
	p.adj <- p.adjust(fit$p.value[,2])		
	sigGene.adj <- probeList[ p.adj < 0.01]
	## without FDR adjustment
	sigGene <- probeList[ fit$p.value[,2] < 0.01]
} 


###################################################
### code chunk number 32: GO analysis (eval = FALSE)
###################################################
## if (require(GOstats) & require(lumiHumanAll.db)) {
## 
## 	## Convert the probe Ids as  Entrez Ids and make them unique
## 	sigLL <- unique(unlist(lookUp(sigGene,'lumiHumanAll.db','ENTREZID')))
## 	sigLL <- as.character(sigLL[!is.na(sigLL)])
## 	params <- new("GOHyperGParams",
##               geneIds= sigLL,
##               annotation="lumiHumanAll.db",
##               ontology="BP",
##               pvalueCutoff= 0.01,
##               conditional=FALSE,
##               testDirection="over")
##          
##          hgOver <- hyperGTest(params)
##          	
## 	## Get the p-values of the test
## 	gGhyp.pv <- pvalues(hgOver)
## 	
## 	## Adjust p-values for multiple test (FDR)
## 	gGhyp.fdr <- p.adjust(gGhyp.pv, 'fdr')
## 
## 	## select the Go terms with adjusted p-value less than 0.01
## 	sigGO.ID <- names(gGhyp.fdr[gGhyp.fdr < 0.01])
## 	
## 	## Here only show the significant GO terms of BP (Molecular Function)
## 	## 	For other categories, just follow the same procedure.
## 	sigGO.Term <- getGOTerm(sigGO.ID)[["BP"]]
## }


###################################################
### code chunk number 33: Create selected GO table (eval = FALSE)
###################################################
## if (require(GOstats) & require(lumiHumanAll.db) & require(xtable)) {
## 
## 	##get gene counts at each GO category
## 	gg.counts <- geneCounts(hgOver)[sigGO.ID]
## 	total.counts <- universeCounts(hgOver)[sigGO.ID]
## 
## 	ggt <- unlist(sigGO.Term)
## 	numCh <- nchar(ggt)
## 	ggt2 <- substr(ggt, 1, 17)
## 	ggt3 <- paste(ggt2, ifelse(numCh > 17, "...", ""), sep="")
## 	
## 	## output the significant GO categories as a table
## 	ggMat <- matrix(c(names(sigGO.Term), ggt3, signif(gGhyp.pv[sigGO.ID],5), gg.counts, total.counts),
##     		byrow=FALSE, nc=5, dimnames=list(1:length(sigGO.Term), c("GO ID",
##    		"Term", "p-value","Significant Genes No.", "Total Genes No.")))
## 	xtable.matrix(ggMat,
##   		caption="GO terms, p-values and counts.", label="ta:GOggterms")
## }


###################################################
### code chunk number 34: sessionInfo
###################################################
toLatex(sessionInfo())


