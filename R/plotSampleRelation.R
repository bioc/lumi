`plotSampleRelation` <-
function(x, selProbe=NULL, cv.Th=0.1, standardize=TRUE, method=c('cluster', 'mds'), dimension=c(1,2), color=NULL) {
	if (is(x, 'ExpressionSet')) {
		dataMatrix <- exprs(x)
	} else if (is.matrix(x)) {
		dataMatrix <- x
	} else {
		stop('The class of "x" should be matrix or LumiBatch!')
	}
	
	if (standardize) dataMatrix <- scale(dataMatrix)

	if (is.null(selProbe)) {
		## Filter the genes with most of the experiments "Absent"
		cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
		probeList <- rownames(dataMatrix)
		selProbe <- probeList[cv.gene > cv.Th]
	}

	dd <- dist(t(dataMatrix[selProbe,]))
	geneNum <- length(selProbe)
	main <- paste('Clusters of the samples based on', geneNum, 'genes with sd/mean >', cv.Th)
	
	method <- match.arg(method)
	if (method == 'cluster') {
		hc = hclust(dd, 'ave')
		plot(hc, xlab='Sample', main=main)
	} else {
		## Multi-Dimension Scaling
		a1 <- cmdscale(dd, k=max(dimension))
		if (is.null(color)) {
			color <- 1
		} else {
			if (!is.numeric(color)) {
				allColor <- colors()
				if (!all(is.element(color, allColor))) {
					color <- as.numeric(factor(color))
				} 
			}
		}
		plot(a1[,dimension[1]],a1[,dimension[2]], type='n', xlab=paste('Principle component', dimension[1]),ylab=paste('Principle component', dimension[2]), main=main)
		text(a1[,dimension[1]],a1[,dimension[2]], col=color, labels=colnames(dataMatrix), cex=1)
	}
	attr(dd, 'geneNum') <- length(selProbe)
	attr(dd, 'threshold') <- cv.Th
	return(invisible(dd))	
}
