`lumiT` <-
function(x.lumi, method=c('vst', 'log2', 'cubicRoot'), ifPlot=FALSE, simpleOutput = TRUE, verbose = TRUE, ...) {
	# if (!is(x.lumi, 'LumiBatch')) stop('The object should be class "LumiBatch"!')
	method <- match.arg(method)
	if (is(x.lumi, 'eSet')) {
		if (method == 'vst') {
			if (is.null(se.exprs(x.lumi)))  stop('Slot se.exprs is required!')
			if (!all(dim(se.exprs(x.lumi)) == dim(exprs(x.lumi))))
				stop('Dimensions of slots exprs and se.exprs do not match!')
		}
	} else {
		stop('The object should be class "eSet" inherited!')
	}
	## the stdCorrection was moved to lumiR function on April 15, 2011.
	# if (!is(x.lumi, 'LumiBatch')) stdCorrection <- FALSE

	## check the negative values
	if (method == 'log2') {
		if (min(exprs(x.lumi), na.rm=TRUE) <= 0) {
			x.lumi <- lumiB(x.lumi, method='forcePositive')
		}
	}
	if (verbose) cat(paste('Perform', method, 'transformation ...\n'))
	history.submitted <- as.character(Sys.time())

    new.lumi <- x.lumi 
	exprs <- exprs(x.lumi) 
	if (max(exprs, na.rm=TRUE) < 50) {
		warning("Looks like the data has been log transformed! \n It should be in the raw scale! \n Please double check!\n")
	}
	if (method == 'log2') {
		exprs(new.lumi) <- log2(exprs)
	} else if (method == 'cubicRoot') {
		exprs(new.lumi) <- sign(exprs) * (abs(exprs))^1/3
	} else {
		se.exprs <- se.exprs(x.lumi)
		## the stdCorrection was moved to lumiR function on April 15, 2011.
#		if (stdCorrection & is(x.lumi, 'LumiBatch')) {
#			bn <- beadNum(x.lumi)
#			if (is.null(bn)) {
#				cat('No Standard Deviation correction was applied because of missing bead number information.\n')
#			} else {
#				se.exprs <- se.exprs * sqrt(bn)
#			}
#		}
		nArray <- ncol(exprs)
		if (is(x.lumi, 'LumiBatch')) {
			detectCall <- detectionCall(x.lumi, Th=0.01, type='matrix')			
		} else {
			detectCall <- NULL
		}
		transExpr <- NULL
		transPara <- NULL
		transFun <- NULL
		for (i in 1:nArray) {
			cat(as.character(Sys.time()), ", processing array ", i, "\n")
			u <- exprs[,i]; std <- se.exprs[,i]
			lowCutoff <- 1/3
			backgroundStd <- NULL
			if (!is.null(detectCall)) {
				ind <- which(detectCall[,i] == 'A')
				if (nrow(se.exprs) - length(ind) < 3000) {
					warning('Too few probes are detectable based on detection p-values!\n 
						Iteration method will be used for VST.')
					lowCutoff <- 1/4
				} else {
					backgroundStd <- mean(se.exprs[ind,i], na.rm=TRUE)
				}
			} else if (is(x.lumi, 'AffyBatch')) {
				lowCutoff <- 1/4
				# afbatch.i <- x.lumi[,i]
				# pmInd <- unlist(indexProbes(afbatch.i, 'pm'))
				# mmInd <- unlist(indexProbes(afbatch.i, 'mm'))
				# pm <- exprs(afbatch.i)[pmInd]
				# mm <- exprs(afbatch.i)[mmInd]
				# bgInd <- pmInd[mm > pm]
				# bgStd <- se.exprs[bgInd,i]
				# hh <- hist(bgStd, 1000, plot=FALSE)
				# Th <- hh$breaks[which.max(hh$counts) + 1] * 2
				# dd <- density(bgStd[bgStd < Th], na.rm=TRUE)
				# backgroundStd <- dd$x[which.max(dd$y)]
				# backgroundStd <- mean(bgStd)
			}
		    x <- vst(u=u, std=std, backgroundStd=backgroundStd, lowCutoff=lowCutoff, ifPlot=ifPlot, ...)
			transExpr <- cbind(transExpr, x)
			
			transPara <- rbind(transPara, attr(x, 'parameter'))
			transFun <- c(transFun, attr(x, 'transformFun'))
		}
         dimnames(transExpr) <- dimnames(exprs)
#        if (!is.null(transExpr)) rownames(transExpr) <- colnames(exprs(x.lumi))
		if (!is.null(transPara))	rownames(transPara) <- colnames(exprs(x.lumi))
		if (!is.null(transFun))	names(transFun) <- colnames(exprs(x.lumi))
		exprs(new.lumi) <- transExpr
	}
	colnames(exprs(new.lumi)) <- colnames(exprs(x.lumi))
	rownames(exprs(new.lumi)) <- rownames(exprs(x.lumi))

	if (simpleOutput) {
		if (is(x.lumi, 'LumiBatch')) {
			se.exprs(x.lumi) <- NULL
			beadNum(x.lumi) <- NULL
			# storage.mode <- storageMode(new.lumi)
			# if ("lockedEnvironment" == storage.mode) {
			# 	aData <- copyEnv(assayData(new.lumi))
			# 	# rm(list=c('se.exprs', 'detection', 'beadNum'), envir=aData)
			# 	rm(list=c('se.exprs', 'beadNum'), envir=aData)
			# 	lockEnvironment(aData, bindings = TRUE)
			# 	assayData(new.lumi) <- aData
			# } else {
			# 	aData <- assayData(new.lumi)
			# 	# rm(list=c('se.exprs', 'detection', 'beadNum'), envir=aData)
			# 	rm(list=c('se.exprs', 'beadNum'), envir=aData)
			# 	assayData(new.lumi) <- aData
			# }
		}
		if (is(x.lumi, 'AffyBatch')) {
			se.exprs(x.lumi) <- new('matrix')
		}
	}

	if (is(x.lumi, 'LumiBatch')) {
		# history tracking
		history.finished <- as.character(Sys.time())
		# history.command <- capture.output(print(match.call(lumiT)))
		history.command <- paste(deparse(match.call(lumiT)), collapse='') 

		if (is.null(new.lumi@history$lumiVersion)) new.lumi@history$lumiVersion <- rep(NA, nrow(new.lumi@history))
		lumiVersion <- packageDescription('lumi')$Version
		new.lumi@history<- rbind(new.lumi@history, data.frame(submitted=history.submitted, 
				finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	}

	if (method == 'vst') {
		attr(new.lumi, 'vstParameter') <- transPara
		attr(new.lumi, 'transformFun') <- transFun
	}

    return(new.lumi)
}

