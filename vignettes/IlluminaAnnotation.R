### R code from vignette source 'IlluminaAnnotation.Rnw'

###################################################
### code chunk number 1: Loading required libraries
###################################################
library(lumi)
library(annotate)
# library(lumiHumanAll.db)


###################################################
### code chunk number 2: IlluminaAnnotation.Rnw:67-72
###################################################
## provide an arbitrary nucleotide sequence as an example 
seq <- 'ACGTAAATTTCAGTTTAAAACCCCCCG'
## create a nuID for it
id <- seq2id(seq)
print(id)


###################################################
### code chunk number 3: IlluminaAnnotation.Rnw:77-78
###################################################
id2seq(id)


###################################################
### code chunk number 4: IlluminaAnnotation.Rnw:83-84
###################################################
is.nuID(id)


###################################################
### code chunk number 5: IlluminaAnnotation.Rnw:88-89
###################################################
is.nuID('adfqeqe')


###################################################
### code chunk number 6: IlluminaAnnotation.Rnw:118-120
###################################################
  if (require(lumiHumanIDMapping))
 	 lumiHumanIDMapping()


###################################################
### code chunk number 7: IlluminaAnnotation.Rnw:124-129
###################################################
  ## Load Illumina example data in lumi package
  data(example.lumi)
  ## Match the chip information of this example data
  if (require(lumiHumanIDMapping))
	getChipInfo(example.lumi, species='Human')


###################################################
### code chunk number 8: IlluminaAnnotation.Rnw:150-152
###################################################
  if (require(lumiHumanIDMapping))
	lumiHumanIDMapping_nuID()


###################################################
### code chunk number 9: IlluminaAnnotation.Rnw:156-160
###################################################
     nuIDs <- featureNames(example.lumi)
      ## return all mapping information
    if (require(lumiHumanIDMapping))
          nuID2RefSeqID(nuIDs[1:10], lib.mapping='lumiHumanIDMapping')


###################################################
### code chunk number 10: IlluminaAnnotation.Rnw:164-166
###################################################
   if (require(lumiHumanIDMapping))
        nuID2EntrezID(nuIDs[1:10], lib.mapping='lumiHumanIDMapping')


###################################################
### code chunk number 11: IlluminaAnnotation.Rnw:170-174
###################################################
   if (require(lumiHumanIDMapping)) {
        mappingInfo <- nuID2RefSeqID(nuIDs[1:10], lib.mapping='lumiHumanIDMapping', returnAllInfo =TRUE)
        head(mappingInfo)
   }


###################################################
### code chunk number 12: IlluminaAnnotation.Rnw:188-192
###################################################
     data(example.lumi)
     nuIDs <- featureNames(example.lumi)
     if (require(lumiHumanAll.db))
     	getSYMBOL(nuIDs[1:3], 'lumiHumanAll.db')


###################################################
### code chunk number 13: IlluminaAnnotation.Rnw:196-198
###################################################
     if (require(lumiHumanAll.db))
    	 getEG(nuIDs[1:3], 'lumiHumanAll.db')


###################################################
### code chunk number 14: IlluminaAnnotation.Rnw:202-206
###################################################
     if (require(lumiHumanAll.db)) {
 	 goInfo <- getGO(nuIDs[1], 'lumiHumanAll.db')
	 goInfo[[1]][[1]]
     }


###################################################
### code chunk number 15: IlluminaAnnotation.Rnw:210-212
###################################################
     if (require(lumiHumanAll.db))
     	lookUp(nuIDs[1:3], "lumiHumanAll.db", what="SYMBOL")


###################################################
### code chunk number 16: IlluminaAnnotation.Rnw:216-218
###################################################
     if (require(lumiHumanAll.db))
	ls('package:lumiHumanAll.db')


