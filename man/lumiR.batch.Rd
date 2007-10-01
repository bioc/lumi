\name{lumiR.batch}
\alias{lumiR.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read BeadStudio output files in batch}
\description{
  
}
\usage{
lumiR.batch(fileList, lib = NULL, sampleInfoFile = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileList}{ a vector of file names or a directory keeping the data files in the format of .csv }
  \item{lib}{ same as \code{\link{lumiR}} parameter lib (optional) }
  \item{sampleInfoFile}{ a Tab-separated text file or a data.frame keeping the sample information (optional) }
  \item{\dots}{ other parameters used by \code{\link{lumiR}} }
}
\details{
The function basically call lumiR for individual files and then combine the returns. The sampleInfoFile parameter is optional. It provides the sample information (for phenoData slot in LumiBatch object), it is a Tab-separated text file. ID column is required. It represents sample ID, which is defined based on the column names of BeadStudio output file. For example, sample ID of column "1881436070\_A\_STA.AVG\_Signal" is "1881436070\_A\_STA". Another Label column is optional, which is sampleNames of LumiBatch object. All information of sampleInfoFile will be directly added in the phenoData slot in LumiBatch object.
}
\value{
A LumiBatch object which combines the individual LumiBatch object corresponding to each file
}
\author{ Pan Du }

\seealso{ \code{\link{lumiR}} }
\examples{
## fileList <- c('file1.csv', 'file2.cvs')
## x.lumi <- lumiR.batch(fileList)

}
\keyword{ methods }