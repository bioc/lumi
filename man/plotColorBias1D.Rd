\name{plotColorBias1D}
\alias{plotColorBias1D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the color bias density plot of Illumina Infinium Methylation data 
}
\description{
Plot the color bias density plot of Illumina Infinium Methylation data in one dimension (comparing with 2D scatter plot)
}
\usage{
plotColorBias1D(methyLumiM, removeGenderProbes = FALSE, logMode = TRUE, channel = c("both", "unmethy", "methy", "sum"), xlim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{removeGenderProbes}{ determine whether exclude probes on X and Y chromosomes if the chromosome information is provided in the methyLumiM object.}
  \item{logMode}{Whether plot the intensities in log-scale }
  \item{channel}{ estimate the intensity in different methods }
  \item{xlim}{ the xlim of density plot}
  \item{\dots}{ other parameters used by \code{\link[stats]{density}} }
}
\details{
Plot the color bias density plot of Illumina Infinium Methylation data. There are four options using "channel" parameter to plot the density plot. "both": estimate the density by pooling together methylated and unmethylated probe intensities. "unmethy" and "methy": plot either unmethylated or methylated probe density. "sum" plot the density of the sum of methylatled and unmethylated probe intensitys.
}
\value{
Invisibly return TRUE if plot successfully.
}
\references{

}
\author{
Pan DU
}

\seealso{
See Also as \code{\link{plotColorBias2D}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{hplot}