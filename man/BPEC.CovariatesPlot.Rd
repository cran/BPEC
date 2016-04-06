\name{BPEC.CovariatesPlot}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{BPEC.CovariatesPlot}
\title{Plot posterior marginal distribution of each covariate for the inferred clusters}
\description{For each covariate, the posterior marginal distribution for each cluster is shown, with colours representing cluster index. The posterior median as well as .10 and .90 quantiles are shown on the same plot in order to provide a measure of uncertainty.  }
\usage{BPEC.CovariatesPlot(CovNames,MCMCout, colorcode=c(7,5,6,3,2,8,4,9,10))}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{CovNames}{The name of each environmental/phenotypic covariate, if a title (header) row was provided in the CoordsLocsFile.}
   \item{MCMCout}{\R object from \code{BPEC.MCMC} run }
  \item{colorcode}{A vector of color codes}
}
\details{Unlike the contour plot, where all posterior samples are shown on the map, in this case quantiles are used pointwise along the axis of each covariate. The plots are shown in the original scale of the covariates and axis labels can be fed into the function to aid interpretation. }
\value{
Produces plots for all additional (environmental or phenotypic) covariate, where posterior medians as well as credible intervals are shown for each cluster. 
}
\author{Ioanna Manolopoulou}
\references{
I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavaré (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.
}

\examples{
#to use example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
RawSeqs = MacrocnemisRawSeqs
CoordsLocs = MacrocnemisCoordsLocs

##to use your own dataset 
#RawSeqs = read.nexus.data('Haplotypes.nex')
#CoordsLocs = read.table(
# 'CoordsLocsFile.txt',header=FALSE,fill=TRUE,col.names=1:max(count.fields('CoordsLocsFile.txt'))
#  )
## to set phenotypic/environmental covariate names, use (as appropriate)
# colnames(CoordsLocs)[1:dims] = c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#for the analysis: 
#check the helpfile of BPEC.MCMC using ?BPEC.MCMC
colorcode = c(7,5,6,3,2,8) #default colour scheme
MCMCout = BPEC.MCMC(RawSeqs,CoordsLocs,MaxMig=2,iter=100,ds=0,PostSamples=5,dims=8)

#if there are also environmental covariates available: 
par(mfrow=c(2,3)) #this splits the plot window into 2x3 to fit all the covariates
BPEC.CovariatesPlot(colnames(CoordsLocs),MCMCout,colorcode) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{Bayesian}
\keyword{covariates}
\keyword{environmental}
\keyword{phenotypic}

