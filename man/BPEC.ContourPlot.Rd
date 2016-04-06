\name{BPEC.ContourPlot}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{BPEC.ContourPlot}
\title{Plot contours of inferred clusters on geographical map}
\description{Plots the contours for a set of population clusters, adding
points for sampling locations, with larger points representing inferred
ancestral locations. Assuming that the data input is in latitude and
longitude form, the contour is superimposed onto the map of the
corresponding region.    }
\usage{BPEC.ContourPlot(MCMCout,CoordsLocs,GoogleEarth=0,colorcode=c(7,5,6,3,2,8,4,9,10))}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MCMCout}{\R object from \code{BPEC.MCMC} run }
  \item{CoordsLocs}{The coordinates of the sampling locations}
\item{GoogleEarth}{If 1, .kml files are produced which can be opened with GoogleEarth.}
 \item{colorcode}{(optional) A vector of color codes}	 
}
\details{Contours with level 0.5 of the posterior mean centres amd covariances of the clusters are shown on the map, with colour representing cluster index. 
Posterior samples are shown with semi-transparent filled contours, so that uncertainty can be represented directly on the map together with the posterior means. The sampling locations are shown as black dots, with the top three (in terms of posterior probability) ancestral locations represented by larger dots. 
If \code{GoogleEarth=1}, this same plot is exported in files that can be directly loaded into Google Earth. 
 }
\value{
Produces a contour plot overlayed on a map. If \code{GoogleEarth=1}, files which may be loaded directly into Google Earth are produced. 
}
\author{Ioanna Manolopoulou, Axel Hille and Steve Brooks}
\references{
I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavaré (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.

S.P. Brooks, I. Manolopoulou, and B.C. Emerson (2007). \emph{Assessing the Effect of Genetic Mutation - A Bayesian Framework for Determining Population History from DNA Sequence Data}. Bayesian Statistics 8. Oxford University Press. 
}

\examples{
#to use example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
RawSeqs = MacrocnemisRawSeqs
CoordsLocs = MacrocnemisCoordsLocs

##to use your own dataset: 
#RawSeqs = read.nexus.data('Haplotypes.nex')
#CoordsLocs = read.table(
#  'CoordsLocsFile.txt',header=FALSE,fill=TRUE,col.names=1:max(count.fields('CoordsLocsFile.txt'))
#    )
## to set phenotypic/environmental covariate names, use (as appropriate)
# colnames(CoordsLocs)[1:dims] = c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#for the analysis: 
#check the helpfile of BPEC.MCMC using ?BPEC.MCMC
colorcode = c(7,5,6,3,2,8) #default colour scheme
MCMCout = BPEC.MCMC(RawSeqs,CoordsLocs,MaxMig=2,iter=20,ds=0,PostSamples=1,dims=8)
BPEC.ContourPlot(MCMCout,CoordsLocs,GoogleEarth=0,colorcode) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{map} 
\keyword{Bayesian}
\keyword{contours}

