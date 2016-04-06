\name{BPEC.GeoTree}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{BPEC.GeoTree}
\title{Phylogeographic GoogleEarth plot}
\description{Given the output of an MCMC run, outputs phylogeographic code to show migration on a map through GoogleEarth.  }
\usage{BPEC.GeoTree(MCMCout,CoordsLocs,file="GoogleEarthTree.kml")}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MCMCout}{\R object from \code{BPEC.MCMC} run }
   \item{CoordsLocs}{A matrix where each row shows latitude, longitude (plus any additional covariates), plus all the haplotype numbers found at each location.  }
    \item{file}{Filename for the .kml file to be created }
}
\details{ Additionally to the output variables, a .kml file is created which can be loaded into GoogleEarth.  }

\value{
\item{Phy}{Set of nodes and branches that can be loaded into \code{Phylo2GE()}}
\item{Geo}{Set of haplotypes and locations that can be loaded into \code{Phylo2GE()}}
}

\references{
I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

G. Valiente (2009). \emph{Combinational Pattern Matching Algorithms in Computational Biology Using Perl and R}. CRC Press.

N. Arrigo, L.P. Albert, P.G. Mickelson and M.S. Barker (2012). Quantitative visualization of biological data in Google Earth using R2G2, an R CRAN package. \emph{Molecular Ecology Resources}, 12(6), 1177-1179. 

}
\author{Axel Hille}
\examples{
\dontrun{
#to use example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
CoordsLocs = MacrocnemisCoordsLocs
RawSeqs = MacrocnemisRawSeqs

##to use your own dataset
#RawSeqs = read.nexus.data('Haplotypes.nex')
#CoordsLocs = read.table(
# 'CoordsLocsFile.txt',header=FALSE,fill=TRUE,col.names=1:max(count.fields('CoordsLocsFile.txt'))
#  )
## to set phenotypic/environmental covariate names, use (as appropriate)
# colnames(CoordsLocs)[1:dims] = c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#to run the MCMC sampler: 
MCMCout = BPEC.MCMC(RawSeqs,CoordsLocs,MaxMig=2,iter=50,ds=0,PostSamples=5,dims=8)
BPEC.Geo = BPEC.GeoTree(MCMCout,CoordsLocs,file="GoogleEarthTree.kml")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree} 
\keyword{phylogenetics}
\keyword{haplotype}

