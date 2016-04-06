\name{BPEC.TreePlot}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{BPEC.TreePlot}
\title{Treeplot indicating cluster membership}
\description{Given a set of haplotype nodes in a tree structure, creates
  a haplotype tree plot, starting from the root at the top. The size of each haplotype node represents the
number of times it was sampled, and color represents population
cluster. Unsampled haplotypes are represented by black dots. Thickness of tree edges represents posterior certainty, with thinner edges corresponding to higher uncertainty. }
\usage{BPEC.TreePlot(MCMCout,colorcode=c(7,5,6,3,2,8,4,9,10))}

\arguments{
  \item{MCMCout}{\R object from \code{BPEC.MCMC} run. }
   \item{colorcode}{A vector of color codes to use, ideally the same ones used in BPEC.ContourPlot.}
 }

\value{

\item{GraphEdgesSub}{Set of nodes and branch that can be loaded into \code{plot.igraph()}.}
\item{GraphEdgesTree}{Phylogenetic tree representation that can be plotted directly using \code{plot()}.}
}

\details{The algorithm will try to plot the tree such that branches
  don't cross. The label which will be shown correspond to the original labels provided by the user. In cases where 2 sequences (with different labels) were collapsed to the same haplotype, the smallest (in number) of the two labels will be shown. }

 \references{
G. Csardi, T. Nepusz (2006). The igraph software package for complex network research. \emph{InterJournal, Complex Systems}, 1695(5), 1-9.

I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

G. Valiente (2009). \emph{Combinational Pattern Matching Algorithms in Computational Biology Using Perl and R}. CRC Press.
}

\author{Ioanna Manolopoulou & Axel Hille}
\examples{
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
 colorcode = c(7,5,6,3,2,8) #default colour scheme 

MCMCout = BPEC.MCMC(RawSeqs,CoordsLocs,MaxMig=2,iter=50,ds=0,PostSamples=1,dims=8)	
BPEC.Tree = BPEC.TreePlot(MCMCout,colorcode)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree} 
\keyword{phylogenetics}
\keyword{haplotype}
