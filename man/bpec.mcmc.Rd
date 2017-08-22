\name{bpec.mcmc}	
\encoding{ISO-8859-1}
\Rdversion{1.1}		
\alias{bpec.mcmc}	
\title{Markov chain Monte Carlo sampler for BPEC}	
\description{Markov chain Monte Carlo for Bayesian Phylogeographic and Ecological Clustering, implemented in C. Given a dataset of DNA sequences (non-recombinant, typically mtDNA) and their respective geographical locations (longitude and latitude), the algorithm simultaneously draws inferences about the genealogy (in the form of a haplotype tree) and the population clustering (with an unknown number of clusters). In addition, the algorithm identifies locations with high posterior probability of being ancestral. In case where additional covariates are available (e.g. climate data), these may be added to the 2-dimensional data and inserted into the analysis. }
\usage{bpec.mcmc(rawSeqs,coordsLocs,maxMig,iter,ds,postSamples,dims=-1)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{rawSeqs}{The input DNA sequences.}	
  \item{coordsLocs}{A matrix where each row shows latitude, longitude (plus any additional covariates), plus all the haplotype numbers found at each location.  }

\item{maxMig}{The maximum number of migration events (this means that the maximum number of clusters will be \code{maxmig+1}). The number you enter here is just an upper bound, so start with \code{maxmig=6} and only increase it if you are really getting 7 clusters in return (which could mean that more clusters are appropriate). If, say, 4 clusters are needed, whether you use \code{maxMig=6} or \code{maxMig=10} (or similar), the number of clusters will collapse down to 4.  }
  \item{iter}{The number of iterations for the MCMC sampler, must be a multiple of 10. You will need quite a large number here, like 100,000. Two MCMC chains will run, after which convergence is checked. If convergence has not been reached, the output will say "NO CONVERGENCE" and you should increase the number of iterations. }
\item{ds}{This represents the parsimony relaxation parameter, with 0 being the minimum. Generally the higher \code{ds}, the more candidate trees are considered, but this comes at a computational cost. Start with \code{ds=0} and increase to \code{ds=1}, etc, observing any changes. }
\item{postSamples}{How many posterior samples per chain to save for use post-processing. A value of at least \code{PostSamples=1000} would provide a reasonable assessment of posterior uncertainty. \code{PostSamples} must not be greater than \code{iter/10}. Also, only up to 20 saved (thinned) samples are used in the \code{bpec.contourPlot} function. }
\item{dims}{The dimension, 2 for purely geographical data, +1 for each covariate (for example if environmental or phenotypic characteristics are also available).}
}
\details{This is the main command for BPEC, the details of which can be found in the provided references. In short, the model is such that, for a given haplotype tree, clusterings are the result of the migration of an individual, such that all the descendants of that individual belong to the new cluster (unless they subsequently migrate themselves). The number of migrating individuals is itself a parameter. Conditional on the clustering allocations for each sequence (note: not all samples of a haplotype need to belong to the same cluster), the geographical (and covariate) distribution for each cluster is assumed to be Gaussian. The distributions of longitude and latitude in each cluster are assumed to have an unknown covariance, whereas additional covariates are assumed to be independent. 

Posterior samples are obtained through Markov chain Monte Carlo (MCMC) sampling. Two Markov chain Monte Carlo runs are carried out. In each, the sampler sweeps through the updates of all the parameters: the haplotype tree, the root of the tree, the number of clusters, the precise clustering (as a result of migrations) the locations and variances of the clusters. Metropolis-Hastings updates are performed for most of the parameters. In the case of the number of clusters, Metropolis-Hastings Reversible Jump updates are performed. 

The space of possible clusterings is vast and highly multi-modal. The migration of an individual to a new cluster implies that all of their descendants will belong to the new cluster. This results in a combinatorial parameter space which is challenging to explore. A number of sophisticated tricks are used in order to overcome this challenge, alternating between local and global MCMC proposal moves. The biggest challenge is to converge to the region of high posterior probability in terms of the number of clusters and the cluster allocation. As such, the first 90 percent of the iterations are discarded as burn-in and only the final 10 percent are used as potential posterior samples. 

\code{bpec.mcmc} requires 2 input files in order to run:

\code{haplotypes.nex} :
   The sequence file in NEXUS format. Sequence labels should either be integers, or contain unique integers which correspond to the labels in the \code{CoordsLocsFile.txt}. For example, '1', '1_label', 'label1_label' will all be treated as haplotype 1. NOTE: BPEC will currently ignore unknown nucleotides in the inference.  
  
\code{coordsLocsFile.txt} :
   For each observation, the coordinates (latitude and longitude, please use a +/- to  indicate W or E), any other environmental or phenotypic covariates (the latitude and longitude MUST come first), plus the ID number of the haplotype (must match the number in the sequence NEXUS file). If more than one haplotype were found at a single location, these can be entered one after the other, eg:

\preformatted{
  36.88	-5.42	 24	 25
  37.00	-3.98	245	251	243	142	143	244	246	247
  43.35  1.48	153		
}
   so, in the first location (lat/long 36.88, -5.42) you have 2 sampled individuals with   haplotypes 24,25, in the second location  eight individuals etc. Sequences don't necessarily need to be collapsed onto haplotypes, the program should take  care of it.

}
\value{\code{bpec} object with the following variables. 
\itemize{
\item{\code{input}: The values and data used as inputs.
\itemize{
\item \code{seqCountOrig}: The number of input sequences.
\item  \code{seqLengthOrig}: The length of the input sequences. 
\item  \code{iterR}: The number of MCMC iterations. 
\item  \code{dsR}: The parsimony relaxation parameter. 
\item  \code{coordsLocsR}: The input coordinate and observation file. 
\item  \code{coordsDimsR}: The input dimension (2 for purely geographical data) .
\item  \code{locNoR}: The number of distinct sampling locations. 
\item  \code{locDataR}: The list of coordinates of each observation. 
}
}

\item{\code{data}: Various post-processing quantities (before the MCMC).
\itemize{
\item \code{seqR}: The output DNA sequences of distinct haplotypes, collapsed to effective nucleotide sites (both sampled and missing sequences which were inferred).
\item \code{seqsFileR}: A vector of the numerical labels of each haplotype.
\item \code{seqLabelsR}: Correspondence vector for each of the processed observations to the original haplotype labels. 
\item \code{seqIndicesR}: Correspondence vector for each of the original observations to the resulting haplotype labels.
\item \code{seqLengthR}: The effective length of the input sequences, given by the number of variable nucleotide sites which are informative. In other words, if two different nucleotide sites are variable in exactly the same haplotypes, then they effectively provide information of a single site.
\item \code{noSamplesR}: The number of times each haplotype was observed in the sample. 
\item \code{countR}: The number of output sequences. 
}
}
\item{\code{clust}: Various posterior sampling summaries about the clustering.
\itemize{
\item \code{sampleMeansR}  A set of posterior samples of the cluster means (i.e. centres).
\item \code{sampleCovsR}: A set of posterior samples of the cluster covariances (i.e. shapes).
\item \code{sampleIndicesR}: A set of posterior samples of the cluster allocations of each observation. 
\item \code{clusterProbsR}: For each haplotype, posterior probabilities that it belongs to each cluster.
}
}

\item{\code{tree}: Various posterior sampling summaries about the tree.
\itemize{
\item \code{cladoR}: The MAP adjacency matrix for the tree in vectorised format: this means that for two haplotypes i,j, the (i,j)th entry of the matrix is 1 if the haplotypes are connected in the network and 0 otherwise.
\item \code{levelsR}: Starting from the root (level 0) all the way to the tips, the discrete depth for the Maximum A Posteriori tree plot.       
\item \code{edgeTotalProbR}: Posterior probabilities of each edge being present, i.e. corresponding to a mutation which occurred.
\item \code{rootProbsR}: The posterior probability per chain that each haplotype was the root of the tree.
\item \code{treeEdgesR}: The set of edges (from and to haplotypes) of the Maximum A Posteriori haplotype tree (could be used in another program if needed).
\item \code{rootLocProbsR}: Vector of posterior probabilities of each sampling location being the ancestral location.
\item \code{migProbsR}: The posterior probability of 0...maxMig migrations.
}
}

\item{\code{mcmc}: Various MCMC tuning parameters and error codes.
\itemize{
\item \code{MCMCparamsR}: Various MCMC tuning parameters, useful for development. 
\item \code{codaInput}: Posterior samples from the two MCMC chains for the cluster means, cluster covariance entries, as well as the root haplotype. Note that, since the number of clusters varies from iteration to iteration, some samples are simply draws from the prior (corresponding to empty clusters). This variable can be loaded directly into the \code{coda} package for convergence analysis.   
\item \code{errorCodeR}: error code, 0 if no problems with running the MCMC, 1 otherwise.     
}
}
}
}
\note{\code{bpec.mcmc} uses  \code{cexcept.h 2.0.1} (an interface for exception-handling in ANSI C) developed by Adam M. Costello and Cosmin Truta.
}

\references{
I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavaré (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.

 Adam M. Costello and Cosmin Truta (2008) \code{cexcept.h} exception handling interface in C http://www.nicemice.net/cexcept/.

}
\author{Ioanna Manolopoulou}
\examples{
#to use the example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
coordsLocs = MacrocnemisCoordsLocs
rawSeqs = MacrocnemisRawSeqs

##to use your own dataset 
#rawSeqs = read.nexus.data('Haplotypes.nex')
#coordsLocs = bpec.loadCoords("coordsLocsFile.txt")

## to set phenotypic/environmental covariate names, use (as appropriate)
# colnames(coordsLocs)[1:dims] = c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#to run the MCMC sampler: 
bpecout = bpec.mcmc(rawSeqs,coordsLocs,maxMig=2,iter=50,ds=0,postSamples=5,dims=8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogeography} 
\keyword{clustering}
\keyword{Bayesian}
\keyword{MCMC}
\keyword{convergence}



