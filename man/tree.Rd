\name{tree}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{tree.bpec}
\alias{tree}
\title{Load location data}
\description{Posterior output for the tree model. }
\usage{
tree(bpecout)
\method{tree}{bpec}(bpecout)
}

\arguments{
\item{bpecout}{\R object from \code{bpec.mcmc} run}		  
 }

\value{
\item{cladoR}{The MAP adjacency matrix for the tree in vectorised format: this means that for two haplotypes i,j, the (i,j)th entry of the matrix is 1 if the haplotypes are connected in the network and 0 otherwise.}
\item{levelsR}{Starting from the root (level 0) all the way to the tips, the discrete depth for the Maximum A Posteriori tree plot.}
\item{edgeTotalProbR}{Posterior probabilities of each edge being present, i.e. corresponding to a mutation which occurred.}
\item{rootProbsR}{The posterior probability per chain that each haplotype was the root of the tree.}
\item{treeEdgesR}{The set of edges (from and to haplotypes) of the Maximum A Posteriori haplotype tree (could be used in another program if needed).}
\item{rootLocProbsR}{Vector of posterior probabilities of each sampling location being the ancestral location.}
\item{migProbsR}{The posterior probability of 0...maxMig migrations.}
}

\author{Ioanna Manolopoulou & Axel Hille}
\examples{
## if you want to load the `mini' example Brown Frog dataset
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocsMini)
rawSeqs <- MacrocnemisRawSeqs
coordsLocs <- MacrocnemisCoordsLocsMini

dims <- 3 #this is 2 if you only have geographical longitude/latitude. 
#(add 1 for each environmental or phenotypic covariate)
maxMig <- 2 #you will need a higher maximum number of migrations, suggest 7
ds <- 0 #start with ds=0 and increase to 1 and then to 2
iter <- 1000 #you will need far more iterations for convergence, start with 100,000
postSamples <- 100 #you will need at least 100 saved posterior samples

#run the Markov chain Monte Carlo sampler
bpecout <- bpec.mcmc(rawSeqs,coordsLocs,maxMig,iter,ds,postSamples,dims)
tree(bpecout)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
