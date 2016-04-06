\name{BPEC-package}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{BPEC-package}
\docType{package}
\title{
Bayesian Phylogeographic and Ecological Clustering (BPEC)
}
\description{
Bayesian Phylogeographic and Ecological Clustering (BPEC) is aimed at identifying genetically, geographically and ecologically distinct population clusters and drawing inferences about ancestral locations. Given a dataset of DNA sequences (non-recombinant, typically mtDNA) and their respective geographical locations (longitude and latitude) along with additional environmental or phenotypic characteristics, the algorithm simultaneously draws inferences about the genealogy (in the form of a haplotype tree) and the population clustering. In addition, the algorithm identifies locations with high posterior probability of being ancestral. 
}


\details{
BPEC requires 2 input files in order to run:

\code{Haplotypes.nex} :
  The sequence file in NEXUS format. Sequence labels should either be integers, or contain unique integers which correspond to the labels in the \code{CoordsLocsFile.txt}. For example, '1', '1_label', 'label1_label' will all be treated as haplotype 1. NOTE: BPEC will currently ignore unknown nucleotides in the inference.  

\code{CoordsLocsFile.txt} :
   For each observation, the coordinates (latitude and longitude, please use a +/- to  indicate W or E), any other environmental or phenotypic covariates (the latitude and longitude MUST come first), plus the ID number of the haplotype (must match the number in the sequence NEXUS file). If more than one haplotype were found at a single location, these can be entered one after the other, eg:

\preformatted{
  36.88	-5.42	 24	 25
  37.00	-3.98	245	251	243	142	143	244	246	247
  43.35  1.48	153		
}
   so, in the first location (lat/long 36.88, -5.42) you have 2 sampled individuals with   haplotypes 24,25, in the second location  eight individuals etc. Sequences don't necessarily need to be collapsed onto haplotypes, the program should take  care of it.
 
 
The main function is \code{BPEC.MCMC} and runs a Markov chain Monte Carlo sampler in order to obtain posterior samples for all the parameters simultaneously: the haplotype tree, the ancestral nodes, the number of population clusters as well as their means and covariances. 
 
 
Three plotting functions are available. \code{BPEC.ContourPlot} shows the inferred population clusters superimposed on a world map. \code{BPEC.TreePlot} shows the Maximum A Posteriori rooted haplotype tree, indicating posterior cluster membership and number of times each haplotype was sampled.  It can also output .kml files which can be loaded into Google Earth.  \code{BPEC.CovariatesPlot} shows the posterior density estimates for each additional covariate of each cluster.  

}
\author{
Ioanna Manolopoulou <ioanna.manolopoulou@gmail.com>, Axel Hille <axel.hille@gmx.net>
}
\references{
I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavaré (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.

S.P. Brooks, I. Manolopoulou, and B.C. Emerson (2007). \emph{Assessing the Effect of Genetic Mutation - A Bayesian Framework for Determining Population History from DNA Sequence Data}. Bayesian Statistics 8. Oxford University Press. 
}

\keyword{phylogeography} 
\keyword{clustering}
\keyword{Bayesian}
\keyword{phylogenetics}

\examples{
\dontrun{
# if you want to load the example burnet moth dataset
data(TransalpinaRawSeqs)
data(TransalpinaCoordsLocs)
RawSeqs = TransalpinaRawSeqs
CoordsLocs = TransalpinaCoordsLocs

##if you want to use your own dataset, use setwd() to enter the correct folder, 
##then run the command below, changing the input parameters if necessary
 #RawSeqs = read.nexus.data('Haplotypes.nex')
 #CoordsLocs = read.table(
 # 'CoordsLocsFile.txt',header=FALSE,fill=TRUE,col.names=1:max(count.fields('CoordsLocsFile.txt'))
 # )
 ## to set phenotypic/environmental covariate names, use (as appropriate)
# colnames(CoordsLocs)[1:dims] = c('lat','long','cov1','cov2','cov3')  
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 
 
dims = 2 #this is 2 if you only have geographical longitude/latitude. 
#(add 1 for each environmental or phenotypic covariate)
MaxMig = 5 #you will need a higher maximum number of migrations, suggest 7
ds = 0 #start with ds=0 and increase to 1 and then to 2
iter = 10000 #you will need far more iterations for convergence, start with 100,000
PostSamples = 2 #you will need at least 100 saved posterior samples
colorcode = c(7,5,6,3,2,8) #default colour scheme


#run the Markov chain Monte Carlo sampler
MCMCout = BPEC.MCMC(RawSeqs,CoordsLocs,MaxMig,iter,ds,PostSamples,dims)

par(mar=c(0,0,0,0),pty="m",mfrow=c(1,2)) #No plot margins. Contours and tree side-by-side
# plot geographical cluster contour map
BPEC.ContourPlot(MCMCout,CoordsLocs,GoogleEarth=0,colorcode) 

# plot tree network with cluster indicators
BPEC.Tree = BPEC.TreePlot(MCMCout,colorcode)  

## if you want to load the example Brown Frog dataset
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
RawSeqs = MacrocnemisRawSeqs
CoordsLocs = MacrocnemisCoordsLocs

dims = 8 #this is 2 if you only have geographical longitude/latitude. 
#(add 1 for each environmental or phenotypic covariate)
MaxMig = 4 #you will need a higher maximum number of migrations, suggest 7
ds = 2 #start with ds=0 and increase to 1 and then to 2
iter = 10000 #you will need far more iterations for convergence, start with 100,000
PostSamples = 2 #you will need at least 100 saved posterior samples
colorcode = c(7,5,6,3,2,8) #default colour scheme

#run the Markov chain Monte Carlo sampler
MCMCout = BPEC.MCMC(RawSeqs,CoordsLocs,MaxMig,iter,ds,PostSamples,dims)

par(mar=c(0,0,0,0),pty="m",mfrow=c(1,2))  #no plot margins, plot contours and tree side-by-side
# plot geographical cluster contour map
BPEC.ContourPlot(MCMCout,CoordsLocs,GoogleEarth=0,colorcode) 

# plot tree network with cluster indicators
BPEC.Tree = BPEC.TreePlot(MCMCout,colorcode)   

# now also plot the environmental covariates
 par(mfrow=c(2,3)) #split the plot window into 2x3 to fit all the covariates
 BPEC.CovariatesPlot(colnames(CoordsLocs),MCMCout,colorcode) 
 
 BPEC.Geo = BPEC.GeoTree(MCMCout,CoordsLocs,file="GoogleEarthTree.kml")
 }
}
