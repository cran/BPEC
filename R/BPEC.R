`bpec` <- function(seqsFile, coordsFile, dims = 2, iter = 100000, postSamples = 100, maxMig = 5, ds = 0, colorCode = c(7,5,6,3,2,8)) {

  rawSeqs = read.nexus.data(seqsFile)
  coordsLocs = bpec.loadCoords(coordsFile)
  
  #run the Markov chain Monte Carlo sampler
  bpecout = bpec.mcmc(rawSeqs, coordsLocs, maxMig, iter, ds, postSamples, dims)
  return(bpecout)
}

`plot.bpec` <-  function(x, GoogleEarth = 0, colorCode = c(7,5,6,3,2,8), ... ) {
  # plot geographical cluster contour map
  bpec.contourPlot(x, GoogleEarth = 0, colorCode = colorCode) 
  
  # plot tree network with cluster indicators
  bpec.Tree = bpec.treePlot(x, colorCode = colorCode)
  
  # now also plot the environmental covariates
  par(mfrow = c(2, 3)) #split the plot window into 2x3 to fit all the covariates
  if (x$input$coordsDimsR > 2) {
    bpec.covariatesPlot(x, colorCode) 
  }
  if (GoogleEarth == 1) {
    bpec.Geo = bpec.geoTree(x,file="GoogleEarthTree.kml")
  }
}

`summary.bpec` <- function(object,...) {
  cat(paste('bpec ran with the following settings:\n'))
  cat(paste('Number of input sequences',object$input$seqCountOrig,'\n'))
  cat(paste('Input sequence length',object$input$seqLengthOrig,'\n'))
  cat(paste('Number of iterations',object$input$iterR,'\n'))
  cat(paste('Number of saved iterations',dim(object$clust$sampleMeansR)[3]/2,'\n'))
  
  cat(paste('Dimensions',object$input$coordsDimsR,'\n'))
  cat(paste('Parsimony relaxation',object$input$dsR,'\n'))
  cat(paste('Maximum number of migrations',length(object$tree$migProbsR) - 1,'\n'))
  
  cat(paste('\nThe results of bpec are: \n'))
  cat(paste('Number of haplotypes (including missing)',dim(object$data$seqR)[1],'\n'))
  cat(paste('of which',dim(object$data$seqR)[1]-object$input$seqCountOrig,'are missing\n'))
  cat(paste('Effective sequence length',object$data$seqLengthR,'\n'))
  cat(paste('The most likely number of migrations is',which.max(object$tree$migProbsR)))
  seqLabels = object$data$seqsFileR[object$data$seqLabelsR]
  maxMig = length(object$tree$migProbsR) - 1
  rootProbMean = (object$tree$rootProbsR[1, ] + object$tree$rootProbsR[2, ]) / 2
 
  if (which.max(rootProbMean) <= length(seqLabels)) {
    root = seqLabels[which.max(rootProbMean)]
    writeLines(paste("\nThe most likely root node is ", root, sep=""))
  }
  if (which.max(rootProbMean) > length(seqLabels)) {
    root = which.max(rootProbMean)
    writeLines(paste("\nThe most likely root node is extinct", sep=""))
  }
  maxVar = 0
  chainMeans = array(0, dim = c(dim(object$clust$sampleMeansR)[1], dim(object$clust$sampleMeansR)[2], 2))
  flag = 0
  postSamples = length(object$clust$sampleMeansR[1, 1, ]) / 2
  
  for(i in 1:(maxMig+1)) {
    for(j in 1:object$input$coordsDimsR) {
      for(l in 1:2) {
        chainMeans[j, i, l] = mean(object$clust$sampleMeansR[j, i, (1 + (l-1) * postSamples):(l * postSamples)], na.rm = TRUE)
      }
      if (sum(is.na(chainMeans[j, i, ])) == 0) {
        if (abs(chainMeans[j, i, 1] - chainMeans[j, i, 2]) > 0.05 * (max(object$input$coordsLocsR[, j]) - min(object$input$coordsLocsR[, j]))) {
          writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations")
          flag = 1
          break
        }
      }
    }
    if (flag == 1) {
      break
    }        
  }
  
  rootProbs1 = object$tree$rootProbsR[1:object$data$countR]
  rootProbs2 = object$tree$rootProbsR[(object$data$countR + 1):(object$data$countR * 2)]
  
  if (max(abs(rootProbs1 / sum(rootProbs1) - rootProbs2 / sum(rootProbs2))) > 0.5 / object$data$countR) {
    writeLines("NO ROOT CONVERGENCE: You need to re-run the sampler with more iterations");
  }   
}


`print.bpec` <- function(x,...) {  
  cat(paste('\nThe results of bpec are: \n'))
  cat(paste('Number of haplotypes (including missing)',dim(x$data$seqR)[1],'\n'))
  cat(paste('of which',dim(x$data$seqR)[1]-x$input$seqCountOrig,'are missing\n'))
  cat(paste('Effective sequence length',x$data$seqLengthR,'\n'))
  cat(paste('The most likely number of migrations is',which.max(x$tree$migProbsR)))
  seqLabels = x$data$seqsFileR[x$data$seqLabelsR]
  maxMig = length(x$tree$migProbsR) - 1
  rootProbMean = (x$tree$rootProbsR[1, ] + x$tree$rootProbsR[2, ]) / 2

  if (which.max(rootProbMean) <= length(seqLabels)) {
    root = seqLabels[which.max(rootProbMean)]
    writeLines(paste("\nThe most likely root node is ", root, sep=""))
  }
  if (which.max(rootProbMean) > length(seqLabels)) {
    root = which.max(rootProbMean)
    writeLines(paste("\nThe most likely root node is extinct", sep=""))
  }
  maxVar = 0
  chainMeans = array(0, dim = c(dim(x$clust$sampleMeansR)[1], dim(x$clust$sampleMeansR)[2], 2))
  flag = 0
  postSamples = length(x$clust$sampleMeansR[1, 1, ]) / 2
  
  for(i in 1:(maxMig+1)) {
    for(j in 1:x$input$coordsDimsR) {
      for(l in 1:2) {
        chainMeans[j, i, l] = mean(x$clust$sampleMeansR[j, i, (1 + (l-1) * postSamples):(l * postSamples)], na.rm = TRUE)
      }
      if (sum(is.na(chainMeans[j, i, ])) == 0) {
        if (abs(chainMeans[j, i, 1] - chainMeans[j, i, 2]) > 0.05 * (max(x$input$coordsLocsR[, j]) - min(x$input$coordsLocsR[, j]))) {
          writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations")
          flag = 1
          break
        }
      }
    }
    if (flag == 1) {
      break
    }        
  }
  
  rootProbs1 = x$tree$rootProbsR[1:x$data$countR]
  rootProbs2 = x$tree$rootProbsR[(x$data$countR + 1):(x$data$countR * 2)]
  
  if (max(abs(rootProbs1 / sum(rootProbs1) - rootProbs2 / sum(rootProbs2))) > 0.5 / x$data$countR) {
    writeLines("NO ROOT CONVERGENCE: You need to re-run the sampler with more iterations");
  }   
}

`mean.bpec` <- function(x, ...) {
    cat('The mean cluster centres are:\n')
    print(apply(x$clust$sampleMeansR, c(1,2), 'mean', na.rm = TRUE))

    cat('The mean cluster covariances are:\n')
    print(apply(x$clust$sampleCovsR, c(1,2,3), 'mean', na.rm = TRUE))

    cat('The mean root posterior probabilities are:\n')
    print(apply(x$tree$rootProbs, 2, 'mean',  na.rm = TRUE))
}

