### R code from vignette source 'BPEC.Rnw'

###################################################
### code chunk number 1: BPEC.Rnw:382-384
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
# setwd("~/Dropbox/Research/BPECpackage/ManolopoulouHilleEmersonJSS_final/Code")


###################################################
### code chunk number 2: BPEC.Rnw:387-389
###################################################
library("BPEC")
rawSeqs <- bpec.loadSeq("haplotypes.nex")


###################################################
### code chunk number 3: BPEC.Rnw:422-423
###################################################
coordsLocs <- bpec.loadCoords("coordsLocsFile.txt", header = TRUE)


###################################################
### code chunk number 4: BPEC.Rnw:427-431
###################################################
data("MacrocnemisRawSeqs")
data("MacrocnemisCoordsLocs")
rawSeqs <- MacrocnemisRawSeqs
coordsLocs <- MacrocnemisCoordsLocs


###################################################
### code chunk number 5: BPEC.Rnw:439-441 (eval = FALSE)
###################################################
## bpecout <- bpec.mcmc(rawSeqs, coordsLocs, maxMig = 3, iter = 1000000,
##                     ds = 3, postSamples = 1000, dims = 8)


###################################################
### code chunk number 6: BPEC.Rnw:618-620
###################################################
#save(bpecout,file="bpecout.RData")
load("bpecout.RData")


###################################################
### code chunk number 7: BPEC.Rnw:626-629 (eval = FALSE)
###################################################
## par(mar = c(0, 0, 0, 0))
## bpec.contourPlot(bpecout, GoogleEarth = 0, mapType = "osm",
##                  colorCode = c(7, 5, 6, 3, 2), mapCentre = NULL, zoom = 7)


###################################################
### code chunk number 8: BPEC.Rnw:664-666 (eval = FALSE)
###################################################
## par(mfrow = c(2, 3))
## bpec.covariatesPlot(bpecout, colorCode = c(7, 5, 6, 3, 2))


###################################################
### code chunk number 9: BPEC.Rnw:683-684 (eval = FALSE)
###################################################
## bpec.tree <- bpec.treePlot(bpecout, colorCode = c(7, 5, 6, 3, 2))


###################################################
### code chunk number 10: BPEC.Rnw:690-692
###################################################
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), ask = T)
bpec.tree <- bpec.treePlot(bpecout, colorCode = c(7, 5, 6, 3, 2))


###################################################
### code chunk number 11: BPEC.Rnw:704-705 (eval = FALSE)
###################################################
## bpec.geo <- bpec.geoTree(bpecout, file = "GoogleEarthTree.kml")


