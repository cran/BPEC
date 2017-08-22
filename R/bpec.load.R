bpec.loadSeq = function(seqsFile) {
  rawSeqs = read.nexus.data(seqsFile)
  return(rawSeqs)
}

bpec.loadCoords = function(coordsFile) {
 coordsLocs = read.table(
  coordsFile, header = FALSE, fill = TRUE, col.names = 1:max(count.fields(coordsFile))
  )
 return(coordsLocs)
}

