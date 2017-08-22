\name{bpec.loadCoords}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec.loadCoords}
\title{Load location data}
\description{Create \R objects of  locations from  corresponding file. }
\usage{bpec.loadCoords(coordsFile)}

\arguments{
\item{coordsFile}{the name of the coordinate and sequence file in full, eg "coordsLocs.txt".}
 }

\value{
\item{coordsLocs}{\R object of locations and corresponding haplotypes.}
}

\details{Given the names of the file, this command will create an \R objects of the coordinates and locations that can be fed into \code{bpec.mcmc}. }


\author{Ioanna Manolopoulou & Axel Hille}
\examples{
\dontrun{
coordsLocs = bpec.loadCoords('coordsLocs.txt')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
