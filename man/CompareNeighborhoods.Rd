% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manifold_to_lowd.R
\name{CompareNeighborhoods}
\alias{CompareNeighborhoods}
\title{Compare Neighborhoods}
\usage{
CompareNeighborhoods(nn1, nn2)
}
\arguments{
\item{nn1}{A matrix of cells by nearest neighborhood, with members being the
identify of a given nearest neighborhood. For example, an ID of 5 would
correspond to the fifth row of the original data matrix.}

\item{nn2}{A matrix of cells by nearest neighborhood with the same format
as nn1, but the KNN were generated in a different way or on a different
piece or representation of the data.}
}
\value{
A vector with each cell corresponding to its position in the
original data matrix. The values are the number of cells within the 2 KNN
that were in agreement, which in turn can be converted into a percent
fidelity
}
\description{
Calculates the sum of the union of each cell's KNN, one of
which comes from one method (eg original manifold) and the other of which
comes from another method (eg lowd embedding)
}
