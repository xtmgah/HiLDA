% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hilda_functions.R
\name{visPMS}
\alias{visPMS}
\title{visualize probabisitic mutaiton signature for the independent model}
\usage{
visPMS(vF, numBases, baseCol = NA, trDir, charSize = 1.2,
  scale = TRUE, alpha = 2, charLimit = 0.25)
}
\arguments{
\item{vF}{a matrix for mutation signature}

\item{numBases}{the number of flanking bases}

\item{baseCol}{the colour of the bases (A, C, G, T, plus strand, minus strand)}

\item{trDir}{the index whether the strand direction is plotted or not}

\item{charSize}{the size of the character}

\item{alpha}{the parameter for the Renyi entropy (applicable only if the isScale is TRUE)}

\item{charLimit}{}

\item{isScale}{the index whether the height of the flanking base is changed or not}
}
\description{
Generate visualization of mutation signatures for the model with
  substitution patterns and flanking bases represented by the indepenent
  representation.
}
