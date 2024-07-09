\name{loe-package}
\alias{loe-package}
\alias{loe}
\docType{package}
\title{
Local Ordinal Embedding
%%~~ package title ~~
}
\description{
Soft Ordinal Embedding (SOE) is a simple new method for ordinal embedding problem.
SOE can deal with a part of ordinal comparisons:
\eqn{dissimilarity(i,j)<dissimilarity(k,l)}{dissimilarity(i,j)<dissimilarity(k,l)} for certain quadruples of indices \eqn{(i,j,k,l)}{(i,j,k,l)}.

The local point of view suggests
ordinal embedding as an interesting alternative to graph drawing
algorithms. If vertex \eqn{i}{i} is connected by an edge to vertex \eqn{j}{j}, but
not to vertex \eqn{k}{k}, we interpret this constellation as a constraint of
the form \eqn{\xi_{ij} < \xi_{ik}}{dissimilarity(i,j)<dissimilarity(i,k)}. 
With this interpretation, 
graph embedding (graph drawing) becomes a special case of ordinal embedding. 
We call SOE for graph embedding problem Local Ordinal Embedding (LOE). 
For the connectivity matrix of given unweighted graph, 
LOE provides a low-dimensional embedding which preserves the original graph structure.

This package provides SOE and LOE. 
%%~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab loe\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-01-28\cr
License: \tab GPL-2\cr
}
The functions for SOE on ordinal comparison matrices and LOE on unweighted adjacency matrices are \code{SOE()} and \code{LOE()}, respectively. 

%%~~ An overview of how to use the package, including the most ~~
%%~~ important functions ~~
}
\author{
Yoshikazu Terada and Ulrike von Luxburg

Maintainer: Yoshikazu Terada \email{terada@sigmath.es.osaka-u.ac.jp}
%~~ The author and/or maintainer of the package ~~
}
%\references{
%~~ Literature or other references for background information ~~
%}
%%~~ Optionally other standard keywords, one per line, from file ~~
%%~~ KEYWORDS in the R documentation directory ~~
\section{Acknowledgement}{
	The authors wishes to express their thanks to Dr. Michio Yamamoto and Mr. Yuichi Saito for their help in creating this package.
	This work was partially supported by Grant-in-Aid for JSPS Fellows Number 24-2466.
	}
%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%}
