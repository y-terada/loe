\name{spec.emb}
\alias{spec.emb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Spectral embedding with the normalized or the unnormalized Laplacian (Laplacian eigenmaps).
%%  ~~function to do ... ~~
}
\description{
	Performs spectral embedding for a given adjacency matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spec.emb(A, p, norm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
  	Given adjacency matrix. If \code{A} is not symmetric, it is automatically symmetrized.
%%     ~~Describe \code{A} here~~
}
  \item{p}{
  	The number of dimensions.
%%     ~~Describe \code{p} here~~
}
  \item{norm}{
  	 If \code{TRUE}, then the normalized Laplacian is used for embedding.
 	 If \code{FASLE}, then the unnormalized Laplacian is used for embedding.
%%     ~~Describe \code{norm} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	The coordinate matrix with \code{p} columns whose rows give the coordinates of the vertexes.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
	Yoshikazu Terada
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(igraph)
ADM <- as.matrix( get.adjacency(graph.famous("Icosahedral")) )

#Apply some graph embedding methods
LE <-spec.emb(A=ADM,2,norm=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
