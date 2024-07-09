\name{GARI}
\alias{GARI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Graph Adjusted Rand Index
%%  ~~function to do ... ~~
}
\description{
	Computes the graph adjusted Rand index measuring a recovery rate of ordinal information in an unweighted graph.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GARI(ADM, EADM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ADM}{
  	The given unweighted adjacency matrix
%%     ~~Describe \code{ADM} here~~
}
  \item{EADM}{
  	An recovered unweighted adjacency matrix from an embedding. The size of EADM should be same as that of \code{ADM}.
%%     ~~Describe \code{EADM} here~~
}
}
\details{
	GARI is bounded from above by 1, and 
	\eqn{\mathrm{GARI}(A_n,\hat{A}_{n})=1\iff A_n=\hat{A}_{n}}{GARI(A_n,\hat{A}_n) iff A_n=\hat{A}_{n}}. 
A high GARI score implies that many of the ordinal constraints have
been satisfied by the solution.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The graph adjusted rand index measuring a recovery rate of ordinal information (a scalar).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
	Yoshikazu Terada
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(igraph)
###########
#Based on the distance matrix of an embedding, 
#this function provides the recovered adjacency matrix for given number of nearest neighbors.
###########
rec.graph <- function(DM, k, symm =FALSE, weight=FALSE ) {
N <- nrow(DM)
ADM <- matrix(0, N, N)
#Search kNN point
if(weight==TRUE){
 for (i in 1:N) {
  nid <- order(DM[i,])
  ADM[ i, nid[2:(k[i]+1)] ] <- DM[ i, nid[2:(k[i]+1)] ]
 }
}else{
 for (i in 1:N) {
  nid <- order(DM[i,])
  ADM[i,nid[2:(k[i]+1)] ] <- 1
 }
}
if(symm==TRUE){
 SADM <- ADM+t(ADM)
 SADM[SADM==2*ADM] <- ADM[SADM==2*ADM]
 ADM <- SADM
}
return(ADM)
}
						

ADM <- as.matrix( get.adjacency(graph.famous("Thomassen")) )
#Apply LOE
result.LOE <- LOE(ADM=ADM,p=2,c=0.1,method="BFGS",report=1000)

#Compute the vector of numbers of nearest neighbors with each verteces
true.nn <- apply(X=ADM,1,sum)

#Reconstracte the adjacency matrix based on the result embedding
EDM <- as.matrix( dist(result.LOE$X) )
EADM <- rec.graph(EDM,k=true.nn )

#Compute GARI between ADM and EADM
GARI(ADM,EADM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
