\name{make.kNNG}
\alias{make.kNNG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Making the kNN Graph
%%  ~~function to do ... ~~
}
\description{
	Making the kNN graph from given distance matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.kNNG(DM, k = as.integer(2 * log(nrow(DM))), symm = FALSE, weight = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DM}{
  	A distance matrix.
%%     ~~Describe \code{DM} here~~
}
  \item{k}{
  	Number of neighbers
%%     ~~Describe \code{k} here~~
}
  \item{symm}{
  	If TRUE, then the connectivity matrix is symmetrized.
%%     ~~Describe \code{symm} here~~
}
  \item{weight}{
  	If TRUE, then the weighted kNN graph is created.
%%     ~~Describe \code{weight} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	The adjacency matrix of a kNN graph.
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
#Create a toy data
x <- seq(-5,5,by=1)
y <- seq(1,6,by=1)
hx1 <- seq(-3.5,-1.5,by=0.5)
hx2 <- seq(1.5,3.5,by=0.5)
hy <- seq(2.5,4.5,by=0.5)
D1 <- matrix(0,66,2)
for(i in 1:11){
	for(j in 1:6){
		D1[i+11*(j-1),] <- c(x[i],y[j])
	}
}
D2n <- matrix(0,25,2)
D2p <- matrix(0,25,2)
for(i in 1:5){
	for(j in 1:5){
		D2n[i+5*(j-1),] <- c(hx1[i],hy[j])
		D2p[i+5*(j-1),] <- c(hx2[i],hy[j])
	}
}
D2n <- D2n[-c(7,9,17,19),]
D2p <- D2p[-c(7,9,17,19),]
Data <- rbind(D1,D2n,D2p)
Data <- scale(Data[order(Data[,1]),], scale=FALSE)

#Visualization of Data
plot(Data,pch=20,xlab="",ylab="",cex=1,col=rainbow(108,start=.7,end=.9),
xlim=c(-7,7),ylim=c(-7,7))

#Creating a k-NN graph based on Data
DM <- as.matrix(dist(Data))
ADM <- make.kNNG(DM,k=25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
