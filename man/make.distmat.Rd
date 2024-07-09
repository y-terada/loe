\name{make.distmat}
\alias{make.distmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Euclidean distance matrix computation
%%  ~~function to do ... ~~
}
\description{
	This function computes and returns the distance matrix computed by using the Euclidean distance between the row of a data matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.distmat(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  	A numeric matrix.
%%     ~~Describe \code{X} here~~
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	The Euclidean distance matrix based on a given corrdinate matrix \code{X}.
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

#Creating a k-NN graph based on Data
DM <- make.distmat(Data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
