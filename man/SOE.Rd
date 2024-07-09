\name{SOE}
\alias{SOE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Soft Ordinal Embedding (SOE)
%%  ~~function to do ... ~~
}
\description{
	Performs SOE for a given part of ordinal informaiton.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SOE(CM, N, p = 2, c = 0.1, maxit = 1000, report = 100, iniX = "rand",rnd=10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CM}{
  	A matrix whose rows represent ordinal information: 
  	A row vector \eqn{(i,j,k,l)}{(i,j,k,l)} means \eqn{\xi_{ij}<\xi_{kl}}{xi_{ij} < xi_{kl}}, where \eqn{\xi_{ij}}{xi_{ij}} is a dissimilarity 
  	between object \eqn{i}{i} and object \eqn{j}{j}. 
  	If the number of given ordinal comparisons is more than 100000, then SOE randomly chooses \code{rnd} ordinal comarisons from given comparisons.
%%     ~~Describe \code{CM} here~~
}
  \item{N}{
  	The number of objects.
%%     ~~Describe \code{N} here~~
}
  \item{p}{
  	The number of dimensions.
%%     ~~Describe \code{p} here~~
}
  \item{c}{
  	Scale parameter which only takes strictly positive value.
%%     ~~Describe \code{c} here~~
}
  \item{maxit}{
  	Maximum number of iteretions.
%%     ~~Describe \code{maxit} here~~
}
  \item{report}{
  	The frequency of reports. Defaults to every 100 iterations.
%%     ~~Describe \code{report} here~~
}
  \item{iniX}{
  	Matrix with starting values for embedding (optional).
%%     ~~Describe \code{iniX} here~~
}
  \item{rnd}{
  	The number of ordinal comarisons which are randomly chosen from given ordinal comparisons \code{CM}.
  	\code{rnd} is used only if the number of given ordinal comparisons is more than 100000.
%%     ~~Describe \code{iniX} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{SOE} returns a list with components:
\item{X}{The best corrdinate matrix with \eqn{p}{p} columns whose rows give the coordinates of objects.}
\item{str}{The value of the stress function of SOE corresponding to \code{X}.}
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
library(MASS)
data(eurodist)
OIF <- get.order(as.matrix(eurodist))
sid <- sample(1:nrow(OIF))
POI <- OIF[sid[1:1000],]
result <- SOE(CM=POI, N=nrow(as.matrix(eurodist)),p=2, c=0.1,maxit =1000,report=100)
plot(result$X,type="n",xlab="Dim. 1", ylab="Dim. 2")
text(result$X,rownames(as.matrix(eurodist)),cex=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
