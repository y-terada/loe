\name{get.order}
\alias{get.order}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Making ordinal information.
%%  ~~function to do ... ~~
}
\description{
	For a given dissimilarity matrix \eqn{\Xi=(\xi_{ij})_{n\times n}}{Xi=(xi_{ij})_{n x n}}, 
	creating the ordinal information matrix whose row vectors are ordinal quadruples of indices: 
	a row vector \eqn{(i,j,k,l)}{(i,j,k,l)} means \eqn{\xi_{ij}<\xi_{kl}}{xi_{ij}<xi_{kl}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.order(DM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DM}{
  	A dissimilarity matrix.
%%     ~~Describe \code{DM} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	The matrix whose row vectors are ordinal quadruples of indices: a row vector (i,j,k,l) means \eqn{\xi_{ij}<\xi_{kl}}{xi_{ij}<xi_{kl}}.
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
data(eurodist)
OIF <- get.order(as.matrix(eurodist))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
