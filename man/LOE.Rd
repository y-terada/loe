\name{LOE}
\alias{LOE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	LOE for a unweighted Graph.
}
\description{
	Performs LOE for a given unweighted grah.
}
\usage{
LOE(ADM, p=2, c=0.1,eps= 1e-5,maxit =1000,method=c("BFGS","SD","MM"),
iniX = "auto",report=100,DEL=1,H=0.5)
}
\arguments{
	\item{ADM}{
	The adjacency matrix.
}
	\item{p}{
	Number of dimensions.
}
	\item{c}{
  	Scale parameter which only takes strictly positive value.
}
	\item{eps}{
	Convergence criterion for the majorization algorithm or the steepest descent algorithm.
}
	\item{maxit}{
	Maximum number of iteretions.
}
	\item{method}{
	If "BFGS", then the BFGS method is used for optimizing the stress function.
	If "SD", then the steepest descent method is used.
	If "MM", then the majorization minimization algortihm is used.
}
	\item{iniX}{
	Matrix with starting values for embedding (optional). If "auto", then Laplacian eigenmaps is used as a starting value.
}
	\item{report}{
	The frequency of reports. Defaults to every 100 iterations.
}
	\item{DEL}{
	The initial step size in the steepest descent algorithm. Defaults to 1.
}
	\item{H}{
The rate parameter of the backtracking line search in the steepest descent algorithm. 
This parameter only takes value in \eqn{(0,1)}{(0,1)}.
Defaults to 0.5.
}
}

\value{
\code{LOE} returns a list with components:
\item{X}{The best corrdinate matrix with p columns whose rows give the coordinates of the vertexes.}
\item{str}{ If \code{method} is "BFGS", then the value of the stress function of LOE corresponding to \code{X} is returned.
 If "SD" or "MM", then the vector of values on each iteration is returned.}
}
%\references{
%}
\author{
Yoshikazu Terada
}
\examples{
#################
#Realizable case#
#################
library(igraph)
###############################
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
N <- nrow(Data)
#Visualization of the original data
plot(Data,pch=20,xlab="",ylab="",cex=1,col=rainbow(N,start=.7,end=.9),
xlim=c(-7,7),ylim=c(-7,7),main="Original data")

#Creating a k-NN graph based on Data
DM <- as.matrix(dist(Data))
ADM <- make.kNNG(DM,k=25)

#plot of the adjacency matrix
AADM <- ADM
diag(AADM) <- NA
image(AADM[N:1,],col=topo.colors(3),ann=FALSE,axes=FALSE)

#Apply some graph embedding methods
LE <-spec.emb(A=ADM,2,norm=FALSE)
result.LOE <- LOE(ADM=ADM,p=2,c=0.1,method="BFGS",report=1000,iniX=LE)

#Procrustes transform
library(vegan)
LOEX <- procrustes(X=Data,Y=result.LOE$X)$Yrot
plot(LOEX,pch=20,xlab="",ylab="",cex=1,col=rainbow(N,start=.7,end=.9),
xlim=c(-7,7),ylim=c(-7,7),main="LOE")
###############################

#############
#This function provide appropriate vectors of xlim and ylim
#for given embedding matrix X.
#############
make.lim <- function(X){
			mima <- matrix(0, 2,2)
			mima[,1] <- apply(X, 2, min)
			mima[,2] <- apply(X, 2, max)
			han <- mima[,2] - mima[,1]
			cent <-  (mima[,2] + mima[,1])/2
			tmpr <- max(han)+max(han)*0.05
			for(s in 1:2){
				mima[s,] <- c(cent[s]-tmpr/2,cent[s]+tmpr/2)
			}
			return(mima)
		}
##############################
#Standered graph-drawing task#
##############################
###############################
ADM <- as.matrix( get.adjacency(graph.famous("Thomassen")) )
TG <- graph.adjacency(ADM)

#Apply some graph embedding methods
LE <-spec.emb(A=ADM,2,norm=FALSE)
KK <- layout.kamada.kawai(TG,maxiter=1000,start=LE)
FR <- layout.fruchterman.reingold(TG,maxiter=1000,start=LE)
result.LOE <- LOE(ADM=ADM,p=2,c=0.1,method="MM",report=1000,maxit=2000)


#Visualization of embeddings
par(mfrow=c(2,3),oma = c(0, 0, 4, 0))
#plot of the adjacency matrix
	AADM <- ADM
	N <- nrow(AADM)
	diag(AADM) <- NA
	image(AADM[N:1,],col=topo.colors(3),ann=FALSE,axes=FALSE)
#plot of Laplacian eigenmaps
	tmplim <- make.lim(LE)
	vsize <- (tmplim[1,2] -tmplim[1,1])*4
	plot(TG, layout=LE, main="Laplacian eigenmaps",
		vertex.size=vsize,vertex.color="blue",vertex.label.color="white",vertex.label=NA,
		edge.arrow.size=0.1,xlim=tmplim[1,],ylim=tmplim[2,],axes=TRUE,rescale=FALSE)
#plot of Kamada and Kawai
	tmplim <- make.lim(KK)
	vsize <- (tmplim[1,2] -tmplim[1,1])*4
	plot(TG, layout=KK, main="Kamada and Kawai",
		vertex.size=vsize,vertex.color="blue",vertex.label.color="white",vertex.label=NA,
		edge.arrow.size=0.1,xlim=tmplim[1,],ylim=tmplim[2,],axes=TRUE,rescale=FALSE)
#plot of Fruchterman Reingold
	tmplim <- make.lim(FR)
	vsize <- (tmplim[1,2] -tmplim[1,1])*4
	plot(TG, layout=FR, main="Fruchterman Reingold",
		vertex.size=vsize,vertex.color="blue",vertex.label.color="white",vertex.label=NA,
		edge.arrow.size=0.1,xlim=tmplim[1,],ylim=tmplim[2,],axes=TRUE,rescale=FALSE)
#plot of LOE
	tmplim <- make.lim(result.LOE$X)
	vsize <- (tmplim[1,2] -tmplim[1,1])*4
	plot(TG, layout=result.LOE$X, main="LOE",
		vertex.size=vsize,vertex.color="blue",vertex.label.color="white",vertex.label=NA,
		edge.arrow.size=0.1,xlim=tmplim[1,],ylim=tmplim[2,],axes=TRUE,rescale=FALSE)
	plot(result.LOE$str,type="l",xlab="Number of iter.", ylab="Stress")
#Make title
	mtext(side = 3, line=1, outer=TRUE, text = "Thomassen", cex=2)
###############################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

