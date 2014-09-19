\name{MST.plot}
\alias{MST.plot}
\title{Plot Multivariate Survival Trees}
\description{Plot trees from the \code{MST} output}
\usage{MST.plot(tree, textDepth = 3, lines = "rectangle")}
\arguments{
  \item{tree}{Tree to be plotted}
  \item{textDepth}{Depth of text}
  \item{lines}{Shape of lines connecting nodes: must be either "rectangle" or "triangle"}
}
\value{Plot of tree structure}
\author{Peter Calhoun, Xiaogang Su, Juanjuan Fan}
\note{This function gives a simple, efficient method of plotting the multivariate survival trees.  True statements in the plot move to the left.  This function does not provide much flexibility when customizing the trees.  User may need to format figures manually}
\seealso{\pkg{tree}, \pkg{rpart}}
\examples{
set.seed(409871)
#Training dataset:
training <- rmultime(beta=c(-1, 2, 2, 0, 0), cutoff=c(.5, .5, 0, 0),
  digits=1, icensor=1, model = "marginal.multivariate.exponential",
  rho=.65, N=200, K=4)$dat

#Test dataset:
test <- rmultime(beta=c(-1, 2, 2, 0, 0), cutoff=c(.5, .5, 0, 0),
	digits=1, icensor=1, model = "marginal.multivariate.exponential", 
	rho=.65, N=100, K=4)$dat

#Construct Multivariate Survival Tree:
fit <- MST(training=training,	test=test, method="marginal",	min.ndsz=30, n0=5,
	col.time=3, col.status=4, col.id=1, col.split.var=5:8, col.ctg=col.ctg,
	max.depth=6, mtry=length(col.split.var), selection.method = "test.sample",
	plot.it=TRUE, filename=NULL, horizontal=TRUE, details=FALSE, cont.split="distinct")


MST.plot(fit$tree0, textDepth=5)          #Plot of initial tree
MST.plot(fit$best.tree.structure[[4]])    #Plot of best-sized tree
}
\keyword{Plot}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
