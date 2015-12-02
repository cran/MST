\name{MST.plot}
\alias{MST.plot}
\title{Plot Multivariate Survival Trees}
\description{Plot trees from the \code{MST} output}
\usage{MST.plot(tree, textDepth = 4, digits = 3, nsmall = 0L, varText = c("vname", "var"),
lines = c("rectangle", "triangle"), ...)}
\arguments{
  \item{tree}{Tree to be plotted}
  \item{textDepth}{Depth of text}
  \item{digits}{Number of significant digits to be displayed}
  \item{nsmall}{The minimum number of digits to the right of the decimal point}
  \item{varText}{Indicates whether variable name or number should be displayed: must be either \code{"vname"} or \code{"var"}}
  \item{lines}{Shape of lines connecting nodes: must be either \code{"rectangle"} or \code{"triangle"}}
  \item{\dots}{Arguments to be passed to methods, such as graphical parameters}
}
\value{Plot of tree structure}
\author{Peter Calhoun, Xiaogang Su, Juanjuan Fan}
\note{This function gives a simple, efficient method of plotting the multivariate survival trees.  True statements in the plot move to the left.  This function does not provide much flexibility when customizing the trees.  User may need to format figures manually}
\seealso{\pkg{tree}, \pkg{rpart}}
\examples{
set.seed(409871)
#Training dataset:
training <- rmultime(N=200, K=4, beta=c(-1, 2, 2, 0, 0), cutoff=c(.5, .5, 0, 0),
  digits=1, icensor=1, model = "marginal.multivariate.exponential", rho=.65)$dat

#Test dataset:
test <- rmultime(N=100, K=4, beta=c(-1, 2, 2, 0, 0), cutoff=c(.5, .5, 0, 0),
  digits=1, icensor=1, model = "marginal.multivariate.exponential", rho=.65)$dat

#Construct Multivariate Survival Tree:
fit <- MST(training=training,	test=test, method="marginal",
  col.time=3, col.status=4, col.id=1, col.split.var=5:8,
  minsplit=30, min.nevents=5, selection.method = "test.sample")

MST.plot(fit$tree0, textDepth=5)          #Plot of initial tree
MST.plot(fit$best.tree.structure$Ga.4)    #Plot of best-sized tree with alpha=4
}
\keyword{Plot}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}