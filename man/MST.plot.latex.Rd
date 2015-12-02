\name{MST.plot.latex}
\alias{MST.plot.latex}
\title{Plot Multivariate Survival Trees in LaTex}
\description{Plot trees in LaTex}
\usage{
MST.plot.latex(tree, file = "tree-code.tex", digits = 5)
}
\arguments{
  \item{tree}{Tree to be plotted}
  \item{file}{Name of the file}
  \item{digits}{Number of significant digits to be displayed}
}
\value{Creates file of LaTex code in working directory}
\author{
Xiaogang Su, Peter Calhoun, and Juanjuan Fan
}
\note{In LaTex, include the command: \verb{\usepackage{pstricks,pst-node,pst-tree}}.  Users familiar with the PSTricks package can make further modifications to the LaTex file.  Use dvips to get proper document (see \url{https://en.wikibooks.org/wiki/LaTeX/PSTricks} to get PDF from PS files)}
\seealso{\code{\link{MST.plot}}}
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

#Generate LaTex files in working directory:
MST.plot.latex(fit$best.tree.structure$Ga.4)
}
\keyword{Plot}
\keyword{LaTex}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
