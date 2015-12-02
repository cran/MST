\name{send.down}
\alias{send.down}
\title{Adds Terminal Node to Dataset}
\description{Adds terminal node to dataset}
\usage{
send.down(data, tree)
}
\arguments{
  \item{data}{Data}
  \item{tree}{Tree}
}
\value{Outputs the data provided with the added node column}
\author{Xiaogang Su, Peter Calhoun, Juanjuan Fan}
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

#Add terminal nodes to training dataset
training<-send.down(training,fit$best.tree.structure$Ga.4)
table(training$node)
}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
