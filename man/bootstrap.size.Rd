\name{bootstrap.size}
\alias{bootstrap.size}
\title{Bootstrap Approach to Select Best-Sized Multivariate Survival Tree}
\description{Selects the best-sized multivariate survival tree using a bootstrap approach}
\usage{
bootstrap.size(bootstrap.grow.prune.results, plot.Ga = TRUE, filename = NULL,
horizontal = TRUE)
}
\arguments{
  \item{bootstrap.grow.prune.results}{Results from \code{bootstrap.grow.prune}}
  \item{plot.Ga}{Logical: Indicates if goodness-of-fit vs. tree size should be plotted}
  \item{filename}{Name of the file plotted}
  \item{horizontal}{Logical: Indicates if plot should be landscape}
}
\value{
  \item{initial.tree}{The initial tree}
  \item{btree}{The best-sized tree structure based on the penalty used}
}
\author{Xiaogang Su, Peter Calhoun, Juanjuan Fan}
\seealso{\code{\link{MST}}}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
\keyword{Bootstrap}