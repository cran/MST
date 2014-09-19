\name{bootstrap.size}
\alias{bootstrap.size}
\title{Bootstrap Approach to Select Best-Sized Multivariate Survival Tree}
\description{Selects the best-sized multivariate survival tree using a bootstrap approach}
\usage{
bootstrap.size(bootstrap.grow.prune.results, plot.it = TRUE, filename = NULL,
horizontal = TRUE)
}
\arguments{
  \item{bootstrap.grow.prune.results}{Results from \code{bootstrap.grow.prune}}
  \item{plot.it}{Logical: Indicates if you want to plot the goodness-of-fit vs. tree size}
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