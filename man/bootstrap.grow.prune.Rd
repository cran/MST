\name{bootstrap.grow.prune}
\alias{bootstrap.grow.prune}
\title{Bootstrap Approach to Grow and Prune Multivariate Survival Trees}
\description{Grows and prunes multivariate survival trees using a bootstrap approach}
\usage{
bootstrap.grow.prune(data, method=c("marginal", "gamma.frailty", "exp.frailty"),
col.time, col.status, col.id, col.split.var, col.ctg=NULL, minsplit=20,
min.nevents=3, max.depth=10, mtry=length(col.split.var), 
cont.split=c("distinct", "percentiles"), delta=0.05, nCutPoints=50,
B=30, LeBlanc=TRUE, min.boot.tree.size=1, details=FALSE)
}
\arguments{
  \item{data}{Dataset}
  \item{method}{Indicates method of handling correlation: must be either \code{"marginal"}, \code{"gamma.frailty"}, or \code{"exp.frailty"}}
  \item{col.time}{Column number for time}
  \item{col.status}{Column number for status}
  \item{col.id}{Column number for id}
  \item{col.split.var}{Column numbers of all variables wanting to split}
  \item{col.ctg}{Column numbers of categorical variables wanting to split; should be a subset of \code{col.split.var}}
  \item{minsplit}{Number: Controls the minimum node size}
  \item{min.nevents}{Number: Controls the minimum number of uncensored event times at either child node}
  \item{max.depth}{Number: Maximum depth of tree}
  \item{mtry}{Number of variables considered at each split}
  \item{LeBlanc}{Logical: Indicates if entire sample used (alternative is out-of-bag sample)}
  \item{min.boot.tree.size}{Number: Minimum size of tree grown at each bootstrap}
  \item{cont.split}{Indicates candidate splits: must be \code{"distinct"} or \code{"percentiles"}}
  \item{delta}{Consider cutpoints from delta to 1 \eqn{-} delta.  Only used when \code{cont.split="percentiles"}}
  \item{nCutPoints}{Number of cutpoints (percentiles) considered.  Only used when \code{cont.split="percentiles"}}
  \item{B}{Number of bootstrap samples}
  \item{details}{Logical: Indicates if detailed information should be printed}
}
\value{
  \item{initial.tree}{The initial tree}
  \item{boot.prune}{Information on trees pruned for each bootstrap sample}
  \item{boot.tree}{Pruned trees for each bootstrap sample}
}
\author{Xiaogang Su, Peter Calhoun, Juanjuan Fan}
\note{This function is computationally intensive and may take some time}
\seealso{\code{\link{MST}}}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
\keyword{Bootstrap}