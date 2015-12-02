\name{grow.MST}
\alias{grow.MST}
\title{Grows Multivariate Survival Trees}
\description{Grows the initial multivariate survival tree}
\usage{
grow.MST(dat, test = NULL, method = c("marginal", "gamma.frailty", "exp.frailty"),
col.time, col.status, col.id, col.split.var, col.ctg = NULL, minsplit = 20,
min.nevents = 3, max.depth = 10, mtry = length(col.split.var), 
cont.split = c("distinct", "percentiles"), delta=0.05, nCutPoints=50,
details = FALSE)
}
\arguments{
  \item{dat}{Data or training sample}
  \item{test}{Test sample if available}
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
  \item{cont.split}{Indicates candidate splits: must be \code{"distinct"} or \code{"percentiles"}}
  \item{delta}{Consider cutpoints from delta to 1 \eqn{-} delta.  Only used when \code{cont.split="percentiles"}}
  \item{nCutPoints}{Number of cutpoints (percentiles) considered.  Only used when \code{cont.split="percentiles"}}
  \item{details}{Logical: Indicates if detailed information should be printed}
}
\value{This function outputs the initial tree}
\author{Xiaogang Su, Peter Calhoun, Juanjuan Fan}
\seealso{\code{\link{MST}}}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
