\name{grow.MST}
\alias{grow.MST}
\title{Grows Multivariate Survival Trees}
\description{Grows the initial multivariate survival tree}
\usage{
grow.MST(dat, test = NULL, method = c("marginal", "gamma.frailty", "exp.frailty"),
min.ndsz = 20, min.nevents = 3, col.time, col.status, col.id, col.split.var,
col.ctg = NULL, max.depth = 10, mtry = length(col.split.var), details = FALSE,
cont.split = c("distinct", "percentiles"), delta=0.05, nCutPoints=50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data or training sample}
  \item{test}{Test sample if available}
  \item{method}{Indicates method of handling correlation: must be either \code{"marginal"}, \code{"gamma.frailty"}, or \code{"exp.frailty"}}
  \item{min.ndsz}{Number: Controls the minimum node size}
  \item{min.nevents}{Number: Controls the minimum number of uncensored event times at either child node}
  \item{col.time}{Column number for time}
  \item{col.status}{Column number for status}
  \item{col.id}{Column number for id}
  \item{col.split.var}{Column numbers of all variables wanting to split}
  \item{col.ctg}{Column numbers of categorical variables wanting to split; should be a subset of \code{col.split.var}}
  \item{max.depth}{Number: Maximum depth of tree}
  \item{mtry}{Number of variables considered at each split}
  \item{details}{Logical: Indicates if detailed information should be printed}
  \item{cont.split}{Indicates candidate splits: must be \code{"distinct"} or \code{"percentiles"}}
  \item{delta}{Consider cutpoints from delta to 1 \eqn{-} delta.  Only used when \code{cont.split="percentiles"}}
  \item{nCutPoints}{Number of cutpoints (percentiles) considered.  Only used when \code{cont.split="percentiles"}}
}
\value{This function outputs the initial tree}
\author{Xiaogang Su, Peter Calhoun, Juanjuan Fan}
\seealso{\code{\link{MST}}}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
