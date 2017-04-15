\name{getTree}
\alias{getTree}
\title{Extract initial or best-sized tree}
\description{
This function extracts the tree based on the split penalty.
}
\usage{
getTree(mstObj, Ga = c("0", "2", "3", "4", "log_n"))
}
\arguments{
  \item{mstObj}{The output from the MST fit}
  \item{Ga}{The split penalty}
}
\value{The tree of object class "constparty"}
\author{
Peter Calhoun <calhoun.peter@gmail.com>
}

\seealso{
\code{\link{MST}}
}
