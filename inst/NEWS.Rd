\name{NEWS}
\title{NEWS file for the MST package}

\section{Changes in Version 2.2 (2020-04-05)}{
  \itemize{
    
    \item{Updated to work for R version 4.0.  Fixed error message when there were few predictors}

  }
}

\section{Changes in Version 2.1 (2018-02-16)}{
  \itemize{
    
    \item{MST package published in Journal of Statistical Software.  Added citation.}

  }
}

\section{Changes in Version 2.0 (2017-04-15)}{
  \itemize{
    
    \item{Changed input names (e.g., changed training to data and max.depth to maxdepth)}

    \item{Changed minevents to represent minimum number of observations needed to split in parent node (instead of in each child node))}
    
    \item{Changed input parameter \code{cont.split} to \code{distinct}}

    \item{Added \code{weights} and \code{subset} optional input parameters to \code{MST()} function}

    \item{Added accessor function \code{getTree()}}
  }
}

\section{Changes in Version 1.3 (2016-05-26)}{
  \itemize{
    
    \item{\code{MST()} now uses a formula}

    \item{Added stratified and independence logrank statistics}
    
    \item{Converted trees to a constparty object}
  }
}

\section{Changes in Version 1.0-1.2 (2014-09-20)}{
  \itemize{
    
    \item{First CRAN release.}
    
    \item{Can construct trees using marginal, gamma frailty, or exponential frailty models}
  }
}