\name{MST}
\alias{MST}
\title{Multivariate Survival Trees}
\description{
Constructs trees for multivariate survival data using marginal and frailty models.  Combines \code{\link{grow.MST}}, \code{\link{prune.size.testsample}}, \code{\link{bootstrap.grow.prune}}, and \code{\link{bootstrap.size}} into one final wrapper function.
}
\usage{
MST(training, test = NULL, method = c("marginal", "gamma.frailty", "exp.frailty"),
min.ndsz = 20, n0 = 3, col.time, col.status, col.id, col.split.var, col.ctg = NULL,
max.depth = 10, mtry = length(col.split.var), 
selection.method = c("test.sample", "bootstrap"), B = 30, LeBlanc = TRUE,
min.boot.tree.size = 1, plot.it = TRUE, filename = NULL, horizontal = TRUE,
details = FALSE, cont.split = c("distinct", "percentiles"), delta=0.05, nCutPoints=50)
}

\arguments{
  \item{training}{Data or training sample}
  \item{test}{Test sample if available}
  \item{method}{Indicates method of estimating correlation: must be either "marginal", "gamma.frailty", or "exp.frailty"}
  \item{min.ndsz}{Number: Controls the minimum node size}
  \item{n0}{Number: Controls the minimum number of uncensored event times at either child node}
  \item{col.time}{Column number for time}
  \item{col.status}{Column number for status}
  \item{col.id}{Column number for id}
  \item{col.split.var}{Column numbers of all variables wanting to split}
  \item{col.ctg}{Column numbers of categorical variables wanting to split; should be a subset of col.split.var}
  \item{max.depth}{Number: Maximum depth of tree}
  \item{mtry}{Number of variables considered at each split}
  \item{selection.method}{Indicates method of selecting the best-sized subtree: "test.sample" or "bootstrap"}
  \item{B}{Number of bootstrap samples.  Only used if selection.method="bootstrap"}
  \item{LeBlanc}{Logical: Indicates if entire sample used (alternative is out-of-bag sample).  Only used if selection.method="bootstrap"}
  \item{min.boot.tree.size}{Number: Minimum size of tree grown at each bootstrap}
  \item{plot.it}{Logical: Indicates if you want to plot the goodness-of-fit vs. tree size}
  \item{filename}{Name of the file plotted}
  \item{horizontal}{Logical: Indicates if plot should be landscape}
  \item{details}{Logical: Indicates if detailed information should be printed}
  \item{cont.split}{Indicates candidate splits: must be "distinct" or "percentiles"}
  \item{delta}{Consider cutpoints from delta to 1 \eqn{-} delta.  Only used when cont.split="percentiles"}
  \item{nCutPoints}{Number of cutpoints (percentiles) considered.  Only used when cont.split="percentiles"}
}
\details{
Marginal and frailty models are the two main ways to analyze correlated failure times.  Let \eqn{Z_{ij}} represent the covariate vector for the \eqn{j}th member in the \eqn{i}th cluster.

\eqn{\vspace{-1mm}}

The marginal model uses the Cox (1972) proportional hazards model:
\deqn{\lambda_{ij}(t|Z_{ij})=\lambda_{0}(t) \exp(\beta \cdot I(Z_{ij} \leq c))}
where \eqn{\lambda_{0}(t)} is an unspecified baseline hazard function and \eqn{I(\cdot)} is the indicator function.

\eqn{\vspace{-3mm}}

The gamma frailty model uses the proportional hazards model:
\deqn{\lambda_{ij}(t|Z_{ij}, w_{i})=\lambda_{0}(t) \exp(\beta \cdot I(Z_{ij} \leq c)) w_{i}}
where \eqn{\lambda_{0}(t)} is an unspecified baseline hazard function, \eqn{I(\cdot)} is the indicator function, and \eqn{w_{i}} is the frailty term for the \eqn{i}th cluster.

\eqn{\vspace{-3mm}}

The exponential frailty model uses the proportional hazards model:
\deqn{\lambda_{ij}(t|Z_{ij}, w_{i})=\exp(\beta_{0} + \beta_{1} \cdot I(Z_{ij} \leq c)) w_{i}}
where \eqn{I(\cdot)} is the indicator function and \eqn{w_{i}} is the frailty term for the \eqn{i}th cluster.

\eqn{\vspace{-1mm}}

For the marginal model, a robust logrank statistic is calculated for each covariate \eqn{Z} and possible cutpoint \eqn{c}.  The estimate of the score function and likelihood of \eqn{\beta} can be obtained assuming independence.  However, the variance-covariance structure adjusts for the dependence using a sandwich-type estimator.  The best split is the one with the largest robust logrank statistic.

For the frailty models, a score test statistic is calculated from the maximum integrated log likelihood for each covariate \eqn{Z} and possible cutpoint \eqn{c}.  The frailty term must follow some known positive distribution; one common choice is \eqn{w_{i} \sim \Gamma(1/\nu, 1/\nu)} where \eqn{\nu} represents an unknown variance.  Note, the exponential frailty model replaces the baseline hazard function with a constant, yielding different score test statistics and typically computationally faster splits.  The best split is the one with the largest score test statistic.

Note, for continuous variables with many distinct cutpoints, the number of cutpoints considered can be reduced to percentiles.  Using percentiles increases efficiency at the expense of less accuracy.

Growing the initial tree is done by splitting nodes (as described above) reiteratively until the maximum depth of the tree is reached or a small number of observations remain at terminal node.  However, as the final tree model can be any subtree of the initial tree, the number of subtrees can become massive.  A goodness-of-fit with an added penalty for the number of internal nodes is used to prune the trees (i.e. reduce the number of subtrees considered).  The best-sized tree is selected by the largest goodness-of-fit with the added penalty using either the test sample or bootstrap samples.
}
\value{
  \item{tree0}{The initial tree}
  \item{prunining.info}{Trees pruned and considered in the best tree selection}
  \item{best.tree.size}{The best tree size based on the penalty used}
  \item{best.tree.structure}{The best tree structure based on the penalty used}
\eqn{\vspace{-1mm}}

Note, G represents the goodness-of-fit statistic with the respective penalties.  Nodes are represented using binary numbers.  Observations with vname \eqn{\leq} cut moves to the left of the tree; otherwise moves to the right.  See \code{\link{MST.plot}}
}
\references{
Cox D.R. (1972) Regression models and life-tables (with discussion).  Journal of the Royal Statistical Society Series B, 34, 187-202

Fan J., Su X., Levine R., Nunn M., LeBlanc M. (2006) Trees for Correlated Survival Data by Goodness of Split, With Applications to Tooth Prognosis. Journal of American Statistical Association, 101, 959-967

Fan J., Nunn M., Su X. (2009) Multivariate exponential survival trees and their application to tooth prognosis.  Computational Statistics and Data Analysis, 53, 1110-1121

Su X., Fan J. (2004) Multivariate Survival Trees: A Maximum Likelihood Approach Based on Frailty Models.  Biometrics, 60, 93-99
}
\author{
Xiaogang Su, Peter Calhoun, and Juanjuan Fan
}
\note{
Code may take awhile to implement large datasets.  To decrease computation time, user should use test sample (selection.method="test.sample"). User can also split continuous variables based on percentiles (cont.split="percentiles") at the expense of slightly less accuracy.  Method="gamma.frailty" is more computationally intensive
}

\section{Warning}{Error messages in the gamma.frailty model sometimes occur when using the bootstrap method.  Increasing min.ndsz may help fix these errors}

\seealso{
\code{\link{MST.plot}}, \pkg{rpart}
}
\examples{
set.seed(409871)
#Training dataset:
training <- rmultime(beta=c(-1, 2, 2, 0, 0), cutoff=c(.5, .5, 0, 0),
  digits=1, icensor=1, model = "marginal.multivariate.exponential",
	rho=.65, N=200, K=4)$dat

#Test dataset:
test <- rmultime(beta=c(-1, 2, 2, 0, 0), cutoff=c(.5, .5, 0, 0),
	digits=1, icensor=1, model = "marginal.multivariate.exponential", 
	rho=.65, N=100, K=4)$dat

#Construct Multivariate Survival Tree:
fit <- MST(training=training,	test=test, method="marginal",	min.ndsz=30, n0=5,
	col.time=3, col.status=4, col.id=1, col.split.var=5:8, col.ctg=col.ctg,
	max.depth=6, mtry=length(col.split.var), selection.method = "test.sample",
	plot.it=TRUE, filename=NULL, horizontal=TRUE, details=FALSE, cont.split="distinct")
}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}