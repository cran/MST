\name{MST}
\alias{MST}
\alias{splitting.stat.MST1}
\alias{splitting.stat.MST2}
\alias{splitting.stat.MST3}
\alias{splitting.stat.MST4}
\title{Multivariate Survival Trees}
\description{
Constructs trees for multivariate survival data using marginal and frailty models.  A wrapper function that grows a large initial tree, prunes the tree, and selects the best sized tree.
}
\usage{
MST(formula, data, test = NULL, weights_data, weights_test, subset,
  method = c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"),
  minsplit = 20, minevents = 3, minbucket = round(minsplit/3), maxdepth = 10,
  mtry = NULL, distinct = TRUE, delta = 0.05, nCutPoints = 50,
  selection.method = c("test.sample", "bootstrap"),
  B = 30, LeBlanc = TRUE, min.boot.tree.size = 1,
  plot.Ga = TRUE, filename = NULL, horizontal = TRUE, details = FALSE, sortTrees = TRUE)
}
\arguments{
  \item{formula}{A linear survival model with the response on the left of a ~ operator and the predictors, separated by + operators, on the right.  Cluster (or id) variable is distinguished by a vertical bar \code{|} (e.g. \code{Surv(time,status) ~ x1 + x2 | id}).  Categorical predictors must be treated as a factor.}
  \item{data}{Data to grow and prune the tree}
  \item{test}{Test sample if available}
  \item{weights_data}{An optional vector of weights to grow the tree}
  \item{weights_test}{An optional vector of weights to select the best-sized tree}
  \item{subset}{An optional vector specifying a subset of observations to be used to grow the tree}
  \item{method}{Indicates method of handling correlation: must be either \code{"marginal"},
\code{"gamma.frailty"}, \code{"exp.frailty"}, \code{"stratified"}, or \code{"independence"}}
  \item{minsplit}{Number: Controls the minimum node size}
  \item{minevents}{Number: Controls the minimum number of uncensored event times}
  \item{minbucket}{Number: Controls the minimum number of observations in any terminal node}
  \item{maxdepth}{Number: Maximum depth of tree}
  \item{mtry}{Number of variables considered at each split.  The default is to consider all variables}
  \item{distinct}{Logical: Indicates if all distinct cutpoints or only percentiles considered}
  \item{delta}{Consider cutpoints from delta to 1 \eqn{-} delta.  Only used when \code{distinct = TRUE}}
  \item{nCutPoints}{Number of cutpoints (percentiles) considered.  Only used when \code{distinct = TRUE}}
  \item{selection.method}{Indicates method of selecting the best-sized subtree: \code{"test.sample"} or \code{"bootstrap"}}
  \item{B}{Number of bootstrap samples.  Only used if \code{selection.method = "bootstrap"}}
  \item{LeBlanc}{Logical: Indicates if entire sample used (alternative is out-of-bag sample).  Only used if \code{selection.method = "bootstrap"}}
  \item{min.boot.tree.size}{Number: Minimum size of tree grown at each bootstrap}
  \item{plot.Ga}{Logical: Indicates if goodness-of-fit vs. tree size should be plotted}
  \item{filename}{Name of the file plotted}
  \item{horizontal}{Logical: Indicates if plot should be landscape}
  \item{details}{Logical: Indicates if detailed information on the construction should be printed}
  \item{sortTrees}{Logical: Indicates if trees should be sorted such that each split to the left has lower risk of failure}
}
\details{
Marginal and frailty models are the two main ways to analyze correlated failure times.  Let \eqn{X_{ij}} represent the covariate vector for the \eqn{j}th member in the \eqn{i}th cluster.

The marginal model uses the Cox (1972) proportional hazards model:
\deqn{\lambda_{ij}(t|X_{ij})=\lambda_{0}(t) \exp(\beta \cdot I(X_{ij} \leq c))}
where \eqn{\lambda_{0}(t)} is an unspecified baseline hazard function and \eqn{I(\cdot)} is the indicator function.

The gamma frailty model uses the proportional hazards model:
\deqn{\lambda_{ij}(t|X_{ij}, w_{i})=\lambda_{0}(t) \exp(\beta \cdot I(X_{ij} \leq c)) w_{i}}
where \eqn{\lambda_{0}(t)} is an unspecified baseline hazard function, \eqn{I(\cdot)} is the indicator function, and \eqn{w_{i}} is the frailty term for the \eqn{i}th cluster.

The exponential frailty model uses the proportional hazards model:
\deqn{\lambda_{ij}(t|X_{ij}, w_{i})=\exp(\beta_{0} + \beta_{1} \cdot I(X_{ij} \leq c)) w_{i}}
where \eqn{I(\cdot)} is the indicator function and \eqn{w_{i}} is the frailty term for the \eqn{i}th cluster.

For the marginal model, a robust logrank statistic is calculated for each covariate \eqn{X} and possible cutpoint \eqn{c}.  The estimate of the score function and likelihood of \eqn{\beta} can be obtained assuming independence.  However, the variance-covariance structure adjusts for the dependence using a sandwich-type estimator.  The best split is the one with the largest robust logrank statistic.

For the frailty models, a score test statistic is calculated from the maximum integrated log likelihood for each covariate \eqn{X} and possible cutpoint \eqn{c}.  The frailty term must follow some known positive distribution; one common choice is \eqn{w_{i} \sim \Gamma(1/\nu, 1/\nu)} where \eqn{\nu} represents an unknown variance.  Note, the exponential frailty model replaces the baseline hazard function with a constant, yielding different score test statistics and typically computationally faster splits.  The best split is the one with the largest score test statistic.

Stratified model grows a tree by minimizing the within-strata variation.  This method should be used with care because the tree will not split on variables with a fixed value within each stratum.  The independence model ignores the dependence and uses the logrank statistic as the splitting rule.

For continuous variables with many distinct cutpoints, the number of cutpoints considered can be reduced to percentiles.  Using percentiles increases efficiency at the expense of less accuracy.

Growing the initial tree is done by splitting nodes (as described above) reiteratively until the maximum depth of the tree is reached or a small number of observations remain at terminal node.  However, as the final tree model can be any subtree of the initial tree, the number of subtrees can become massive.  A goodness-of-fit with an added penalty for the number of internal nodes is used to prune the trees (i.e. reduce the number of subtrees considered).  The best-sized tree is selected by the largest goodness-of-fit with the added penalty using either the test sample or bootstrap samples.
}

\value{
  \item{tree0}{The initial tree.  Tree listed as constparty object}
  \item{prunining.info}{Trees pruned and considered in the best tree selection}
  \item{best.tree.size}{The best tree size based on the penalty used}
  \item{best.tree.structure}{The best tree structure based on the penalty used.  Tree listed as constparty object}

Note, the constparty object requires a constant fit from each terminal node.  Thus, the \code{predict} and \code{plot} functions ignore the dependence, so users are recommended to fit their own model when making predictions (see example)
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
Code may take awhile to implement large datasets.  To decrease computation time, user should use test sample (\code{selection.method = "test.sample"}). User can also split continuous variables based on percentiles (\code{distinct = FALSE}) at the expense of slightly less accuracy.  Gamma frailty models are more computationally intensive
}
\section{Warning}{Error messages in the gamma frailty models sometimes occur when using the bootstrap method.  Increasing \code{minsplit} may help fix these errors.  The exponential frailty model can have problems for large, extremely unbalanced designs.  Currently weights can only be applied to marginal and gamma frailty models.}

\seealso{
\pkg{rpart}
}

\examples{
set.seed(186117)
data <- rmultime(N = 200, K = 4, beta = c(-1, 0.8, 0.8, 0, 0), cutoff = c(0.5, 0.3, 0, 0),
    model = "marginal.multivariate.exponential", rho = 0.65)$dat
test <- rmultime(N = 100, K = 4, beta = c(-1, 0.8, 0.8, 0, 0), cutoff = c(0.5, 0.3, 0, 0),
    model = "marginal.multivariate.exponential", rho = 0.65)$dat

#Construct Multivariate Survival Tree:
fit <- MST(formula = Surv(time, status) ~ x1 + x2 + x3 + x4 | id, data, test,
    method = "marginal", minsplit = 100, minevents = 20, selection.method = "test.sample")

(tree_final <- getTree(fit, 4))
plot(tree_final)

#Fit a model from the final tree
data$term_nodes <- as.factor(predict(tree_final, newdata = data, type = 'node'))
coxph(Surv(time, status) ~ term_nodes + cluster(id), data = data)
}
\keyword{Trees}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
