\name{rmultime}
\alias{rmultime}
\title{Random Multivariate Survival Data}
\description{Generates multivariate survival data}
\usage{
rmultime(N = 100, K = 4, beta = c(-1, 2, 1, 0, 0), cutoff = c(0.5, 0.5, 0, 0),
  digits = 1, icensor = 1, model = c("gamma.frailty", "log.normal.frailty",
  "marginal.multivariate.exponential", "marginal.nonabsolutely.continuous",
  "nonPH.weibull"), v = 1, rho = 0.65, a = 1.5, lambda = 0.1)
}
\arguments{
  \item{N}{Number of clusters (ids)}
  \item{K}{Number of units per cluster}
  \item{beta}{Vector of beta coefficients (first number is baseline hazard coefficient (\eqn{\beta_{0}}), remaining numbers are slope coefficients for covariates (\eqn{\beta_{1}}))}
  \item{cutoff}{Cutoff values for each covariate}
  \item{digits}{Rounding digits}
  \item{icensor}{Control for censoring rate: 1 - 50\%}
  \item{model}{Model for simulating data: must be either \code{"gamma.frailty"}, \code{"log.normal.frailty"}, \code{"marginal.multivariate.exponential"}, \code{"marginal.nonabsolutely.continuous"},
  or \code{"nonPH.weibull"}}
  \item{v}{Scale parameter for \code{"gamma.frailty"} and \code{"nonPH.weibull"} or variance parameter for \code{"log.normal.frailty"} models.  Not used in marginal models}
  \item{rho}{Correlation for marginal models.  Not used in other models}
  \item{a}{Parameter for \code{"nonPH.weibull"} model.  Not used in other models}
  \item{lambda}{Parameter for \code{"nonPH.weibull"} model.  Not used in other models}
}
\details{
This function generates multivariate survival data.  Letting \eqn{i=1,...,N} number of clusters, \eqn{j=1,...,K} number of units per cluster, and \eqn{X_{ij}} be a candidate covariate, the following multivariate survival models can be used:

gamma.frailty: \eqn{\hspace{2mm}} \eqn{\lambda_{ij}(t)=\exp(\beta_{0}+\beta_{1} \cdot I(X_{ij} \leq c)) w_{i}} with \eqn{w_{i} \sim \Gamma(1/v, 1/v)}

log.normal.frailty: \eqn{\hspace{2mm}}  \eqn{\lambda_{ij}(t)=\exp(\beta_{0}+\beta_{1} \cdot I(X_{ij} \leq c) + w_{i})} with \eqn{w_{i} \sim N(0, v)}

marginal.multivariate.exponential: \eqn{\hspace{2mm}}  \eqn{\lambda_{ij}(t)=\exp(\beta_{0}+\beta_{1} \cdot I(X_{ij} \leq c))} absolutely continuous

marginal.nonabsolutely.continuous: \eqn{\hspace{2mm}}  \eqn{\lambda_{ij}(t)=\exp(\beta_{0}+\beta_{1} \cdot I(X_{ij} \leq c))} not absolutely continuous

nonPH.weibull: \eqn{\hspace{2mm}} \eqn{\lambda_{ij}(t)=\lambda_{0}(t) \exp(\beta_{0}+\beta_{1} \cdot I(X_{ij} \leq c)) w_{i}} with \eqn{w_{i} \sim \Gamma(1/v ,1/v)} and

\eqn{\hspace{96mm}} \eqn{\lambda_{0}(t)=\alpha \lambda t^{\alpha-1}}

The user specifies the coefficients (\eqn{\beta_{0}} and \eqn{\beta_{1}}), the cutoff values, the censoring rate, and the model with the respective parameters.
}
\value{
  \item{dat}{The simulated data}
  \item{model}{The model used}
}
\references{
Fan J., Nunn M., Su X. (2009) Multivariate exponential survival trees and their application to tooth prognosis.  \emph{Computational Statistics and Data Analysis}, \bold{53}(4), 1110--1121.

Su X., Fan J., Wang A., Johnson M. (2006) On Simulating Multivariate Failure Times.  \emph{International Journal of Applied Mathematics & Statistics}, \bold{5}, 8--18
}
\author{Xiaogang Su, Peter Calhoun, Juanjuan Fan}
\seealso{\pkg{genSurv}, \pkg{complex.surv.dat.sim}, \pkg{survsim}}
\examples{
randMarginalExp <- rmultime(N = 200, K = 4, beta = c(-1, 2, 2, 0, 0), cutoff = c(0.5, 0.5, 0, 0),
    digits = 1, icensor = 1, model = "marginal.multivariate.exponential", rho = .65)$dat

randFrailtyGamma <- rmultime(N = 200, K = 4, beta = c(-1, 1, 3, 0), cutoff = c(0.4, 0.6, 0),
    digits = 1, icensor = 1, model = "gamma.frailty", v = 1)$dat
}
\keyword{Simulation}
\keyword{Survival}
\keyword{Multivariate}
\keyword{Correlated}
