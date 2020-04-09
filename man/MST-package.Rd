\name{MST-package}
\alias{MST-package}
\alias{is.odd}
\alias{as.numeric.factor}
\alias{ordinalizeFunc}
\alias{power.set}
\alias{loglik0}
\alias{gr0}
\alias{partition.MST}
\alias{grow.MST}
\alias{de}
\alias{prune.size.testsample}
\alias{prune.size}
\alias{bootstrap.grow.prune}
\alias{bootstrap.size}
\alias{obtain.btree}
\alias{send.down}
\alias{sortTree}
\alias{MST.plot}
\alias{MST.plot.latex}
\alias{listIntoParty}
\alias{listIntoTree}
\docType{package}
\title{Multivariate Survival Trees Package}
\description{
This package constructs trees for multivariate survival data using marginal and frailty models
}
\details{
\tabular{ll}{
Package: \tab MST\cr
Type: \tab Package\cr
Version: \tab 2.2\cr
Date: \tab 2020-04-05\cr
License: \tab GPL-2 \cr
}
Decision trees require few statistical assumptions, handle a variety of data structures, and provide meaningful interpretations.  There are several R packages that provide functions to construct survival trees (see \pkg{rpart}, \pkg{partykit}, and \pkg{DStree}); this package extends the implementation to multivariate survival data. There are two main approaches to analyzing correlated failure times.  One is the marginal approach studied by authors Wei et al. (1989) and Liang et al. (1993).  In the marginal model, the correlation is modeled implicitly using generalized estimating equations on the marginal distribution formulated by the Cox (1972) proportional hazards model.  The other approach is the frailty model studied by Clayton (1978) and Clayton and Cuzick (1985).  In the frailty model, the correlation is modeled explicitly by a multiplicative random effect called frailty, which corresponds to some common unobserved characteristics shared by all correlated times.


The construction of the tree adopts a modified CART procedure controlling for the correlated failure times.  The procedure consists of three stages: growing the initial tree, pruning the tree, and selecting the best-sized subtree; details of these steps are described elsewhere (Fan et al. [2006], Su and Fan [2004], and Fan et al. [2009]).  There are two methods for selecting the best-sized subtree.  When the dataset is large, one may divide the dataset into a training sample to grow and prune the initial tree and a test sample to select the best-sized tree.  When the dataset is small, one can resample the dataset to choose the best-sized subtree.
}
\author{
Xiaogang Su, Peter Calhoun, & Juanjuan Fan

Maintainer: Peter Calhoun <calhoun.peter@gmail.com>
}
\references{

Calhoun P., Su X., Nunn M., Fan J. (2018) Constructing Multivariate Survival Trees: The MST Package for R. \emph{Journal of Statistical Software}, \bold{83}(12), 1--21.

Clayton D.G. (1978) A model for association in bivariate life tables and its application in epidemiologic studies of familial tendency in chronic disease incidence. \emph{Biometrika}, \bold{65}(1), 141--151

Clayton D.G. and Cuzick J. (1985) Multivariate generalization of the proportional hazards model.  \emph{Journal of the Royal Statistical Society Series A}, \bold{148}(2), 82--108

Cox D.R. (1972) Regression models and life-tables (with discussion).  \emph{Journal of the Royal Statistical Society Series B}, \bold{34}(2), 187--220.

Fan J., Su X., Levine R., Nunn M., LeBlanc M. (2006) Trees for Correlated Survival Data by Goodness of Split, With Applications to Tooth Prognosis. \emph{Journal of American Statistical Association}, \bold{101}(475), 959--967.

Fan J., Nunn M., Su X. (2009) Multivariate exponential survival trees and their application to tooth prognosis.  \emph{Computational Statistics and Data Analysis}, \bold{53}(4), 1110--1121.

Liang K.Y., Self S.G., Chang Y. (1993) Modeling marginal hazards in multivariate failure time data.  \emph{Journal of the Royal Statistical Society Series B}, \bold{55}(2), 441--453

Su X., Fan J. (2004) Multivariate Survival Trees: A Maximum Likelihood Approach Based on Frailty Models.  \emph{Biometrics}, \bold{60}(1), 93--99.

Su X., Fan J., Wang A., Johnson M. (2006) On Simulating Multivariate Failure Times.  \emph{International Journal of Applied Mathematics & Statistics}, \bold{5}, 8--18

Wei L.J., Lin D.Y., Weissfeld L. (1989) Regression analysis of multivariate incomplete failure time data by modeling marginal distributions.  \emph{Journal of the American Statistical Association}, \bold{84}(408), 1065--1073
}
