\name{MST-package}
\alias{MST-package}
\alias{as.numeric.factor}
\alias{de}
\alias{gr0}
\alias{is.odd}
\alias{loglik0}
\alias{obtain.btree}
\alias{ordinalize}
\alias{partition.MST}
\alias{power.set}
\alias{prune.size}
\alias{send.down}
\alias{splitting.stat.MST1}
\alias{splitting.stat.MST2}
\docType{package}
\title{Multivariate Survival Trees Package}
\description{
This package constructs trees for multivariate survival data using marginal and frailty models
}
\details{
\tabular{ll}{
Package: \tab MST\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-19\cr
License: \tab GPL-2 \cr
}
Regression trees require few statistical assumptions, handle a variety of data structures, and provide meaningful interpretations.  There are several functions to construct regression trees (see \pkg{rpart}, \pkg{tree}, and \pkg{party}); this package extends the implementation to multivariate survival data. There are two main approaches to analyzing correlated failure times.  One is the marginal approach studied by authors Wei et al. (1989) and Liang et al. (1993).  In the marginal model, the correlation is modeled implicitly using estimating equations on the marginal distribution formulated by the Cox (1972) proportional hazards model.  The other approach is the frailty model studied by Clayton (1978) and Clayton and Cuzick (1985).  In the frailty model, the correlation is modeled explicitly by a multiplicative random effect called frailty, which corresponds to some common unobserved characteristics shared by all correlated times.


Construction of the tree consists of three stages: growing the initial tree, pruning the tree, and selecting the best-sized subtree; details of these steps are described elsewhere (Fan et al. [2006], Su and Fan [2004], and Fan et al. [2009]).  There are two methods for selecting the best-sized subtree.  When the dataset is large, one may divide the dataset into a training sample to grow and prune the initial tree and a test sample to select the best-sized tree.  When the dataset is small, one can resample the dataset to choose the best-sized subtree. 
}
\author{
Xiaogang Su, Peter Calhoun, & Juanjuan Fan

Maintainer: Peter Calhoun <calhoun.peter@gmail.com>
}
\references{

Clayton D.G. (1978) A model for association in bivariate life tables and its application in epidemiologic studies of familial tendency in chronic disease incidence. Biometrika, 65, 141-151

Clayton D.G. and Cuzick J. (1985) Multivariate generalization of the proportional hazards model.  Journal of the Royal Statistical Society Series A, 148, 82-108

Cox D.R. (1972) Regression models and life-tables (with discussion).  Journal of the Royal Statistical Society Series B, 34, 187-202

Fan J., Su X., Levine R., Nunn M., LeBlanc M. (2006) Trees for Correlated Survival Data by Goodness of Split, With Applications to Tooth Prognosis. Journal of American Statistical Association, 101, 959-967

Fan J., Nunn M., Su X. (2009) Multivariate exponential survival trees and their application to tooth prognosis.  Computational Statistics and Data Analysis, 53, 1110-1121

Liang K.Y., Self S.G., Chang Y. (1993) Modeling marginal hazards in multivariate failure time data.  Journal of the Royal Statistical Society Series B, 55, 441-453

Su X., Fan J. (2004) Multivariate Survival Trees: A Maximum Likelihood Approach Based on Frailty Models.  Biometrics, 60, 93-99

Su X., Fan J., Wang A., Johnson M. (2006) On Simulating Multivariate Failure Times.  International Journal of Applied Mathematics & Statistics, 5, 8-18

Wei L.J., Lin D.Y., Weissfeld L. (1989) Regression analysis of multivariate incomplete failure time data by modeling marginal distributions.  Journal of the American Statistical Association, 84, 1065-1073
}