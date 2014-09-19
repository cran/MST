rmultime <- function(beta=c(-1, 2, 1, 0, 0),      						# beta vector values with a beta0 (dim=p+1) 
                     cutoff=c(.5, .5, 0, 0),								# cutoff values for each x (dim <= p)
                     digits=1,										# rounding digits, suggested to be 1 or 2
                     icensor=1, 										# control for censoring rate: 1 - 50%;
                     model = c("gamma.frailty", "log.normal.frailty", 				# Models A-E in Fan, Nunn, and Su (2009)
                               "marginal.multivariate.exponential",
                               "marginal.nonabsolutely.continuous", "nonPH.weibull"),
                     v=1, 											# v for the frailty model - the scale(or 1/rate) parameter in gamma or the variance parameter in normal  
                     rho=.65, 										# the correlation used in marginal model with multivariate exponential (Model C)
                     a=1.5, lambda=0.1, 								# two parameters used in the non-PH model (Model E)
                     N=100, K=4)										# N- number of units; K - number of correlated failures for each unit
{
  if(all(model==c("gamma.frailty", "log.normal.frailty","marginal.multivariate.exponential",
                  "marginal.nonabsolutely.continuous", "nonPH.weibull"))){model="gamma.frailty"
  } else if (!(model %in% c("gamma.frailty", "log.normal.frailty","marginal.multivariate.exponential",
                            "marginal.nonabsolutely.continuous", "nonPH.weibull"))){stop("Wrong specification of model= argument")}
  
  #### Generate Covariates
  p <- length(beta)-1 
  X <- matrix(0, nrow=N*K, ncol=p)
  for (j in 1:p) {
    if (is.odd(j)) X[,j] <- rep(runif(N, 0, 1), rep(K,N))
    else X[, j] <- runif(N*K)
  }
  X <- round(X, digits=digits)
  Z <- X
  for (j in 1:length(cutoff)) {
    if (cutoff[j] > 0 && cutoff[j] <1) Z[,j] <- sign(X[,j] <= cutoff[j])
  }
  Z <- cbind(1, Z)
  eta <- Z%*%beta	
  
  #### Generate Multivariate Survival and Censoring Time Data from Different Models
  if (model=="gamma.frailty") {
    if (v==0){w <- 1
    } else {w <- rep(rgamma(N, 1/v, 1/v), rep(K, N))}
    rate <- exp(eta)*w 
    xobs <- rexp(N*K, rate)
    dind <- 1
    if (icensor!=0){     
      c <- rexp(N*K, rate)
      dind <- sign(xobs <= c*icensor)
      xobs <- pmin(xobs, c*icensor) 
    } 
  } else if (model=="log.normal.frailty") {
    w <- 1
    if ( v!=0) w <- rep(rnorm(N, 0, v), rep(K, N))
    rate <- exp(eta + w)
    xobs <- rexp(N*K, rate)
    dind <- 1
    if (icensor!=0){     
      c <- rexp(N*K, rate)
      dind <- sign(xobs <= c*icensor)
      xobs <- pmin(xobs, c*icensor) 
    }
  } else if (model=="marginal.multivariate.exponential") {
    cor <- matrix(rho, K, K)
    diag(cor) <- 1
    x <- -log(1-pnorm(as.vector(t(mvrnorm(N,mu=rep(0,K),Sigma=cor)))))
    c <- -log(1-pnorm(as.vector(t(mvrnorm(N,mu=rep(0,K),Sigma=cor)))))
    if (icensor == 0) {
      xobs <- x* exp(-eta)
      dind <- 1
    } else {
      dind <- sign(x <= c*icensor)
      xobs <- pmin(x, c*icensor)*exp(-eta)
    }    
  } else if (model =="marginal.nonabsolutely.continuous"){
    lambda0 <- 2 - 2/(rho+1)    
    x <- pmin(rep(rexp(N, lambda0), rep(K,N)),  rexp(N*K, 1-lambda0))
    c <- pmin(rep(rexp(N, lambda0), rep(K,N)),  rexp(N*K, 1-lambda0))
    if (icensor==0) {
      xobs <- x* exp(-eta)
      dind <- 1
    } else {     
      xobs <- pmin(x, c*icensor)* exp(-eta)
      dind <- sign(x <= c*icensor)
    }
  } else if (model=="nonPH.weibull"){
    if (v==0) w <- 1
    else {w <- rep(rgamma(N, 1/v, 1/v), rep(K, N))}
    rate <- exp(eta)*w
    lambda1 <- lambda*rate
    b <- (lambda1)^-(1/a) 
    xobs <- rweibull(N*K, shape=a, scale=b) 
    dind <- 1
    if (icensor!=0) {
      c <- rweibull(N*K, shape=a, scale=b) 
      dind <- sign(xobs <= c*icensor)
      xobs <- pmin(xobs, c*icensor) 
    }
  }
  
  ##### Output
  colnames(X) <- paste("x", 1:p, sep="")
  dat <- data.frame(id=rep(1:N, rep(K,N)),rep=rep(1:K, N), 
                    time=xobs, status=dind, X)
  return(list(dat=dat, model=model))			  
}