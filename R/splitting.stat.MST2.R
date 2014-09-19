splitting.stat.MST2 <- function(time, status, id, z, n0, 
                                method=c("c.loglik", "loglik", "wald.test")){
  if(all(method==c("c.loglik", "loglik", "wald.test"))){method="c.loglik"
  } else if (!(method %in% c("c.loglik", "loglik", "wald.test"))){stop("Wrong specification of method= argument")}
  
  n1 <- sum(z==1&status==1); n2 <- sum(z==0&status==1); 
  score <- NA
  # USE SMALLER NUMBER OF ITERATION TO HELP SPEED IT UP.
  options(warn=-1)
  control0 <- coxph.control(eps = 1e-04, toler.chol = .Machine$double.eps^0.75,
                            iter.max = 10, toler.inf = sqrt(1e-09), outer.max = 5)
  if (all(!is.na(c(n1, n2))) && min(n1, n2)>=n0){
    # print(cbind(length(time), length(z), length(id)))    		######################
    fit <- coxph(Surv(time, status) ~ z + frailty.gamma(id, method ="em"), control=control0)
    if (method=="c.loglik") score <- fit$history$`frailty.gamma(id, method = "aic")`$c.loglik
    else if (method=="loglik") score <- fit$loglik[[2]]
    else if (method=="wald.test") score <- as.vector(fit$wald.test)
  }
  options(warn=0)
  score
}
