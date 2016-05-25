splitting.stat.MST3 <-
function(time, status, id, z, min.nevents){  
  n1 <- sum(z==1&status==1); n2 <- sum(z==0&status==1); 
  score <- NA
  if (all(!is.na(c(n1, n2))) && min(n1, n2)>=min.nevents){
    options(warn=-1)
    #Unfortunately, survdiff spits out error when singular.  Create new function
    fit<-try(survdiff(Surv(time, status) ~ z + strata(id)), silent=TRUE) 
    if(is(fit, "try-error")){score <- NA
    } else {score <- fit$chisq}
    options(warn=0)
  }
  score
}
