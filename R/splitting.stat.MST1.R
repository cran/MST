splitting.stat.MST1 <-
function(time, status, id, z, n0){  
  n1 <- sum(z==1&status==1); n2 <- sum(z==0&status==1); 
  score <- NA
  if (all(!is.na(c(n1, n2))) && min(n1, n2)>=n0){
    options(warn=-1)
    fit <- coxph(Surv(time, status) ~ z + cluster(id))
    options(warn=0)
    score <- fit$rscore
  }
  score
}
