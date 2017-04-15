splitting.stat.MST1 <-
function(surv, id, z, weights) {
  score <- NA
  options(warn = -1)
  fit <- try(coxph(surv ~ z + cluster(id), weights=weights), silent = TRUE)
  options(warn = 0)
  if (!inherits(fit, "try-error")) { score <- fit$rscore }
  score
}
