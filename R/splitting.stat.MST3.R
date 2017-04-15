splitting.stat.MST3 <-
function(surv, id, z) {
  score <- NA
  options(warn = -1)
  fit <- try(survdiff(surv ~ z + strata(id)), silent = TRUE)
  options(warn = 0)
  if (!inherits(fit, "try-error")) { score <- fit$chisq }
  score
}
