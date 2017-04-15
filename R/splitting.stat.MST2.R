splitting.stat.MST2 <-
function(surv, id, z, weights,
                                splitstat = c("c.loglik", "loglik", "wald.test")) {
  splitstat <- match.arg(splitstat, c("c.loglik", "loglik", "wald.test"))

  score <- NA
  # USE SMALLER NUMBER OF ITERATION TO HELP SPEED IT UP.
  control0 <- coxph.control(eps = 1e-04, toler.chol = .Machine$double.eps^0.75,
                            iter.max = 10, toler.inf = sqrt(1e-09), outer.max = 5)
  options(warn = -1)
  fit <- try(coxph(surv ~ z + frailty.gamma(id, method ="em"), weights = weights, control = control0), silent = TRUE)
  options(warn = 0)
  if (!inherits(fit, "try-error")) {
    if (splitstat == "c.loglik") { score <- fit$history$`frailty.gamma(id, method = "aic")`$c.loglik
    } else if (splitstat == "loglik") { score <- fit$loglik[[2]]
    } else if (splitstat == "wald.test") { score <- as.vector(fit$wald.test) }
    if (length(score) == 0) { score <- NA }
  }
  score
}
