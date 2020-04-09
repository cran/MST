partition.MST <-
function(dat, test = NULL, weights_data, weights_test, name = "1",
                          method = c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"),
                          col.surv, col.id, col.split.var, col.ctg = NULL,
                          minsplit = 20, minevents = 3, minbucket = round(minsplit/3), maxdepth = 10, mtry = length(col.split.var),
                          distinct = TRUE, delta = 0.05, nCutPoints = 50, details = FALSE) {
  method <- match.arg(method, c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"))

  call <- match.call(); out <- match.call(expand.dots = FALSE)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  n <- nrow(dat);  vnames <- colnames(dat); var <- NA; cut <- NA;
  max.score <- -1e50; score.test <- NA;
  depth <- nchar(name) # CONTROL THE MAX TREE DEPTH
  surv <- dat[, col.surv]; id <- dat[, col.id]
  n.event <- sum(surv[ , 2])  # COMPUTE NUMBER OF EVENT TIMES
  if (details) print(paste("Tree ", name, ": the sample size: ", n, " and num of events: ", n.event))
  if (details) print(depth)
  if (sum(!is.element(col.ctg, col.split.var)) > 0) warning("col.ctg should be a subset of col.split.var.")
  if (depth <= maxdepth && n >= minsplit && n.event >= minevents) {
    if (method == "exp.frailty") {
      # FIT THE NULL MODEL
      X.i <- aggregate(x = surv[ , 1], by = list(id), FUN = sum)$x
      Delta.i <- aggregate(x = surv[ , 2], by = list(id), FUN = sum)$x
      dat1 <- data.frame(X.i = X.i, Delta.i = Delta.i, stringsAsFactors=FALSE)
      # print(dat1)
      x <- optim(par = c(1, 1), fn = loglik0, gr = gr0, method = "SANN", control = list(maxit = 800), hessian = FALSE, dat = dat1)    	### GLOBAL OPTIMIZATION (SIMULATED ANNEALING)
      # print(x$par); print(x$value)
      x <- optim(par = x$par, fn = loglik0, gr = gr0, method = "BFGS", hessian = TRUE, control = list(maxit = 50), dat = dat1)
      # x <- optim(par = x$par, fn = loglik0, method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 50), dat = dat1)  							### NUMERIC DERIVATIVES
      if (!is.null(x$message)) print(paste("At the ", i, "-th run, there are ", x$message, " in fitting the null model."))
      beta0 <- x$par[1]; v0 <- x$par[2]
      I0.inv <- try(solve(x$hessian), silent = TRUE)
      Ai <- 1/v0 + Delta.i; Bi <- 1 / v0 + exp(beta0) * X.i
    }

    # SEARCH FOR BEST SPLIT
    #If one col.split.var, sample function does not work correctly.  Specifically handle 1 split variable
    if (length(col.split.var) == 1) { col.split.varTemp <- rep(col.split.var, 2)
    } else { col.split.varTemp <- col.split.var }

    for (i in sort(sample(col.split.varTemp, size = mtry, replace = FALSE))) {
      z <- dat[ , i]; v.name <- vnames[i]; temp <- sort(unique(na.omit(z)));
      if (length(temp) > 1) {
        if (is.element(i, col.ctg)) { zcut <- power.set(temp)
        } else {
          if (distinct || length(temp) <= 50) { zcut <- temp[-length(temp)]
          } else if (delta > 0.2 || delta < 0.01 || nCutPoints < 5) { stop("Choice of percentile cutpoints too small")
          } else { zcut <- quantile(temp, probs = seq(delta, 1 - delta, length = nCutPoints)) }
        }
        for (j in zcut) {
          score <- NA
          if (is.element(i, col.ctg)) {
            grp <- sign(is.element(z, j)); cut1 <- paste(j, collapse = " ")
            nbucket <- min(table(grp))
          } else {
            grp <- sign(z <= j); cut1 <- as.character(j)
            nbucket <- min(table(grp))
          }
          if (nbucket >= minbucket) {
            if (method == "marginal") { score <- splitting.stat.MST1(surv, id, z = grp, weights = weights_data)
            } else if (method == "gamma.frailty") { score <- splitting.stat.MST2(surv, id, z = grp, weights = weights_data, splitstat = "wald.test")
            } else if (method == "exp.frailty") {
              mi <- aggregate(x = surv[ , 1] * grp, by = list(id), FUN = sum)$x
              zdi <- aggregate(x = surv[ , 2] * grp, by = list(id), FUN = sum)$x
              U <- sum(zdi) - exp(beta0) * sum(mi * Ai / Bi)
              I11 <- exp(beta0) * sum(mi * Ai * (Bi - exp(beta0) * mi) / Bi^2)
              I1theta <- c(exp(beta0) * sum(mi * Ai * (Bi - exp(beta0) * X.i) / Bi^2), exp(beta0) / v0^2 * sum(mi * (Ai - Bi) / Bi^2))
              if (!inherits(I0.inv, "try-error")) {
                score <- U^2/(I11 - t(I1theta) %*% I0.inv %*% I1theta)
                score <- max(0, score)
              }
            } else if (method == "stratified") { score <- splitting.stat.MST3(surv, id, z = grp)
            } else if (method == "independence") { score <- splitting.stat.MST4(surv, id, z = grp)
            } else { stop("Unexpected method") }
          }

          if (identical(score, numeric(0))) score <- NA		# TO DEAL WITH THE PROBLEM THAT THE WALD TEST COULD RETURN numeric(0).
          if (!is.na(score) && score >= max.score) {
            max.score <- score; var <- i; vname <- v.name
            cut <- cut1; best.cut<-j; grp.best <- grp
          }
          if (details) { print(cbind(var = i, v.name = v.name, cut = j, score = score, max.score = max.score)) } # print(is.null(score)); print(score)}
        }
      }
    }
  }

  # THE TEST SAMPLE
  if (!(is.null(test))) {
    n.test <- nrow(test);
    if (!(is.na(var)) && max.score != 0) {
      surv.test <- test[ , col.surv]; id.test <- test[ , col.id]
      # COMPUTE THE SCORE STAT BASED ON THE TEST SAMPLE
      if (is.element(var,col.ctg)) { grp.test <- sign(is.element(test[,var], best.cut))
      } else { grp.test <- sign(test[,var] <= best.cut) }
      neventsL <- sum(grp.test == 1 & surv.test[ , 2] == 1); neventsR <- sum(grp.test == 0 & surv.test[ , 2] == 1)
      if (method == "marginal") { score.test <- splitting.stat.MST1(surv.test, id.test, z = grp.test, weights = weights_test)
      } else if (method == "gamma.frailty") { score.test <- splitting.stat.MST2(surv.test, id.test, z = grp.test, weights = weights_test, splitstat = "wald.test")
      } else if (method == "exp.frailty") {
        X.i <- aggregate(x = surv.test[ , 1], by = list(id.test), FUN = sum)$x
        Delta.i <- aggregate(x = surv.test[ , 2], by = list(id.test), FUN = sum)$x
        dat2 <- data.frame(X.i = X.i, Delta.i = Delta.i, stringsAsFactors=FALSE)
        x1 <- optim(c(2, 1.5), fn = loglik0, gr = gr0, method = "Nelder-Mead", hessian = TRUE, dat = dat2)
        beta0 <- x1$par[1]; v0 <- x1$par[2]
        I0.inv <- try(solve(x1$hessian), silent = TRUE)
        if (!inherits(I0.inv, "try-error")) {
          Ai <- 1 / v0 + Delta.i; Bi <- 1 / v0 + exp(beta0) * X.i
          # compute the score test
          mi <- aggregate(x = surv.test[ , 1] * grp.test, by = list(id.test), FUN = sum)$x
          zdi <- aggregate(x = surv.test[ , 2] * grp.test, by = list(id.test), FUN = sum)$x
          U <- sum(zdi) - exp(beta0) * sum(mi * Ai / Bi)
          I11 <- exp(beta0) * sum(mi * Ai * (Bi - exp(beta0) * mi) / Bi^2)
          I1theta <- c(exp(beta0) * sum(mi * Ai * (Bi - exp(beta0) * X.i) / Bi^2), exp(beta0) / v0^2 * sum(mi * (Ai - Bi) / Bi^2))
          score.test <- U^2 / (I11 - t(I1theta) %*% I0.inv %*% I1theta)
          score.test <- ifelse(score.test < 0, NA, score.test)
        }
      } else if (method == "stratified") { score.test <- splitting.stat.MST3(surv.test, id.test, z = grp.test)
      } else if (method=="independence") { score.test <- splitting.stat.MST4(surv.test, id.test, z = grp.test)
      } else { stop("Unexpected method") }
    }
  } else { n.test <- n; score.test <- max.score; neventsL <- neventsR <- 1}

  if (!is.na(score.test) && !is.na(var) && min(c(neventsL, neventsR)) >= 1) {
    out$name.l <- paste(name, 1, sep = ""); out$name.r <- paste(name, 2, sep = "")
    if (!is.null(test)) { out$left.test <- test[grp.test == 1, ]; out$right.test <- test[grp.test == 0, ]
    } else { out$left.test <- out$right.test <- NULL }
    out$left  <- dat[grp.best == 1, ];  out$right <- dat[grp.best == 0, ]
    
    out$left.weights.data <- weights_data[grp.best == 1]; out$right.weights.data <- weights_data[grp.best == 0]
    if (!is.null(test)) { out$left.weights.test <- weights_test[grp.test == 1]; out$right.weights.test <- weights_test[grp.test == 0] }
    
  } else { var <- NA; vname <- NA; cut <- NA;  max.score <- score.test <- NA }
  out$info <- data.frame(node = name, size = n, var = var, vname = vname, cut = cut,
                         score = ifelse(max.score == -1e10, NA, max.score), size.test = n.test, score.test, stringsAsFactors=FALSE)
  out
}
