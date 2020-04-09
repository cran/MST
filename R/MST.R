MST <-
function(formula, data, test = NULL, weights_data, weights_test, subset,
                method = c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"),
                minsplit = 20, minevents = 3, minbucket = round(minsplit/3), maxdepth = 10, mtry = NULL,
                distinct = TRUE, delta = 0.05, nCutPoints = 50,
                selection.method = c("test.sample", "bootstrap"),
                B = 30, LeBlanc = TRUE, min.boot.tree.size = 1,
                plot.Ga = TRUE, filename = NULL, horizontal = TRUE, details = FALSE, sortTrees = TRUE) {
  method <- match.arg(method, c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"))
  selection.method <- match.arg(selection.method, c("test.sample", "bootstrap"))
  stopifnot(is.logical(LeBlanc), is.logical(plot.Ga), is.logical(horizontal), is.logical(details), is.logical(sortTrees))
  #Set minsplit to be 2*minbucket if too small
  minsplit <- max(minsplit, 2 * minbucket)
  if (is.null(test) & selection.method == "test.sample") {
    print("No test sample supplied, changed selection.method = 'bootstrap'")
    selection.method <- "bootstrap"
  }
  #Convert character variables into factors
  data[ , names(data)] <- lapply(data[ , names(data)], function(x){ if (is.character(x)) { factor(x) } else { x } })
  if (!is.null(test)) { test[ , names(test)] <- lapply(test[ , names(test)], function(x){ if (is.character(x)) { factor(x) } else { x } }) }

  mf <- match.call(expand.dots = FALSE)
  m_data <- match(c("formula", "data", "subset", "weights_data"), names(mf), 0)
  m_test <- match(c("formula", "weights_test"), names(mf), 0)
  mf_data <- mf[c(1, m_data)]; mf_test <- mf[c(1, m_test)]
  mf_test$data = test  #Change test to data to get formula to work
  if (!is.null(test)) mf_test$subset = 1:nrow(test)  #Only apply subset to data
  names(mf_data)[grepl("weights", names(mf_data))] = "weights"
  names(mf_test)[grepl("weights", names(mf_test))] = "weights"
  
  formula <- Formula::Formula(formula)
  mf_data$formula <- mf_test$formula <- formula
  mf_data$drop.unused.levels <- mf_test$drop.unused.levels <- FALSE
  mf_data[[1]] <- mf_test[[1]] <- quote(stats::model.frame); 

  mf_data <- eval(mf_data)
  #Extract weights and remove from data
  weights_data <- model.weights(mf_data)
  if (is.null(weights_data)) weights_data <- rep(1, nrow(mf_data))
  mf_data <- mf_data[ , colnames(mf_data) != "(weights)"]
  
  if (!is.null(test)) {
    mf_test <- eval(mf_test, parent.frame())
    weights_test <- model.weights(mf_test)
    if (is.null(weights_test)) weights_test <- rep(1, nrow(mf_test))
    mf_test <- mf_test[ , colnames(mf_test) != "(weights)"]
  } else { weights_test = NULL }
  
  if(method %in% c("exp.frailty", "stratified", "independence") && 
     (length(unique(weights_data)) > 1 || length(unique(weights_test)) > 1)){
    print("Note: weights ignored.  Currently weights only work for marginal or gamma frailty model")
  }
  
  # OBTAIN THE COLUMNS
  vnames <- colnames(mf_data)
  Vs <- all.vars(formula); n.V <- length(Vs)
  col.surv <- 1
  col.id <- ncol(mf_data)
  col.split.var <- 2:(ncol(mf_data)-1)
  if (is.null(mtry)) { mtry <- length(col.split.var) }
  # EXTRACT THE COLS FOR CATEGORICAL VARIABLES
  cat.var <- vapply(mf_data[col.split.var], is.factor, logical(1))
  col.ctg <- col.split.var[cat.var]
  if (length(col.ctg) == 0) col.ctg <- NULL
  col.ctg.ordinal <- col.ctg  ### DELETE THIS LATER.

  #Will convert factors to numeric.  Save data dataset beforehand
  dataTemp <- mf_data

  if (!is.null(col.ctg.ordinal)) {
    temp <- ordinalizeFunc(mf_data, col.surv = col.surv, col.id = col.id,
                           col.ctg = col.ctg.ordinal, min.levels = 4, details = details)
    data <- temp$dat
    if (!is.null(test)) { test <- ordinalizeFunc(mf_test, col.surv = col.surv, col.id = col.id,
                                                 col.ctg = col.ctg.ordinal, min.levels = 4, details = details)$dat
    }
  } else {
    data <- mf_data
    test <- mf_test
  }

  # THE TEST SAMPLE METHOD
  # ========================
  if (selection.method == "test.sample") {
    # GROW A LARGE INITIAL TREE
    tree0 <- grow.MST(dat = data, test = test, weights_data = weights_data, weights_test = weights_test, method = method,
                      col.surv = col.surv, col.id = col.id, col.split.var = col.split.var, col.ctg = col.ctg,
                      minsplit = minsplit, minevents = minevents, minbucket = minbucket, maxdepth = maxdepth, mtry = mtry,
                      distinct = distinct, delta = delta, nCutPoints = nCutPoints, details = details)
    # PRUNING AND TREE SIZE SELECTION
    prn <- prune.size.testsample(tree0)
    # PRUNING INFORMATION, LOOK FOR THE MAXIMUM Ga.2, Ga.3, Ga.4, OR Ga.log_n
    pruning.info <- tmp <- prn$result

    col.Ga <- 7:10
    for (j in c(4, col.Ga)) tmp[,j] <- as.numeric(tmp[ , j])

    # PLOT THE G.a WITH DIFFERENT CHOICES OF a
    if (plot.Ga) {
      if (!is.null(filename)) {
        if (substr(filename, nchar(filename) - 2, nchar(filename)) == 'pdf') {
          pdf(file = filename, width = 7 + 3 * horizontal, height = 7 + 3 * (1 - horizontal))
        } else { postscript(file = filename, horizontal = horizontal) }
      }
      par(mfrow=c(1, 1), mar = rep(4, 4))   ### SET THE PLOTTING PARAMETERS
      y.min <- min(tmp[ , col.Ga]); y.max <- max(tmp[ , 7:10])
      x.min <- min(tmp[ , 4]); x.max <- max(tmp[ , 4])
      plot(c(0, x.max), c(y.min, y.max), type = "n", xlab = "# of Terminal Nodes", ylab = "G.a Values")
      for (j in col.Ga) lines(tmp$size.tmnl, tmp[ , j], col = j, lty = 1, lwd = 2)
      legend(x.min + 4, (y.min + y.max) / 2, col = col.Ga, lty = 1, lwd = 2, legend = c("Ga.2", "Ga.3", "Ga.4", "Ga.ln(n)"))
      if (!is.null(filename)) dev.off()
    }

    # OBTAIN THE BEST TREE SIZE AND BEST TREE STRUCTURE
    best.tree.size <- best.tree.structure <- as.list(NULL)
    Ga.cols <- c("Ga.2", "Ga.3", "Ga.4", "Ga.log_n")
    for (j in Ga.cols) {
      best.size <- tmp$size.tmnl[which.max(tmp[ , j])]
      best.tree.size[[j]] <- best.size
      best.tree.structure[[j]] <- obtain.btree(tree0, bsize = best.size)
    }

    # BOOTSTRAP METHOD
    # =================
  } else if (selection.method == "bootstrap") {
    # GROW AND PRUNE B BOOTSTRAP TREES
    boot.result <- bootstrap.grow.prune(B = B, data = data, method = method, weights_data = weights_data, 
                                        col.surv = col.surv, col.id = col.id, col.split.var = col.split.var, col.ctg = col.ctg,
                                        minsplit = minsplit, minevents = minevents, minbucket = minbucket, maxdepth = maxdepth, mtry = mtry,
                                        LeBlanc = LeBlanc, min.boot.tree.size = min.boot.tree.size, details = details)
    OUT <- bootstrap.size(boot.result, plot.Ga = plot.Ga, filename = filename, horizontal = horizontal)
    # THE INITIAL LARGE TREE CAN BE FOUND FROM THE RESULTS
    tree0 <- OUT$initial.tree
    pruning.info <- OUT$G.a
    best.tree.size <- OUT$bsize	# BEST TREE SIZES
    best.tree.structure <- OUT$btree  	# BEST TREE MODELS
  }

  if (sortTrees == TRUE) {
    tree0 <- sortTree(tree0, data, col.surv, col.id, col.ctg, col.ctg.ordinal)
    for (name in names(best.tree.structure)) {
      best.tree.structure[[name]] <- sortTree(best.tree.structure[[name]], data, col.surv, col.id, col.ctg, col.ctg.ordinal)
    }
  } else {
    operator <- factor(rep("<=", length(tree0$var)), levels = c("<=", ">", "in", "not in"))
    operator[is.na(tree0$var)] <- NA
    tree0 <- cbind(tree0, operator = operator)[ , c(1, 2, 3, 4, 9, 5, 6, 7, 8)]
    for (name in names(best.tree.structure)) {
      operator <- factor(rep("<=", length(best.tree.structure[[name]]$var)), levels = c("<=", ">", "in", "not in"))
      operator[is.na(best.tree.structure[[name]]$var)] <- NA
      best.tree.structure[[name]] <- cbind(best.tree.structure[[name]], operator = operator)[ , c(1, 2, 3, 4, 9, 5, 6, 7, 8)]
    }
  }

  if (!is.null(col.ctg.ordinal)) {
    treeCut <- as.character(tree0$cut)
    for (i in 1:NROW(tree0)) {
      if (tree0$var[i] %in% col.ctg.ordinal) {
        info <- get(as.character(tree0$vname[i]), temp$info)
        if (tree0$operator[i] == "<=") { treeCut[i] <- paste(info[as.numeric(info[ , 3]) <= as.numeric(as.character(tree0$cut[i])), 1], collapse = ",")
        } else if (tree0$operator[i] == ">") { treeCut[i]<-paste(info[as.numeric(info[ , 3]) > as.numeric(as.character(tree0$cut[i])), 1], collapse = ",")
        } else { stop("Unexpected operator") }
      }
    }
    tree0$cut <- treeCut
    tree0$operator[tree0$var %in% col.ctg] <- "in"

    for (name in names(best.tree.structure)) {
      tree <- best.tree.structure[[name]]
      treeCut <- as.character(tree$cut)
      for (i in 1:NROW(tree)) {
        if (tree$var[i] %in% col.ctg.ordinal) {
          info <- get(as.character(tree$vname[i]), temp$info)
          if (tree$operator[i] == "<=") { treeCut[i] <- paste(info[as.numeric(info[ , 3]) <= as.numeric(as.character(tree$cut[i])), 1], collapse = ",")
          } else if (tree$operator[i] == ">") { treeCut[i]<- paste(info[as.numeric(info[,3]) > as.numeric(as.character(tree$cut[i])), 1], collapse = ",")
          } else { stop("Unexpected operator") }
        }
      }
      best.tree.structure[[name]]$cut <- treeCut
      best.tree.structure[[name]]$operator[best.tree.structure[[name]]$var %in% col.ctg] <- "in"
    }
  }

  tree0 <- listIntoTree(tree = tree0, data = dataTemp, formula = formula, weights = weights_data)
  for (name in names(best.tree.structure)) {
    best.tree.structure[[name]] <- listIntoTree(tree = best.tree.structure[[name]], data = dataTemp, formula = formula, weights = weights_data)
  }

  return(list(tree0 = tree0, pruning.info = pruning.info,
              best.tree.size = best.tree.size, best.tree.structure = best.tree.structure))
}
