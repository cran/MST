grow.MST <-
function(dat, test = NULL, weights_data, weights_test,
                     method = c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"),
                     col.surv, col.id, col.split.var, col.ctg = NULL,
                     minsplit = 20, minevents = 3, minbucket = round(minsplit/3), maxdepth = 10, mtry = length(col.split.var),
                     distinct = TRUE, delta = 0.05, nCutPoints = 50, details = FALSE) {
  method <- match.arg(method, c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"))

  out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- temp.weights.data <- temp.weights.test <- NULL
  list.nd <- list(dat);
  if (!is.null(test)) list.test <- list(test)
  list.weights.data <- list(weights_data)
  list.weights.test <- list(weights_test)
  name <- 1
  
  while (length(list.nd) != 0) {
    for (i in 1:length(list.nd)) {
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 0) {
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        split <- partition.MST(dat = list.nd[[i]], test = test0,
                               weights_data = list.weights.data[[i]], weights_test = list.weights.test[[i]],
                               name = name[i], method = method,
                               col.surv = col.surv, col.id = col.id, col.split.var = col.split.var, col.ctg = NULL,
                               minsplit = minsplit, minevents = minevents, minbucket = minbucket, maxdepth = maxdepth, mtry = mtry,
                               distinct = distinct, delta = delta, nCutPoints = nCutPoints, details = details)
        # print(split$info)
        out <- rbind(out, split$info)
        
        if (!is.null(split$left) && is.null(test)) {
          temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.weights.data <- c(temp.weights.data, list(split$left.weights.data, split$right.weights.data))
          temp.weights.test <- c(temp.weights.test, list(split$left.weights.test, split$right.weights.test))
        } else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.weights.data <- c(temp.weights.data, list(split$left.weights.data, split$right.weights.data))
          temp.weights.test <- c(temp.weights.test, list(split$left.weights.test, split$right.weights.test))
        }
      }
    }
    list.nd <- temp.list; list.test <- temp.test; name <- temp.name
    list.weights.data <- temp.weights.data; list.weights.test <- temp.weights.test
    temp.list <- temp.test <- temp.name <- temp.weights.data <- temp.weights.test <- NULL
  }
  out$node <- as.character(out$node)
  out <- out[order(out$node), ]
  out
}
