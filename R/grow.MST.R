grow.MST <- function(
  dat,       						# the training data	
  test=NULL, 								# the test data; If it is NULL, the splitting info is replaced with the training data.
  method=c("marginal", "gamma.frailty", "exp.frailty"),			# Choose among the three MST methods, with choices "marginal", "gamma.frailty", and "exp.frailty" 
  min.ndsz=20, 							# min.ndsz controls the minimum node size 
  n0=3, 								# n0=5 controls the minimum number of uncensored event times at either child node
  col.time, col.status, col.id,					# columns for time, status, and id, respectively
  col.split.var, 							# columns of splitting variables
  col.ctg=NULL, 							# columns of categorical variables; this should be a subset of col.split.var
  max.depth=10, 							# maximum depth of tree
  mtry=length(col.split.var),					# the number of variables considered at each split
  details=FALSE,								# whether or not to print detailed information
  cont.split=c("distinct","percentiles"),  							# whether candidate splits are every distinct value or percentiles
  delta=0.05,                # when using percentiles split, only split from delta to 1-delta
  nCutPoints=50)             # when using percentiles split, specify the number cutpoints (percentiles) considered
{
  if(all(method==c("marginal", "gamma.frailty", "exp.frailty"))){method="marginal"
  } else if (!(method %in% c("marginal", "gamma.frailty", "exp.frailty"))){stop("Wrong specification of method= argument")}
  if(all(cont.split==c("distinct","percentiles"))){cont.split="distinct"
  } else if (!(cont.split %in% c("distinct","percentiles"))){stop("Wrong specification of cont.split= argument")}
  
  out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
  list.nd <- list(dat); 
  if (!is.null(test)) list.test <- list(test)
  name <- 1
  while (length(list.nd)!=0) {
    for (i in 1:length(list.nd)){
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        if (method=="marginal") split <- partition.MST(dat=list.nd[[i]], test=test0, name=name[i],
                                                       method="marginal", min.ndsz=min.ndsz, n0=n0, 
                                                       col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=NULL, 
                                                       max.depth=max.depth, details=details,cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
        else if (method=="gamma.frailty") split <- partition.MST(dat=list.nd[[i]], test=test0, name=name[i],
                                                                 method="gamma.frailty", min.ndsz=min.ndsz, n0=n0, 
                                                                 col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=NULL, 
                                                                 max.depth=max.depth, details=details, cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
        else if (method=="exp.frailty") split <- partition.MST(dat=list.nd[[i]], test=test0, name=name[i],
                                                               method="exp.frailty", min.ndsz=min.ndsz, n0=n0, 
                                                               col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=NULL, 
                                                               max.depth=max.depth, details=details, cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
        # print(split$info)
        out <- rbind(out, split$info)
        if (!is.null(split$left) && is.null(test)) {
          temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        } else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
        }		
      }
    }
    list.nd <- temp.list; list.test <- temp.test; name <- temp.name
    temp.list <- temp.test <- temp.name <- NULL
  }
  # print(out)
  out$node <- as.character(out$node)
  out <- out[order(out$node), ]
  out
}
