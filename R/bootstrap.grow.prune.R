bootstrap.grow.prune <-
function(data, 									# data frame containing the training set
                                 method=c("marginal", "gamma.frailty", "exp.frailty"),	# Choose among the three MST methods, with choices "marginal", "gamma.frailty", and "exp.frailty" 
                                 min.ndsz=20, 								# min.ndsz controls the minimum node size 
                                 min.nevents=3, 									# min.nevents controls the minimum number of uncensored event times at either child node
                                 col.time, col.status, col.id,					# columns for time, status, and id, respectively
                                 col.split.var, 							# columns of splitting variables
                                 col.ctg=NULL, 								# columns of categorical variables; this should be a subset of col.split.var
                                 max.depth=10, 								# maximum depth of tree
                                 mtry=length(col.split.var),					# the number of variables considered at each split
                                 B=30,     								# number of bootstrap samples used
                                 LeBlanc=TRUE, 								# OPTION LeBlanc IS TO APPLY THE WHOLE SAMPLE (TRUE) OR THE OUT-OF-BAD SAMPLE (FALSE) IN THE BOOTSTRAP PROCEDURE
                                 min.boot.tree.size=1,						# OPTION min.boot.tree.size IS TO MAKE SURE A NON-NULL TREE IS GROWN AT EACH BOOTSTRAP
                                 details=FALSE,  					# whether or not to print detailed information
                                 cont.split=c("distinct","percentiles"),    						# whether candidate splits are every distinct value or percentiles
                                 delta=0.05,                # when using percentiles split, only split from delta to 1-delta
                                 nCutPoints=50)             # when using percentiles split, specify the number cutpoints (percentiles) considered
{
  if(all(method==c("marginal", "gamma.frailty", "exp.frailty"))){method="marginal"
  } else if (!(method %in% c("marginal", "gamma.frailty", "exp.frailty"))){stop("Wrong specification of method= argument")}
  if(all(cont.split==c("distinct","percentiles"))){cont.split="distinct"
  } else if (!(cont.split %in% c("distinct","percentiles"))){stop("Wrong specification of cont.split= argument")}
  
  call <- match.call(); out <- match.call(expand.dots = FALSE)
  out$boot.tree <- out$boot.prune <- out$... <- NULL
  time.start <- date()
  tree0 <- grow.MST(dat=data, test=NULL, method=method, min.ndsz=min.ndsz, min.nevents=min.nevents, 
                    col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                    max.depth=max.depth, mtry=mtry, details=details, cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
  if (NROW(tree0)==1) {print(tree0); stop("Your initial tree is a NULL tree! There is no need for tree size selection.")}
  # print(tree0);
  prune0 <- prune.size(tree0); 
  boot.tree <- list(tree0); boot.prune <- list(prune0) 
  b <- 1
  while (b <= B) {
    # SAMPLING OBSERVATION
    samp <- sample(1:nrow(data), size=nrow(data), replace=TRUE) 
    dat <- data[samp, ]
    dat.oob <- data[-unique(samp),]
    n.oob <- nrow(dat.oob); # print(n.oob)
    if (LeBlanc) tree <- grow.MST(dat=dat, test=data, method=method, min.ndsz=min.ndsz, min.nevents=min.nevents,
                                  col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                                  max.depth=max.depth, mtry=mtry, details=details, cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
    else  tree <- grow.MST(dat=dat, test=dat.oob, method=method, min.ndsz=min.ndsz, min.nevents=min.nevents,
                           col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                           max.depth=max.depth, mtry=mtry, details=details, cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
    if (details) print(tree)
    if (nrow(tree)> min.boot.tree.size) {
      print(paste("###################### b = ", b, " ###########################", sep=""))
      boot.tree <- c(boot.tree, list(tree))
      prune <- prune.size(tree); # print(prune)
      boot.prune <- c(boot.prune, list(prune))
      b <- b+1
    }
  }
  time.end <- date()
  print(paste("The Start and End time for ", B, "bootstrap runs is:"))
  print(rbind(time.start, time.end))
  out$boot.tree <- boot.tree
  out$boot.prune <- boot.prune
  # THE INITIAL LARGE TREE
  out$initial.tree <- tree0
  out
}
