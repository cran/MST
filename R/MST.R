MST <- function(training,   							# data frame containing the training set
                test=NULL,							# The test sample if available
                method=c("marginal", "gamma.frailty", "exp.frailty"),			# Choose among the three MST methods, with choices "marginal", "gamma.frailty", and "exp.frailty" 
                # Arguments related to tree constrution
                min.ndsz=20, 							# min.ndsz controls the minimum node size 
                min.nevents=3, 								# min.nevents controls the minimum number of uncensored event times at either child node
                col.time, col.status, col.id,					# columns for time, status, and id, respectively
                col.split.var, 							# columns of splitting variables
                col.ctg=NULL, 							# columns of categorical variables; this should be a subset of col.split.var
                col.ctg.ordinal=col.ctg,           # columns of categorical variables that should be ordinalized; should be a subset of col.ctg
                max.depth=10, 							# maximum depth of tree
                mtry=length(col.split.var),					# the number of variables considered at each split
                # Argumetns related to tree size selection
                selection.method=c("test.sample", "bootstrap"),	# two tree size selection method "test.sample" or "bootstrap"
                B = 30, 								# number of bootstrap samples used
                LeBlanc=TRUE, 							# OPTION LeBlanc IS TO APPLY THE WHOLE SAMPLE (TRUE) OR THE OUT-OF-BAD SAMPLE (FALSE) IN THE BOOTSTRAP PROCEDURE
                min.boot.tree.size=1,						# OPTION min.boot.tree.size IS TO MAKE SURE A NON-NULL TREE IS GROWN AT EACH BOOTSTRAP
                plot.Ga=TRUE, filename=NULL, horizontal=TRUE,		# Whether or not you want to plot the G.a vs. tree size
                details=FALSE,   						# whether or not to print detailed information
                cont.split=c("distinct","percentiles"),     					# whether candidate splits are every distinct value or percentiles
                delta=0.05,                # when using percentiles split, only split from delta to 1-delta
                nCutPoints=50,           # when using percentiles split, specify the number cutpoints (percentiles) considered
                sortTrees=TRUE)           # sort trees so that splits with less risk of failure go to the left
{
  if(all(method==c("marginal", "gamma.frailty", "exp.frailty"))){method="marginal"
  } else if (!(method %in% c("marginal", "gamma.frailty", "exp.frailty"))){stop("Wrong specification of method= argument")}
  if(all(selection.method==c("test.sample", "bootstrap"))){selection.method=ifelse(is.null(test),"bootstrap","test.sample")
  } else if (!(selection.method %in% c("test.sample", "bootstrap"))){stop("Wrong specification of selection.method= argument")}
  if(all(cont.split==c("distinct","percentiles"))){cont.split="distinct"
  } else if (!(cont.split %in% c("distinct","percentiles"))){stop("Wrong specification of cont.split= argument")}
  
  if(any(is.na(training[,c(col.time,col.status,col.id,col.split.var)])) | 
     (ifelse(!is.null(test),any(is.na(test[,c(col.time,col.status,col.id,col.split.var)])),FALSE))){
    print("Note: Data with missing values deleted")
    training<-training[complete.cases(training[,c(col.time,col.status,col.id,col.split.var)]),]
    test<-test[complete.cases(test[,c(col.time,col.status,col.id,col.split.var)]),]
  }
  
  if(!is.null(col.ctg.ordinal)){
    temp<-ordinalizeFunc(training,col.time=col.time,col.status=col.status,col.id=col.id,
                         col.ctg=col.ctg.ordinal,min.levels=4,details=details)
    training<-temp$dat
    test<-ordinalizeFunc(test,col.time=col.time,col.status=col.status,col.id=col.id,
                         col.ctg=col.ctg.ordinal,min.levels=4,details=details)$dat
  }
  
  # THE TEST SAMPLE METHOD
  # ========================
  if (selection.method == "test.sample") {
    if (is.null(test)) stop("You cannot use the test.sample method without a test sample.")
    # GROW A LARGE INITIAL TREE
    if (method=="marginal") {
      tree0 <- grow.MST(dat=training, test=test, 
                        method="marginal", min.ndsz=min.ndsz, min.nevents=min.nevents, 
                        col.time=col.time, col.status=col.status, col.id=col.id,
                        col.split.var=col.split.var, col.ctg=col.ctg,
                        max.depth=max.depth, mtry=length(col.split.var), details=details,
                        cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
    } else if (method=="gamma.frailty") {
      tree0 <- grow.MST(dat=training, test=test, 
                        method="gamma.frailty", min.ndsz=min.ndsz, min.nevents=min.nevents,
                        col.time=col.time, col.status=col.status, col.id=col.id,
                        col.split.var=col.split.var, col.ctg=col.ctg,
                        max.depth=max.depth, mtry=length(col.split.var), details=details,
                        cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
    } else if (method=="exp.frailty") {
      tree0 <- grow.MST(dat=training, test=test, 
                        method="exp.frailty", min.ndsz=min.ndsz, min.nevents=min.nevents, 
                        col.time=col.time, col.status=col.status, col.id=col.id,
                        col.split.var=col.split.var, col.ctg=col.ctg,
                        max.depth=max.depth, mtry=length(col.split.var), details=details,
                        cont.split=cont.split,delta=delta,nCutPoints=nCutPoints)
    }
    
    # PRUNING AND TREE SIZE SELECTION
    prn <- prune.size.testsample(tree0)
    # PRUNING INFORMATION, LOOK FOR THE MAXIMUM Ga.2, Ga.3, Ga.4, OR Ga.BIC
    pruning.info <- tmp <- prn$result
    
    # PLOT THE G.a WITH DIFFERENT CHOICES OF a
    if (plot.Ga) {
      if (!is.null(filename)){
        if (substr(filename,nchar(filename)-2,nchar(filename))=='pdf'){
          pdf(file=filename, width=7+3*horizontal, height=7+3*(1-horizontal))
        } else {postscript(file=filename, horizontal=horizontal)}
      }
      par(mfrow=c(1, 1), mar=rep(4, 4))   ##################### SET THE PLOTTING PARAMETERS 
      col.Ga <- 7:10
      for (j in c(4, col.Ga)) tmp[,j] <- as.numeric.factor(tmp[,j])
      y.min <- min(tmp[, col.Ga]); y.max <- max(tmp[, 7:10])
      x.min <- min(tmp[, 4]); x.max <- max(tmp[, 4])
      plot(c(0, x.max), c(y.min, y.max), type="n", xlab="# of Terminal Nodes", ylab="G.a Values")
      for (j in col.Ga) lines(tmp$size.tmnl, tmp[,j], col=j, lty=1, lwd=2)
      legend((x.min+x.max)/2, (y.min+y.max)/2, col=col.Ga, lty=1, lwd=2, legend=colnames(tmp)[col.Ga])
      if (!is.null(filename)) dev.off()
    }
    
    # OBTAIN THE BEST TREE SIZE AND BEST TREE STRUCTURE
    best.tree.size <- best.tree.structure <- as.list(NULL)
    Ga.cols <- c("Ga.2", "Ga.3", "Ga.4", "Ga.BIC")
    for (j in Ga.cols) {
      best.size <- tmp$size.tmnl[which.max(tmp[,j])]
      best.tree.size[[j]] <- best.size
      best.tree.structure[[j]] <- obtain.btree(tree0, bsize=best.size)
    }
    
    # BOOTSTRAP METHOD
    # =================
  } else if (selection.method == "bootstrap") {
    # GROW AND PRUNE B BOOTSTRAP TREES
    boot.result <- bootstrap.grow.prune(B=B, data=training, method=method, 
                                        min.ndsz=min.ndsz, min.nevents=min.nevents, 
                                        col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                                        max.depth=max.depth, mtry=mtry, LeBlanc=LeBlanc, min.boot.tree.size=min.boot.tree.size, details=details)
    OUT <- bootstrap.size(boot.result, plot.Ga=plot.Ga, filename=filename, horizontal=horizontal)
    # THE INITIAL LARGE TREE CAN BE FOUND FROM THE RESULTS
    tree0 <- OUT$initial.tree    
    pruning.info <- OUT$G.a
    best.tree.size <- OUT$bsize	# BEST TREE SIZES
    best.tree.structure <- OUT$btree  	# BEST TREE MODELS
  }
  
  if(sortTrees==TRUE){
    tree0<-sortTree(tree0,training,col.time,col.status,col.id,col.ctg,col.ctg.ordinal)
    for(name in names(best.tree.structure)){
      best.tree.structure[[name]]=sortTree(best.tree.structure[[name]],training,col.time,col.status,col.id,col.ctg,col.ctg.ordinal)
    }
  } else {
    operator=ifelse(tree0$var %in% col.ctg,"in","<="); operator[is.na(tree0$var)]=NA
    tree0<-cbind(tree0,operator=operator)[,c(1,2,3,4,9,5,6,7,8)]
    for(name in names(best.tree.structure)){
      operator=ifelse(best.tree.structure[[name]]$var %in% col.ctg,"in","<="); operator[is.na(best.tree.structure[[name]]$var)]=NA
      best.tree.structure[[name]]=cbind(best.tree.structure[[name]],operator=operator)[,c(1,2,3,4,9,5,6,7,8)]
    }
  }
  
  if(!is.null(col.ctg.ordinal)){
    treeCut<-as.character(tree0$cut)
    for(i in 1:NROW(tree0)){
      if(tree0$var[i] %in% col.ctg.ordinal){
        info<-get(as.character(tree0$vname[i]), temp$info)
        if(tree0$operator[i] == "in"){treeCut[i]<-paste(info[as.numeric(info[,3]) <= as.numeric(as.character(tree0$cut[i])),1],collapse=",")
        } else if (tree0$operator[i] == "not in"){treeCut[i]<-paste(info[as.numeric(info[,3]) > as.numeric(as.character(tree0$cut[i])),1],collapse=",")
        } else {stop("Unexpected operator")}
      }
    }
    tree0$cut<-treeCut
    for(name in names(best.tree.structure)){
      tree<-best.tree.structure[[name]]
      treeCut<-as.character(tree$cut)
      for(i in 1:NROW(tree)){
        if(tree$var[i] %in% col.ctg.ordinal){
          info<-get(as.character(tree$vname[i]), temp$info)
          if(tree$operator[i] == "in"){treeCut[i]<- paste(info[as.numeric(info[,3]) <= as.numeric(as.character(tree$cut[i])),1],collapse=",")
          } else if (tree$operator[i] == "not in"){treeCut[i]<- paste(info[as.numeric(info[,3]) <= as.numeric(as.character(tree$cut[i])),1],collapse=",")
          } else {stop("Unexpected operator")}
        }
      }
      best.tree.structure[[name]]$cut<-treeCut
    }
  }
  
  return(list(tree0=tree0, pruning.info=pruning.info, 
              best.tree.size=best.tree.size, best.tree.structure=best.tree.structure))
}
