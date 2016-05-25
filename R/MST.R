MST <-
function(formula, training, test=NULL, method=c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"), 
                minsplit=20, min.nevents=3, max.depth=10, mtry=NULL,
                cont.split=c("distinct","percentiles"), delta=0.05, nCutPoints=50,
                selection.method=c("test.sample", "bootstrap"),
                B = 30, LeBlanc=TRUE, min.boot.tree.size=1,
                plot.Ga=TRUE, filename=NULL, horizontal=TRUE, details=FALSE,sortTrees=TRUE)
{
  method<-match.arg(method,c("marginal", "gamma.frailty", "exp.frailty", "stratified", "independence"))
	cont.split<-match.arg(cont.split,c("distinct","percentiles"))
	selection.method<-match.arg(selection.method,c("test.sample", "bootstrap"))
	if (is.null(test) & selection.method=="test.sample"){
  		print("No test sample supplied, changed selection.method='bootstrap'")
  		selection.method<-"bootstrap"
	}
	#Convert character variables into factors
	training[,names(training)] <- lapply(training[,names(training)] , function(x){if(is.character(x)){factor(x)} else {x}})
	if(!is.null(test)){test[,names(test)] <- lapply(test[,names(test)] , function(x){if(is.character(x)){factor(x)} else {x}})}
	
	# OBTAIN THE COLUMNS 
	form0 <- formula; vnames <- colnames(training)
	Vs <- all.vars(form0); n.V <- length(Vs)
	col.time <- which(vnames == (Vs[1]))
	col.status <- which(vnames == Vs[2])
	col.id <- which(vnames == Vs[n.V])
	col.split.var <- which(is.element(vnames, Vs[3: (n.V-1)]))
	if(is.null(mtry)){mtry<-length(col.split.var)}
	# EXTRACT THE COLS FOR CATEGORICAL VARIABLES
	xs <- gsub(" ", "", unlist(strsplit(unlist(strsplit(as.character(form0)[3], split="[|]"))[1], split="[+]")),fixed = TRUE)
	cat.var <- sapply(xs,function(x){is.factor(training[[x]])})
	#cat.var <- grep("factor", xs)
	col.ctg <- col.split.var[cat.var]
	if (length(col.ctg)==0) col.ctg=NULL
	col.ctg.ordinal <- col.ctg  ############ I THINK WE MIGHT NOT NEED THIS ONE. 

  if(any(is.na(training[,c(col.time,col.status,col.id,col.split.var)])) | 
    (ifelse(!is.null(test),any(is.na(test[,c(col.time,col.status,col.id,col.split.var)])),FALSE))){
    print("Note: Data with missing values deleted")
    training<-training[complete.cases(training[,c(col.time,col.status,col.id,col.split.var)]),]
    test<-test[complete.cases(test[,c(col.time,col.status,col.id,col.split.var)]),]
  }
	
	#Will convert factors to numeric.  Save training dataset beforehand
	trainingTemp<-training
  
  if(!is.null(col.ctg.ordinal)){
    temp<-ordinalizeFunc(training,col.time=col.time,col.status=col.status,col.id=col.id,
                         col.ctg=col.ctg.ordinal,min.levels=4,details=details)
    training<-temp$dat
    if(!is.null(test)){test<-ordinalizeFunc(test,col.time=col.time,col.status=col.status,col.id=col.id,
                         col.ctg=col.ctg.ordinal,min.levels=4,details=details)$dat
    }
  }

  # THE TEST SAMPLE METHOD
  # ========================
  if (selection.method == "test.sample"){
    # GROW A LARGE INITIAL TREE
    tree0 <- grow.MST(dat=training, test=test, method=method,
                      col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                      minsplit=minsplit, min.nevents=min.nevents, max.depth=max.depth, mtry=mtry,
                      cont.split=cont.split, delta=delta, nCutPoints=nCutPoints, details=details)
    # PRUNING AND TREE SIZE SELECTION
    prn <- prune.size.testsample(tree0)
    # PRUNING INFORMATION, LOOK FOR THE MAXIMUM Ga.2, Ga.3, Ga.4, OR Ga.log_n
    pruning.info <- tmp <- prn$result

    col.Ga <- 7:10
    for (j in c(4, col.Ga)) tmp[,j] <- as.numeric.factor(tmp[,j])
    
    # PLOT THE G.a WITH DIFFERENT CHOICES OF a
    if (plot.Ga) {
      if (!is.null(filename)){
        if (substr(filename,nchar(filename)-2,nchar(filename))=='pdf'){
          pdf(file=filename, width=7+3*horizontal, height=7+3*(1-horizontal))
        } else {postscript(file=filename, horizontal=horizontal)}
      }
      par(mfrow=c(1, 1), mar=rep(4, 4))   ##################### SET THE PLOTTING PARAMETERS 
      y.min <- min(tmp[, col.Ga]); y.max <- max(tmp[, 7:10])
      x.min <- min(tmp[, 4]); x.max <- max(tmp[, 4])
      plot(c(0, x.max), c(y.min, y.max), type="n", xlab="# of Terminal Nodes", ylab="G.a Values")
      for (j in col.Ga) lines(tmp$size.tmnl, tmp[,j], col=j, lty=1, lwd=2)
      legend(x.min+4, (y.min+y.max)/2, col=col.Ga, lty=1, lwd=2, legend=c("Ga.2", "Ga.3", "Ga.4", "Ga.ln(n)"))
      if (!is.null(filename)) dev.off()
    }
    
    # OBTAIN THE BEST TREE SIZE AND BEST TREE STRUCTURE
    best.tree.size <- best.tree.structure <- as.list(NULL)
    Ga.cols <- c("Ga.2", "Ga.3", "Ga.4", "Ga.log_n")
    for (j in Ga.cols){
      best.size <- tmp$size.tmnl[which.max(tmp[,j])]
      best.tree.size[[j]] <- best.size
      best.tree.structure[[j]] <- obtain.btree(tree0, bsize=best.size)
    }
    
    # BOOTSTRAP METHOD
    # =================
  } else if (selection.method == "bootstrap") {
    # GROW AND PRUNE B BOOTSTRAP TREES
    boot.result <- bootstrap.grow.prune(B=B, data=training, method=method,
                                        col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                                        minsplit=minsplit, min.nevents=min.nevents, max.depth=max.depth, mtry=mtry,
                                        LeBlanc=LeBlanc, min.boot.tree.size=min.boot.tree.size, details=details)
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
        if(tree0$operator[i] == "<="){treeCut[i]<-paste(info[as.numeric(info[,3]) <= as.numeric(as.character(tree0$cut[i])),1],collapse=",")
        } else if (tree0$operator[i] == ">"){treeCut[i]<-paste(info[as.numeric(info[,3]) > as.numeric(as.character(tree0$cut[i])),1],collapse=",")
        } else {stop("Unexpected operator")}
      }
    }
    tree0$cut<-treeCut
    tree0$operator[tree0$var %in% col.ctg]="in"
    
    for(name in names(best.tree.structure)){
      tree<-best.tree.structure[[name]]
      treeCut<-as.character(tree$cut)
      for(i in 1:NROW(tree)){
        if(tree$var[i] %in% col.ctg.ordinal){
          info<-get(as.character(tree$vname[i]), temp$info)
          if(tree$operator[i] == "<="){treeCut[i]<- paste(info[as.numeric(info[,3]) <= as.numeric(as.character(tree$cut[i])),1],collapse=",")
          } else if (tree$operator[i] == ">"){treeCut[i]<- paste(info[as.numeric(info[,3]) > as.numeric(as.character(tree$cut[i])),1],collapse=",")
          } else {stop("Unexpected operator")}
        }
      }
      best.tree.structure[[name]]$cut<-treeCut
      best.tree.structure[[name]]$operator[best.tree.structure[[name]]$var %in% col.ctg]="in"
    }
  }

	tree0<-listIntoTree(tree=tree0, data=trainingTemp, formula=formula)
	for(name in names(best.tree.structure)){
  	best.tree.structure[[name]]<-listIntoTree(tree=best.tree.structure[[name]], data=trainingTemp, formula=formula)
	}
	
  return(list(tree0=tree0, pruning.info=pruning.info, 
              best.tree.size=best.tree.size, best.tree.structure=best.tree.structure))
}
