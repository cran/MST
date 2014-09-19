ordinalize <-
function(dat, 		# A data frame
                       col.time, 				# The column number in data for storing the 'time' variable
                       col.status, 			# The column number in data for storing the 'status' variable
                       col.id, 				# The column number in data for storing the 'id' variable
                       cols.cat, 				# The column numbers in data for storing the categorical variables that need to be ordinalized.
                       details=TRUE,				# If true, will print out some details. 	
                       n0=3)					# The minimum number for suggesting level-collaping before 'ordinalization'.
{
  vnames <- colnames(dat)
  time <- dat[, col.time]; status <- dat[,col.status];  id <- dat[, col.id]; 
  n <- NROW(dat)	
  p <- length(cols.cat)
  OUT <- as.list(1:p); names(OUT) <- vnames[cols.cat]
  for (j in 1:p){
    col.cat <- cols.cat[j]
    vname <- vnames[col.cat]
    x <- as.character(dat[, col.cat])
    x.level <- sort(unique(x))
    
    # SHOULD WE MERGE SOME LEVELS FIRST?
    tab <- table(x, status)
    if (min(tab) <= n0) print(paste("You might want to merge some levels of ", x, "first before ordinalization.", sep="")) 
    if (details) {print(vname); print(x.level); print(tab)}
    
    # FITTING MARGINAL MODEL
    options(warn=-1)
    fit <- coxph(Surv(time, status) ~ x + cluster(id), data=dat)
    options(warn=0)
    betas <- c(baseline=0, coef(fit))
    
    # THE INFORMATION REGARDING "ORDINALZIATION"
    OUT[[j]] <- cbind(x.level, betas, rank=rank(betas))
    if (details) print(OUT[[j]])
    
    # REPLACE THE COLUMN IN DATA
    x.level.ordered <- x.level[order(betas)]
    x1 <- ordered(x, levels = x.level.ordered)
    dat[, col.cat] <- as.numeric(x1)
  }
  return(list(dat=dat, info=OUT))
}
