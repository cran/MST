bootstrap.size <-
function(bootstrap.grow.prune.results, plot.it=TRUE, filename=NULL, horizontal=TRUE){
  # EXTRACT NEEDED COMPONENTS
  boot.prune <- bootstrap.grow.prune.results$boot.prune
  tree0 <- bootstrap.grow.prune.results$initial.tree
  n <- tree0$size[1]
  
  OUT <- as.list(NULL)
  #  COMPUTE THE ALPHA PRIME'S
  prune0 <- boot.prune[[1]] 
  n.subtree <- nrow(prune0)
  alpha <- as.numeric(as.vector(prune0$alpha))
  # temp <- c(alpha[1], alpha[-length(alpha)])    	##############
  temp <- c(0, alpha[-length(alpha)])  			############## CHANGE FIRST VALUE OF ALPHA TO 0
  alpha.prime <- sqrt(alpha*temp)
  # cbind(alpha,  alpha.prime=prune0$alpha.prime)
  b <- length(boot.prune)
  G <- as.numeric(as.vector(prune0$G))
  size.tmnl <- as.numeric(as.vector(prune0$size.tmnl))
  subtree <- as.numeric(as.vector(prune0$subtree))
  # tree.penalty <- log(nrow(teeth))
  G.a <- matrix(0, n.subtree, 5)
  penalty <- c(0, 2:4, log(n))
  for (i in 1:n.subtree) {
    a.i <- alpha.prime[i]
    bias <- 0
    for (j in 2:b){
      prune.bs <- boot.prune[[j]]
      alpha.bs <- as.numeric(as.vector(prune.bs$alpha))
      g <- as.numeric(as.vector(prune.bs$G))
      g.test <- as.numeric(as.vector(prune.bs$G.test))
      indx <- 1
      if (sum(alpha.bs <= a.i)>0) {
        temp1 <- which.max(which(alpha.bs<=a.i))
        indx <- ifelse(is.null(temp1), 1, temp1)
      }
      temp2 <- (g-g.test)[indx]
      bias <- bias + temp2
      # print(cbind(i, a.i, j, bias, indx, temp2))
    }
    G.honest <- G[i] - bias/(b-1) 
    G.a[i,] <- G.honest - penalty*(size.tmnl[i]-1)
  }
  out <- data.frame(cbind(size.tmnl, G.a))
  colnames(out) <- c("tmnl.size", "G", "G.2", "G.3", "G.4", "G.log(n)")
  G.a <- out
  
  n.subtrees <- nrow(G.a)
  subtree.size <- G.a[,1]   
  # PLOT THE G.a WITH DIFFERENT CHOICES OF a 	
  if (plot.it) {
    if (!is.null(filename)) postscript(file=filename, horizontal=horizontal)
    par(mfrow=c(1, 1), mar=rep(4, 4))   ##################### SET THE PLOTTING PARAMETERS
    
    min.x <- min(subtree.size); max.x <- max(subtree.size)
    min.y <- min(G.a$"G.log(n)"); max.y <- max(G.a$G.2)
    plot(x=c(min.x, max.x), y=c(min.y, max.y), type="n", xlab="tree size", ylab="G(a)")
    for (j in 3:6) lines(subtree.size, G.a[,j], lty=j-1, col=j-1, lwd=2)
    legend(x=min.x, y=(max.y+min.y)/2, lty=2:5, col=2:5, legend=c("G(2)", "G(3)", "G(4)", "G(ln(n))")) 
    if (!is.null(filename)) dev.off()
  }
  # OBTAIN THE BEST TREE SIZE 
  bsize <- btree <- as.list(NULL)
  Ga.cols <- c("G.2", "G.3", "G.4", "G.log(n)")
  for (j in Ga.cols) {
    best.size <- subtree.size[which.max(G.a[,j])]
    bsize[[j]] <- best.size
    btree[[j]] <- obtain.btree(tree0, bsize=best.size)
  }
  OUT$initial.tree <- tree0; OUT$G.a <- G.a; OUT$bsize <- bsize; OUT$btree <- btree
  return(OUT)
}
