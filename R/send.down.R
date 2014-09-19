send.down <-
function(data, tree, char.var=1000){
  call <- match.call(); out <- match.call(expand.dots = FALSE)
  out$tree <- out$data <- out$... <- NULL
  dat <- cbind(data, node=1); tree <- cbind(tree, n.test=NA)
  cut.point <- as.vector(tree$cut)
  split.var <- as.numeric(as.vector(tree$var))
  for (i in 1:nrow(tree)){
    in.node <- (dat$node)==(tree$node[i])
    tree$n.test[i] <- sum(in.node)
    if (!is.na(split.var[i])){
      # print(cbind(i, var=tree$var[i], cut=tree$cut[i]))
      var.split <- dat[,split.var[i]]
      cut <- cut.point[i]
      if (!is.element(split.var[i], char.var)) { 
        cut1 <- as.numeric(cut)
        l.nd <- dat$node[in.node & var.split <= cut1]
        r.nd <- dat$node[in.node & var.split > cut1]
        dat$node[in.node & var.split <= cut1] <- paste(l.nd, 1, sep="")
        dat$node[in.node & var.split >  cut1] <- paste(r.nd, 2, sep="")  
      }
      else {
        var.split <- as.character(var.split)
        cut1 <- unlist(strsplit(as.character(cut), split=" ")) #####################
        l.nd <- dat$node[in.node & is.element(var.split, cut1)] 
        r.nd <- dat$node[in.node & !is.element(var.split, cut1)]                  
        dat$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 1, sep="")  
        dat$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 2, sep="")}                   
    }}
  # print(tree)
  out$data <- dat
  out$tree <- tree
  out 
}
