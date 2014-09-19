MST.plot <-
function(tree, textDepth=3, lines="rectangle"){
  depth<-max(nchar(tree[,1]))
  par(xaxs='i')
  par(mar=c(1,1,1,1))
  par(xpd=TRUE)
  plot(1, type="n", xlab="",ylab="",xlim=c(0,1),ylim=c(0,1), axes=FALSE,xaxs="i",yaxs="i")
  nodes<-tree$node
  nodesBin<-gsub("1", "0", nodes)
  nodesBin<-gsub("2", "1", nodesBin)
  lastObs<-nchar(nodesBin)
  nodesBin<-substr(nodesBin,2,lastObs)
  var <- tree$vname
  cut <- as.character(tree$cut)
  size <- tree$size
  
  for(i in 1:length(nodesBin)){
    nChar<-nchar(nodesBin[i])
    if(!is.na(var[i])){
      if(lines=="rectangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
      } else if(lines=="triangle"){
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+1)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
        
        lines(c((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),
                (1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+3)),
              c((depth-nChar)/(depth+1),(depth-nChar-1)/(depth+1)))
      }         
      
      if(nChar <= textDepth){
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),(depth-nChar)/(depth+1)+1/(depth+20),bquote(.(as.character(var[i]))<=.(cut[i])),cex=1)
      }
    } else {
      if(nChar <= textDepth){
        text((1/(2^(nChar+2)))*(4*strtoi(nodesBin[i], base = 2)+2),(depth-nChar)/(depth+1),paste("N=",size[i],sep=""),cex=1,offset=1)
      }
    }
  }
}
