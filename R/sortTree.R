sortTree <-
function(tree, data, col.surv, col.id, col.ctg, col.ctg.ordinal) {
  if (nrow(tree) == 1) { return(cbind(tree, operator = NA)[ , c(1, 2, 3, 4, 9, 5, 6, 7, 8)]) }
  leq <- rep(NA, NROW(tree))
  nodesTemp <- tree$node
  #Initialize for root node
  if (any(tree$var[1] %in% col.ctg) & !any(tree$var[1] %in% col.ctg.ordinal)) { split <- (data[ , tree$var[1]] %in% as.character(tree$cut[1]))
  } else { split <- (data[,tree$var[1]] <= as.numeric(as.character(tree$cut[1]))) }

  fitCoxph<-suppressWarnings(coef(coxph(data[ , col.surv] ~ split + cluster(data[ , col.id]))))
  if (!is.na(fitCoxph)) { leq[1] <- (fitCoxph < 0)
  } else {
    #If previous coxph is NA, then try sorting by coxph without cluster()
    fitCoxph <- suppressWarnings(coef(coxph(data[ , col.surv] ~ split)))
    if (!is.na(fitCoxph)) { leq[1] <- (fitCoxph < 0)
    #If all coxph is NA, then assign left split as <=
    } else { leq[1] <- TRUE }
  }

  for (j in 2:NROW(tree)) {
    if (!is.na(tree$score.test[j])) {
      if (any(tree$var[j] %in% col.ctg) & !any(tree$var[j] %in% col.ctg.ordinal)) { split <- (data[ , tree$var[j]] %in% as.character(tree$cut[j]))
      } else { split <- (data[ , tree$var[j]] <= as.numeric(as.character(tree$cut[j]))) }
      parentNode <- rep(TRUE, NROW(data))
      for (i in 1:(nchar(tree$node[j]) - 1)) {
        #Direction of split depends on next node (1 or 2)
        lastNode <- substr(tree$node[j], nchar(tree$node[j]) - i + 1, nchar(tree$node[j]) - i + 1)
        varTemp <- tree$var[substr(tree$node[j], 1, nchar(tree$node[j]) - i) == tree$node]
        cutTemp <- tree$cut[substr(tree$node[j], 1, nchar(tree$node[j]) - i) == tree$node]
        if (lastNode == 1) {
          if (any(varTemp %in% col.ctg) & !any(varTemp %in% col.ctg.ordinal)) {
            parentNode[!(data[ , varTemp] %in% as.character(cutTemp))] <- FALSE
          } else { parentNode[data[,varTemp] > as.numeric(as.character(cutTemp)) + 1e-8] <- FALSE}
        } else {
          if (any(varTemp %in% col.ctg) & !any(varTemp %in% col.ctg.ordinal)) {
            parentNode[data[ , varTemp] %in% as.character(cutTemp)] <- FALSE
          } else { parentNode[data[,varTemp] <= as.numeric(as.character(cutTemp)) + 1e-8] <- FALSE}
        }
      }
      #if(sum(parentNode) != tree$size[j]){print(j)}
      #print(c(tree$size[j],sum(parentNode)))
      #First try sorting by coxph with cluster()
      fitCoxph <- suppressWarnings(coef(coxph(data[ , col.surv][parentNode] ~ split[parentNode] + cluster(data[ , col.id][parentNode]))))
      if (!is.na(fitCoxph)) { leq[j] <- (fitCoxph < 0)
      } else {
        #If previous coxph is NA, then try sorting by coxph without cluster()
        fitCoxph <- suppressWarnings(coef(coxph(data[ , col.surv][parentNode] ~ split[parentNode])))
        if (!is.na(fitCoxph)) { leq[j] <- (fitCoxph < 0)
        #If all coxph is NA, then assign left split as <=
        } else { leq[j] <- TRUE }
      }
    }
    lastTempNode <- nodesTemp[substr(tree$node[j], 1, nchar(tree$node[j]) - 1) == tree$node]
    endNode <- substr(nodesTemp[j], nchar(tree$node[j]), nchar(tree$node[j]))
    lastLEQ <- ifelse(leq[substr(tree$node[j], 1, nchar(tree$node[j]) - 1) == tree$node], endNode, 3 - as.numeric(endNode))
    nodesTemp[j] <- paste0(lastTempNode, lastLEQ)
  }
  operator <- ifelse(leq, "<=", ">")
  #operator[tree$var %in% col.ctg]=ifelse(leq[tree$var %in% col.ctg],"in","not in")
  tree$operator <- operator
  tree$node <- nodesTemp
  tree <- tree[order(tree$node), c(1, 2, 3, 4, 9, 5, 6, 7, 8)]
  return(tree)
}
