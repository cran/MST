listIntoParty <-
function(tree, data, subset = NULL, id = 1) {

  if (is.null(subset)) { subset <- rep(TRUE, nrow(data)) }

  if (is.na(tree$var[id])) { return(partynode(id = id))
  } else {
    if (tree$operator[id] == "in") {
      #Index is used for categorical variable splits
      indices <- rep(2L, length(levels(factor(data[[tree$var[id]]]))))
      indices[levels(factor(data[[tree$var[id]]])) %in% strsplit(tree$cut[id], split=',')[[1]]] <- 1L
      indices[!(levels(factor(data[[tree$var[id]]])) %in% levels(factor(data[subset, ][[tree$var[id]]])))] <- NA
      sp <- partysplit(varid = tree$var[id],
                       index = indices,
                       info = list(stats = tree$score[id]))
    } else if (tree$operator[id] == "not in") {
      #Index is used for categorical variable splits
      indices <- rep(1L, length(levels(factor(data[[tree$var[id]]]))))
      indices[levels(factor(data[[tree$var[id]]])) %in% strsplit(tree$cut[id],split=',')[[1]]] <- 2L
      indices[!(levels(factor(data[[tree$var[id]]])) %in% levels(factor(data[subset, ][[tree$var[id]]])))] <- NA
      sp <- partysplit(varid = tree$var[id],
                       index = indices,
                       info = list(stats = tree$score[id]))
    } else {
      #Breaks is used for continuous variable splits
      if (is.factor(tree$cut[id])) { breaks <- as.numeric.factor(tree$cut[id])
      } else { breaks <- as.numeric(tree$cut[id]) }
      indexTemp <- sort(1:2, decreasing = (tree$operator[id] == ">"))
      sp <- partysplit(varid = tree$var[id],
                       breaks = breaks,
                       index = indexTemp,
                       info = list(stats = tree$score[id]))
    }
  }
  kidids <- kidids_split(sp, data = data)
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))

  for (kidid in 1:length(kids)) {
    s <- subset
    s[kidids != kidid] <- FALSE
    if (kidid > 1) { myid <- max(nodeids(kids[[kidid - 1]]))
    } else { myid <- id }
    # Start recursion on this daugther node
    kids[[kidid]] <- listIntoParty(tree, subset = s, data = data, id = as.integer(myid + 1))
  }

  return(partynode(id = as.integer(id), split = sp, kids = kids,
                   info = list(stats = info_split(sp)$stats)))
}
