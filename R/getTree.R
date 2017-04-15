getTree <-
function(mstObj, Ga = c("0", "2", "3", "4", "log_n")) {
  Ga <- match.arg(as.character(substitute(Ga)), c("0", "2", "3", "4", "log_n"))
  if (Ga == 0) {
    if (is.null(mstObj$tree0)) { stop("No initial tree in ", deparse(substitute(mstObj)))
    } else { return(mstObj$tree0) }
  } else {
    if (is.null(mstObj$best.tree.structure)) { stop("No best-sized tree in ", deparse(substitute(mstObj)))
    } else { return(mstObj$best.tree.structure[[paste0("Ga.", Ga)]]) }
  }
}
