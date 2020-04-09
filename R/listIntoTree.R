listIntoTree <-
function(tree, data, formula, weights) {

  nodes <- listIntoParty(tree, data = data)
  Vs <- all.vars(formula)

  fitted <- fitted_node(nodes, data = data)
  response <- data[ , 1]
  ret <- party(nodes, data = data,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = response,
                                   "(weights)" = weights,
                                   check.names = FALSE, stringsAsFactors=FALSE),
               terms = terms(formula))
  return(as.constparty(ret))
}
