listIntoTree <-
function(tree, data, formula){
  
  nodes<-listIntoParty(tree, data=data)
  Vs <- all.vars(formula)
  
  fitted<-fitted_node(nodes, data=data)
  response<-Surv(data[[Vs[1]]],data[[Vs[2]]])
  ret <- party(nodes, data = data,
               fitted = data.frame("(fitted)" = fitted,
                                   "(response)" = response,
                                   check.names = FALSE),
               terms = terms(formula))
  return(as.constparty(ret))
}
