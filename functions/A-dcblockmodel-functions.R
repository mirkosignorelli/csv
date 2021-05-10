weights_dcblockm = function(group_assignment) {
  v = length(group_assignment)
  p = length(unique(group_assignment))
  weights = numeric(v)
  for (r in 1:p) {
    units = which(group_assignment == r)
    n_r = length(units)
    raw_weights = runif(n_r, min = 0.05, max = 1)
    weights[units] = n_r * raw_weights / sum(raw_weights)
  }
  return(weights)
}

dcblockmodel = function(group_assignment, degree_weights, blockprob, selfloops = F) {
  v = length(group_assignment)
  p = length(unique(group_assignment))
  if (dim(blockprob)[1] != p) warning('Different number of groups in group_assignment and blockprob')
  adj = matrix(NA, v, v)
  for (i in 1:v) {
    for (j in i:v) {
      pij = degree_weights[i]*degree_weights[j]*blockprob[group_assignment[i],group_assignment[j]]
      adj[i,j] = rbinom(1, 1, min(pij, 1))
      adj[j,i] = adj[i,j]
    }
  }
  if (selfloops == F) diag(adj) = 0
  out = list(group_assignment, degree_weights, blockprob, adj)
  names(out) = c('group_assignment', 'degree_weights',
                 'blockprob', 'adj')
  return(out)
}

library(compiler)
weights_dcblockm = cmpfun(weights_dcblockm)
dcblockmodel = cmpfun(dcblockmodel)
