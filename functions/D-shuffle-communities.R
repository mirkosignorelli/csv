shuffle_groups = function(groups, proportion) {
  if (proportion < 0 | proportion > 1) {
    warning('Proportion should be in [0,1]')
  }
  v = length(groups)
  nchanges = round(v*proportion, 0)
  if (nchanges == 0) return(groups)
  else if (nchanges > 0) {
    p = length(unique(groups))
    units = sample(1:v, nchanges)
    out = groups
    for (i in 1:nchanges) {
      possible_groups = (1:p)[(1:p)!=groups[units[i]]]
      out[units[i]] = sample(possible_groups, 1)
    }
  }
  return(out)
}

library(compiler)
shuffle_groups = cmpfun(shuffle_groups)