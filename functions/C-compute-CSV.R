compute_indexes = function(communities, sizes, adj, v, alpha, 
                           correction = c('bonferroni', 'BH')) {
  keep = which(sizes>=5)
  communities = communities[keep]
  test_between = neat_1tailed(alist = communities, 
                              network = adj, nodes = 1:v,
                              nettype = 'undirected', type = 'under') 
  test_within = neat_within(alist = communities, 
                             network = adj, nodes = 1:v,
                             nettype = 'undirected', type = 'over')
  pvalues = c(test_between$pvalue, test_within$pvalue)
  pcorr = p.adjust(pvalues, method = correction)
  pbetween = pcorr[1:(dim(test_between)[1])]
  pwithin = pcorr[(dim(test_between)[1]+1):length(pcorr)]
  dummyb = pbetween <= alpha
  dummyw = pwithin <= alpha
  numerator1 = sum(dummyb) + sum(dummyw)
  numerator2 = sum(dummyb*(1-pbetween/alpha)) +
    sum(dummyw*(1-pwithin/alpha))
  denominator = dim(test_between)[1]+dim(test_within)[1]
  unw_index = numerator1/denominator
  w_index = numerator2/denominator
  return(list('unw_index'=unw_index, 'w_index'=w_index))
}

neat_within = function(alist, network, nettype, nodes, type = c('over','under')) {
  q = length(alist)
  out = neat_1tailed(alist = alist[1], blist = alist[1], network, nodes, nettype = nettype, type = type)
  for (i in 2:q) {
    out = rbind(out, neat_1tailed(alist = alist[i], blist = alist[i], network, nodes, nettype = nettype, type = type))
  }
  return(out)
}

library(compiler)
compute_indexes = cmpfun(compute_indexes)