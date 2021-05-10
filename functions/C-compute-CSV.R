compute_csv = function(communities, adj, v, alpha, 
                           mtc = c('Heyse', 'BH', 'bonferroni')) {
  if (mtc == 'Heyse') require(discreteMTP)
  # remove communities with less than 5 nodes:
  keep = which(lengths(communities) >= 5)
  communities = communities[keep]
  # testing procedure:
  test_between = neat_1tailed(alist = communities, 
                              network = adj, nodes = 1:v,
                              nettype = 'undirected', type = 'under') 
  test_within = neat_within(alist = communities, 
                             network = adj, nodes = 1:v,
                             nettype = 'undirected', type = 'over')
  if (mtc == 'Heyse') {
    pvalues = c(test_between$pvalue, test_within$pvalue)
    # compute the cdfs:
    n_betw = nrow(test_between)
    n_within = nrow(test_within)
    pCDF = vector('list', n_betw + n_within)
    k = 1
    for (i in 1:n_betw) {
      # support, from 0 to max:
      supp = test_between$min[i]:test_between$max[i]
      # possible p-values:
      temp = pvalue_1tailed(x = supp, 
                            ss = test_between$ss[i],
                            K = test_between$K[i],
                            N = test_between$N[i],
                            type = 'under')
      # sorted:
      pCDF[[k]] = sort(temp)
      k = k+1
    }
    for (i in 1:n_within) {
      # support, from 0 to max:
      supp = test_within$min[i]:test_within$max[i]
      # possible p-values:
      temp = pvalue_1tailed(x = supp, 
                            ss = test_within$ss[i],
                            K = test_within$K[i],
                            N = test_within$N[i],
                            type = 'over')
      # sorted:
      pCDF[[k]] = sort(temp)
      k = k+1
    }
    # include here computation with p.discrete.adjust
    pcorr = p.discrete.adjust(p = pvalues, pCDF = pCDF,
                              method = 'DBH')
    # method = 'DBH' is the step-up procedure of Heyse (2011)
  }
  if (mtc %in% c('BH', 'bonferroni')) {
    pvalues = c(test_between$pvalue, test_within$pvalue)
    pcorr = p.adjust(pvalues, method = mtc)
  }
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
  return(list('unw_index' = unw_index, 'w_index' = w_index,
              'pcorr' = pcorr))
}

neat_within = function(alist, network, nettype, nodes, type = c('over','under')) {
  q = length(alist)
  out = neat_1tailed(alist = alist[1], blist = alist[1], network, nodes, nettype = nettype, type = type)
  for (i in 2:q) {
    out = rbind(out, neat_1tailed(alist = alist[i], blist = alist[i], network, nodes, nettype = nettype, type = type))
  }
  return(out)
}

cdf_midp = function(n, K, N) {
  
}

library(compiler)
compute_csv = cmpfun(compute_csv)
