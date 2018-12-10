clustering_methods = function(adj) {
  graph = graph.adjacency(adj, 'undirected')
  membership = vector('list', 4)
  communities = vector('list', 4)
  sizes = vector('list', 4)
  modularity = rep(NA, 4)
  fg = cluster_fast_greedy(graph)
  le = cluster_leading_eigen(graph)
  lo = cluster_louvain(graph)
  wt = cluster_walktrap(graph)
  sizes[[1]] = sizes(fg)
  sizes[[2]] = sizes(le)
  sizes[[3]] = sizes(lo)
  sizes[[4]] = sizes(wt)
  membership[[1]] = fg$membership
  membership[[2]] = le$membership
  membership[[3]] = lo$membership
  membership[[4]] = wt$membership
  communities[[1]] = communities(fg)
  communities[[2]] = communities(le)
  communities[[3]] = communities(lo)
  communities[[4]] = communities(wt)
  modularity[1] = modularity(graph, membership[[1]])
  modularity[2] = modularity(graph, membership[[2]])
  modularity[3] = modularity(graph, membership[[3]])
  modularity[4] = modularity(graph, membership[[4]])
  return(list('communities'=communities, 'membership'=membership,
              'modularity'=modularity, 'sizes' = sizes))
}

library(compiler)
clustering_methods = cmpfun(clustering_methods)
