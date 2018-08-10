library(igraph)
library(Matrix)
source('B-one-tailed-neat.R')
source('C-compute-CSV.R')
load('adjlist-Gambardella.RData')

compute_CSV = function(adj1, adj2, mtc = c('BH', 'none')) {
  keep1 = which(rowSums(adj1) != 0)
  keep2 = which(rowSums(adj2) != 0)
  keepboth = intersect(keep1, keep2)
  adj1 = adj1[keepboth, keepboth]
  adj2 = adj2[keepboth, keepboth]
  graph1 = graph.adjacency(adj1, 'undirected')
  graph2 = graph.adjacency(adj2, 'undirected')
  cl1 = cluster_louvain(graph1)
  cl2 = cluster_louvain(graph2)
  nc1 = max(cl1$membership)
  setlist1 = vector('list', nc1)
  for (i in 1:nc1) {
    id = as.numeric(which(membership(cl1) == i))
    setlist1[[i]] = keepboth[id]
  }
  names(setlist1) = paste('community', 1:nc1)
    nc2 = max(cl2$membership)
  setlist2 = vector('list', nc2)
  for (i in 1:nc2) {
    id = as.numeric(which(membership(cl2) == i))
    setlist2[[i]] = keepboth[id]
  }
  names(setlist2) = paste('community', 1:nc2)
  csv_within_g1 = compute_indexes(communities = setlist1,
                                  adj = adj1, nodes = keepboth, alpha = 0.05, 
                                  correction = mtc)
  csv_within_g2 = compute_indexes(communities = setlist2,
                                  adj = adj2, nodes = keepboth, alpha = 0.05, 
                                  correction = mtc)
  csv_c1_g2 = compute_indexes(communities = setlist1,
                              adj = graph2, nodes = keepboth, alpha = 0.05, 
                              correction = mtc)
  csv_c2_g1 = compute_indexes(communities = setlist2,
                              adj = graph1, nodes = keepboth, alpha = 0.05, 
                              correction = mtc)
  out = list('csv_c1_g2' = csv_c1_g2, 'csv_c2_g1' = csv_c2_g1, 
             'csv_within_g1' = csv_within_g1, 'csv_within_g2' = csv_within_g2)
  return(out)
}

absolute_csv = matrix(NA, 30, 30)
relative_csv = matrix(NA, 30, 30)

for (i in 1:29) {
  for (j in (i+1):30) {
    temp <- compute_CSV(adjlist[[i]], adjlist[[j]], mtc = 'BH')
    absolute_csv[i,j] = temp$csv_c1_g2$unw_index
    absolute_csv[j,i] = temp$csv_c2_g1$unw_index
    relative_csv[i,j] = temp$csv_c1_g2$unw_index/temp$csv_within_g1$unw_index
    relative_csv[j,i] = temp$csv_c2_g1$unw_index/temp$csv_within_g2$unw_index
  }
  #print(i)
}

save.image(file = 'results.RData')

