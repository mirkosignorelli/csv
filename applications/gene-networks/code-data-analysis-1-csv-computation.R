if (on.shark) load('data/adjlist-Gambardella.RData')
if (!on.shark) load('Gambardella_v1_2018/adjlist-Gambardella.RData')
# all adj matrices have the same size:
sapply(adjlist, dim)

library(igraph)
library(Matrix)
source('functions/B-one-tailed-neat.R')
source('functions/C-compute-CSV.R')

CSV_two_graphs = function(adj1, adj2, mtc = 'Heyse') {
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
    setlist1[[i]] = which(membership(cl1) == i)
  }
  names(setlist1) = paste('community', 1:nc1)
  nc2 = max(cl2$membership)
  setlist2 = vector('list', nc2)
  for (i in 1:nc2) {
    setlist2[[i]] = which(membership(cl2) == i)
  }
  names(setlist2) = paste('community', 1:nc2)
  csv_within_g1 = compute_csv(communities = setlist1,
                     adj = adj1, v = length(keepboth), 
                     alpha = 0.05, mtc = mtc)
  csv_within_g2 = compute_csv(communities = setlist2,
                     adj = adj2, v = length(keepboth), 
                     alpha = 0.05, mtc = mtc)
  csv_c1_g2 = compute_csv(communities = setlist1,
                     adj = graph2, v = length(keepboth),
                     alpha = 0.05, mtc = mtc)
  csv_c2_g1 = compute_csv(communities = setlist2,
                     adj = graph1, v = length(keepboth),
                     alpha = 0.05, mtc = mtc)
  out = list('csv_c1_g2' = csv_c1_g2, 
             'csv_c2_g1' = csv_c2_g1, 
             'csv_within_g1' = csv_within_g1, 
             'csv_within_g2' = csv_within_g2)
  return(out)
}

absolute_csv = matrix(NA, 30, 30)
relative_csv = matrix(NA, 30, 30)

for (i in 1:29) {
  for (j in (i+1):30) {
    temp <- CSV_two_graphs(adjlist[[i]], adjlist[[j]], mtc = 'BH')
    absolute_csv[i,j] = temp$csv_c1_g2$unw_index
    absolute_csv[j,i] = temp$csv_c2_g1$unw_index
    relative_csv[i,j] = temp$csv_c1_g2$unw_index/temp$csv_within_g1$unw_index
    relative_csv[j,i] = temp$csv_c2_g1$unw_index/temp$csv_within_g2$unw_index
    cat(j)
  }
  print(paste('i =', i, 'done'))
}

save(relative_csv, absolute_csv,
  file = 'results/1_CSV_values.RData')
