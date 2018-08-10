library(igraph)
library(neat)
library(mclust)
source('A-dcblockmodel-functions.R')
source("B-one-tailed-neat.R")
source("C-compute-CSV.R")
source("D-shuffle-communities.R")
source("E-clustering-methods.R")

case = 1 # from 1 to 6
p_out = c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10)[case]

set.seed(0202)

v = 1000 # number of vertices
p = 8 # number of communities

group_assignment = sample(1:p, v, T)
group_sizes = rep(NA, p)
group_list = vector('list', p)
names(group_list) = 1:p
for (i in 1:p) {
  group_list[[i]] = which(group_assignment == i)
  group_sizes[i] = length(group_list[[i]])
}

blockprob1 = matrix(p_out, p, p)
diag(blockprob1) = 0.3

nrep = 100
truemodul = rep(NA, nrep)
ucsv_bh = matrix(NA, nrep, 4)
wcsv_bh = matrix(NA, nrep, 4)
modul = matrix(NA, nrep, 4)
ari = matrix(NA, nrep, 4)

for (i in 1:nrep) {
  degree_weights = weights_dcblockm(group_assignment)
  # graph generation
  graph = dcblockmodel(group_assignment, degree_weights, 
                               blockprob1, selfloops = F)
  # compute modularity
  truemodul[i] = modularity(graph.adjacency(graph$adj,
                  'undirected'), group_assignment)
  # perform clustering and corresponding modularity
  clu = clustering_methods(graph$adj)
  modul[i,] = clu$modularity
  # compute CSV index for each clustering method
  for (j in 1:4) {
    ari[i,j] = adjustedRandIndex(group_assignment, clu$membership[[j]])
    if (length(clu$communities[[j]]) == 1) index[i,j] = NA
    else {
      indexbh = compute_indexes(clu$communities[[j]], 
                                  sizes = clu$sizes[[j]],
                                  adj = graph$adj, v, alpha=0.05, 
                                  correction = 'BH')
      ucsv_bh[i,j] = indexbh$unw_index
      wcsv_bh[i,j] = indexbh$w_index
    }
  }
  if (i%%10 == 0) cat(i,'')
}

output.file = paste('results-sim-3-', case, '.RData', sep = '')

save.image(output.file)