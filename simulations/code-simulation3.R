library(igraph)
library(mclust)
source('functions/A-dcblockmodel-functions.R')
source("functions/B-one-tailed-neat.R")
source("functions/C-compute-CSV.R")
source("functions/D-shuffle-communities.R")
source("functions/E-clustering-methods.R")

# case = 1 to 6
p_out = c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10)[case]
print(paste('Case ', case, ': p_out = ', p_out, sep = ''))

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

n.repl = 100
truemodul = rep(NA, n.repl)
ucsv_heyse = matrix(NA, n.repl, 4)
wcsv_heyse = matrix(NA, n.repl, 4)
modul = matrix(NA, n.repl, 4)
ari = matrix(NA, n.repl, 4)

for (i in 1:n.repl) {
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
      csv_val = compute_csv(clu$communities[[j]], 
                            adj = graph$adj, v, alpha=0.05, 
                            mtc = 'Heyse')
      ucsv_heyse[i,j] = csv_val$unw_index
      wcsv_heyse[i,j] = csv_val$w_index
    }
  }
  if (i%%10 == 0) cat(i,'')
}

output.file = paste('results/results-sim-3-', case, '.RData', sep = '')

save.image(output.file)