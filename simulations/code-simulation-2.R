library(igraph)
library(neat)
source('A-dcblockmodel-functions.R')
source("B-one-tailed-neat.R")
source("C-compute-CSV.R")
source("D-shuffle-communities.R")

case = 1 # from 1 to 6
p_out = c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10)[case]

set.seed(1002)

v = 1000 # number of vertices
p = 8 # number of communities

# DEFINITION OF THE REFERENCE GRAPH:
group_assignment = sample(1:p, v, T)
degree_weights = weights_dcblockm(group_assignment)
blockprob = matrix(p_out, p, p)
diag(blockprob) = 0.3

reference_graph = dcblockmodel(group_assignment, degree_weights, 
                               blockprob, selfloops = F)

cases = 21; trials = 100
ucsv_bh = matrix(NA, trials, cases)
wcsv_bh = matrix(NA, trials, cases)

reference_setlist = vector('list', p)
names(reference_setlist) = paste('group',1:p)
for (i in 1:v) reference_setlist[[group_assignment[i]]] = c(reference_setlist[[group_assignment[i]]], i)
group_sizes = rep(NA, p)
for (i in 1:p) {
  group_sizes[i] = length(reference_setlist[[i]])
}

for (k in 1:cases) {
  # reshuffle 0, 5, 10, ..., 100% of nodes' communities
  proportion = 0.05*(k-1)
  for (i in 1:trials) {
    new_communities = shuffle_groups(group_assignment, proportion)
    new_degree_weights = weights_dcblockm(new_communities)
    # graph generation:
    newgraph = dcblockmodel(new_communities, new_degree_weights, 
                            blockprob, selfloops = F)
    # computation of the index
    indexbh = compute_indexes(reference_setlist, 
                              sizes = group_sizes,
                              adj = newgraph$adj, v, alpha=0.05,
                              correction = 'BH')
    ucsv_bh[i,k] = indexbh$unw_index
    wcsv_bh[i,k] = indexbh$w_index
  }
  cat(paste(k,''))
}

output.file = paste('results-sim-2-', case, '.RData', sep = '')

save.image(output.file)