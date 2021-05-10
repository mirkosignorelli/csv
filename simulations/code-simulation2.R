library(igraph)
source('functions/A-dcblockmodel-functions.R')
source("functions/B-one-tailed-neat.R")
source("functions/C-compute-CSV.R")
source("functions/D-shuffle-communities.R")

# case = 1 to 6
p_out = c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10)[case]
print(paste('Case ', case, ': p_out = ', p_out, sep = ''))

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

n.subcases = 21; n.repl = 100
ucsv_heyse = matrix(NA, n.repl, n.subcases)
wcsv_heyse = matrix(NA, n.repl, n.subcases)

reference_setlist = vector('list', p)
names(reference_setlist) = paste('group',1:p)
for (i in 1:v) reference_setlist[[group_assignment[i]]] = c(reference_setlist[[group_assignment[i]]], i)
group_sizes = rep(NA, p)
for (i in 1:p) {
  group_sizes[i] = length(reference_setlist[[i]])
}

for (k in 1:n.subcases) {
  # reshuffle 0, 5, 10, ..., 100% of nodes' communities
  proportion = 0.05*(k-1)
  for (i in 1:n.repl) {
    new_communities = shuffle_groups(group_assignment, proportion)
    new_degree_weights = weights_dcblockm(new_communities)
    # graph generation:
    newgraph = dcblockmodel(new_communities, new_degree_weights, 
                            blockprob, selfloops = F)
    # computation of the index
    csv_val = compute_csv(reference_setlist, 
                          adj = newgraph$adj, v, alpha=0.05,
                          mtc = 'Heyse')
    ucsv_heyse[i,k] = csv_val$unw_index
    wcsv_heyse[i,k] = csv_val$w_index
  }
  cat(paste(k,''))
}

output.file = paste('results/results-sim-2-', case, '.RData', sep = '')
save.image(output.file)