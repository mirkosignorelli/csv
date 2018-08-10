rm(list = ls()) 

library(igraph)
library(neat)
library(mclust)
source('A-dcblockmodel-functions.R')
source("B-one-tailed-neat.R")
source("C-compute-CSV.R")

set.seed(0202)

case = 1 # from 1 to 4

if (case == 1) v = 100
if (case == 2) v = 500
if (case == 3) v = 1000
if (case == 4) v = 5000

p = 8 # number of communities

# FIRST SEQUENCE OF GRAPHS:
group_assignment = sample(1:p, v, T)
group_sizes = rep(NA, p)
group_list = vector('list', p)
names(group_list) = 1:p
for (i in 1:p) {
  group_list[[i]] = which(group_assignment == i)
  group_sizes[i] = length(group_list[[i]])
}

pout = seq(0,100)*0.3/100
pin = c(0.25, 0.26, 0.28, 0.30, 0.30, 0.32, 0.34, 0.35)

ncases = length(pout)
truemodul = rep(NA, ncases)
ucsv_bh = rep(NA, ncases)
wcsv_bh = rep(NA, ncases)

for (i in 1:ncases) {
  degree_weights = weights_dcblockm(group_assignment)
  blockprob = matrix(pout[i], p, p)
  diag(blockprob) = pin
  # graph generation
  graph = dcblockmodel(group_assignment, degree_weights, 
                               blockprob, selfloops = F)
  # compute modularity
  truemodul[i] = modularity(graph.adjacency(graph$adj,
                  'undirected'), group_assignment)
  # compute CSV index
  indexbh = compute_indexes(group_list, sizes = group_sizes,
                              adj = graph$adj, v, alpha=0.05, 
                              correction = 'BH')
  ucsv_bh[i] = indexbh$unw_index
  wcsv_bh[i] = indexbh$w_index
  if (i%%10 == 0) cat(i,'')
}

output.file = paste('results-sim-1-', case, '.RData', sep = '')

save.image(output.file)

