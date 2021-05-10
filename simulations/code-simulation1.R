library(igraph)
library(mclust)
source('functions/A-dcblockmodel-functions.R')
source("functions/B-one-tailed-neat.R")
source("functions/C-compute-CSV.R")

set.seed(0202)

case = 4 # from 1 to 4

if (case == 1) v = 100
if (case == 2) v = 500
if (case == 3) v = 1000
if (case == 4) v = 5000

p = 6 # number of communities

# FIRST SEQUENCE OF GRAPHS:
group_assignment = sample(1:p, v, T)
# group_sizes = rep(NA, p) not needed anymore
group_list = vector('list', p)
names(group_list) = 1:p
for (i in 1:p) {
  group_list[[i]] = which(group_assignment == i)
  #group_sizes[i] = length(group_list[[i]])
}

pout = seq(0,100)*0.3/100
pin = c(0.22, 0.25, 0.28, 0.32, 0.35, 0.38)

ncases = length(pout)
truemodul = rep(NA, ncases)
ucsv_bh = ucsv_heyse = rep(NA, ncases)
wcsv_bh = wcsv_heyse = rep(NA, ncases)

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
  csv.bh = compute_csv(communities = group_list, 
                              adj = graph$adj, v = v,
                              alpha = 0.05, mtc = 'BH')
  csv.heyse = compute_csv(communities = group_list, 
                            adj = graph$adj, v = v,
                            alpha = 0.05, mtc = 'Heyse')
  ucsv_bh[i] = csv.bh$unw_index
  wcsv_bh[i] = csv.bh$w_index
  ucsv_heyse[i] = csv.heyse$unw_index
  wcsv_heyse[i] = csv.heyse$w_index
  if (i%%10 == 0) cat(i,'')
}

output.file = paste('results/res_sim1_', case, '.RData', sep = '')

save(truemodul, v, ucsv_bh, ucsv_heyse, wcsv_bh, wcsv_heyse,
  file = output.file)
