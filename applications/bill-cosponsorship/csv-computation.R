library(igraph)
source('functions/B-one-tailed-neat.R')
source('functions/C-compute-CSV.R')

load('results/2-results-no-metadata.RData')

alpha = 0.05

csv_lou = compute_csv(communities = communities(lou), 
                  adj = adjacency, v = 663, alpha = alpha, 
                  mtc = 'Heyse')

csv_wt = compute_csv(communities = communities(wt), 
                              adj = adjacency, v = 663, alpha = alpha, 
                              mtc = 'Heyse')

csv_le = compute_csv(communities = communities(le), 
                         adj = adjacency, v = 663, alpha = alpha, 
                         mtc = 'Heyse')
csv_fg = compute_csv(communities = communities(fg), 
                         adj = adjacency, v = 663, alpha = alpha, 
                         mtc = 'Heyse')

load('results/3-results-metadata.RData')

names(comm.party.3) = names(comm.gender.3) = LETTERS[1:3]
names(comm.party.4) = names(comm.gender.4) = LETTERS[1:4]

csv_newman_party_3 = compute_csv(communities = comm.party.3, 
                         adj = adjacency, v = 663, alpha = alpha, 
                         mtc = 'Heyse')

csv_newman_gender_3 = compute_csv(communities = comm.gender.3, 
                                   adj = adjacency, v = 663, alpha = alpha, 
                                   mtc = 'Heyse')

csv_newman_party_4 = compute_csv(communities = comm.party.4, 
                                     adj = adjacency, v = 663, alpha = alpha, 
                                     mtc = 'Heyse')

csv_newman_gender_4 = compute_csv(communities = comm.gender.4, 
                                   adj = adjacency, v = 663, alpha = alpha, 
                                   mtc = 'Heyse')

# CSV for partitions with 4 clusters
csv_wt
csv_le
csv_lou
csv_newman_party_4
csv_newman_gender_4

# CSV for partitions with 3 clusters
csv_fg
csv_newman_party_3
csv_newman_gender_3

