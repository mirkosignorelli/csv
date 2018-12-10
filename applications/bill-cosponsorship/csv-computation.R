rm(list = ls())
library(igraph)
source('B-one-tailed-neat.R')
source('C-compute-CSV.R')

load('clustering-without-metadata.RData')

alpha = 0.05

csv_lou = compute_indexes(communities = communities(lou), 
                  sizes = table(lou$membership), 
                  adj = adjacency, v = 663, alpha = alpha, 
                  correction = 'BH')

csv_wt = compute_indexes(communities = communities(wt), 
                              sizes = table(wt$membership), 
                              adj = adjacency, v = 663, alpha = alpha, 
                              correction = 'BH')

csv_le = compute_indexes(communities = communities(le), 
                         sizes = table(wt$membership), 
                         adj = adjacency, v = 663, alpha = alpha, 
                         correction = 'BH')
csv_fg = compute_indexes(communities = communities(fg), 
                         sizes = table(fg$membership), 
                         adj = adjacency, v = 663, alpha = alpha, 
                         correction = 'BH')

load('clustering-with-metadata.RData')

names(comm.party.3) = names(comm.gender.3) = LETTERS[1:3]
names(comm.party.4) = names(comm.gender.4) = LETTERS[1:4]

csv_newman_party_3 = compute_indexes(communities = comm.party.3, 
                         sizes = table(newman.party.3), 
                         adj = adjacency, v = 663, alpha = alpha, 
                         correction = 'BH')

csv_newman_gender_3 = compute_indexes(communities = comm.gender.3, 
                                   sizes = table(newman.gender.3), 
                                   adj = adjacency, v = 663, alpha = alpha, 
                                   correction = 'BH')

csv_newman_party_4 = compute_indexes(communities = comm.party.4, 
                                     sizes = table(newman.party.4), 
                                     adj = adjacency, v = 663, alpha = alpha, 
                                     correction = 'BH')

csv_newman_gender_4 = compute_indexes(communities = comm.gender.4, 
                                   sizes = table(newman.gender.4), 
                                   adj = adjacency, v = 663, alpha = alpha, 
                                   correction = 'BH')

csv_wt
csv_le
csv_lou
csv_newman_party_4
csv_newman_gender_4

csv_fg
csv_newman_party_3
csv_newman_gender_3

