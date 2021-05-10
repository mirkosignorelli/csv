load('results/1_CSV_values.RData')
load('Gambardella_v1_2018/adjlist-Gambardella.RData')

rownames(absolute_csv) = names(adjlist)
colnames(absolute_csv) = names(adjlist)
rownames(relative_csv) = names(adjlist)
colnames(relative_csv) = names(adjlist)

names(adjlist)
# mGland = ghiandola mammaria
# adrenal = ghiandole adrenali, sopra il fegato
# skm = sketal muscle
# lymph = linfonodo
# umbell = cordone ombellicale

group1 = c('cerebellum', 'cerebrum', 'midBrain')
group2 = c('bronchus', 'lung')
group3 = c('colon', 'duodenum', 'intestine')
group4 = c('ovary', 'prostate', 'testis', 'uterus')
group5 = c('kidney', 'liver')

summary(as.vector(absolute_csv), na.rm=T)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.02105 0.33333 0.40000 0.40938 0.50000 0.85714      30 

absolute_csv[group1, group1]
absolute_csv[group2, group2]
absolute_csv[group3, group3]
absolute_csv[group4, group4]
absolute_csv[group5, group5]

summary(as.vector(relative_csv), na.rm=T)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.02105 0.33333 0.40000 0.40938 0.50000 0.85714      30 

relative_csv[group1, group1]
relative_csv[group2, group2]
relative_csv[group3, group3]
relative_csv[group4, group4]
relative_csv[group5, group5]

q = as.vector(relative_csv[!is.na(relative_csv)])
ordered_relative_csv = q[order(q)]; rm(q)

relative_csv_list = matrix(NA, 30*29, 3)
k = 1
for (i in 1:30) {
  for (j in ((1:30)[-i])) {
    relative_csv_list[k,1] = names(adjlist)[i]
    relative_csv_list[k,2] = names(adjlist)[j]
    relative_csv_list[k,3] = relative_csv[i,j]
    k = k+1
  }
}

row_ids = which(relative_csv_list[,3] %in% tail(ordered_relative_csv, 30))

relative_csv_list[row_ids,]

row_ids = which(relative_csv_list[,3] %in% head(ordered_relative_csv, 30))

relative_csv_list[row_ids,]

rownames(relative_csv)[2] = 'adrenal gland'
rownames(relative_csv)[4] = 'bone marrow'
rownames(relative_csv)[16] = 'lymph node'
rownames(relative_csv)[17] = 'mammary gland'
rownames(relative_csv)[25] = 'skeletal muscle'
rownames(relative_csv)[26] = 'brain stem'
rownames(relative_csv)[29] = 'umbell. cord'


symm_relative = (relative_csv + t(relative_csv)) / 2
diag(symm_relative) = 1
isSymmetric(symm_relative)

dist_relative = 1-symm_relative
hc = hclust(as.dist(dist_relative), method = 'complete')

#install.packages('dendextend')
#install.packages('circlize')

library(dendextend)
library(circlize)
dend = as.dendrogram(hc)
dend <- color_branches(dend, h = 0.46)

pdf('figs/2_dendrogram.pdf', width = 9)
plot(dend, horiz = T)
circlize_dendrogram(dend, cex = 0.8, labels_track_height = 0.3)
dev.off()

