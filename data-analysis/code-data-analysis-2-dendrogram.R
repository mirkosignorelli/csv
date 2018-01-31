
load('results.RData')

rownames(absolute_csv) = names(adjlist)
colnames(absolute_csv) = names(adjlist)
rownames(relative_csv) = names(adjlist)
colnames(relative_csv) = names(adjlist)

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

pdf('dendrogram.pdf', width = 9)
#plot(dend, horiz = T)
circlize_dendrogram(dend, cex = 0.8, labels_track_height = 0.3)
dev.off()

