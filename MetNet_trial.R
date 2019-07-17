# molecular network with MetNet

## load packages
library("MetNet")
library("mass2adduct")

## read intensity matrix
msi <- read.csv("MetaMix_SDHBp_m2a.csv",
                header = FALSE,
                sep = ",")
row.names(msi) <- msi[,1]
msi <- msi[,2:ncol(msi)]
msi <- t(msi)
colnames(msi)[1] <- "mz"
row.names(msi) <- msi[,"mz"]


## create transformations
transformations <- adducts
colnames(transformations) <- c("group", "formula", "mass")

## create structural adjacency matrix
struct_adj <- createStructuralAdjacency(x = msi,
                                        transformation = transformations,
                                        ppm = 10)
## create statistical adjacency matrix
msi_int <- msi[,2:dim(msi)[2]]
stat_adj <- createStatisticalAdjacency(msi_int,
                                       model = c("pearson", "spearman"),
                                       correlation_adjust = "BH")
## comnbine structural and statistical matrices
cons_adj <- combineStructuralStatistical(structure = struct_adj[[1]],
                                         statistical = stat_adj)

g <- igraph::graph_from_adjacency_matrix(cons_adj,
                                         mode = "undirected")
plot(g,
     edge.width=5,
     vertex.label.cex=0.5,
     edge.color="grey")
