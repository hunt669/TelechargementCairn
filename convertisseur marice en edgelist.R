
g <- read.csv("C:\\Users\\rayanneBENMAHMOUD\\Desktop\\Cairn amyas\\R_test_\\r_test__sirencsv.csv", sep=";")
g <- as.matrix(g)

colnames(g) <- gsub("^X", "", colnames(g))
rownames(g) <- gsub("^X", "", rownames(g))

library(igraph)
g <- graph_from_adjacency_matrix(g)


edgelist <- as_edgelist(g)


head(edgelist)
edgelist_df <- as.data.frame(edgelist)
colnames(edgelist_df) <- c("Source", "Id")


chemin_sortie <- "C:/Users/rayanneBENMAHMOUD/Desktop/edgelist_gephisiren2.csv"
write.csv(edgelist_df, 
          file = chemin_sortie,
          row.names = FALSE,
          quote = FALSE,
          fileEncoding = "UTF-8")

