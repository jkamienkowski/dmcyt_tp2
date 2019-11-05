
# Si igraph está instalado, lo carga; si no, lo instala y luego lo carga
if(!require(igraph)) install.packages("igraph"); require(igraph)
if(!require(ggplot2)) install.packages("ggplot2"); require(ggplot2)

download.file("http://moreno.ss.uci.edu/beach.dat", destfile = "windsurfers.dat")
ws <- read.table("windsurfers.dat", skip = 7)
dim(ws)

ws.obs <- as.matrix(ws[1:43, ])
ws.per <- as.matrix(ws[44:86, ])
# image(ws.obs)
# image(ws.per)

ws.obs.red <- graph.adjacency(ws.obs, mode="undirected", diag=FALSE, weighted = T)
plot(ws.obs.red)

hist(ws.per[lower.tri(ws.per)], main = "histograma de las interacciones percibidas")

umbral <- 0.5
ws.per.2 <- ws.per
ws.per.2[which(ws.per.2 <= umbral)] <- 0
ws.per.red <- graph.adjacency(ws.per.2, mode="undirected", diag=FALSE, weighted = T)
image(ws.per.2)
plot(ws.per.red)

