## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(cache=TRUE)

## ------------------------------------------------------------------------
# Si igraph está instalado, lo carga; si no, lo instala y luego lo carga
if(!require(igraph)) install.packages("igraph"); require(igraph)

if(!require(ggplot2)) install.packages("ggplot2"); require(ggplot2)


## ------------------------------------------------------------------------
download.file("http://moreno.ss.uci.edu/beach.dat", destfile = "windsurfers.dat") 
ws <- read.table("windsurfers.dat", skip = 7)
dim(ws)

## ------------------------------------------------------------------------
ws.obs <- as.matrix(ws[1:43, ])
ws.per <- as.matrix(ws[44:86, ])

## ------------------------------------------------------------------------
ws.obs.red <- graph.adjacency(ws.obs, mode="undirected", diag=FALSE, weighted = T)
plot(ws.obs.red)

## ------------------------------------------------------------------------
hist(ws.per[lower.tri(ws.per)], main = "histograma de las interacciones percibidas")
umbral <- 0.5
ws.per.2 <- ws.per
ws.per.2[which(ws.per.2 <= umbral)] <- 0
ws.per.red <- graph.adjacency(ws.per.2, mode="undirected", diag=FALSE, weighted = T)
plot(ws.per.red)

## ------------------------------------------------------------------------
summary(ws.obs.red)
summary(ws.per.red)

# Para información más detallada
# str(ws.obs.red)

# Si necesito sólo los recuentos
vcount(ws.obs.red)
ecount(ws.obs.red)

# Para ver los nodos y aristas
V(ws.obs.red)
E(ws.obs.red)

V(ws.per.red)
E(ws.per.red)


## ------------------------------------------------------------------------
is.simple(ws.obs.red)
is.simple(ws.per.red)

is.connected(ws.obs.red)
is.connected(ws.per.red)

## ------------------------------------------------------------------------
# Red de interacciones observadas
diameter(ws.obs.red)
get.diameter(ws.obs.red)

# Red de interacciones percibidas
diameter(ws.per.red)
get.diameter(ws.per.red)


## ------------------------------------------------------------------------
graph.density(ws.obs.red)
graph.density(ws.per.red)


## ------------------------------------------------------------------------
head (transitivity(ws.obs.red, type = "local"))
head (transitivity(ws.per.red, type = "local"))

# red de interacciones observadas
transitivity(ws.obs.red, type = "global")
# red de interacciones percibidas
transitivity(ws.per.red, type = "global")


## ------------------------------------------------------------------------
par(mfrow = c(1,2))
hist(transitivity(ws.obs.red, type = "local"), main = "observada", 
     breaks = seq(0.2, 1, 0.1), xlab = "coefs. de clustering")
hist(transitivity(ws.per.red, type = "local"), main = "percibida", 
     breaks = seq(0.2, 1, 0.1), xlab = "coefs. de clustering")
# El display gráfico vuelve a la configuración de un gráfico por panel
par(mfrow = c(1,1))


## ------------------------------------------------------------------------
degree(ws.obs.red)
sort(degree(ws.obs.red), decreasing = T)

## ------------------------------------------------------------------------
qplot(degree(ws.per.red), degree(ws.obs.red))
cor(degree(ws.per.red), degree(ws.obs.red))

## ------------------------------------------------------------------------
head( degree.distribution(ws.obs.red ), 15)
head( degree.distribution(ws.obs.red, cumulative = T ))

par( mfrow = c(1,2) )
plot( degree.distribution(ws.obs.red), 
      xlab = "grados", ylab = "proporción de nodos", type = "h", main = "observadas")

plot( degree.distribution(ws.per.red), 
      xlab = "grados", ylab = "proporción de nodos", type = "h", main = "percibidas")

# El display gráfico vuelve a la configuración de un gráfico por panel
par( mfrow = c(1,2) )


## ------------------------------------------------------------------------
par( mfrow = c(1,2) )

plot( degree.distribution(ws.obs.red, cumulative = T), type = "l", xlab = "grado", ylab = "proporción de nodos con grado > x", main = "observadas")

plot( degree.distribution(ws.per.red, cumulative = T), type = "l", xlab = "grado", ylab = "proporción de nodos con grado > x", main = "percibidas")

# El display gráfico vuelve a la configuración de un gráfico por panel
par( mfrow = c(1,1))


## ------------------------------------------------------------------------
ws.obs.red.plf <- power.law.fit(degree(ws.obs.red))
ws.per.red.plf <- power.law.fit(degree(ws.per.red))

## ------------------------------------------------------------------------
ws.obs.red.plf$alpha 
ws.per.red.plf$alpha 


## ------------------------------------------------------------------------
ws.obs.red.plf$KS.p 
ws.per.red.plf$KS.p 


## ------------------------------------------------------------------------
ws.obs.red.plf$xmin
ws.per.red.plf$xmin 


## ------------------------------------------------------------------------
ws.obs.red.plf.2 <- power.law.fit(degree(ws.obs.red), xmin = 9)
ws.obs.red.plf.2$alpha
ws.obs.red.plf.2$KS.p

ws.per.red.plf.2 <- power.law.fit(degree(ws.per.red), xmin = 17)
# Elegimos 17 porque los grados mínimos en ws.per.red son mayores
ws.per.red.plf.2$alpha
ws.per.red.plf.2$KS.p


## ------------------------------------------------------------------------
rg.transitivity.barabasi <- array()
rg.transitivity.erdos <- array()

for(i in 1:1000){
  rg.1 <- barabasi.game(43, power = 2.4, m=8, directed = F)
  rg.2 <- sample_gnm(43, 336)
  rg.transitivity.barabasi[i] <-  mean(transitivity(rg.1, "local", isolates="zero"))
  rg.transitivity.erdos[i] <- mean(transitivity(rg.2, "local", isolates="zero"))
}
red.transitivity <- mean(transitivity(ws.obs.red, "local", isolates="zero"))

## ------------------------------------------------------------------------
table(red.transitivity > rg.transitivity.barabasi)
hist(rg.transitivity.barabasi, main = "coef. clustering, grafos de Barabasi-Albert")
abline( v = red.transitivity, col ="red", lwd= 2)


## ------------------------------------------------------------------------
table(red.transitivity > rg.transitivity.erdos)
hist(rg.transitivity.erdos, xlim=c(0.34, 0.7), main = "coef. clustering, grafos al azar")
abline( v = red.transitivity, col ="red", lwd= 2)

## ------------------------------------------------------------------------
assortativity.degree(ws.obs.red)
assortativity.degree(ws.per.red)


## ------------------------------------------------------------------------
# Intermediación
head( sort( betweenness(ws.obs.red), decreasing = T) )
head( sort( betweenness(ws.per.red), decreasing = T) )

# Cercania
head( sort( closeness(ws.obs.red), decreasing = T) )
head( sort( closeness(ws.per.red), decreasing = T) )

# Centralidad de autovectores
head( sort( eigen_centrality(ws.obs.red)$vector, decreasing = T) )
head( sort( eigen_centrality(ws.per.red)$vector, decreasing = T) )


## ------------------------------------------------------------------------
ws.obs.red.cl.eb <- cluster_edge_betweenness(ws.obs.red, directed = F, merges = T)
plot(ws.obs.red, vertex.color = ws.obs.red.cl.eb$membership)

ws.per.red.cl.eb <- cluster_edge_betweenness(ws.per.red, directed = F, merges = T)
plot(ws.per.red, vertex.color = ws.per.red.cl.eb$membership)

table(ws.per.red.cl.eb$membership, ws.obs.red.cl.eb$membership)
plot(ws.obs.red, vertex.color = ws.per.red.cl.eb$membership)

plot(ws.obs.red, vertex.color = ws.obs.red.cl.eb$membership, vertex.shape = ifelse(ws.per.red.cl.eb$membership == 1, "circle", "square") )


## ------------------------------------------------------------------------
modularity(ws.obs.red, ws.per.red.cl.eb$membership)

## ------------------------------------------------------------------------
# Un array de modularidades al azar
random.membership <- array()

# Valores de modularidad para 1000 agrupamientos al azar
for(i in 1:1000) random.membership[i] <- modularity( ws.obs.red, sample(1:2, 43, replace = T) )

# test del clustering basado en las modularidades
table( modularity(ws.obs.red, ws.per.red.cl.eb$membership) > random.membership )

## ------------------------------------------------------------------------
# Comparar comunidades

adjusted.rand <- compare(ws.obs.red.cl.eb,ws.per.red.cl.eb, method = "adjusted.rand")

random.rand <- array()
tmp.comm <- ws.per.red.cl.eb
for(i in 1:1000) {
  tmp.comm$membership <- as.numeric(sample(1:2, 43, replace = T))
  random.rand[i] <- compare( ws.obs.red.cl.eb, tmp.comm, method = "adjusted.rand" )
}

table( adjusted.rand > random.rand )






