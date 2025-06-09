library(igraph)

motter_lai_degree <- function(network, p, alpha){
  N <- vcount(network)
  V(network)$name <- as.character(1:N)  # Ensure vertex names

  L_i_initial <- betweenness(network)
  C_i_initial <- (1 + alpha) * L_i_initial
  names(C_i_initial) <- V(network)$name

  degr <- degree(network)
  top_nodes <- order(degr, decreasing = TRUE)[1:floor(p * N)]
  net <- delete_vertices(network, top_nodes)

  repeat {
    if (vcount(net) == 0) return(0)

    L_i <- betweenness(net)
    node_names <- V(net)$name

    C_i <- C_i_initial[node_names]

    overloaded <- node_names[L_i >= C_i]

    if (length(overloaded) == 0) break

    net <- delete_vertices(net, overloaded)
  }

  comps <- components(net)
  lcc_size <- if (comps$no > 0) max(comps$csize) / N else 0
  return(lcc_size)
}



motter_lai_random <- function(network, p, alpha){
  N <- vcount(network)
  V(network)$name <- as.character(1:N) 
  
  L_i_initial <- betweenness(network)
  C_i_initial <- (1 + alpha) * L_i_initial
  names(C_i_initial) <- V(network)$name
  
  remove<-round(p*N)

  node_chose<-sample(V(network), remove)
  net<-delete_vertices(network, node_chose)
  
  repeat {
    if (vcount(net) == 0) return(0)
    
    L_i <- betweenness(net)
    node_names <- V(net)$name
    
    C_i <- C_i_initial[node_names]
    
    overloaded <- node_names[L_i >= C_i]
    
    if (length(overloaded) == 0) break
    
    net <- delete_vertices(net, overloaded)
  }
  
  comps <- components(net)
  lcc_size <- if (comps$no > 0) max(comps$csize) / N else 0
  return(lcc_size)
}


motter_lai_btw <- function(network, p, alpha) {
  N <- vcount(network)
  V(network)$name <- as.character(1:N)
  
  L_i_initial <- betweenness(network)
  C_i_initial <- (1 + alpha) * L_i_initial
  names(C_i_initial) <- V(network)$name
  
  remove <- floor(p * N)
  sort_btw <- order(L_i_initial, decreasing = TRUE)
  node_chose <- sort_btw[1:remove]
  
  net <- delete_vertices(network, node_chose)
  
  repeat {
    if (vcount(net) == 0) return(0)
    
    L_i <- betweenness(net)
    node_names <- V(net)$name
    C_i <- C_i_initial[node_names]
    
    overloaded <- node_names[L_i >= C_i]
    
    if (length(overloaded) == 0) break
    
    net <- delete_vertices(net, overloaded)
  }
  
  comps <- components(net)
  lcc_size <- if (comps$no > 0) max(comps$csize) / N else 0
  return(lcc_size)
}




# Initial Parameters
alphas <- 10^seq(-5, 10, length.out = 20)
n_graphs <- 20
N_nodes <- 1000
avg_k <- 5   #RE and BA have the same degree 
p_edge_er <- avg_k/(N_nodes - 1) 
m_edge_ba <- avg_k/2 
p_attack <- 0.01


#################################

######### RANDOM ###############

#################################


#Erdos-Renyi
lcc_matrix_random_er <- matrix(NA, nrow = n_graphs, ncol = length(alphas))

for (i in 1:n_graphs) {
  g <- sample_gnp(N_nodes, p_edge_er)
  lcc_matrix_random_er[i, ] <- sapply(alphas, function(a) {
    motter_lai_random(g, p = p_attack, alpha = a)
  })
}

lcc_means_random_er <- colMeans(lcc_matrix_random_er)


#Barabasi model
lcc_matrix_random_ba <- matrix(NA, nrow = n_graphs, ncol = length(alphas))

for (i in 1:n_graphs) {
  g_ba <- sample_pa(N_nodes, m = m_edge_ba, directed = FALSE)
  lcc_matrix_random_ba[i, ] <- sapply(alphas, function(a) {
    motter_lai_random(g_ba, p = p_attack, alpha = a)
  })
}

lcc_means_random_ba <- colMeans(lcc_matrix_random_ba, na.rm = TRUE)

################################
plot(alphas, lcc_means_random_er, type = "b", log = "x", col = "blue",
     xlab = expression("Tolerance Parameter " * alpha * " (log scale)"), 
     ylab = "Relative Size of LCC",
     main = "Random Attacks on RE and BA Graphs")

lines(alphas, lcc_means_random_ba, type = "b", col = "red")
grid()
legend("bottomleft", legend = c("RE Graph", "BA Graph"),
       col = c("blue", "red"), lty = 1, pch = 1)
max_val_er <- max(lcc_means_random_er)
max_index_er <- which.max(lcc_means_random_er)
max_alpha_er <- alphas[max_index_er]

max_val_ba <- max(lcc_means_random_ba)
max_index_ba <- which.max(lcc_means_random_ba)
max_alpha_ba <- alphas[max_index_ba]
text(max_alpha_er, max_val_er, labels = paste0("Max: ", round(max_val_er, 3)),
     pos = 3, offset = 0.5, col = "blue")
text(max_alpha_ba, max_val_ba, labels = paste0("Max: ", round(max_val_ba, 3)),
     pos = 3, offset = 0.5, col = "red")
################################




#################################

######### DEGREE ###############

#################################
#Erdos-Renyi
lcc_matrix_degree_er <- matrix(NA, nrow = n_graphs, ncol = length(alphas))

for (i in 1:n_graphs) {
  g <- sample_gnp(N_nodes, p_edge_er)
  lcc_matrix_degree_er[i, ] <- sapply(alphas, function(a) {
    motter_lai_degree(g, p = p_attack, alpha = a)
  })
}

lcc_means_degree_er <- colMeans(lcc_matrix_degree_er)


#Barabasi model
lcc_matrix_degree_ba <- matrix(NA, nrow = n_graphs, ncol = length(alphas))

for (i in 1:n_graphs) {
  g_ba <- sample_pa(N_nodes, m = m_edge_ba, directed = FALSE)
  lcc_matrix_degree_ba[i, ] <- sapply(alphas, function(a) {
    motter_lai_degree(g_ba, p = p_attack, alpha = a)
  })
}

lcc_means_degree_ba <- colMeans(lcc_matrix_degree_ba, na.rm = TRUE)

################################
plot(alphas, lcc_means_degree_er, type = "b", log = "x", col = "blue",
     xlab = expression("Tolerance Parameter " * alpha * " (log scale)"), 
     ylab = "Relative Size of LCC",
     main = "Degree-Rank Attacks on RE and BA Graphs")

lines(alphas, lcc_means_degree_ba, type = "b", col = "red")
grid()
legend("bottomleft", legend = c("RE Graph", "BA Graph"),
       col = c("blue", "red"), lty = 1, pch = 1)
max_val_er <- max(lcc_means_degree_er)
max_index_er <- which.max(lcc_means_degree_er)
max_alpha_er <- alphas[max_index_er]

max_val_ba <- max(lcc_means_degree_ba)
max_index_ba <- which.max(lcc_means_degree_ba)
max_alpha_ba <- alphas[max_index_ba]
text(max_alpha_er, max_val_er, labels = paste0("Max: ", round(max_val_er, 3)),
     pos = 3, offset = 0.5, col = "blue")
text(max_alpha_ba, max_val_ba, labels = paste0("Max: ", round(max_val_ba, 3)),
     pos = 3, offset = 0.5, col = "red")
################################



#################################

######### BETWEENESS ############

#################################


#Erdos-Renyi
lcc_matrix_btw_er <- matrix(NA, nrow = n_graphs, ncol = length(alphas))

for (i in 1:n_graphs) {
  g <- sample_gnp(N_nodes, p_edge_er)
  lcc_matrix_btw_er[i, ] <- sapply(alphas, function(a) {
    motter_lai_btw(g, p = p_attack, alpha = a)
  })
}

lcc_means_btw_er <- colMeans(lcc_matrix_btw_er)


#Barabasi model
lcc_matrix_btw_ba <- matrix(NA, nrow = n_graphs, ncol = length(alphas))

for (i in 1:n_graphs) {
  g_ba <- sample_pa(N_nodes, m = m_edge_ba, directed = FALSE)
  lcc_matrix_btw_ba[i, ] <- sapply(alphas, function(a) {
    motter_lai_btw(g_ba, p = p_attack, alpha = a)
  })
}

lcc_means_btw_ba <- colMeans(lcc_matrix_btw_ba, na.rm = TRUE)

################################
plot(alphas, lcc_means_btw_er, type = "b", log = "x", col = "blue",
     xlab = expression("Tolerance Parameter " * alpha * " (log scale)"), 
     ylab = "Relative Size of LCC",
     main = "Load-Rank Attacks on RE and BA Graphs")

lines(alphas, lcc_means_btw_ba, type = "b", col = "red")
grid()
legend("bottomleft", legend = c("RE Graph", "BA Graph"),
       col = c("blue", "red"), lty = 1, pch = 1)
max_val_er <- max(lcc_means_btw_er)
max_index_er <- which.max(lcc_means_btw_er)
max_alpha_er <- alphas[max_index_er]

max_val_ba <- max(lcc_means_btw_ba)
max_index_ba <- which.max(lcc_means_btw_ba)
max_alpha_ba <- alphas[max_index_ba]
text(max_alpha_er, max_val_er, labels = paste0("Max: ", round(max_val_er, 3)),
     pos = 3, offset = 0.5, col = "blue")
text(max_alpha_ba, max_val_ba, labels = paste0("Max: ", round(max_val_ba, 3)),
     pos = 3, offset = 0.5, col = "red")
################################





#################################

######### POWER GRID ############

#################################


USpower_data<-read.table('~/Documentos/semester2/Complex Systems/Session VIII - IX - Robustness and resilience-20250520/USpowergrid.edges', header=F)
USpower_data<-USpower_data[,-3] #Drop last column
USpower_g<- graph_from_edgelist(as.matrix(USpower_data), directed=F)

vcount(USpower_g)
mean(degree(USpower_g))

degs <- degree(USpower_g)
sizes <- sqrt(degs/max(degs))/degs
node_colors <- 'red4'

E(USpower_g)$weight<-1
lay <- layout_with_drl(USpower_g, options = list(liquid.iteration=100, expasion.iterations=100))
plot(USpower_g, layout=lay, vertex.label=NA, vertex.color=node_colors, vertex.frame.color='white', vertex.size=10*sizes, edge.color="gray80", edge.width=0.5*E(USpower_g)$weight/max(E(USpower_g)$weight))

lcc_matrix_random_USpower<- sapply(alphas, function(a) {
  motter_lai_random(USpower_g, p = p_attack, alpha = a)
})

lcc_matrix_degree_USpower<- sapply(alphas, function(a) {
  motter_lai_degree(USpower_g, p = p_attack, alpha = a)
})

lcc_matrix_btw_USpower<- sapply(alphas, function(a) {
  motter_lai_btw(USpower_g, p = p_attack, alpha = a)
})

plot(alphas, lcc_matrix_btw_USpower, type = "b", log = "x", col = "blue",
     xlab = expression("Tolerance Parameter " * alpha * " (log scale)"), 
     ylab = "Relative Size of LCC",
     main = "Us Power Grid Cascade Failure")
max_val_US <- max(lcc_matrix_btw_USpower)
max_index_US <- which.max(lcc_matrix_btw_USpower)
max_alpha_US <- alphas[max_index_US]
lines(alphas, lcc_matrix_degree_USpower, col = 'red', type = 'b')
lines(alphas, lcc_matrix_random_USpower, col='darkgreen', type = 'b')
legend("bottomright", legend = c("Load", "Degree", 'Random'),
       col = c("blue", "red", 'darkgreen'), lty = 1, pch = 1)
grid()



################################


#INTERNET
internet_data<-read.table('~/Documentos/semester2/Complex Systems/Session VIII - IX - Robustness and resilience-20250520/internet_AS_20000102.edges', header=F)

internet_g<- graph_from_edgelist(as.matrix(internet_data), directed=F)
mean(degree(internet_g))
vcount(internet_g)
degs <- degree(internet_g)
sizes <- 5*sqrt(degs/max(degs))/degs
node_colors <- 'red4'

E(internet_g)$weight<-1
lay <- layout_with_drl(internet_g, options = list(liquid.iteration=100, expasion.iterations=100))
plot(internet_g, layout=lay, vertex.label=NA, vertex.color=node_colors, vertex.frame.color='white', vertex.size=10*sizes, edge.color="gray80", edge.width=0.5*E(internet_g)$weight/max(E(internet_g)$weight))


lcc_matrix_random_Internet<- sapply(alphas, function(a) {
  motter_lai_random(internet_g, p = p_attack, alpha = a)
})

lcc_matrix_degree_Internet<- sapply(alphas, function(a) {
  motter_lai_degree(internet_g, p = p_attack, alpha = a)
})

lcc_matrix_btw_Internet<- sapply(alphas, function(a) {
  motter_lai_btw(internet_g, p = p_attack, alpha = a)
})



plot(alphas, lcc_matrix_btw_Internet, type = "b", log = "x", col = "blue2",
     xlab = expression("Tolerance Parameter " * alpha * " (log scale)"), 
     ylab = "Relative Size of LCC",
     ylim = c(0, 0.45),
     main = "Internet AS Network Cascade Failure")
max_val_I <- max(lcc_matrix_btw_Internet)
max_index_I <- which.max(lcc_matrix_btw_Internet)
max_alpha_I <- alphas[max_index_I]
lines(alphas, lcc_matrix_degree_Internet, col = 'red', type = 'b')
lines(alphas, lcc_matrix_random_Internet, col='darkgreen', type = 'b')
legend("bottomright", legend = c("Load", "Degree", 'Random'),
       col = c("blue2", "red", 'darkgreen'), lty = 1, pch = 1)
grid()






