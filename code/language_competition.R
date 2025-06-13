library(deSolve)
library(igraph)
library(ggplot2)



KURAMOTO_L <- function(network, size, dt=5e-2, sd.meas.noise=0.1, sd.dyn.noise=0., sigma=1, a=1.3, s = 0.5){
  if(is.null(igraph::E(network)$weight)) stop("No edge weights specified.")
  
  S <- c()    
  M <- 1    
  idxs <- list()
  
  W <- igraph::get.adjacency(network, attr="weight", sparse=F)
  Nodes <- length(igraph::V(network))
  
  omegas <- rnorm(Nodes, 0, 1/Nodes)
  

  for(m in 1:Nodes){
    #state <- c(X = 0)
    #n-state equivalent:
    S <- c(S, runif(1, 0, 1))
  }
  
  for(m in 1:M){
    idxs[[m]] <- (1:Nodes-1)*M + m
  }
  
  #   MODEL
  #KURAMOTO:  d x/dt = peer_influence + coupling*sum_j A_ij*(x_j-x_i) 
  model <- function(t, S, parameters) {
    dS <- c()
    
    ###################################
    # self term
    ###################################
    
    dS[idxs[[1]]] <- omegas[idxs[[1]]]
    
    ###################################
    # interaction term
    ###################################
    
    a <- parameters$a
    s <- parameters$s
    sigma <- parameters$sigma
    
    #strength <- rowSums(W)
    coupling <- sigma/degree(network)
    #coupling <- sigma/Nodes
    
    x <- S[idxs[[1]]]  
    peer_influence <- W %*% x - rowSums(W) * x
    dS[idxs[[1]]] <- omegas + coupling * peer_influence
    
    
    P_yx <- function(x, s, a){(1 - x)^a * (1 - s)}
    P_xy <- function(x, s, a){x^a * s}
    
  
    nonlinear_shift <- P_yx(x, s, a) - P_xy(x, s, a)
    dS[idxs[[1]]] <- dS[idxs[[1]]] + nonlinear_shift
    
    ###################################
    # dynamic noise term
    ###################################
    
    if(sd.dyn.noise>0){
      dS <- dS + rnorm(length(dS), 0, sd.dyn.noise)
    }
    list(dS)
  }
  
  times <- seq(0, size*dt, by = dt)
  

  parameters <- list(a=a, s=s, sigma=sigma)
  multi <- ode(y = S, times = times, func = model, parms = parameters)

  x <- list()
  for(m in 1:Nodes){
    x[[m]] <- multi[2:nrow(multi), 1 + (m-1)*M + 1]
  }
  
  if(sd.meas.noise>0){
    for(m in 1:Nodes){
      x[[m]] <- x[[m]] + rnorm(size,0,sd.meas.noise)
    }
  }
  
  return(x)
}


plot.MultiTS_L <- function(MultiTS){
  library(ggplot2)
  
  dat <- data.frame()
  
  tseq <- 1:length(MultiTS[[1]])
  
  for(m in 1:length(MultiTS)){
    dat <- rbind(dat, data.frame(node=m, time=tseq, value=MultiTS[[m]]))
  }
  
  return( 
    ggplot(dat, aes(time, sin(value), group=node, color=node)) + theme_bw() + theme(panel.grid=element_blank()) +
      geom_line(alpha=0.3) + 
      scale_color_viridis_c() +
      xlab("Time") + ylab("Relative population of Y")
  )
}




N<-100
M<-100

g<- sample_smallworld(dim = 1, size = N, nei=2, p=0.3)
E(g)$weight <- 1

layout<-layout_in_circle(g)
plot(g, layout = layout,
     vertex.label = NA, 
     vertex.frame.color = NA,
     vertex.size=1,
     vertex.color = 'red', edge.color = "gray80")

MultiTS <- KURAMOTO_L(g, size=M, sd.meas.noise=0., sigma=0.8, a=1.3,s=0.8)
plot.MultiTS_L(MultiTS)

#Barabasi vs small-world
avg_k <- 3
m<- avg_k/2
nei<-  round(avg_k / 2)

g_ba<- sample_pa(N, m)
E(g_ba)$weight <- 1

g_sw<- sample_smallworld(dim = 1, size = N, nei=nei, p=0.3)
E(g_sw)$weight <- 1

MultiTS_ba <- KURAMOTO_L(g_ba, size=M, sd.meas.noise=0., sigma=0.8, a=1.3,s=0.8)
plot.MultiTS_L(MultiTS_ba)

MultiTS_sw <- KURAMOTO_L(g_sw, size=M, sd.meas.noise=0., sigma=0.8, a=1.3,s=0.8)
plot.MultiTS_L(MultiTS_sw)



