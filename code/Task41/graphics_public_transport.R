library(igraph)
library(leaflet)
library(poweRlaw)



setwd('~/Documentos/semester2/Complex Systems/Project/')



##########################

########  TOKYO  #########

##########################

nodes_tokyo <- read.csv('nodes/city_nodes_Tokyo.csv')
edges_tokyo <- read.csv('edges/city_edges_Tokyo.csv')

print(unique(edges_tokyo$mode))


nodes_tok <- nodes_tokyo[,c("nodeID", "nodeLabel", "longitud", "latitud", "mode")]
edges_tok <- data.frame(from=edges_tokyo$nodeID_from, to=edges_tokyo$nodeID_to)

edges_tok<-na.omit(edges_tok)

g <- graph_from_data_frame(edges_tok, directed = FALSE, vertices = nodes_tok)
g <- simplify(g)
  
layout <- cbind(nodes_tok$longitud, nodes_tok$latitud)
  
plot(g, layout = layout, 
       vertex.size = 3,
       vertex.color = 'blue',
       vertex.frame.color = 'red4', 
       edge.color = "gray80",
       edge.width = 0.5,
       vertex.label = NA) 

  

map_tok <- leaflet(data = nodes_tok) %>%
    addTiles() %>%
    addCircleMarkers(
      lng = ~longitud, lat = ~latitud,
      label = ~nodeLabel,
      radius = 4,
      color = "blue",
      stroke = FALSE, fillOpacity = 0.7
    )

map_tok


##########################

########  CDMX  ##########

##########################


nodes_cdmx <- read.csv('nodes/city_nodes_Mexico City.csv')
edges_cdmx <- read.csv('edges/city_edges_Mexico City.csv')


print(unique(edges_cdmx$mode))

nodes_mex <- nodes_cdmx[,c("nodeID", "nodeLabel", "longitud", "latitud", "mode")]
edges_mex <- data.frame(from=edges_cdmx$nodeID_from, to=edges_cdmx$nodeID_to)



g <- graph_from_data_frame(edges_mex, directed = FALSE, vertices = nodes_mex)
g <- simplify(g)
  
layout <- cbind(nodes_mex$longitud, nodes_mex$latitud)
  
plot(g, layout = layout, 
     vertex.size = 3,
     vertex.color = 'blue',
     vertex.frame.color = 'red4', 
     edge.color = "gray80",
     edge.width = 0.5,
     vertex.label = NA) 



map_mex <- leaflet(data = nodes_mex) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~longitud, lat = ~latitud,
    label = ~nodeLabel,
    radius = 4,
    color = "blue",
    stroke = FALSE, fillOpacity = 0.7
  )

map_mex



##########################

########  ROME  ##########

##########################

nodes_rome <- read.csv('nodes/city_nodes_Rome.csv')
edges_rome <- read.csv('edges/city_edges_Rome.csv')

print(unique(edges_rome$mode))

nodes_ita <- nodes_rome[,c("nodeID", "nodeLabel", "longitud", "latitud", "mode")]
edges_ita <- data.frame(from=edges_rome$nodeID_from, to=edges_rome$nodeID_to)



g <- graph_from_data_frame(edges_ita, directed = FALSE, vertices = nodes_ita)
g <- simplify(g)

layout <- cbind(nodes_ita$longitud, nodes_ita$latitud)

plot(g, layout = layout, 
     vertex.size = 3,
     vertex.color = 'blue',
     vertex.frame.color = 'red4', 
     edge.color = "gray80",
     edge.width = 0.5,
     vertex.label = NA) 



map_rome <- leaflet(data = nodes_ita) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~longitud, lat = ~latitud,
    label = ~nodeLabel,
    radius = 4,
    color = "blue",
    stroke = FALSE, fillOpacity = 0.7
  )

map_rome



nodes_rome <- read.csv('nodes_trams/city_nodes_Rome.csv')
edges_rome <- read.csv('edges_trams/city_edges_Rome.csv')

print(unique(edges_rome$mode))

nodes_ita <- nodes_rome[,c("nodeID", "nodeLabel", "longitud", "latitud", "mode")]
edges_ita <- data.frame(from=edges_rome$nodeID_from, to=edges_rome$nodeID_to)



g <- graph_from_data_frame(edges_ita, directed = FALSE, vertices = nodes_ita)
g <- simplify(g)

layout <- cbind(nodes_ita$longitud, nodes_ita$latitud)

plot(g, layout = layout, 
     vertex.size = 3,
     vertex.color = 'blue',
     vertex.frame.color = 'red4', 
     edge.color = "gray80",
     edge.width = 0.5,
     vertex.label = NA) 



#We can create a table that contains the topological information for each city. From this we can study, 
#in general the topological aspects of transport networks as it can be the degree, average path length, mixing states, etc.

cities <- as.data.frame(read.csv('Data /cities.csv'))   

cities_ids <- unique(cities$id)

cities_stats <- list()


for (cid in cities_ids) {
  city_name <- cities$name[cities$id == cid]
  
  nodes_path <- paste0("nodes/city_nodes_", city_name, ".csv")
  edges_path <- paste0("edges/city_edges_", city_name, ".csv")
  
  if (!file.exists(nodes_path) || !file.exists(edges_path)) next
  
  nodes_df <- read.csv(nodes_path)
  edges_df <- read.csv(edges_path)
  
  
  nodes <- nodes_df[,c("nodeID", "nodeLabel", "longitud", "latitud", "mode")]
  edges <- data.frame(from=edges_df$nodeID_from, to=edges_df$nodeID_to)
  edges<-na.omit(edges)
  
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  comps <- components(g)
  giant <- induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
  
  
  city_stats <- tibble(
    city_id = cid,
    city_name = city_name,
    n_nodes = vcount(giant),
    n_edges = ecount(giant),
    avg_degree = mean(degree(giant)),
    clustering = transitivity(giant, type = "average"),
    assortativity = assortativity_degree(giant, directed = FALSE),
    avg_path_len = mean_distance(giant, unconnected = TRUE),
    diameter = diameter(giant, unconnected = TRUE),
    density = edge_density(giant)
  )
  
  cities_stats[[city_name]] <- city_stats
}


all_stats_df <- bind_rows(cities_stats)


mean(na.omit(all_stats_df$density))

mean(all_stats_df$avg_degree)

hist(all_stats_df$avg_degree, main = 'Histogram of average degree for all the cities', xlab = expression('<k>'), ylab = 'Frequency', col = 'navy', density = 40)
hist(all_stats_df$assortativity, main = 'Histogram of assortative values for all the cities', xlab = 'Assortative coefficient', ylab = 'Frequency', col = 'navy', density = 40)





