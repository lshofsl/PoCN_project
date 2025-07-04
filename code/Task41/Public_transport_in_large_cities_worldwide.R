library(dplyr)
library(stringr)

setwd('~/Documentos/semester2/Complex Systems/Project/')

cities <- as.data.frame(read.csv('Data /cities.csv'))   
lines <- as.data.frame(read.csv('Data /lines.csv'))
stations <- as.data.frame(read.csv('Data /stations.csv'))
trams <- as.data.frame(read.csv('Data /trams.csv'))
lines_stations <- as.data.frame(read.csv('Data /lines_stations.csv'))
lines_trams <- as.data.frame(read.csv('Data /lines_trams.csv'))
means_of_transport <- as.data.frame(read.csv('Data /means_transport.csv'))
#systems <- read.csv('Data /systems.csv') #Dont cotribute 


#We can check if the columns of id in trams and line_trams refer to the same object
length(intersect(trams$tid, lines_trams$id)) / nrow(lines_trams)

#Because is cero, they are 2 completely different objects 

#First we clean the data 
cities <- cities %>% select(-coords, -country_state, -length) %>% rename( city_id = id)
lines <- lines %>% select(-url_name, -color, -system_id) %>% rename(line_id = id, id = transport_mode_id) %>%
  left_join(means_of_transport, by = 'id') %>%rename(mode = name.y)
stations <- stations %>% select(-buildstart)
if (all(stations$closure != 99999)) {
  stations <- stations %>% select(-closure)
}
stations <- stations %>% rename(station_id = id)
trams <- trams %>% select(-buildstart, -opening, -closure) %>% rename( trams_id = id)
lines_stations <- lines_stations %>% select(-deprecated_line_group, -toyear, -fromyear,-created_at)


all_city_ids <- unique(lines$city_id)  #Identify each city 

for (cid in all_city_ids) {

  city_stations <- stations %>% filter(city_id == cid)
  city_lines <- lines %>% filter(city_id == cid)
  city_ls <- lines_stations %>% filter(city_id == cid)
  #city_trams <- trams %>% filter(city_id == cid)
  
  if (nrow(city_stations) == 0 || nrow(city_ls) == 0 ) next
  
  nodes_df <- city_ls %>%
    left_join(city_stations, by = "station_id") %>%
    left_join(city_lines, by = "line_id") %>%
    filter(!is.na(name.x))  
  

  
  # From the geometry column of station we can obtain the long and lat of each station 
  coords_list <- city_stations %>%
    mutate(
      coords = str_remove_all(geometry, "POINT\\(|\\)"), 
      longitud = as.numeric(str_split_fixed(coords, " ", 2)[,1]),
      latitud  = as.numeric(str_split_fixed(coords, " ", 2)[,2])
    )
  
  #New index
  n <- min(nrow(nodes_df), nrow(coords_list))
  
  # We assign the values for each column 
  final_nodes <- nodes_df[1:n, ] %>%
    mutate(
      nodeID = 1:n,
      nodeLabel = name,
      longitud  = coords_list[1:n, "longitud"],
      latitud   = coords_list[1:n, "latitud"],
      mode = mode,
      year = opening
    ) %>%
    select(nodeID, nodeLabel, longitud,
           latitud, mode, year)
  
  # Save as CSV
  city_name <- cities$name[cities$city_id == cid]
  write.csv(final_nodes, paste0("city_nodes_", city_name, ".csv"), row.names = FALSE)
  
  
  #For the edges
  final_nodes_to_edges <- nodes_df[1:n, ] %>%
    mutate(
      nodeID = 1:n,
      nodeLabel = name,
      longitud  = coords_list[1:n, "longitud"],
      latitud   = coords_list[1:n, "latitud"],
      mode = mode,
      year = opening
    ) %>%
    select(nodeID, nodeLabel, longitud, latitud, mode, year, station_id)
  
  station_id_to_nodeID <- final_nodes_to_edges %>%
    distinct(station_id, .keep_all = TRUE) %>%
    select(station_id, nodeID, year)
  
  city_lines_stations <- lines_stations %>% filter(city_id == cid)
  
  #Join  of the stations with the lines
  city_edges_base <- city_lines_stations %>%
    left_join(station_id_to_nodeID, by = "station_id", relationship = "many-to-many") %>%
    left_join(lines, by = c("line_id", "city_id")) %>%
    group_by(line_id) %>%
    mutate(order = row_number()) %>%  #Order of one station with the other (sequentially)
    arrange(line_id, order) %>%
    mutate(
      next_station_id = lead(station_id),
      next_nodeID     = lead(nodeID)
    ) %>%
    ungroup() %>%
    filter(!is.na(next_nodeID))
  
  edges_df <- city_edges_base %>%
    transmute(
      nodeID_from = nodeID,
      nodeID_to   = next_nodeID,
      mode        = mode,   
      line        = name.x,
      year        = year
    )
  
  write.csv(edges_df, paste0("city_edges_", city_name, ".csv"), row.names = FALSE)
  
}
