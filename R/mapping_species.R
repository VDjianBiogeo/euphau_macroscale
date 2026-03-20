
mapping_species <- function(data){
  
  # Mapping species occurrence
  data(World, package = "tmap")
  
  ## Worldwide
  bio.map <- ggplot() +
    geom_sf(data = World) +
    geom_point(aes(x=x, y=y), 
               fill="green", 
               color="black", 
               pch = 21, 
               size = 2,
               data = sp_df[,c(1:3)]) +
    coord_sf(crs = sf::st_crs(4326), expand = F) +
    labs(title = paste(sp, 
                       " world distribution",
                       sep = "")) +
    theme_bw(base_size = 18) +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  return(bio.map)
  
}