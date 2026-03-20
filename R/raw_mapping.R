#'@description
#'This function maps the spatial distribution of raw occurrence data from a raw
#'database.
#'
#'@return
#'This function returns a map as a ggplot object.
#'
#'@export

raw_mapping <- function(db_to_map){
  
  # Importing background geographic info
  data(World, package = "tmap")
  
  # Mapping
  all_map <- ggplot() +
    geom_sf(data = World) +
    geom_point(aes(x = lon,
                   y = lat,
                   color = source),
               size = 1, 
               data = na.omit(db_to_map[, c("lon", "lat", "source")])) +
    coord_sf(expand = F, crs = 4326, default_crs = sf::st_crs(4326)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_bw(base_size = 18) +
    theme(axis.title = element_blank())
  
  if(dir.exists(here::here("figures",
                           "data_treatment")) == F){
    dir.create(here::here("figures",
                          "data_treatment"))
  }
  ggsave(here::here("figures",
                    "data_treatment",
                    "raw_occurrence_map.png"),
         all_map,
         dpi = "retina",
         width = 3660,
         height = 2940,
         units = "px")
  
  return(all_map)
}