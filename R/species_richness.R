#'@description
#'This function takes a spatial occurrence dataframe and computes species 
#'richness for each geographic site.
#'
#'@param data a list of occurrence dataframes.
#'@param list_sp several species lists to keep for analysis.
#'@export

species_richness <- function(data, list_sp){
  library(ggplot2)
  
  for(i in 1:length(data)){
    df <- data[[i]]
    sp_to_keep <- list_sp[[i]]
    
    df$richness <- rowSums(x = df[, sp_to_keep], na.rm = T)
    
    # Plotting species richness
    data("World", package = "tmap")
    
    richness_map <- ggplot() +
      geom_sf(data = World) +
      geom_raster(aes(x = x,
                      y = y, 
                      fill = richness),
                  data = df) +
      scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 7,
                                                              name = "Spectral"))) +
      coord_sf(expand = F, 
               crs = 4326) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank())
    
    # Exporting map
    ggsave(filename = here::here("figures",
                                 "explor_figures",
                                 paste0("sp_richn_map_",
                                        i,
                                        ".png")), 
           plot = richness_map,
           dpi = "retina",
           height = 2000,
           width = 3660, 
           units = "px")
  }
  
  return("done")
}