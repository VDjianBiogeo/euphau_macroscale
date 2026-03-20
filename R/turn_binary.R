#' @description
#' This function turns probability rasters into binary rasters.
#' 
#' @param sp_stack A raster stack of species distribution.
#' @param raw_occ A dataframe with raw occurrences data.
#' 
#' @export
#'

library(tidyverse)
library(terra)
turn_binary <- function(sp_stack, raw_occ, sp_list){
  
  binary <- terra::rast(xmin = -180,
                        xmax = 180,
                        ymin = -90,
                        ymax = 90)
  
  for(sp in sp_list){
    # Transforming rasters ----
    
    # Turning raster into a dataframe
    rast_df <- sp_stack[,c("x","y",paste(sp, "median", sep = "_"))]
    rast <- terra::rast(rast_df, type = "xyz")
    
    # Joining raster with raw occurrence data
    joined_df <- inner_join(x = rast_df,
                            y = raw_occ[raw_occ$species == sp,],
                            by = c("x", "y")) %>%
      na.omit() %>%
      as.data.frame()
    
    # Computing 5% quantile
    for_binary <- quantile(joined_df[,3],
                           probs = c(0.05))
    
    # Binary transformation
    binary_proj <- biomod2::bm_BinaryTransformation(data = rast, 
                                                    threshold = for_binary[[1]])

        # Changing CRS
    terra::crs(binary_proj) <- "epsg:4326"
    
    # Changing layer name
    names(binary_proj) <- paste(sp, "binary", sep = "_")
    
    # Writing species binary file
    terra::writeRaster(x = binary_proj,
                       filename = here::here("outputs",
                                             "post_mod_outputs",
                                             sp,
                                             "parsim",
                                             paste(sp,
                                                   "binary_proj.tif",
                                                   sep = "_")),
                       overwrite = T)
    
    # Joining into a raster stack
    binary <- c(binary, binary_proj)
    
    # Mapping ----
    data(World, package = "tmap")
    mapping_df <- as.data.frame(binary_proj, xy = T)
    
    ## Mapping world map
    world_map <- ggplot() +
      geom_sf(data = World) +
      geom_raster(aes(x = x,
                      y = y), 
                  fill = "green",
                  data = mapping_df[mapping_df[,3] == 1,]) +
      geom_point(aes(x = x,
                     y = y),
                 data = raw_occ[raw_occ$species == sp,],
                 pch = 23, 
                 colour = "black", 
                 fill = "red", 
                 size=1) +
      coord_sf(crs = 4326, expand = F) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank())
    
    ## Saving species world map 
    ggsave(plot = world_map,
           filename = here::here("figures",
                                 "post_mod_figures",
                                 sp,
                                 "parsim",
                                 paste(sp,
                                       "_parsim_world_binary.png",
                                       sep = "")), 
           dpi = "retina", 
           width = 3200,
           height = 1500,
           units = "px")
    
    ## Mapping subantarctic map
    subant_map <- ggplot() +
      geom_sf(data = World) +
      geom_raster(aes(x = x,
                      y = y), 
                  fill = "green",
                  data = mapping_df[mapping_df[,3] == 1,]) +
      geom_point(aes(x = x,
                     y = y),
                 data = raw_occ[raw_occ$species == sp,],
                 pch = 23, 
                 colour = "black", 
                 fill = "red", 
                 size=1) +
      coord_sf(crs = 4326, expand = F, 
               xlim = c(-20,160), ylim = c(-20, -60)) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank())
    
    ## Saving species subantarctic map 
    ggsave(plot = subant_map,
           filename = here::here("figures",
                                 "post_mod_figures",
                                 sp,
                                 "parsim",
                                 paste(sp,
                                       "_parsim_subant_binary.png",
                                       sep = "")), 
           dpi = "retina", 
           width = 3660,
           height = 1500,
           units = "px")
  }
  
  # Writing binary stack into a tif file
  terra::writeRaster(x = binary,
                     filename = here::here("outputs",
                                           "post_mod_outputs",
                                           "euphau_parsim_binary.tif"),
                     overwrite = T)
  
  binary <- as.data.frame(binary, xy = T)
  
  return(binary)
  
}