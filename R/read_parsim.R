#' @description
#' This function takes a list of species and reads raster files. Each raster is 
#' formated to be in WGS84 and joined in a raster stack.
#' 
#' @param sp_list a list of species names.
#' @param grouped_occ a dataframe containing raw occurrences for mapping.
#' 

library(tidyverse)
library(terra)

read_parsim <- function(sp_list, grouped_occ){
  
  dir.create(here::here("outputs",
                        "post_mod_outputs"))
  
  dir.create(here::here("figures",
                        "post_mod_figures"))
  
  parsim_rasters <- terra::rast(xmin = -180,
                                xmax = 180,
                                ymin = -90,
                                ymax = 90)
  
  for(sp in sp_list){
    
    dir.create(here::here("figures",
                          "post_mod_figures",
                          sp))
    
    dir.create(here::here("outputs",
                          "post_mod_outputs",
                          sp))
    
    # Data treatment ----
    dir.create(here::here("outputs",
                          "post_mod_outputs",
                          sp,
                          "parsim"))
    
    # Importing species raster 
    parsim <- terra::rast(here::here("outputs",
                                     "mod_outputs",
                                     sp,
                                     "parsim",
                                     "ensemble",
                                     paste(sp,
                                           "_ensemble_rasts.tif",
                                           sep = "")))[[1:2]]
    
    # Changing crs
    terra::crs(parsim) <- "epsg:4326"
    
    # Changing layer name
    names(parsim) <- c(paste(sp, "median", sep = "_"), 
                       paste(sp, "cv", sep = "_"))
    
    # Writing species raster file
    terra::writeRaster(x = parsim,
                       filename = here::here("outputs",
                                             "post_mod_outputs",
                                             sp,
                                             "parsim",
                                             paste(sp,
                                                   "_parsim_ensemble_rasts.tif",
                                                   sep = "")),
                       overwrite = T)
    
    # Joining median rasters
    parsim_rasters <- c(parsim_rasters,
                               parsim[[1]])
    
    # Filtering species raw occurrences
    df <- grouped_occ[grouped_occ$species == sp,]
    
    # Turning raster into dataframe
    sp_prob <- as.data.frame(parsim, xy = T)
    
    # Mapping ----
    dir.create(here::here("figures",
                          "post_mod_figures",
                          sp,
                          "parsim"))
    
    ## Plotting species raster world map ----
    data(World, package = "tmap")
    
    world_map <- ggplot() + 
      geom_sf(data = World) +
      geom_raster(aes(x = x, 
                      y = y,
                      fill = sp_prob[,3]/1000),
                  data = sp_prob) +
      scale_fill_viridis_c(name = "Probability of\nocccurence") +
      geom_point(aes(x = x,
                     y = y),
                 data = df,
                 pch = 23, 
                 colour = "black", 
                 fill = "red", 
                 size=1) +
      coord_sf(crs = 4326, expand = F) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank())
    
    ## Saving species world map ----

    ggsave(plot = world_map,
           filename = here::here("figures",
                                 "post_mod_figures",
                                 sp,
                                 "parsim",
                                 paste(sp,
                                       "_parsim_world_median.png",
                                       sep = "")), 
           dpi = "retina", 
           width = 3660,
           height = 1500,
           units = "px")
    
    ## Plotting species subantarctic map ----
    subant_map <- ggplot() + 
      geom_sf(data = World) +
      geom_raster(aes(x = x, 
                      y = y,
                      fill = sp_prob[,3]/1000),
                  data = sp_prob) +
      scale_fill_viridis_c(name = "Probability of\nocccurence") +
      geom_point(aes(x = x,
                     y = y),
                 data = df,
                 pch = 23, 
                 colour = "black", 
                 fill = "red", 
                 size=1) +
      coord_sf(crs = 4326, expand = F, 
               xlim = c(-20,160), ylim = c(-20, -60)) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank())
    
    ## Saving subantarctic map ----
    ggsave(plot = subant_map,
           filename = here::here("figures",
                                 "post_mod_figures",
                                 sp,
                                 "parsim",
                                 paste(sp,
                                       "_parsim_subant_median.png",
                                       sep = "")), 
           dpi = "retina", 
           width = 3660,
           height = 1500,
           units = "px")
  }
  
  # Writing all species raster file
  terra::writeRaster(x = parsim_rasters,
                     filename = here::here("outputs",
                                           "post_mod_outputs",
                                           "euphau_parsim_median.tif"),
                     overwrite = T)
  
  parsim_df <- as.data.frame(parsim_rasters, xy = T)
  return(parsim_df)
}
  