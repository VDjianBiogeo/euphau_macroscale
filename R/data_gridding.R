#'@description
#'This function turns a biodiversity dataframe into a regular grid for a 
#'chosen spatial resolution.
#'
#'@return
#'This function returns either one dataframe/raster or a list of dataframes 
#'if several resolutions were chosen.
#'
#'@param db_to_grid a species occurrence dataframe to be gridded resulting from
#'the previous steps.
#'@param resolution Either one or several spatial resolutions.
#'@export

data_gridding <- function(db_to_grid, resolution){
  
  if(length(resolution) == 1){
    
    # Formating into a raster format
    grid <- terra::rast(xmin = -180, 
                        xmax = 180, 
                        ymin = -90, 
                        ymax = 90, 
                        res = resolution)
    
    sp.stack <- grid
    
    for(sp in unique(db_to_grid$Species)){
      # Filtering to keep only one species at a time
      sp.df <- dplyr::filter(db_to_grid[,c("lon",
                                         "lat", 
                                         "Species")], 
                     Species == sp) %>%
        as.data.frame()
      
      # Turning dataframe into raster of matching resolution
      sp.rast <- terra::rasterize(x = as.matrix(sp.df[,1:2]), 
                                  y = grid, 
                                  fun = function (x, ...) 1)
      
      # Changing layer name to keep species identity
      names(sp.rast) <- str_replace_all(string = sp, 
                                        pattern = " ", 
                                        replacement = "_")
      
      # Changing CRS to WGS84
      terra::crs(sp.rast) <- terra::set.crs("epsg:4326")
      
      # Adding layer to raster stack
      sp.stack <- c(sp.stack, sp.rast)
      
      
    }
    
    # Exporting output in a GeoTIFF format
    terra::writeRaster(x = sp.stack, 
                       filename = here::here("data",
                                             "treated",
                                             "bio",
                                             paste("stacked_euphau",
                                                   resolution,
                                                   "_occurrences.tif", 
                                                   sep = "")))
    
    # Same but in a CSV format
    write.csv(x = as.data.frame(sp.stack, xy = T),
              filename = here::here("data",
                                    "treated",
                                    "bio",
                                    paste("stacked_euphau",
                                          resolution,
                                          "_occurrences.csv", 
                                          sep = "")), 
              row.names = F)
    
    return(sp.stack)
    
  } else {
    
    # Creating a list to store all raster stacks
    list_rast <- list()
    
    for(i in 1:length(resolution)){
      
      # Formating into a raster format
      grid <- terra::rast(xmin = -180, 
                          xmax = 180, 
                          ymin = -90, 
                          ymax = 90, 
                          res = resolution[i])
      
      sp.stack <- grid
      
      for(sp in unique(db_to_grid$Species)){
        # Filtering to keep only one species at a time
        sp.df <- dplyr::filter(db_to_grid[,c("lon",
                                                   "lat", 
                                                   "Species")], 
                               Species == sp) %>%
          as.data.frame()
        
        # Turning dataframe into raster of matching resolution
        sp.rast <- terra::rasterize(x = as.matrix(sp.df[,1:2]), 
                                    y = grid, 
                                    fun = function (x, ...) 1)
        
        # Changing layer name to keep species identity
        names(sp.rast) <- str_replace_all(string = sp, 
                                          pattern = " ", 
                                          replacement = "_")
        names(sp.rast) <- paste(str_sub(string = names(sp.rast),
                                        start = 1, 
                                        end = 1), 
                                str_split(names(sp.rast),
                                          pattern = "_",
                                          simplify = T)[,2],
                                sep = ".")
        
        # Changing CRS to WGS84
        terra::crs(sp.rast) <- terra::crs("epsg:4326")
        
        # Adding layer to raster stack
        sp.stack <- c(sp.stack, sp.rast)
        
      }
      
      # Adding spatial resolution to layers' names
      names(sp.stack) <- paste(names(sp.stack),
                                     resolution[i],
                                     sep = "_")
      
      # Export stack for each res as GeoTIFF format
      terra::writeRaster(x = sp.stack,
                         filename = here::here("data",
                                               "treated",
                                               "bio",
                                               paste("euphau",
                                                     resolution[i],
                                                     "stack.tif", 
                                                     sep = "_")),
                         overwrite = T)
      
      # Same but in CSV format
      write.csv(x = as.data.frame(sp.stack, xy = T),
                file = here::here("data",
                                  "treated",
                                  "bio",
                                  paste("euphau",
                                        resolution[i],
                                        "stack.csv", 
                                        sep = "_")),
                row.names = F)
      
      # Adding raster stack to the list
      list_rast[[i]] <- as.data.frame(sp.stack, xy = T)
    }
    
    # Export list to a RDS format
    saveRDS(list_rast, 
            file = here::here("data",
                              "treated",
                              "bio",
                              "euphau_stacks_several_res.rds"))
    
    return(list_rast)
  }
  
}