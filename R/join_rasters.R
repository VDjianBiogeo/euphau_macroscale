#' @description
#' This function joins a raster stack of species distribution with an environmental
#' raster stack.
#' @param sp_stack A raster stack of species distribution.
#' @param env A raster stack of environmental values.
#' 
#' @export
#' 

library(terra)
library(tidyverse)

join_rasters <- function(sp_stack, env){
  
  # env <- terra::rast(env, type = "xyz")
  # terra::crs(env) <- "epsg:4326"
  # 
  # stack <- c(sp_stack, env)
  
  stack <- inner_join(sp_stack,
                      env[,-grep(pattern = "northv|eastv|po4|no3", 
                                      colnames(env))],
                      by = c("x", "y")) |>
    na.omit() %>%
    as.data.frame()
  
  return(stack)
}