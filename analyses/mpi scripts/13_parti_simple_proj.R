# Setup ----
rm(list = ls())
gc()

set.seed(4951)

library(doMPI)
library(tidyverse)
library(biomod2)
library(terra)

# Importing species list ----
list_sp <- str_split(readRDS("./outputs/explor_outputs/species_lists.rds")[[1]], 
                     pattern = "_",
                     simplify = T)[,1]
list_sp <- list_sp[order(list_sp)]

cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

env_simple_proj <- function(sp_name){
  
  
  # Component class for each variable ----
  variable_classes <- read.csv("./data/treated/env/env_classes.csv",
                               sep = ";") 
  
  components <- unique(variable_classes$env_comp)[2:4]
  
  for(comp in components){
    if(file.exists(paste("./outputs/pre_mod_outputs/",
                         sp_name,
                         "/",
                         sp_name,
                         "_",
                         comp,
                         "_biomod_format.rds",
                         sep = "")) == T){
      
      # Reading needed data for modelling ----
      bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                               sp_name,
                               "/",
                               sp_name,
                               "_",
                               comp,
                               "_biomod_format.rds",
                               sep = ""))
      
      # Importing environmental data ----
      env <- terra::rast(read.csv(paste("./outputs/pre_mod_outputs/",
                                        sp_name,
                                        "/",
                                        sp_name,
                                        "_selected_env.csv",
                                        sep = "")),
                         type = "xyz")[[colnames(bm_data@data.env.var)]]
      
      # Get simple models ----
      sdm_models <- readRDS(paste("./outputs/mod_outputs/",
                                  sp_name,
                                  "/env_part/",
                                  comp,
                                  "/",
                                  sp_name,
                                  "_",
                                  comp,
                                  "_simple_model.rds",
                                  sep = ""))
      
      # Simple projection ----
      mySDMProj <- BIOMOD_Projection(bm.mod = sdm_models,
                                     new.env = env,
                                     proj.name = 'pres',
                                     models.chosen = 'all',
                                     compress = 'gzip',
                                     build.clamping.mask = T,
                                     output.format = '.tif',
                                     do.stack=T)
      
      saveRDS(mySDMProj,
              paste("./outputs/mod_outputs/",
                    sp_name,
                    "/env_part/",
                    comp,
                    "/",
                    sp_name,
                    "_",
                    comp,
                    "_simple_proj.rds",
                    sep = ""))
      
      # Save projections to a raster format ----
      Proj <- biomod2::get_predictions(obj = mySDMProj)
      
      # Changing crs
      terra::crs(Proj) <- "epsg:4326"
      
      terra::writeRaster(x = Proj, 
                         filename = paste("./outputs/mod_outputs/",
                                          sp_name,
                                          "/env_part/",
                                          comp,
                                          "/",
                                          sp_name,
                                          "_",
                                          comp,
                                          "_simple_rast.tif",
                                          sep = ""),
                         overwrite = T)
    } else {
      print(paste("no data available for ", comp, " for ", sp_name, sep = ""))
    }
  }
  
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  env_simple_proj(sp_name)
  NULL
  
}

# script termination
closeCluster(cl)
mpi.quit()