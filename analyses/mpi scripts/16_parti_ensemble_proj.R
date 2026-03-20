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

env_ensemble_proj <- function(sp_name){
  
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
                         sep = "")) == F){
      print(paste("No data for ",
                  comp,
                  sep = ""))
    } else {
      # Reading needed data for modelling ----
      bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                               sp_name,
                               "/",
                               sp_name,
                               "_",
                               comp,
                               "_biomod_format.rds",
                               sep = ""))
      
      # Get ensemble models ----
      myEnsemble <- readRDS(paste("./outputs/mod_outputs/",
                                  sp_name,
                                  "/env_part/",
                                  comp,
                                  "/",
                                  sp_name,
                                  "_",
                                  comp,
                                  "_ensemble_model.rds",
                                  sep = ""))
      
      # Importing environmental data ----
      env <- terra::rast(read.csv(paste("./outputs/pre_mod_outputs/",
                                        sp_name,
                                        "/",
                                        sp_name,
                                        "_selected_env.csv",
                                        sep = "")),
                         type = "xyz")[[colnames(bm_data@data.env.var)]]
      
      # Ensemble projection ----
      myEnsProj <- BIOMOD_EnsembleForecasting(bm.em = myEnsemble,
                                              new.env = env,
                                              proj.name = "ensemble",
                                              models.chosen = "all")
      
      saveRDS(myEnsProj,
              paste("./outputs/mod_outputs/",
                    sp_name,
                    "/env_part/",
                    comp,
                    "/",
                    sp_name,
                    "_",
                    comp,
                    "_ensemble_proj.rds",
                    sep = ""))
      
      # Save projections into a raster format ----
      EMProj <- biomod2::get_predictions(obj = myEnsProj)
      
      # Changing crs
      terra::crs(EMProj) <- "epsg:4326"
      
      terra::writeRaster(EMProj,
                         paste("./outputs/mod_outputs/",
                               sp_name,
                               "/env_part/",
                               comp,
                               "/",
                               sp_name,
                               "_",
                               comp,
                               "_ensemble_rasts.tif",
                               sep = ""),
                         overwrite = T)
    }

  }
  
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  env_ensemble_proj(sp_name)
  NULL
  
}

# script termination
closeCluster(cl)
mpi.quit()