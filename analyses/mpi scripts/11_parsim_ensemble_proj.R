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

# Launching clusters ----
cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

parsim_en_proj <- function(sp_name){
  
  # Getting bm_data object ----
  bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                           sp_name,
                           "/",
                           sp_name,
                           "_parsim_biomod_format.rds",
                           sep = ""))
  
  # Get ensemble models ----
  myEnsemble <- readRDS(paste("./outputs/mod_outputs/",
                              sp_name,
                              "/parsim/ensemble/",
                              sp_name,
                              "_parsim_ensemble_model.rds",
                              sep = ""))
  
  # Get environmental data ----
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
                "/parsim/ensemble/",
                sp_name,
                "_ensemble_proj.rds",
                sep = ""))
  
  # Save projections into a raster format ----
  EMProj <- biomod2::get_predictions(obj = myEnsProj)
  
  # Changing crs
  terra::crs(EMProj) <- "epsg:4326"
  
  terra::writeRaster(EMProj,
                     paste("./outputs/mod_outputs/",
                           sp_name,
                           "/parsim/ensemble/",
                           sp_name,
                           "_ensemble_rasts.tif",
                           sep = ""),
                     overwrite = T)
  
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  parsim_en_proj(sp_name)
  
}

# script termination
closeCluster(cl)
mpi.quit()