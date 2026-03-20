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

simple_projection <- function(sp_name){
  
  # Getting simple models ----
  sdm_models <- readRDS(paste("./outputs/mod_outputs/",
                                  sp_name,
                                  "/simple/",
                                  sp_name,
                                  "_simple_model.rds",
                                  sep = ""))
  
  # Getting environmental data ----
  env <- terra::rast(read.csv(paste("./outputs/pre_mod_outputs/",
                                    sp_name,
                                    "/",
                                    sp_name,
                                    "_selected_env.csv",
                                    sep = "")),
                     type = "xyz")
  
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
                "/simple/",
                sp_name,
                "_simple_proj.rds",
                sep = ""))
  
  # Save projections to a raster format ----
  Proj <- biomod2::get_predictions(obj = mySDMProj)
  
  # Changing crs
  terra::crs(Proj) <- "epsg:4326"
  
  terra::writeRaster(Proj,
                     paste("./outputs/mod_outputs/",
                           sp_name,
                           "/simple/",
                           sp_name,
                           "_simple_rasts.tif",
                           sep = ""),
                     overwrite = T)
  
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  simple_projection(sp_name)
  NULL
  
}

# script termination
closeCluster(cl)
mpi.quit()