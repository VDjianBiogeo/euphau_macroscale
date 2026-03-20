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

parsimonious_sproj <- function(sp_name){
  
  bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                           sp_name,
                           "/",
                           sp_name,
                           "_parsim_biomod_format.rds",
                           sep = ""))
  
  # Getting simple models ----
  simple_models <- readRDS(paste("./outputs/mod_outputs/",
                                 sp_name,
                                 "/parsim/simple/",
                                 sp_name,
                                 "_parsim_simple_model.rds",
                                 sep = ""))
  
  # Getting environmental data ----
  env <- terra::rast(read.csv(paste("./outputs/pre_mod_outputs/",
                                    sp_name,
                                    "/",
                                    sp_name,
                                    "_selected_env.csv",
                                    sep = "")),
                     type = "xyz")[[colnames(bm_data@data.env.var)]]
  
  # Simple projection ----
  mySDMProj <- BIOMOD_Projection(bm.mod = simple_models,
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
                "/parsim/simple/",
                sp_name,
                "_parsim_simple_proj.rds",
                sep = ""))
  
  # Save projections to a raster format ----
  Proj <- biomod2::get_predictions(obj = mySDMProj)
  
  # Changing crs
  terra::crs(Proj) <- "epsg:4326"
  
  terra::writeRaster(Proj,
                     paste("./outputs/mod_outputs/",
                           sp_name,
                           "/parsim/simple/",
                           sp_name,
                           "_parsim_simple_rasts.tif",
                           sep = ""),
                     overwrite = T)
  
  return("done")
}


# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  parsimonious_sproj(sp_name)
  
}

# script termination
closeCluster(cl)
mpi.quit()