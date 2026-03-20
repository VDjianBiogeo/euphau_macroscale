# Setup ----
rm(list = ls())
gc()

set.seed(4951)

library(doMPI)
library(tidyverse)
library(biomod2)
library(terra)
library(h2o)

# # Creating output subdirectory ----
dir.create("./outputs/mod_outputs")

# Creating figures subdirectory ----
dir.create("./figures/mod_figures")

# Importing species list ----
list_sp <- str_split(readRDS("./outputs/explor_outputs/species_lists.rds")[[1]], 
                     pattern = "_",
                     simplify = T)[,1]
list_sp <- list_sp[order(list_sp)]

cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

# Directory creation ----
dir_create <- function(sp_name){
  
  # Creating species output directory if needed
  if(dir.exists(paste("./outputs/mod_outputs/",
                      sp_name,
                      sep = "")) == F){
    dir.create(paste("./outputs/mod_outputs/",
                     sp_name,
                     sep = ""))
  }
  
  # Creating species figure directory if needed
  if(dir.exists(paste("./figures/mod_figures/",
                      sp_name,
                      sep = "")) == F){
    dir.create(paste("./figures/mod_figures/",
                     sp_name,
                     sep = ""))
  }
  
  # Checking if simple directory exists
  if(dir.exists(paste("./outputs/mod_outputs/",
                      sp_name,
                      "/simple",
                      sep = "")) == F){
    dir.create(paste("./outputs/mod_outputs/",
                     sp_name,
                     "/simple",
                     sep = ""))
    
    if(dir.exists(paste("./outputs/mod_outputs/",
                        sp_name,
                        "/ensemble",
                        sep = "")) == F){
      dir.create(paste("./outputs/mod_outputs/",
                       sp_name,
                       "/ensemble",
                       sep = ""))
    }
  }
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  dir_create(sp_name)
  NULL
  
}

# # script termination
# closeCluster(cl)
# mpi.quit()

# Simple species modelling ----
# source("./R/simple_modelling.R")
simple_modelling <- function(sp_name){
  
  # Getting modelling option ----
  myOpt <- readRDS(paste("./outputs/pre_mod_outputs/",
                         sp_name,
                         "/",
                         sp_name,
                         "_mod_options.rds",
                         sep = ""))
  
  # Getting bm_data object ----
  bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                           sp_name,
                           "/",
                           sp_name,
                           "_biomod_format.rds",
                           sep = ""))
  
  # Simple modelling ----
  sdm_models <- BIOMOD_Modeling(bm.format = bm_data,
                                models = c("RF", "GBM"),
                                OPT.strategy = "user.defined",
                                OPT.user = myOpt,
                                CV.strategy = "kfold",
                                CV.k = 5,
                                CV.nb.rep = 5,
                                CV.do.full.models = F,
                                metric.eval = c("ROC", "ACCURACY", "TSS"),
                                var.import = 10,
                                modeling.id = "full1")
  
  saveRDS(sdm_models,
          paste("./outputs/mod_outputs/",
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

# cl <- doMPI::startMPIcluster()
# registerDoMPI(cl)

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  simple_modelling(sp_name)
  NULL
  
}

# script termination
closeCluster(cl)
mpi.quit()