# Setup ----
rm(list = ls())
gc()

set.seed(4951)

library(doMPI)
library(tidyverse)
library(biomod2)
library(terra)
library(h2o)

# Importing species list ----
list_sp <- str_split(readRDS("./outputs/explor_outputs/species_lists.rds")[[1]], 
                     pattern = "_",
                     simplify = T)[,1]
list_sp <- list_sp[order(list_sp)]

# Launching clusters ----
cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

parsimonious_simple <- function(sp_name){
  
  # Creating parsimonious output directory
  dir.create(paste("./outputs/mod_outputs/",
                   sp_name,
                   "/parsim",
                   sep = ""))
  
  # Parsimonious simple modelling ----
  
  # Getting modelling option ----
  myOpt <- readRDS(paste("./outputs/pre_mod_outputs/",
                         sp_name,
                         "/",
                         sp_name,
                         "_parsim_mod_options.rds",
                         sep = ""))
  
  # Getting bm_data object ----
  bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                           sp_name,
                           "/",
                           sp_name,
                           "_parsim_biomod_format.rds",
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
                                var.import = 5,
                                modeling.id = "parsim1")
  
  # Saving output
  dir.create(paste("./outputs/mod_outputs/",
                   sp_name,
                   "/parsim/simple",
                   sep = ""))
  
  saveRDS(sdm_models,
          paste("./outputs/mod_outputs/",
                sp_name,
                "/parsim/simple/",
                sp_name,
                "_parsim_simple_model.rds",
                sep = ""))
  
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  parsimonious_simple(sp_name)
  
}

# script termination
closeCluster(cl)
mpi.quit()