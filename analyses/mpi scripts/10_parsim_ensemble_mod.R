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


parsim_en_mod <- function(sp_name){
  
  # Get simple models ----
  sdm_models <- readRDS(paste("./outputs/mod_outputs/",
                              sp_name,
                              "/parsim/simple/",
                              sp_name,
                              "_parsim_simple_model.rds",
                              sep = ""))
  
  # Get simple model metrics ----
  sdm_metrics <- read.csv(paste("./outputs/mod_outputs/",
                                sp_name,
                                "/parsim/simple/",
                                sp_name,
                                "_parsim_auc_metrics.csv",
                                sep = ""))
  
  # Models to choose ----
  models_to_take <- sdm_metrics$model_name[sdm_metrics$auc >= 0.7]
  
  # Ensemble modelling ----
  myEnsemble <- BIOMOD_EnsembleModeling(bm.mod = sdm_models, 
                                        models.chosen = models_to_take,
                                        em.by = 'all',
                                        # metric.select = "user.defined", 
                                        # metric.select.thresh = 0.7, 
                                        # metric.select.table = sdm_metrics[,1:2],
                                        # EMwmean.decay = "proportional",
                                        metric.eval = c('ROC'),
                                        em.algo = c('EMmedian','EMcv'),
                                        var.import = 5,
                                        metric.select.dataset = "validation")
  
  # Saving output
  dir.create(paste("./outputs/mod_outputs/",
                   sp_name,
                   "/parsim/ensemble",
                   sep = ""))
  
  saveRDS(myEnsemble,
          paste("./outputs/mod_outputs/",
                sp_name,
                "/parsim/ensemble/",
                sp_name,
                "_parsim_ensemble_model.rds",
                sep = ""))
  
  return("done")
}


# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  parsim_en_mod(sp_name)
  
}

# script termination
closeCluster(cl)
mpi.quit()