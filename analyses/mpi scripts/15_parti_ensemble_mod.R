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

env_ensemble_mod <- function(sp_name){
  
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
                  comp, sep = ""))
    } else{
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
      
      # Get simple model metrics ----
      sdm_metrics <- read.csv(paste("./outputs/mod_outputs/",
                                    sp_name,
                                    "/env_part/",
                                    comp,
                                    "/",
                                    sp_name,
                                    "_",
                                    comp,
                                    "_auc_metrics.csv",
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
                                            metric.eval = c("TSS",'ACCURACY','ROC'),
                                            em.algo = c('EMmean', 'EMmedian','EMcv'),
                                            var.import = 5,
                                            metric.select.dataset = "validation")
      
      saveRDS(myEnsemble,
              paste("./outputs/mod_outputs/",
                    sp_name,
                    "/env_part/",
                    comp,
                    "/",
                    sp_name,
                    "_",
                    comp,
                    "_ensemble_model.rds",
                    sep = ""))
    }
    
  }
  
  
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  env_ensemble_mod(sp_name)
  NULL
  
}

# script termination
closeCluster(cl)
mpi.quit()