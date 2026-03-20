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

env_simple_mod <- function(sp_name){
  # Component class for each variable ----
  variable_classes <- read.csv("./data/treated/env/env_classes.csv",
                               sep = ";") 
  
  components <- unique(variable_classes$env_comp)[2:4]
  
  # Creating output directory ----
  dir.create(paste("./outputs/mod_outputs/",
                   sp_name,
                   "/env_part",
                   sep = ""), showWarnings = F)
  
  for(comp in components){
    if(file.exists(paste("./outputs/pre_mod_outputs/",
                         sp_name,
                         "/",
                         sp_name,
                         "_",
                         comp,
                         "_biomod_format.rds",
                         sep = "")) == T){
      
      # Creating component output directory ----
      dir.create(paste("./outputs/mod_outputs/",
                       sp_name,
                       "/env_part/",
                       comp,
                       sep = ""))
      
      # Reading needed data for modelling ----
      bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                               sp_name,
                               "/",
                               sp_name,
                               "_",
                               comp,
                               "_biomod_format.rds",
                               sep = ""))
      
      # Getting modelling option ----
      myOpt <- readRDS(paste("./outputs/pre_mod_outputs/",
                             sp_name,
                             "/",
                             sp_name,
                             "_",
                             comp,
                             "_mod_options.rds",
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
                                    modeling.id = paste(comp, "1", sep = ""))
      
      saveRDS(sdm_models,
              paste("./outputs/mod_outputs/",
                    sp_name,
                    "/env_part/",
                    comp,
                    "/",
                    sp_name,
                    "_",
                    comp,
                    "_simple_model.rds",
                    sep = ""))
      
      
    } else {
      print(paste("no data available for ", comp, " for ", sp_name, sep = ""))
    }
  }
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  env_simple_mod(sp_name)
  NULL
  
}

# script termination
closeCluster(cl)
mpi.quit()