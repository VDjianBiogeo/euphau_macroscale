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

cl <- doMPI::startMPIcluster()
registerDoMPI(cl)

# Simple evaluation ----
simple_evaluation <- function(sp_name){
  # Model evaluation ----
  bm_data <- readRDS(paste("./outputs/pre_mod_outputs/",
                           sp_name,
                           "/",
                           sp_name,
                           "_biomod_format.rds",
                           sep = ""))
  
  Proj <- terra::rast(paste("./outputs/mod_outputs/",
                            sp_name,
                            "/simple/",
                            sp_name,
                            "_simple_rasts.tif",
                            sep = ""))
  
  # Extracting PA datasets with spatial coordinates
  PA_set <- cbind(bm_data@coord, bm_data@PA.table)
  
  # Changing PA into a 0/1 format
  PA_set[,3:5] <- sapply(PA_set[,3:5],as.numeric)
  
  ## Turning projections into binary ----
  
  ## Reading raw occurrences 
  occurrences <- read.csv("./outputs/pre_mod_outputs/euphau_longer_occ.csv")
  occurrences <- occurrences[occurrences$species == sp_name,]
  
  ## Joining probability rasters with raw occurrences
  joined_df <- inner_join(x = occurrences,
                          y = as.data.frame(Proj, xy = T),
                          by = c("x", "y")) |>
    na.omit() |>
    as.data.frame()
  
  ## Computing quantiles and turn probability into binary
  quantile_models <- data.frame(model_name = NA,
                                quantile_005 = NA)
  
  ## Turning projections into a dataframe
  simple_df <- as.data.frame(Proj, xy = T)
  
  binary_rast <- terra::rast(xmin = -180,
                             xmax = 180,
                             ymin = -90,
                             ymax = 90,
                             resolution = 0.25,
                             crs = "epsg:4326")
  
  for(i in colnames(joined_df[, 6:ncol(joined_df)])){
    
    # Computing 5% quantile
    for_binary <- quantile(joined_df[,i],
                           probs = c(0.05))[[1]]
    
    # Binary transformation
    binary_proj <- biomod2::bm_BinaryTransformation(data = terra::rast(simple_df[,c("x","y",i)], 
                                                                       type = "xyz"), 
                                                    threshold = for_binary)
    
    # Stocking binary output
    binary_rast <- c(binary_rast, 
                     binary_proj)
    
    
  }
  
  ## Create the dataframe to use to compute the AUC metric for all models
  eval_dataset <- inner_join(PA_set, 
                             as.data.frame(binary_rast, xy = T),
                             by = c("x", "y"))
  
  metric_df <- data.frame(model_name = NA, auc = NA)
  
  # Computing AUC metric
  for(i in names(Proj)){
    
    if(length(grep(pattern = "PA1", x = i)) == 1){
      auc <- pROC::roc(eval_dataset$PA1, 
                       eval_dataset[, i])$auc[1]
      
      metric_df <- rbind(metric_df,
                         data.frame(model_name = i,
                                    auc = auc))
    } else {
      if(length(grep(pattern = "PA2", x = i)) == 1){
        auc <- pROC::roc(eval_dataset$PA2, 
                         eval_dataset[, i])$auc[1]
        
        metric_df <- rbind(metric_df,
                           data.frame(model_name = i,
                                      auc = auc))
      } else {
        if(length(grep(pattern = "PA3", x = i)) == 1){
          auc <- pROC::roc(eval_dataset$PA3, 
                           eval_dataset[, i])$auc[1]
          
          metric_df <- rbind(metric_df,
                             data.frame(model_name = i,
                                        auc = auc)) |>
            na.omit()
        }
      }
    }
    
  }
  
  metric_df$species <- sp_name
  
  # Writing metric dataframe as a csv file
  write.csv(metric_df, 
            paste("./outputs/mod_outputs/",
                  sp_name,
                  "/simple/",
                  sp_name,
                  "_auc_metrics.csv",
                  sep = ""),
            row.names = F)
  return("done")
}

# Parallel computation:
foreach(sp_name = list_sp) %dopar% {
  simple_evaluation(sp_name)
  
}

# script termination
closeCluster(cl)
mpi.quit()