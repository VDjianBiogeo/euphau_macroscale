#' @description
#' This function extracts species projections for each environmental component
#' and full model to compute partial semi-correlations.
#' 
#' @param sp_list a list of species on which to compute semi-correlations.
#' 
#' @export
#' 

semi_corr <- function(sp_list){
  
  # Component class for each variable ----
  variable_classes <- read.csv("./data/treated/env/env_classes.csv",
                               sep = ";") 
  
  components <- unique(variable_classes$env_comp)[2:4]
  
  
  all_corr <- data.frame("water_mass" = NA, 
                        "ocean_dynamics" = NA, 
                        "biologic_production" = NA)
  
  all_dist <- data.frame("water_mass" = NA, 
                         "ocean_dynamics" = NA, 
                         "biologic_production" = NA)
  
  
  for(sp in sp_list){
    
    comp_predictions <- terra::rast(xmin = -180, 
                                    xmax = 180, 
                                    ymin = -90, 
                                    ymax = 90, 
                                    res = 0.25)
    
    if(length(list.dirs(path = here::here("outputs", 
                                          "mod_outputs",
                                          sp, "env_part"))) >= 4){
      
      # full_pred <- readRDS(here::here("outputs",
      #                                 "mod_outputs",
      #                                 sp,
      #                                 "ensemble",
      #                                 paste(sp,
      #                                       "_ensemble_proj.rds",
      #                                       sep = "")))
      # 
      # full_proj <- biomod2::get_predictions(full_pred)[[2]]
      # 
      full_proj <- terra::rast(here::here("outputs",
                                          "mod_outputs",
                                          sp,
                                          "ensemble",
                                          paste(sp,
                                                "_ensemble_rasts.tif",
                                                sep = "")))[[2]]
      
      terra::crs(full_proj) <- "epsg:4326"
      names(full_proj) <- paste(sp, "full", sep = "_")
      
      comp_predictions <- c(comp_predictions, full_proj)
      
      for(comp in components){
        
        pred <- terra::rast(here::here("outputs",
                                       "mod_outputs",
                                       sp,
                                       "env_part",
                                       comp,
                                       paste(sp,
                                             comp,
                                             "ensemble_rasts.tif",
                                             sep = "_")))[[2]]
        terra::crs(pred) <- "epsg:4326"
        names(pred) <- paste(sp, comp, sep = "_")
        comp_predictions <- c(comp_predictions, pred)
      }
      
      # Saving species predictions
      dir.create(here::here("outputs",
                            "post_mod_outputs",
                            sp,
                            "env_part"))
      
      terra::writeRaster(x = comp_predictions,
                         filename = here::here("outputs",
                                               "post_mod_outputs",
                                               sp,
                                               "env_part",
                                               paste(sp,
                                                     "env_part_pred.tif",
                                                     sep = "_")),
                         overwrite = T)
      
      write.csv(x = as.data.frame(comp_predictions, xy = T),
                file = here::here("outputs",
                                  "post_mod_outputs",
                                  sp,
                                  "env_part",
                                  paste(sp,
                                        "env_part_pred.tif",
                                        sep = "_")),
                row.names = F)
      
      # Semi-partial correlations between predictions
      spcor_results <- ppcor::spcor(na.omit(as.data.frame(comp_predictions)), 
                                    "spearman")
      corr_matrix <- spcor_results$estimate
      
      corr_matrix <- as.data.frame(t(corr_matrix))[1,-1]
      
      colnames(corr_matrix) <- colnames(all_corr)
      
      all_corr <- rbind(all_corr, corr_matrix)
      
      # Saving semi-partial correlations
      write.csv(x = corr_matrix,
                file = here::here("outputs",
                                  "post_mod_outputs",
                                  sp,
                                  "env_part",
                                  paste(sp,
                                        "semi_corr.csv",
                                        sep = "_")),
                row.names = F)
      
      # Computing distance matrix
      distMatrix <- 1 - abs(corr_matrix) %>%
        as.data.frame()
      
      all_dist <- rbind(all_dist, distMatrix)
      
      # Saving distance matrix output
      write.csv(x = distMatrix,
                file = here::here("outputs",
                                  "post_mod_outputs",
                                  sp,
                                  "env_part",
                                  paste(sp,
                                        "env_part_dist.csv",
                                        sep = "_")),
                row.names = F)
      
    } else {
      print(paste("No data available for all components for", sp))
    } 
    
  }
  
  # Saving all semi-correlations
  write.csv(x = na.omit(all_corr),
            file = here::here("outputs",
                              "post_mod_outputs",
                              "euphau_all_corr.csv"),
            row.names = F)
  # Saving all distances
  write.csv(x = na.omit(all_dist),
            file = here::here("outputs",
                              "post_mod_outputs",
                              "euphau_all_dist.csv"),
            row.names = F)
  
  return(list(na.omit(all_corr), na.omit(all_dist)))
  
}