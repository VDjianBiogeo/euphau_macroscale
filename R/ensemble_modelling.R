#'@description
#'This function launches ensemble models using biomod2 for a defined species.
#'@param sp_name a species name from a list.
#'
#'@export

ensemble_modelling <- function(sp_name){
  
  # Get simple models ----
  sdm_models <- readRDS(paste("./outputs/mod_outputs/", 
                              sp_name,
                              "/simple/",
                              sp_name,
                              "_simple_model.rds", 
                              sep = ""))
  
  # Ensemble modelling ----
  myEnsemble <- BIOMOD_EnsembleModeling(bm.mod = sdm_models,
                                        models.chosen  = 'all',
                                        em.by = 'all',
                                        metric.select = c('ROC'),
                                        metric.select.thresh = c(0.7),
                                        metric.eval = c('ACCURACY','ROC', "TSS"),
                                        em.algo = c('EMmean', 'EMmedian',
                                                    'EMcv'),
                                        var.import = 5)
  
  saveRDS(myEnsemble,
          paste("./outputs/mod_outputs/", 
                sp_name,
                "/ensemble/",
                sp_name,
                "_ensemble_model.rds", 
                sep = ""))
  
  # Get environmental data ----
  env <- terra::rast(read.csv(paste("./outputs/pre_mod_outputs/",
                                    sp_name,
                                    "/",
                                    sp_name,
                                    "_selected_env.csv",
                                    sep = "")), 
                     type = "xyz")
  
  # Ensemble projection ----
  myEnsProj <- BIOMOD_EnsembleForecasting(bm.em = myEnsemble,
                                          new.env = env,
                                          proj.name = "ensemble",
                                          models.chosen = "all")
  saveRDS(myEnsProj,
          paste("./outputs/mod_outputs/",
                sp_name,
                "/ensemble/",
                sp_name,
                "_ensemble_proj.rds",
                sep = ""))
  
  # Save projections into a raster format ----
  EMProj <- biomod2::get_predictions(obj = mySDMProj)
  
  # Changing crs 
  terra::crs(Proj) <- "epsg:4326"
  
  saveRDS(Proj,
          paste("./outputs/mod_outputs/", 
                sp_name,
                "/ensemble/",
                sp_name,
                "_ensemble_rasts.rds", 
                sep = ""))
  
  return("done")
  
}