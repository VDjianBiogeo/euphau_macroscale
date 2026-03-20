#'@description
#'This function launches single models using biomod2 for a defined species.
#'@param sp_name a species name from a list.
#'
#'@export

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
                                modeling.id = "demo1")
  
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
  
  saveRDS(Proj,
          paste("./outputs/mod_outputs/", 
                sp_name,
                "/simple/",
                sp_name,
                "_simple_rasts.rds", 
                sep = ""))
  
  return("done")
}