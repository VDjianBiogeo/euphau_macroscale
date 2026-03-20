#' @description
#' This function takes a dataset containing both species occurrences and 
#' environmental values associated in order to tune RandomForest and GBM 
#' models.
#' @param dataframe A dataset containing both species occurrences and 
#' environmental values for geographical positions.
#' 
#' @export
#' 

partition_tuning <- function(dataframe){
  
  # Getting species name as response name
  myRespName <- unique(dataframe$species)
  sp <- myRespName
  
  # Getting environmental data as raster
  env_rast <- terra::rast(read.csv(here::here("outputs", 
                                              "pre_mod_outputs", 
                                              sp,
                                              paste(sp,
                                                    "_selected_env.csv",
                                                    sep = ""))), 
                          type = "xyz")
  
  # Component class for each variable
  variable_classes <- read.csv(here::here("data",
                                          "treated",
                                          "env",
                                          "env_classes.csv"),
                               sep = ";", dec = ",")
  
  components <- unique(variable_classes$env_comp)[2:4]
  
  # For loop to tune models for each component
  for(comp in components){
    print(comp)
    
    if(sum(names(env_rast) %in% 
           variable_classes$Variable[variable_classes$env_comp == comp]) == 0){
      print(paste("No enviromental parameters for ",
            comp, " concerning ", sp,
            sep = ""))
    } else {
      
      if(file.exists(here::here("outputs",
                                "pre_mod_outputs",
                                sp,
                                paste(sp,
                                      "_",
                                      comp,
                                      "_mod_options.rds",
                                      sep = ""))) == T){
        print(paste(comp, "already treated for", sp))
      } else {
        
      
      # Generating pseudo-absences using biomod2 r package
      n_PA <- 3
      bm_data <- BIOMOD_FormatingData(resp.name = myRespName, 
                                      resp.var = dataframe[,4], 
                                      expl.var = env_rast[[names(env_rast) %in% 
                                                             variable_classes$Variable[variable_classes$env_comp == comp]]],
                                      resp.xy = dataframe[,1:2], 
                                      PA.nb.rep = n_PA, 
                                      PA.nb.absences = as.integer(nrow(dataframe)),
                                      PA.strategy = "sre",
                                      PA.sre.quant = 0.025,
                                      na.rm = T)
      
      # Saving biomod formating output
      saveRDS(bm_data, 
              here::here("outputs",
                         "pre_mod_outputs",
                         sp,
                         paste(sp,
                               "_",
                               comp,
                               "_biomod_format.rds",
                               sep = "")))
      
      # Extracting all pseudo-absences datasets generated
      list_PA <- list()
      for(i in 1:n_PA){
        # Combining spatial coordinates with PA
        pa <- cbind(bm_data@coord[bm_data@PA.table[,i],],
                    bm_data@data.species[bm_data@PA.table[,i]])
        
        # Changing NA to 0
        pa[,3][is.na(pa[,3])] <- 0
        
        # Changing colname for simplicity
        colnames(pa)[3] <- paste("PA", i, sep="")
        
        list_PA[[i]] <- pa
      }
      names(list_PA) <- paste("PA", 1:n_PA, sep="")
      
      # Finding hyperparameters for each PA dataset
      # Importing h2o package for grid search
      
      library(h2o)
      
      # Initiating h2o 
      h2o.init()
      
      gbm.param.list <- list()
      rf.param.list <- list()
      
      # For loop in order to define best hyperparameters possible
      for(i in 1:n_PA){
        df <- list_PA[[i]]
        
        # Changing response variable to a categorical one
        df[,3] <- as.factor(df[,3]) 
        
        # Adding environmental variables by joining with coordinates
        df <- inner_join(df, as.data.frame(env_rast, xy = T),
                         by = c("x", "y"))
        
        # Converting to h2o format
        h2o.df <- as.h2o(df)
        
        # Creating hyperparameter grid for RF and GBM
        hyper.grid.rf <- list(ntrees = seq(100, 5000, by = 200),
                              mtries = 1:length(env_rast),
                              max_depth = seq(5, 30, by = 5),
                              min_rows = seq(1, 5, by = 1),
                              nbins = seq(10, 30, by = 5))
        #sample_rate = seq(0.5, 0.8, by = 0.05))
        
        hyper.grid.gbm <- list(ntrees = seq(500, 3000, by = 200),
                               max_depth = seq(5, 30, by = 5),
                               min_rows = seq(1, 5, by = 1),
                               nbins = seq(10, 30, by = 5),
                               #sample_rate = seq(0.5, 0.8, by = 0.05),
                               learn_rate = c(0.01, 0.05, 0.1))
        #learn_rate_annealing = 0.99)
        
        # Defining search criteria from both models
        search_criteria <- list(strategy = "RandomDiscrete",
                                max_models = 400,
                                max_runtime_secs = 5*60,
                                stopping_rounds = 2,                
                                stopping_metric = "AUC",
                                stopping_tolerance = 0.01)
        
        # Grid search to find best hyperparameters for RF and GBM
        rf.grid <- h2o.grid(algorithm = "randomForest", 
                            y = 3, 
                            x = 4:ncol(h2o.df), 
                            training_frame = h2o.df, 
                            search_criteria = search_criteria, 
                            hyper_params = hyper.grid.rf, 
                            seed = 4951)
        
        gbm.grid <- h2o.grid(algorithm = "gbm", 
                             y = 3, 
                             x = 4:ncol(h2o.df), 
                             training_frame = h2o.df, 
                             search_criteria = search_criteria, 
                             hyper_params = hyper.grid.gbm, 
                             seed = 4951)
        
        # Extracting hyperparameter grid
        grid_perf_rf <- h2o.getGrid(grid_id = rf.grid@grid_id,
                                    sort_by = "AUC",
                                    decreasing = F)
        
        grid_perf_gbm <- h2o.getGrid(grid_id = gbm.grid@grid_id,
                                     sort_by = "AUC",
                                     decreasing = F)
        
        # Extracting best model evaluated
        best_rf <- grid_perf_rf@summary_table[1,]
        best_rf$PA <- names(list_PA[i])
        
        best_gbm <- grid_perf_gbm@summary_table[1,]
        best_gbm$PA <- names(list_PA[i])
        
        # Adding best hyperparameters to the list
        rf.param.list[[i]] <- as.data.frame(best_rf)
        gbm.param.list[[i]] <- as.data.frame(best_gbm)
      }
      
      # Defining model parameters for each PA dataset
      user.val <- list(RF.binary.randomForest.randomForest = list("_PA1_allRun" = list("mtry"=rf.param.list[[1]]$mtries,
                                                                                       "ntree"=rf.param.list[[1]]$ntrees,
                                                                                       "nodesize"=rf.param.list[[1]]$min_rows),
                                                                  "_PA2_allRun" = list("mtry"=rf.param.list[[2]]$mtries,
                                                                                       "ntree"=rf.param.list[[2]]$ntrees,
                                                                                       "nodesize"=rf.param.list[[2]]$min_rows),
                                                                  "_PA3_allRun" = list("mtry"=rf.param.list[[3]]$mtries,
                                                                                       "ntree"=rf.param.list[[3]]$ntrees,
                                                                                       "nodesize"=rf.param.list[[3]]$min_rows)),
                       GBM.binary.gbm.gbm = list("_PA1_allRun" = list("n.trees" = gbm.param.list[[1]]$ntrees,
                                                                      "interaction.depth" = gbm.param.list[[1]]$max_depth,
                                                                      "n.minobsinnode" = gbm.param.list[[1]]$min_rows,
                                                                      "bag.fraction" = gbm.param.list[[1]]$sample_rate,
                                                                      "shrinkage" = gbm.param.list[[1]]$learn_rate),
                                                 "_PA2_allRun" = list("n.trees" = gbm.param.list[[2]]$ntrees,
                                                                      "interaction.depth" = gbm.param.list[[2]]$max_depth,
                                                                      "n.minobsinnode" = gbm.param.list[[2]]$min_rows,
                                                                      "bag.fraction" = gbm.param.list[[2]]$sample_rate,
                                                                      "shrinkage" = gbm.param.list[[2]]$learn_rate),
                                                 "_PA3_allRun" = list("n.trees" = gbm.param.list[[3]]$ntrees,
                                                                      "interaction.depth" = gbm.param.list[[3]]$max_depth,
                                                                      "n.minobsinnode" = gbm.param.list[[3]]$min_rows,
                                                                      "bag.fraction" = gbm.param.list[[3]]$sample_rate,
                                                                      "shrinkage" = gbm.param.list[[3]]$learn_rate)))
      # Save modelling user values
      saveRDS(user.val,
              here::here("outputs",
                         "pre_mod_outputs",
                         sp,
                         paste(sp,
                               "_",
                               comp,
                               "_user_options.rds",
                               sep = "")))
      
      # Creating biomod modelling options
      myOpt <- bm_ModelingOptions(data.type = 'binary',
                                  models = c("RF", "GBM"),
                                  strategy = "user.defined",
                                  user.val = user.val, 
                                  user.base = "bigboss",
                                  bm.format = bm_data)
      
      # Saving biomod modelling options
      saveRDS(myOpt,
              here::here("outputs",
                         "pre_mod_outputs",
                         sp,
                         paste(sp,
                               "_",
                               comp,
                               "_mod_options.rds",
                               sep = "")))
      
      h2o.removeAll()
      h2o.shutdown(prompt = F)
      }
    }
    

  }
  
  return("done")
}