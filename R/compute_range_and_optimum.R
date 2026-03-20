#' @description
#' This function joins binary species predictions and environmental 
#' descriptors in order to compute environmental ranges and optimums.
#' 
#' @param binary_and_env binary species distributions and environmental descriptors
#' spatially joined as a dataframe.
#' 
#' @export
#' 

compute_range_and_optimum <- function(binary_and_env){
  
  pivot_df <- binary_and_env %>%
    pivot_longer(cols = 3:31, names_to = "species", values_to = "occurrences")
  
  pivot_df$species <- str_split(pivot_df$species, "_", simplify = T)[,1]
  
  pivot_df <- pivot_df[pivot_df$occurrences == 1,]
  
  all_sp_ranges <- data.frame(species = NA,
                              chl_range = NA,
                              chl_optimum = NA,
                              KE_range = NA,
                              KE_optimum = NA,
                              MLD_range = NA,
                              MLD_optimum = NA,
                              o2_range = NA,
                              o2_optimum = NA,
                              sal_range = NA,
                              sal_optimum = NA,
                              SSH_range = NA,
                              SSH_optimum = NA,
                              temp_range = NA,
                              temp_optimum = NA)
  
  for(sp in unique(pivot_df$species)){
    
    sp_df <- data.frame(species = sp)
    
    pivot_quantile <- pivot_df[pivot_df$species == sp,]
    
    for(env_var in colnames(binary_and_env)[grep("median", 
                                                 colnames(binary_and_env))]){
      
      # Computing environmental range as interquantile difference
      range <- (quantile(dplyr::pull(.data = pivot_quantile, 
                                     var = env_var), 
                         probs = c(0.95))[[1]] - quantile(dplyr::pull(.data = pivot_quantile, 
                                                                      var = env_var), 
                                                          probs = c(0.05))[[1]])
      
      # Computing optimum via computing density
      density <- stats::density(x = dplyr::pull(.data = pivot_quantile, 
                                                var = env_var))
      density_df <- data.frame(x = density$x,
                               y = density$y)
      optimum <- density_df$x[density_df$y == max(density_df$y)]
      
      # optimum <- quantile(dplyr::pull(.data = pivot_quantile,
      #                                 var = env_var),
      #                     probs = c(0.5))[[1]]
      
      sp_df <- cbind(sp_df,
                     data.frame(tempo_range = range,
                                tempo_optimum = optimum))
      
      var <- str_split(env_var, "_")[[1]][1]
      colnames(sp_df)[grep("tempo_range", colnames(sp_df))] <- paste(var, 
                                                               "range",
                                                               sep = "_")
      
      colnames(sp_df)[grep("tempo_optimum", colnames(sp_df))] <- paste(var, 
                                                                     "optimum",
                                                                     sep = "_")
      
    }
    
    all_sp_ranges <- rbind(all_sp_ranges,
                           sp_df) %>%
      na.omit()
  }
  
  return(all_sp_ranges)
}