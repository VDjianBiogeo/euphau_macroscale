#' @description
#' This function computes species frequency and occurrences wolrdwide and at
#' the scale of the indian subantarctic from a species contingence dataframe.
#' 
#' @param data a list of species dataframes.
#' @export


species_freq <- function(data){
  
  list_freq <- list()
  for(i in 1:length(data)){
    
    df <- as.data.frame(data[[i]], xy = T)
    
    # Filtering on study area 
    df_subant <- df[df$x < 160 & df$x > -20 & df$y < -20 & df$y > -60,]
    
    # Number of cells where biophiids were seen worldwide
    ncell_world <- nrow(unique(df[,c("x","y")]))
    
    # Number of cells where biophiids were seen in the indian subantarctic
    ncell_subant <- nrow(unique(df_subant[,c("x","y")]))
    
    # Preparing future dataframe columns
    df_pivot <- tidyr::pivot_longer(data = df, cols = 3:ncol(df), 
                             names_to = "Species", values_to = "Presence") |>
      na.omit()
    
    df_subant <- df_pivot[df_pivot$x < 160 & 
                            df_pivot$x > -20 & 
                            df_pivot$y < -20 & 
                            df_pivot$y > -60,]
    
    list.species <- unique(df_pivot$Species[order(df_pivot$Species)])
    list.total_occ <- list()
    list.occ_subant <- list()
    list.freq_world <- list()
    list.freq_subant <- list()
    
    # Computing frequency for each species 
    
    for(sp in list.species){
      total_occ <- length(df_pivot$Presence[df_pivot$Species == sp])
      freq_world <- total_occ/ncell_world * 100
      
      occ_subant <- length(df_subant$Presence[df_subant$Species == sp])
      freq_subant <- occ_subant/ncell_subant * 100
      
      list.total_occ <- c(list.total_occ, total_occ)
      list.occ_subant <- c(list.occ_subant, occ_subant)
      list.freq_world <- c(list.freq_world, freq_world)
      list.freq_subant <- c(list.freq_subant, freq_subant)
      
    }
    
    freq_df <- data.frame(Species = list.species, 
                          Total_occ_world = unlist(list.total_occ),
                          Occ_subant = unlist(list.occ_subant),
                          Frequency_world = unlist(list.freq_world), 
                          Frequency_subant = unlist(list.freq_subant))
    
    # Exporting frequency data in a csv file
    write.csv(x = freq_df,
              file = here::here("outputs",
                                "explor_outputs",
                                paste0("species_freq_",
                                       i,
                                       ".csv")), 
              row.names = F)
    
    list.sp <- intersect(x = freq_df$Species[order(freq_df$Frequency_subant, 
                                                   decreasing = T)],
                         y = freq_df$Species[freq_df$Total_occ_world > 100 
                                             & freq_df$Occ_subant > 5])
    
    list_freq[[i]] <- list.sp
    
    
  
  }
  saveRDS(object = list_freq,
          file = here::here("outputs",
                            "explor_outputs",
                            "species_lists.rds"))
  return(list_freq)
}
