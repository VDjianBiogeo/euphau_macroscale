#' @description
#' This function computes the frequency of species occurrences by Longhurst 
#' provinces.
#' 
#' @param data a sf dataframe containing species and province information.
#' 
#' @export

compute_freq <- function(data){
  
  freq_df <- data.frame(Species = NA, 
                        ProvCode = NA,
                        Count = NA,
                        Freq = NA)
  
  for(sp in unique(data$Species)){
    
    data_filtered <- data[data$Species == sp,]
    for(prov in unique(data_filtered$ProvCode)){
      
      total_count <- nrow(data_filtered)
      
      count <- nrow(data_filtered[data_filtered$ProvCode == prov,])
      
      df <- data.frame(Species = unique(data_filtered$Species),
                       ProvCode = prov,
                       Count = count,
                       Freq = (count/total_count) * 100)
      
      freq_df <- rbind(freq_df, df)
      
    }
  }
  return(na.omit(freq_df))
}