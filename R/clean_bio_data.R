clean_bio_data <- function(raw_data){
  
  sp_list <- readRDS(here::here("outputs",
                                "explor_outputs",
                                "species_lists.rds"))
  
  filtered_df <- species_filtering(dataset = raw_data, 
                    sp_list = sp_list, 
                    xy = T, 
                    xy_name = c("x","y")) 
  
  filtered_df <- na.omit(tidyr::pivot_longer(data = filtered_df, 
                                             cols = 3:ncol(filtered_df), 
                                             names_to = "species_name",
                                             values_to = "occurrence"))
  
  return(filtered_df)
}