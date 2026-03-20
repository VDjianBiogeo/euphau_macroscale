#'@description
#'This function maps species occurrences while specifying which occurrences are 
#'in the native range and which are not.
#'
#'@param data a sf dataframe with species occurrences and information about 
#'occurrences' position compared to species native range.
#'
#'@export

native_range_map <- function(data){
  
  for(sp in unique(data$Species)){
    data("World", package = "tmap")
    
    # Filtering given dataframe to only keep sp occurrences
    map_df <- data[data$Species == sp,]
    
    # Mapping
    sp_map <- ggplot() +
      geom_sf(data = World) +
      geom_sf(aes(color = factor(to_keep),
                  shape = factor(to_keep)), 
              data = map_df) +
      scale_color_manual(values = c("0" = "red", "1" = "green")) +
      scale_shape_manual(values = c("0" = 15, "1" = 17)) +
      coord_sf(expand = F, crs = 4326) +
      theme_bw(base_size = 18) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      labs(color = "To keep", shape = "To keep", title = sp) +
      theme(plot.title = element_text(hjust = 0.5, face = "italic"))
    
    # Saving maps as png
    if(dir.exists(here::here("figures",
                             "data_treatment",
                             "occurrence_filtering")) == F){
      
      dir.create(here::here("figures",
                            "data_treatment",
                            "occurrence_filtering"))
    }
    
    ggsave(filename = here::here("figures",
                                 "data_treatment",
                                 "occurrence_filtering",
                                 paste(paste(str_sub(string = str_split(sp, 
                                                                        " ", 
                                                                        simplify = T)[,1], 
                                                     start = 1, 
                                                     end = 1), 
                                             str_split(sp, 
                                                       " ", 
                                                       simplify = T)[,2], 
                                             sep = "_"),
                                       "in_out_range_map.png",
                                       sep = "_")),
           plot = sp_map,
           dpi = "retina",
           height = 2940,
           width = 3660,
           units = "px")
  }
}