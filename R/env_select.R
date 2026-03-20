#' @description
#' This function takes an occurrence dataset and filter it based on a list of 
#' species to keep.
#' @param filtered_occ A list of species to select in the given dataset.
#' @param env_data A dataset containing species occurrences as columns.
#' @param coord_names names of the columns matching spatial coordinates.
#' 
#' @export
#' 
library(tidyverse)
env_select <- function(filtered_occ, env_data, coord_names){

  # Creating output directory
  if(dir.exists(paths = here::here("outputs",
                                   "pre_mod_outputs")) == F){
    dir.create(here::here("outputs",
                          "pre_mod_outputs"))
  }
  
  # Creating figures directory
  if(dir.exists(paths = here::here("figures",
                                   "pre_mod_figures")) == F){
    dir.create(here::here("figures",
                          "pre_mod_figures"))
  }
  
  # Changing coordinates names to x and y
  colnames(filtered_occ[, coord_names]) <- c("x","y")
  
  
  # Getting species names
  list_sp <- colnames(subset(filtered_occ, 
                             select = -c(x,y)))
    
  list_sp <- list_sp[order(list_sp)]
  
  # Creating an empty raster for output
  env_list <- list()
  
  for(sp in list_sp){
    
    dir.create(path = here::here("figures",
                                 "pre_mod_figures",
                                 sp))
    
    dir.create(path = here::here("outputs",
                                 "pre_mod_outputs",
                                 sp))
    
    # Joining occurrence and environment data
    sp_df <- inner_join(filtered_occ[,c("x","y", sp)],
                        env_data,
                        by = c("x", "y")) |>
      na.omit()
    
    # Mapping species occurrence
    data(World, package = "tmap")
    
    ## Worldwide
    bio.map <- ggplot() +
      geom_sf(data = World) +
      geom_point(aes(x=x, y=y), 
                 fill="green", 
                 color="black", 
                 pch = 21, 
                 size = 2,
                 data = sp_df[,c(1:3)]) +
      coord_sf(crs = sf::st_crs(4326), expand = F) +
      labs(title = paste(sp, 
                         " world distribution",
                         sep = "")) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(plot = bio.map, 
           filename = here::here("figures",
                                 "pre_mod_figures",
                                 sp,
                                 paste(sp, 
                                       "_worldwide_map.png", 
                                       sep="")), 
           height = 2200, width = 3660, units = "px",
           dpi = "retina")
    
    ## Regional scale
    bio.map <- ggplot() +
      geom_sf(data = World) +
      geom_point(aes(x=x, y=y), 
                 fill="green", 
                 color="black", 
                 pch = 21, 
                 size = 2,
                 data = sp_df[,c(1:3)]) +
      coord_sf(crs = sf::st_crs(4326), expand = F, 
               xlim = c(0,110), ylim = c(-60, -20)) +
      labs(title = paste(sp, 
                         " world distribution",
                         sep = "")) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(plot = bio.map, 
           filename = here::here("figures",
                                 "pre_mod_figures",
                                 sp,
                                 paste(sp, 
                                       "_region_map.png", 
                                       sep="")), 
           height = 2200, width = 3660, units = "px",
           dpi = "retina")
    
    # Selecting environmental variables based on multicollinearity and simple correlation 
    ## Computing correlation between environmental variables for future interpretations
    corMatrix <- Hmisc::rcorr(x = as.matrix(sp_df[,-c(1:3)]), type = "pearson")
    distMatrix <- 1 - abs(corMatrix[["r"]])
    
    # Clustering variables based on correlation distance
    library(stats)
    clust <- hclust(d = as.dist(distMatrix), 
                           method = "average")
    dendro <- as.dendrogram(clust)
    
    # Cutting at a correlation coefficient of 0.7
    cut <- as.data.frame(dendextend::cutree(tree = dendro, h = 0.3))
    cut$Variables <- rownames(cut)
    
    # Plotting dendrogram trees
    dendro_graph <- ggdendro::ggdendrogram(dendro) +
      geom_hline(yintercept = 0.3, 
                 color = "red", 
                 linewidth = 0.7, 
                 linetype = 2) +
      labs(y = "Distance (1 - spearman)",
           title = paste("Environmental variables' correlation for ",
                         sp, 
                         " occurrences", 
                         sep="")) +
      theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, size = 18),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size=18),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = here::here("figures",
                                 "pre_mod_figures",
                                 sp,
                                 paste(sp,
                                       "_corr_dendro.png",
                                       sep = "")),
           plot = dendro_graph, 
           width = 3660, 
           height = 3660, 
           units = "px", 
           dpi = "retina")
    
    # Computing vifstep to exclude collinear variables
    vifstep <- usdm::vifstep(sp_df[,4:ncol(sp_df)], 
                             method = "spearman", 
                             th = 5)
    
    vifstep_graph <-  ggplot() +
      geom_bar(data = vifstep@results,
               aes(x = VIF, 
                   y = reorder(Variables, VIF)),
               fill = "steelblue", color = "black", stat = "identity") +
      geom_vline(xintercept = 5, linewidth = 1, linetype = 2, color = "red") +
      theme_bw(base_size = 18) +
      scale_x_continuous(limits = c(0,10)) +
      labs(x = "VIF values", y = "Variables", 
           title = paste("VIFstep results for ",
                         sp,
                         " occurrences", sep="")) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    
    ggsave(vifstep_graph, 
           filename = here::here("figures",
                                 "pre_mod_figures",
                                 sp,
                                 paste(sp, 
                                       "_vifstep_plot.png")),
           width = 3660, height = 2940, units = "px", dpi = "retina")
    
    # Computing vifcor to exclude correlated variables left
    vifcor <- usdm::vifcor(sp_df[,vifstep@results$Variables], 
                           method = "spearman", 
                           th = 0.7)
    
    vifcor_graph <-  ggplot() +
      geom_bar(data = vifcor@results,
               aes(x =VIF, 
                   y = reorder(Variables, VIF)),
               fill = "steelblue", color = "black", stat = "identity") +
      geom_vline(xintercept = 5, linewidth = 1, linetype = 2, color = "red") +
      theme_bw(base_size = 18) +
      scale_x_continuous(limits = c(0,10)) +
      labs(x = "VIF values", y = "Variables", 
           title = paste("VIFcor results for ",
                         sp,
                         " occurrences", sep="")) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    
    ggsave(vifcor_graph, 
           filename = here::here("figures",
                                 "pre_mod_figures",
                                 sp,
                                 paste(sp, 
                                       "_vifcor_plot.png")), 
           width = 3660, height = 2940, units = "px", dpi = "retina")
    
    # Sending selected environmental variables as output
    env_grid <- env_data[,vifcor@results$Variables]
    
    write.csv(x = env_data[,vifcor@results$Variables],
              file = here::here("outputs",
                           "pre_mod_outputs",
                           sp,
                           paste(sp,
                                 "_selected_env.csv",
                                 sep = "")),
                row.names = F)
    
    # Stocking data into a list
    env_grid <- env_data[,vifcor@results$Variables]
    names(env_grid) <- paste(names(env_grid), sp, sep = "_")
    
    env_list[[grep(pattern = sp, x = list_sp)]] <- env_grid

  }
  
  return(env_list)
}
