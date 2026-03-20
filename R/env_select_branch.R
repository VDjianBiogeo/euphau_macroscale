#' @description
#' This function takes a filtered occurrence dataset, joins it with
#' an environmental dataframe based on spatial coordinates and computes
#' a series of operations in order to select which environmental variables will
#' be kept. 
#' @param filtered_occ A dataframe containing species occurrences.
#' @param env_data A dataset containing enviromental values.
#' @param coord_names names of the columns matching spatial coordinates.
#' 
#' @export
#' 

env_select_branch <- function(filtered_occ, env_data, coord_names){
  
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
  
  # Creating a repository by species
  sp <- unique(filtered_occ$species)
  
  dir.create(path = here::here("figures",
                               "pre_mod_figures",
                               sp))
  
  dir.create(path = here::here("outputs",
                               "pre_mod_outputs",
                               sp))
  
  # Joining occurrence and environment data
  sp_df <- inner_join(filtered_occ,
                      env_data[,-grep(pattern = "northv|eastv|po4|no3", 
                                      colnames(env_data))],
                      by = c("x", "y")) |>
    na.omit() %>%
    as.data.frame()
  
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
  corMatrix <- Hmisc::rcorr(x = as.matrix(sp_df[,-c(1:5)]), type = "pearson")
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
  vifstep <- usdm::vifstep(as.data.frame(sp_df[,colnames(distMatrix)]), 
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
                                     "_vifstep_plot.png",
                                     sep = "")),
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
                                     "_vifcor_plot.png",
                                     sep = "")), 
         width = 3660, height = 2940, units = "px", dpi = "retina")
  
  # Sending selected environmental variables as output
  env_grid <- env_data[,c("x","y",vifcor@results$Variables)]
  
  write.csv(x = env_grid,
            file = here::here("outputs",
                              "pre_mod_outputs",
                              sp,
                              paste(sp,
                                    "_selected_env.csv",
                                    sep = "")),
            row.names = F)
  
  final_df <- dplyr::inner_join(filtered_occ,
                                env_grid,
                                by = c("x","y"))
  
  write.csv(final_df,
            here::here("outputs",
                       "pre_mod_outputs",
                       sp,
                       paste(sp,
                             "_occ_and_env.csv",
                             sep = "")),
            row.names = F)
  
  return(final_df)
}