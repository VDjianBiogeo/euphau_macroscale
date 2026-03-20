#' @description
#' This function uses data on predicted occurrences and environmental layers
#' to plot density ridgelines of species occurrences vs all environmental
#' descriptors.
#'
#' @param corr_pca a dataframe species coordinates in the PCA space with their 
#' associated cluster.
#' @param binary_and_env a dataframe containing all environmental descriptors 
#' and predicted occurrences for each species
#' 
#' @export


plotting_env_ranges <- function(binary_and_env, corr_pca) {
  
  # Grouping all species in one column
  pivot_df <- binary_and_env |>
    pivot_longer(cols = 3:31, names_to = "species", values_to = "occurrences")
  
  # Only keeping occurrences and not predicted absences
  pivot_df <- pivot_df[pivot_df$occurrences == 1,]
  
  # Taking out the "_binary" at the end of species names
  pivot_df$species <- str_split(pivot_df$species,
                                "_",
                                simplify = T)[,1]
  
  # Filtering E.superba out 
  pivot_df <- pivot_df[pivot_df$species != "E.superba",]
  
  # Same with H.atlantica
  pivot_df <- pivot_df[pivot_df$species != "H.atlantica",]
  
  # Extracting all environmental descriptors
  env_descriptors <- colnames(pivot_df)[3:20]
  
  # Add cluster values
  pivot_df$cluster <- corr_pca$hac[match(pivot_df$species,
                                         corr_pca$Species)]
  library(ggridges)
  
  # Plotting loop
  dir.create(here::here("figures",
                        "post_mod_figures",
                        "ridgelines"))
  
  for(descr in env_descriptors){
    
    plot_df <- pivot_df[, c("x","y","species", "occurrences", "cluster", descr)] |>
      as.data.frame() |>
      na.omit()
    
    # Plotting ridgelines
    ridgelines <- ggplot() +
      geom_density_ridges2(aes(x = plot_df[,6], y = reorder(species, y, median),
                               
                               fill = cluster),scale = 5, alpha = 0.7,
                           data = na.omit(plot_df), show.legend = F) +
      theme_bw(base_size = 18) + 
      theme(axis.text.y = element_text(face = "italic"),
            axis.title.y = element_blank()) +
      labs(x = descr) +
      scale_fill_manual(name = "HAC cluster", 
                        values = RColorBrewer::brewer.pal(n = 3, 
                                                          name = "Set2"))
    
    ggsave(here::here("figures",
                      "post_mod_figures",
                      "ridgelines",
                      paste(descr,
                            "ridgelines.png",
                            sep = "_")), ridgelines, 
           height = 3660,
           width = 5660,
           units = "px",
           dpi = "retina")
    
    saveRDS(object = ridgelines,
            here::here("figures",
                       "post_mod_figures",
                       "ridgelines",
                       paste(descr,
                             "ridgelines.rds",
                             sep = "_")))
  }
  return(pivot_df)
}