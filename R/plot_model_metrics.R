#' @description
#' This function uses a species list to extract model evaluation metrics and raw
#' occurrences in order to plot modelling evaluation metrics for all species modelled.
#'
#' @param sp_list a list of the modelled species.
#' 
#' @export

plot_model_metrics <- function(sp_list){
  
  # Creating a dataframe to stock AUC values for each species
  metrics <- data.frame(model_name = NA,
                        auc = NA,
                        species = NA)
  for(sp in sp_list){
    sp_metrics <- read.csv(here::here("outputs",
                                      "mod_outputs",
                                      sp,
                                      "simple",
                                      paste(sp,
                                            "_auc_metrics.csv",
                                            sep = "")))
    metrics <- rbind(metrics, sp_metrics) %>%
      na.omit()
    
  }
  
  # Extracting species occurrences
  occurrences <- read.csv(here::here("data",
                                     "treated",
                                     "bio",
                                     "euphau_0.25_stack.csv"))
  
  # Computing total occurrences per species
  occ_sums <- colSums(occurrences[,3:ncol(occurrences)], na.rm = T)
  
  # Formating species names
  names(occ_sums) <- str_split(names(occ_sums), "_", simplify = T)[,1]
  
  occ_sums <- as.data.frame(occ_sums)
  
  # Adding species raw occurrences to the metric dataframe
  metrics$occurrences <- occ_sums$occ_sums[match(metrics$species, 
                                                 rownames(occ_sums))]
  
  # Plotting evaluation metrics for all species
  plot <- ggplot() +
    geom_boxplot(aes(x = auc, y = reorder(species, auc, median), fill = occurrences),
                 data = metrics) +
    geom_vline(aes(xintercept = 0.7), linetype = 2, linewidth = 1) +
    scale_fill_distiller(palette ="Spectral", direction = -1, 
                         name = "Number of\noccurrences (log10)", 
                         transform = "log10") +
    # xlim(0.6, 0.75) + 
    labs(x = "AUC") +
    theme_bw(base_size = 18) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"),
          legend.direction = "horizontal", 
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position="top", 
                                  title.hjust = 0.5, 
                                  theme = theme(legend.key.height  = unit(0.5, "cm"),
                                                legend.key.width = unit(10, "cm"))))
  
  ggsave(filename = "./figures/post_mod_figures/euphau metrics with occurrence number.png", 
         plot = plot,
         height = 3660, 
         width = 4660, 
         units = "px", 
         dpi = "retina")
  
  saveRDS(plot, 
          "./figures/post_mod_figures/euphau metrics with occurrence number.rds")
  
  return(metrics)
}