#' @description
#' This function computes biogeographic metrics for resulting ecoregions.
#' @param ecoregions the output of the ecoregion() function.
#' @param binary_sp a raster stack containing all species' binary distributions.
#' 
#' @export
#' 

biogeographic_metrics <- function(ecoregions, binary_sp){
  
  # Getting regionalisation scales
  scales <- names(ecoregions)
  
  ecoreg_metrics <- list()
  
  for(scale in scales){
    # Importing binary distributions
    sp_df <- binary_sp
    
    # Changing species names in binary stack
    colnames(sp_df)[3:ncol(sp_df)] <- str_split(string = colnames(sp_df)[3:ncol(sp_df)],
                                                pattern = "_",
                                                simplify = T)[,1]
    
    # Importing sites data
    sites <- read.csv(here::here("outputs",
                                 "post_mod_outputs",
                                 "ecoregions",
                                 scale,
                                 paste(scale,
                                       "reg_cells.csv", 
                                       sep = "_")))
    colnames(sites)[2] <- "clusters"
    
    # Importing species data
    species <- read.csv(here::here("outputs",
                                   "post_mod_outputs",
                                   "ecoregions",
                                   scale,
                                   paste(scale,
                                         "reg_species.csv", 
                                         sep = "_")))
    
    colnames(species)[2] <- "clusters"
    
    # Joining binary layers to ecoregions' sites
    sp_df <- inner_join(sp_df, sites, by = c("x", "y"))
    
    # Preparing metrics dataframe
    metrics_all_sp <- data.frame(species = NA,
                                 affinity = NA,
                                 fidelity = NA,
                                 indval = NA,
                                 dilval = NA,
                                 cluster = NA)
    
    for(sp in unique(species$ID)){
      # Getting total species predicted occurrences
      total_occ <- length(sp_df[sp_df[,sp] == 1, sp])
      
      for(clust in unique(sites$clusters)){
        # Total number of cells from the cluster
        total_cells <- nrow(sp_df[sp_df$clusters == clust,])
        
        # Number of cluster's cells where the species is found compared to total number of cluster cells
        affinity <- nrow(sp_df[sp_df[,sp] == 1 & 
                                 sp_df$clusters == clust,])/total_cells
        
        # Number of cluster's cells where the species is found compared to total number of soecies occurrences
        fidelity <- nrow(sp_df[sp_df[,sp] == 1 & 
                                 sp_df$clusters == clust,])/total_occ
        
        # IndVal according to Legendre
        indval <- affinity * fidelity
        dilval <- affinity * (1 - fidelity)
        
        sp_metrics <- data.frame(species = sp,
                                 affinity = affinity,
                                 fidelity = fidelity,
                                 indval = indval,
                                 dilval = dilval,
                                 cluster = clust)
        
        metrics_all_sp <- rbind(metrics_all_sp,
                                sp_metrics) %>%
          na.omit()
      }
    }
    
    # Associating colors to clusters
    metrics_all_sp$color <- sites$colors[match(metrics_all_sp$cluster,
                                               sites$clusters)]
    ecoreg_metrics[[scale]] <- metrics_all_sp
    
    write.csv(metrics_all_sp,
              here::here("outputs",
                         "post_mod_outputs",
                         "ecoregions",
                         scale,
                         paste(scale,
                               "biogeo_metrics.csv", 
                               sep = "_")),
              row.names = F)
    
    # Plotting metrics
    ## Indval
    indval_plot <- ggplot() +
      geom_bar(aes(x = indval,
                   y = species,
                   fill = factor(cluster)), 
               color = "black", alpha = 0.8,
               metrics_all_sp, stat = "identity") +
      scale_fill_manual(name = "Clusters",
                        values = unique(metrics_all_sp$color[order(metrics_all_sp$cluster)])) +
      theme_bw(base_size = 18) +
      labs(y = "Species", x = "IndVal") +
      facet_wrap(~cluster) +
      scale_y_discrete(limits = rev) +
      theme(axis.text.y = element_text(size = 11, face = "italic"))
    
    saveRDS(object = indval_plot,
            file = here::here("figures",
                              "post_mod_figures",
                              "ecoregions",
                              scale,
                              paste(scale,
                                    "indval_plot.rds",
                                    sep = "_")))
    
    ggsave(here::here("figures",
                      "post_mod_figures",
                      "ecoregions",
                      scale,
                      paste(scale,
                            "indval_plot.jpeg",
                            sep = "_")), 
           indval_plot,
           width = 3660, height = 2940, units = "px", dpi = "retina")
    
    ## DilVal
    dilval_plot <- ggplot() +
      geom_bar(aes(x = dilval,
                   y = species,
                   fill = factor(cluster)), 
               color = "black", alpha = 0.8,
               metrics_all_sp, stat = "identity") +
      scale_fill_manual(name = "Clusters",
                        values = unique(metrics_all_sp$color[order(metrics_all_sp$cluster)])) +
      theme_bw(base_size = 18) +
      labs(y = "Species", x = "DilVal") +
      facet_wrap(~cluster) +
      scale_y_discrete(limits = rev) +
      theme(axis.text.y = element_text(size = 11, face = "italic"))
    
    saveRDS(object = dilval_plot,
            file = here::here("figures",
                              "post_mod_figures",
                              "ecoregions",
                              scale,
                              paste(scale,
                                    "dilval_plot.rds",
                                    sep = "_")))
    
    ggsave(here::here("figures",
                      "post_mod_figures",
                      "ecoregions",
                      scale,
                      paste(scale,
                            "dilval_plot.jpeg",
                            sep = "_")), 
           dilval_plot,
           width = 3660, height = 2940, units = "px", dpi = "retina")
    
  }
  return(ecoreg_metrics)
}