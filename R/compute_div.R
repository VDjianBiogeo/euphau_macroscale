#' @description
#' This function uses a list of species to extract predicted distributions and 
#' compute several biodiversity indices, such as species richness and species
#' diversity.
#'
#' @param binary_proj a dataframe containing species distribution in a binary
#' format (0/1).
#' @param prob_proj a dataframe containing species distribution as probability
#' of occurrence.
#' 
#' @export
#' 

compute_div <- function(binary_proj, prob_proj){
  
 # Computing species richness for each cell
  sp_rich <- rowSums(x = binary_proj[,3:ncol(binary_proj)]) %>%
    as.data.frame()
  
  sp_rich <- data.frame(x = binary_proj$x, 
                        y = binary_proj$y,
                        sp_rich = sp_rich$.)
  
  terra::writeRaster(x = terra::rast(sp_rich,
                                     type = "xyz",
                                     crs = "epsg:4326"),
                     here::here("outputs",
                                "post_mod_outputs",
                                "euphau_sp_rich.tif"),
                     overwrite = T)
  
  write.csv(sp_rich,
            here::here("outputs",
                       "post_mod_outputs",
                       "euphau_sp_richness.csv"),
            row.names = F)
  
  # Plotting predicted species richness
  data(World, package = "tmap")
  
  ## Worldwide
  sp_rich_world <- ggplot() +
    geom_sf(data = World) +
    geom_raster(aes(x = x, y = y, fill = sp_rich),
                data = sp_rich) +
    scale_fill_viridis_c(name = "Predicted species richness") + 
    coord_sf(expand = F, crs = 4326) +
    theme_bw(base_size = 18) +
    theme(axis.title = element_blank(), 
          legend.direction = "horizontal", 
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position="top", 
                                  title.hjust = 0.5, 
                                  theme = theme(legend.key.height  = unit(0.5, "cm"),
                                                legend.key.width = unit(10, "cm"))))
  
  # Saving species richness world map
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_sp_rich_world_map.png"),
         plot = sp_rich_world,
         width = 3800,
         height = 2400,
         units = "px",
         dpi = "retina")
  
  ## At subantarctic scale
  sp_rich_subant <- ggplot() +
    geom_sf(data = World) +
    geom_raster(aes(x = x, y = y, fill = sp_rich),
                data = sp_rich) +
    scale_fill_viridis_c(name = "Predicted species richness") + 
    coord_sf(expand = F, crs = 4326,
             xlim = c(-20,160), ylim = c(-20, -60)) +
    theme_bw(base_size = 18) +
    theme(axis.title = element_blank(), 
          legend.direction = "horizontal", 
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position="top", 
                                  title.hjust = 0.5, 
                                  theme = theme(legend.key.height  = unit(0.5, "cm"),
                                                legend.key.width = unit(10, "cm"))))
  
  saveRDS(object = sp_rich_subant, 
          file = here::here("figures",
                            "post_mod_figures",
                            "euphau_sp_rich_subant_map.rds"))
  
  # Saving species richness world map
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_sp_rich_subant_map.png"),
         plot = sp_rich_subant,
         width = 3800,
         height = 2400,
         units = "px",
         dpi = "retina")
  
  # Computing species diversity
  sp_shannon <- vegan::diversity(x = prob_proj[, 3:ncol(prob_proj)],
                   index = "shannon") %>%
    as.data.frame()
  
  sp_shannon <- data.frame(x = binary_proj$x, 
                        y = binary_proj$y,
                        sp_shannon = sp_shannon$.)
  
  write.csv(sp_shannon,
            here::here("outputs",
                       "post_mod_outputs",
                       "amphi_sp_shannon.csv"),
            row.names = F)
  
  # Plotting predicted species diversity (Shannon index)
  
  ## Worldwide
  sp_shannon_world <- ggplot() +
    geom_sf(data = World) +
    geom_raster(aes(x = x, y = y, fill = sp_shannon),
                data = sp_shannon) +
    scale_fill_viridis_c(name = "Predicted Shannon's index") + 
    coord_sf(expand = F, crs = 4326) +
    theme_bw(base_size = 18) +
    theme(axis.title = element_blank(), 
          legend.direction = "horizontal", 
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position="top", 
                                  title.hjust = 0.5, 
                                  theme = theme(legend.key.height  = unit(0.5, "cm"),
                                                legend.key.width = unit(10, "cm"))))
  
  # Saving species shannon world map
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_sp_shannon_world_map.png"),
         plot = sp_shannon_world,
         width = 3800,
         height = 2400,
         units = "px",
         dpi = "retina")
  
  ## At subantarctic scale
  sp_shannon_subant <- ggplot() +
    geom_sf(data = World) +
    geom_raster(aes(x = x, y = y, fill = sp_shannon),
                data = sp_shannon) +
    scale_fill_viridis_c(name = "Predicted Shannon's index") + 
    coord_sf(expand = F, crs = 4326,
             xlim = c(-20,160), ylim = c(-20, -60)) +
    theme_bw(base_size = 18) +
    theme(axis.title = element_blank(), 
          legend.direction = "horizontal", 
          legend.position = "bottom") +
    guides(fill = guide_colourbar(title.position="top", 
                                  title.hjust = 0.5, 
                                  theme = theme(legend.key.height  = unit(0.5, "cm"),
                                                legend.key.width = unit(10, "cm"))))
  
  # Saving species richness world map
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_sp_shannon_subant_map.png"),
         plot = sp_shannon_subant,
         width = 3800,
         height = 2400,
         units = "px",
         dpi = "retina")
  
  return(list(sp_rich, sp_shannon))
  
}