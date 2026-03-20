#' @description
#' This function takes the output of the turn_binary() function and computes
#' a network clustering algorithm in order to determine ecoregions based on
#' species' predicted distribution. 
#' @param sp.stack Either a raster stack of the species' predicted distributions
#' or a dataframe resulting from a raster stack.
#' 
#' @export
#' 

ecoregionalisation <- function(sp.stack){
  
  # Preparing return of the function
  returns <- list()
  
  # Creating output directory ----
  dir.create(here::here("outputs", "post_mod_outputs", "ecoregions"))
  
  # Creating figures directory ----
  dir.create(here::here("figures", "post_mod_figures", "ecoregions"))
  
  # Changing colnames to take out the "_binary" part
  colnames(sp.stack)[3:ncol(sp.stack)] <- str_split(colnames(sp.stack)[3:ncol(sp.stack)],
                                                    "_",
                                                    simplify = T)[,1]
  
  # Worldwide scale regionalisation ----
  ## Creating worldwide scale output directory ----
  dir.create(here::here("outputs", 
                        "post_mod_outputs",
                        "ecoregions", 
                        "worldwide"))
  
  ## Creating worldwide scale figures directory ----
  dir.create(here::here("figures", 
                        "post_mod_figures", 
                        "ecoregions", 
                        "worldwide"))
  
  ## Formating data ----
  data_to_reg <- pivot_longer(sp.stack, cols = 3:ncol(sp.stack), 
                              names_to = "species", 
                              values_to = "occurrence") %>%
    filter(occurrence == 1)
  
  ## Defining cell id ----
  sp.xy <- unique(data_to_reg[, c("x", "y")])
  sp.xy$cells <- rownames(sp.xy)
  
  ## Adding cell id to data to regionalise ----
  sp.biogeo <- inner_join(sp.xy, data_to_reg, by = c("x", "y"))
  sp.biogeo <- unique(sp.biogeo[, c("cells", "species")])
  
  ## Ecoregionalisation ----
  bioregion::install_binaries()
  if(file.exists(here::here("outputs",
                            "post_mod_outputs",
                            "ecoregions",
                            "worldwide",
                            "worldwide_network.rds")) == F){
    
    network <- bioregion::netclu_infomap(net = as.data.frame(sp.biogeo), 
                                         numtrials = 1000, 
                                         bipartite = T,
                                         site_col = "cells",
                                         species_col = "species", 
                                         return_node_type = "both",
                                         weight = F, 
                                         show_hierarchy = T)
    
    # Saving output
    saveRDS(network, 
            here::here("outputs",
                       "post_mod_outputs",
                       "ecoregions",
                       "worldwide",
                       "worldwide_network.rds"))
  }
  network <- readRDS(here::here("outputs",
                                "post_mod_outputs",
                                "ecoregions",
                                "worldwide",
                                "worldwide_network.rds"))
  
  returns[["worldwide"]] <- network 
  
  ## Extracting informations from the network ----
  ### Species informations ----
  # Extracting species
  species <- network$clusters[which(network$clusters$ID %in% 
                                      sp.biogeo$species),]
  
  # Saving species information
  write.csv(species, 
            here::here("outputs",
                       "post_mod_outputs",
                       "ecoregions",
                       "worldwide",
                       "worldwide_reg_species.csv"),
            row.names = F)
  
  ### Cells' informations ----
  # Extracting cells
  sites <- network$clusters[which(network$clusters$ID %in% 
                                    sp.biogeo$cells),]
  
  # Adding spatial coordinates
  sites <- inner_join(sites, 
                      sp.xy, by = c("ID" = "cells"))
  
  # Keeping only clusters tied to species (eliminate really small clusters)
  sites <- sites[sites[,2] %in% unique(species[,2]),]
  
  # Adding colors
  colors <- data.frame(clusters = unique(sites[,2]),
                       color = RColorBrewer::brewer.pal(n = length(unique(sites[,2])), 
                                          name = "Paired"))
  sites$colors <- colors$color[match(sites[,2], colors$clusters)]
  
  # Saving cells informations
  write.csv(sites, 
            here::here("outputs", 
                       "post_mod_outputs",
                       "ecoregions",
                       "worldwide",
                       "worldwide_reg_cells.csv"), 
            row.names = F)
  
  # ## Mapping resulting ecoregions ----
  # # Importing continents
  # data("World", package = "tmap")
  # 
  # ### Worldwide scale ----
  # world_map <- ggplot() +
  #   geom_raster(aes(x = x, y = y, fill = factor(sites[,2])), data = sites) +
  #   geom_sf(data = World) +
  #   scale_fill_manual(values = unique(sites$colors[order(sites[,2])])) +
  #   coord_sf(expand = F, crs = 4326) +
  #   labs(fill = "Clusters") +
  #   theme_bw(base_size = 18) +
  #   theme(axis.title = element_blank())
  # 
  # ggsave(plot = world_map,
  #        filename = here::here("figures",
  #                              "post_mod_figures",
  #                              "ecoregions",
  #                              "worldwide",
  #                              "worldwide_ecoreg.jpeg"),
  #        width = 3660, height = 2000, units = "px", dpi = "retina")
  # 
  # ### Regional scale ----
  # reg_map <- ggplot() +
  #   geom_raster(aes(x = x, y = y, fill = factor(sites[,2])), data = sites) +
  #   geom_sf(data = World) +
  #   scale_fill_manual(values = unique(sites$colors[order(sites[,2])])) +
  #   coord_sf(expand = F, crs = 4326) +
  #   labs(fill = "Clusters") +
  #   theme_bw(base_size = 18) +
  #   theme(axis.title = element_blank())
  # 
  # ggsave(plot = world_map,
  #        filename = here::here("figures",
  #                              "post_mod_figures",
  #                              "ecoregions",
  #                              "worldwide",
  #                              "worldwide_ecoreg_cropped.jpeg"),
  #        width = 3660, height = 2000, units = "px", dpi = "retina")
  
  # Regional scale regionalisation ----
  
  ## Creating regional scale output directory ----
  dir.create(here::here("outputs", 
                        "post_mod_outputs", 
                        "ecoregions", 
                        "regional"))
  
  ## Creating regional scale figure directory ----
  dir.create(here::here("figures", 
                        "post_mod_figures", 
                        "ecoregions", 
                        "regional"))
  
  ## Cropping data to regional extent ----
  data_to_reg <- terra::crop(x = terra::rast(sp.stack, type = "xyz"),
                             y = terra::rast(xmin = -20,
                                             xmax = 160,
                                             ymin = -60,
                                             ymax = -20,
                                             res = 0.25), snap = "near") %>%
    as.data.frame(xy = T)
  
  ## Formating data ----
  data_to_reg <- pivot_longer(data_to_reg, cols = 3:ncol(data_to_reg), 
                              names_to = "species", 
                              values_to = "occurrence") %>%
    filter(occurrence == 1)
  
  ## Defining cell id ----
  sp.xy <- unique(data_to_reg[, c("x", "y")])
  sp.xy$cells <- rownames(sp.xy)
  
  ## Adding cell id to data to regionalise ----
  sp.biogeo <- inner_join(sp.xy, data_to_reg, by = c("x", "y"))
  sp.biogeo <- unique(sp.biogeo[, c("cells", "species")])
  
  ## Ecoregionalisation ----
  if(file.exists(here::here("outputs",
                            "post_mod_outputs",
                            "ecoregions",
                            "regional",
                            "regional_network.rds")) == F){
    
    network <- bioregion::netclu_infomap(net = as.data.frame(sp.biogeo), 
                                         numtrials = 1000, 
                                         bipartite = T,
                                         site_col = "cells",
                                         species_col = "species", 
                                         return_node_type = "both",
                                         weight = F,
                                         show_hierarchy = T)
    
    # Saving output
    saveRDS(network, 
            here::here("outputs",
                       "post_mod_outputs",
                       "ecoregions",
                       "regional",
                       "regional_network.rds"))
  }
  network <- readRDS(here::here("outputs",
                                "post_mod_outputs",
                                "ecoregions",
                                "regional",
                                "regional_network.rds"))
  returns[["regional"]] <- network
  
  ## Extracting informations from the network ----
  ### Species informations ----
  # Extracting species
  species <- network$clusters[which(network$clusters$ID %in% 
                                      sp.biogeo$species),]
  
  # Saving species information
  write.csv(species, 
            here::here("outputs",
                       "post_mod_outputs",
                       "ecoregions",
                       "regional",
                       "regional_reg_species.csv"),
            row.names = F)
  
  ### Cells' informations ----
  # Extracting cells
  sites <- network$clusters[which(network$clusters$ID %in% 
                                    sp.biogeo$cells),]
  
  # Adding spatial coordinates
  sites <- inner_join(sites, 
                      sp.xy, by = c("ID" = "cells"))
  
  # Keeping only clusters tied to species (eliminate really small clusters)
  sites <- sites[sites[,2] %in% unique(species[,2]),]
  
  # Adding colors
  colors <- data.frame(clusters = unique(sites[,2]),
                       color = RColorBrewer::brewer.pal(n = length(unique(sites[,2])), 
                                          name = "Paired"))
  sites$colors <- colors$color[match(sites[,2], colors$clusters)]
  
  # Saving cells informations
  write.csv(sites, 
            here::here("outputs", 
                       "post_mod_outputs",
                       "ecoregions",
                       "regional",
                       "regional_reg_cells.csv"), 
            row.names = F)
  
  ## Mapping resulting ecoregions ----
  ### Regional scale ----
  # world_map <- ggplot() +
  #   geom_raster(aes(x = x, y = y, fill = factor(sites[,2])), data = sites) +
  #   geom_sf(data = World) +
  #   scale_fill_manual(values = unique(sites$colors[order(sites[,2])])) +
  #   coord_sf(expand = F, crs = 4326) +
  #   labs(fill = "Clusters") +
  #   theme_bw(base_size = 18) +
  #   theme(axis.title = element_blank())
  # 
  # ggsave(plot = world_map,
  #        filename = here::here("figures",
  #                              "post_mod_figures",
  #                              "ecoregions",
  #                              "regional",
  #                              "regional_ecoreg.jpeg"),
  #        width = 3660, height = 2000, units = "px", dpi = "retina")

  
  return(returns)

}
  