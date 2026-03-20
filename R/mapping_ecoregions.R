#' @description
#' This function maps the resulting ecoregions from the ecoregionalisation
#' function.
#' @param scales a vector indicating which (or both) of worldwide and regional scales
#' should be mapped. Scales possible are "regional" and "worldwide".
#' 
#' @export

mapping_ecoregions <- function(ecoregions){
  data("World", package = "tmap")
  
  scales <- names(ecoregions)
  
  if(length(scales) > 1){
    for(scale in scales){
      sites <- read.csv(here::here("outputs", 
                                   "post_mod_outputs",
                                   "ecoregions",
                                   scale,
                                   paste(scale, 
                                         "_reg_cells.csv", 
                                         sep = "")))
      if(scale == "worldwide"){
        ### Worldwide scale ----
        world_map <- ggplot() +
          geom_raster(aes(x = x, y = y, fill = factor(sites[,2])), data = sites) +
          geom_sf(data = World) +
          scale_fill_manual(values = unique(sites$colors[order(sites[,2])])) +
          coord_sf(expand = F, crs = 4326) +
          labs(fill = "Clusters") +
          theme_bw(base_size = 18) +
          theme(axis.title = element_blank())
        
        saveRDS(object = world_map,
                file = here::here("figures",
                                  "post_mod_figures",
                                  "ecoregions",
                                  "worldwide",
                                  "worldwide_ecoreg.rds"))
        
        ggsave(plot = world_map,
               filename = here::here("figures",
                                     "post_mod_figures",
                                     "ecoregions",
                                     "worldwide",
                                     "worldwide_ecoreg.jpeg"),
               width = 3660, height = 2000, units = "px", dpi = "retina")
        
        ### Regional scale ----
        reg_map <- ggplot() +
          geom_raster(aes(x = x, y = y, fill = factor(sites[,2])), data = sites) +
          geom_sf(data = World) +
          scale_fill_manual(values = unique(sites$colors[order(sites[,2])])) +
          coord_sf(expand = F, crs = 4326, xlim = c(-20,160), ylim = c(-20,-60)) +
          labs(fill = "Clusters") +
          theme_bw(base_size = 18) +
          theme(axis.title = element_blank())
        
        saveRDS(object = reg_map,
                file = here::here("figures",
                                  "post_mod_figures",
                                  "ecoregions",
                                  "worldwide",
                                  "worldwide_ecoreg_cropped.rds"))
        
        ggsave(plot = reg_map,
               filename = here::here("figures",
                                     "post_mod_figures",
                                     "ecoregions",
                                     "worldwide",
                                     "worldwide_ecoreg_cropped.jpeg"),
               width = 3660, height = 2000, units = "px", dpi = "retina")
      }
      
      ### Regional scale ----
      world_map <- ggplot() +
        geom_raster(aes(x = x, y = y, fill = factor(sites[,2])), data = sites) +
        geom_sf(data = World) +
        scale_fill_manual(values = unique(sites$colors[order(sites[,2])])) +
        coord_sf(expand = F, crs = 4326, xlim = c(-20,160), ylim = c(-20,-60)) +
        labs(fill = "Clusters") +
        theme_bw(base_size = 18) +
        theme(axis.title = element_blank())
      
      saveRDS(object = world_map,
              file = here::here("figures",
                                "post_mod_figures",
                                "ecoregions",
                                "worldwide",
                                "regional_ecoreg.rds"))
      
      ggsave(plot = world_map,
             filename = here::here("figures",
                                   "post_mod_figures",
                                   "ecoregions",
                                   "regional",
                                   "regional_ecoreg.jpeg"),
             width = 3660, height = 2000, units = "px", dpi = "retina")
      
    }
  }
  
  return("done")
}