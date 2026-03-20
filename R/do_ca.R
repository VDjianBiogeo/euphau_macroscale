#' @description
#' This function computes canonical analysis on a list of dataframes containing
#' occurrences data depending on a list of selected species.
#' 
#' @param data a list of species dataframes.
#' @param list_sp a list of species to keep for the CA in the given dataframe.
#' @export

do_ca <- function(data, list_sp){
  library(ggplot2)
  
  for(i in 1:length(data)){
    
    # Importing species list
    sp_to_keep <- list_sp[[i]]
    
    #Importing species dataframe
    df <- data[[i]]
  
    # Changing NA to 0
    df[is.na(df)] <- 0
    
    # Changing rownames in order to "name" sites
    rownames(df) <- paste(df$x,
                          df$y,sep = ";")
    df$cell_names <- rownames(df)
    
    # Computing CA only on species to keep data
    ca <- FactoMineR::CA(X = df[,sp_to_keep], 
                         graph = F)
    
    # Extracting species and sites coordinates
    coord_row <- ca$row$coord |> as.data.frame()
    coord_row$cell_names <- rownames(coord_row)
    
    coord_sp <- ca$col$coord |> as.data.frame()
    
    # Graph CA
    ca_plot <- ggplot() +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_point(aes(x = `Dim 1`, 
                     y = `Dim 2`),
                 data = coord_row, pch = 23, 
                 fill = "red") +
      geom_label(aes(x = `Dim 1`,
                     y = `Dim 2`,
                     label = rownames(coord_sp)),
                 data = coord_sp) +
      labs(x = paste0("Axe 1 (",
                     round(ca$eig[1,2], digits = 2),
                     "%)"),
           y = paste0("Axe 2 (",
                          round(ca$eig[2,2], digits = 2),
                          "%)")) +
      theme_bw(base_size = 18)
    
    # Checking if output directory exists
    
    if(dir.exists(here::here("figures",
                             "explor_figures")) == F){
      dir.create(here::here("figures",
                            "explor_figures"))
    }
    
    # Exporting graph
    ggsave(filename = here::here("figures",
                                 "explor_figures",
                                 paste0("ca_graph_",
                                        i,
                                        ".png")), 
           plot = ca_plot,
           units = "px", 
           height = 2000, 
           width = 3660, 
           dpi = "retina")
    
    # Creating a spatial CA coord dataframe
    map_df <- df[, c("x","y","cell_names")] |> 
      na.omit()
    
    map_df <- dplyr::inner_join(x = map_df, 
                                y = coord_row, 
                                by = "cell_names")
    
    # Mapping in space CA coordinates
    
    data("World", package = "tmap")
    
    # Map for first axis
    ca_map <- ggplot() +
      geom_sf(data = World) +
      geom_raster(aes(x = x,
                      y = y,
                      fill = `Dim 1`),
                  data = map_df) +
      scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 4,
                                                             name = "RdYlBu"))) +
      coord_sf(expand = F,
               crs = 4326) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank())
    
    # Export CA map of first axis
    ggsave(filename = here::here("figures",
                                 "explor_figures",
                                 paste0("first_axis_map_ca_",
                                        i,
                                        ".png")), 
           plot = ca_map, 
           dpi = "retina",
           height = 2000,
           width = 3660,
           units = "px")
    
    # Map for second axis
    ca_map <- ggplot() +
      geom_sf(data = World) +
      geom_raster(aes(x = x,
                      y = y,
                      fill = `Dim 2`),
                  data = map_df) +
      scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 4,
                                                              name = "RdYlBu"))) +
      coord_sf(expand = F,
               crs = 4326) +
      theme_bw(base_size = 18) +
      theme(axis.title = element_blank())
    
    # Export CA map of second axis
    ggsave(filename = here::here("figures",
                                 "explor_figures",
                                 paste0("scnd_axis_map_ca_",
                                        i,
                                        ".png")), 
           plot = ca_map, 
           dpi = "retina",
           height = 2000,
           width = 3660,
           units = "px")

  }
  return("done")
}