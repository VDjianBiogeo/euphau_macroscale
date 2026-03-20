#' @description
#' This function uses a dataframe containing correlations and runs a CA. It then 
#' produces several plots such as CA graphs and boxplots for mean semi-partial
#' correlation coefficients.
#'
#' @param corr_table a dataframe containing correlations.
#' @param sp_range_optimums a dataframe containing all environmental ranges and
#' optimums for selected species.
#' 
#' @export
#' 

make_pca_on_corr <- function(corr_table, sp_range_optimums){
  
  # Selecting one of the two dataframes
  corr <- corr_table[[1]]
  
  # Taking out E.superba
  corr <- corr[-grep("E.superba", 
                     rownames(corr)),]
  corr <- corr[-grep("E.crystallorophias",
                     rownames(corr)),]
  corr$species <- str_split(rownames(corr), "_", simplify = T)[,1]
  
  # Idem
  sp_range_optimums <- sp_range_optimums[sp_range_optimums$species != "E.superba",]
  
  # Joining with environmental ranges and optimums
  corr <- inner_join(x = corr,
                     y = sp_range_optimums,
                     by = c("species"))
  
  rownames(corr) <- corr$species
  
  # dist <- corr_table[[2]]
  
  # Running CA on semi-partial correlation data
  pca_corr <- FactoMineR::PCA(X = corr[,-4], 
                              quanti.sup = 4:ncol(corr[,-4]))
  
  # Extracting species coordinates
  coord_sp <- as.data.frame(pca_corr$ind$coord)
  coord_sp$Species <- str_split(rownames(coord_sp), 
                                "_", 
                                simplify = T)[,1]
  coord_model <- as.data.frame(pca_corr$var$coord)
  
  # Extracting variable coordinates
  coord_var <- as.data.frame(pca_corr$var$coord)
  
  # Extracting supplementary variables coordinates
  coord_supp <- as.data.frame(pca_corr$quanti.sup$coord)
  
  # PCA correlation circle
  pca_circle <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    geom_segment(aes(x = 0, y = 0, xend = `Dim.1`, yend = `Dim.2`),
                 linetype = 2, data = coord_supp,
                 arrow = arrow(length = unit(0.3, "cm"))) +
    ggrepel::geom_label_repel(aes(x = `Dim.1`, 
                                  y = `Dim.2`, 
                                  label = rownames(coord_supp)),
                              size = 6,
                              data = coord_supp, fontface = "italic") +
    geom_segment(aes(x = 0,
                     y = 0, 
                     xend = `Dim.1`, 
                     yend = `Dim.2`), 
                 linewidth = 0.7,
                 data = coord_var, 
                 arrow = arrow(length = unit(0.3, "cm"))) +
    geom_label(aes(x = `Dim.1`, 
                   y = `Dim.2`, 
                   label = rownames(coord_var)),
               size = 6,
               data = coord_var, fontface = "bold") +
    xlab(paste("Dim 1 (", 
               round(pca_corr$eig[1,2], 
                     digits = 2), 
               "%)", 
               sep="")) +
    ylab(paste("Dim 2 (", 
               round(pca_corr$eig[2,2], 
                     digits = 2), 
               "%)", 
               sep="")) +
    theme_bw(base_size = 18) +
    xlim(c(-1,1)) +
    ylim(c(-1,1))
  
  # Save as RDS for later use
  saveRDS(pca_circle, 
          file = here::here("figures",
                            "post_mod_figures",
                            "euphau_pca_circle.rds"))
  # PCA Circle plot
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_pca_circle.png"),
         plot = pca_circle,
         width = 4660,
         height = 3660,
         units = "px",
         dpi = "retina")
  
  # PCA graphic without HAC
  pca_graph <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    geom_point(aes(x = `Dim.1`, y = `Dim.2`), 
               fill = "#fec44f",
               pch = 23, color = "black",
               data = coord_sp, size = 5) +
    ggrepel::geom_text_repel(aes(x = `Dim.1`, y = `Dim.2`, label = Species, 
                                 fontface = 3), size = 5,
                             data = coord_sp) +
    # geom_point(aes(x = `Dim 1`, y = `Dim 2`), 
    #            pch = 23, color = "black", fill = "red",
    #            data = coord_model, size = 3) +
    # geom_label(aes(x = `Dim 1`, y = `Dim 2`, label = rownames(coord_model),
    #                fontface = 2), size = 4.5, 
    #            data = coord_model) +
    # scale_x_continuous(limits = c(-1, 0.6)) +
    xlab(paste("Dim 1 (", round(pca_corr$eig[1,2], digits = 2), "%)", sep="")) +
    ylab(paste("Dim 2 (", round(pca_corr$eig[2,2], digits = 2), "%)", sep="")) +
    theme_bw(base_size = 18)
  
  # Saving PCA plot
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_pca_corr.png"),
         plot = pca_graph,
         width = 4660,
         height = 3660,
         units = "px",
         dpi = "retina")
  
  # Clustering on PCA coordinates
  cluster <- FactoMineR::HCPC(res = pca_corr, nb.clust = -1)
  
  # Joining cluster groups to PCA coordinates
  coord_sp$hac <- cluster$data.clust$clust[match(rownames(coord_sp), 
                                                 rownames(cluster$data.clust))]
  
  # PCA graphic with HAC
  pca_graph <- ggplot() +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") +
    geom_hline(yintercept = 0, linetype = 2, colour = "black") +
    geom_point(aes(x = `Dim.1`, y = `Dim.2`, fill = factor(hac)), 
               pch = 23, color = "black",
               data = coord_sp, size = 6) +
    scale_fill_manual(name = "HAC cluster", 
                      values = RColorBrewer::brewer.pal(n = 3, 
                                                        name = "Set2")) +
    ggrepel::geom_text_repel(aes(x = `Dim.1`, y = `Dim.2`, 
                                 label = Species, 
                                 fontface = 3), size = 6,
                             data = coord_sp) +
    # geom_point(aes(x = `Dim 1`, y = `Dim 2`), 
    #            pch = 23, color = "black", fill = "red",
    #            data = coord_model, size = 3) +
    # geom_label(aes(x = `Dim 1`, y = `Dim 2`, label = rownames(coord_model),
    #                fontface = 2), size = 4.5, 
    #            data = coord_model) +
    # scale_x_continuous(limits = c(-1, 0.6)) +
    xlab(paste("Dim 1 (", round(pca_corr$eig[1,2], digits = 2), "%)", sep="")) +
    ylab(paste("Dim 2 (", round(pca_corr$eig[2,2], digits = 2), "%)", sep="")) +
    theme_bw(base_size = 18)
  
  # Save as RDS for later use
  saveRDS(pca_graph, 
          file = here::here("figures",
                            "post_mod_figures",
                            "euphau_pca_corr_hac.rds"))
  
  # Saving PCA plot
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_pca_corr_hac.png"),
         plot = pca_graph,
         width = 4660,
         height = 3660,
         units = "px",
         dpi = "retina")
  
  # Comparing clusters in terms of mean correlations
  
  library(rstatix)
  corr$hac <- coord_sp$hac
  
  ## For water mass
  # Kruskal test
  kruskal_results <- kruskal.test(corr$water_mass ~ corr$hac)
  
  # Psot-hoc Dunn test
  stat.test <- corr %>%
    dunn_test(water_mass ~ hac, p.adjust.method = "bonferroni") %>%
    add_significance() %>% 
    add_xy_position(x = "hac")
  
  # Boxplots with added p-values
  water_bxp <- ggplot() +
    geom_boxplot(aes(x = factor(hac), 
                     y = water_mass, 
                     fill = factor(hac)), 
                 data = corr) +
    ggpubr::stat_pvalue_manual(stat.test, 
                               label = "p.adj.signif", size = 8) +
    scale_fill_manual(name = "HAC cluster", 
                      values = RColorBrewer::brewer.pal(n = length(unique(coord_sp$hac)), 
                                                        name = "Set2")) +
    theme_bw(base_size = 24) +
    theme(axis.title.x = element_blank()) +
    labs(y = "Mean semi-partial correlation coefficient")
  
  # Saving graph 
  ggsave(filename = here::here("figures", 
                               "post_mod_figures",
                               "water_mass_mean_boxplot.png"),
         plot = water_bxp,
         width = 4000,
         height = 3660,
         units = "px",
         dpi = "retina")
  
  ## For ocean dynamics
  # Kruskal test
  kruskal_results <- kruskal.test(corr$ocean_dynamics ~ corr$hac)
  
  # Psot-hoc Dunn test
  stat.test <- corr %>%
    dunn_test(ocean_dynamics ~ hac, p.adjust.method = "bonferroni") %>%
    add_significance() %>% 
    add_xy_position(x = "hac")
  
  # Boxplots with added p-values
  dynamics_bxp <- ggplot() +
    geom_boxplot(aes(x = factor(hac), 
                     y = ocean_dynamics, 
                     fill = factor(hac)), 
                 data = corr) +
    ggpubr::stat_pvalue_manual(stat.test, 
                               label = "p.adj.signif", size = 8) +
    scale_fill_manual(name = "HAC cluster", 
                      values = RColorBrewer::brewer.pal(n = length(unique(coord_sp$hac)), 
                                                        name = "Set2")) +
    theme_bw(base_size = 24) +
    theme(axis.title.x = element_blank()) +
    labs(y = "Mean semi-partial correlation coefficient")
  
  # Saving graph 
  ggsave(filename = here::here("figures", 
                               "post_mod_figures",
                               "ocean_dynamics_mean_boxplot.png"),
         plot = dynamics_bxp,
         width = 4000,
         height = 3660,
         units = "px",
         dpi = "retina")
  
  ## For biological production
  # Kruskal test
  kruskal_results <- kruskal.test(corr$biologic_production ~ corr$hac)
  
  # Psot-hoc Dunn test
  stat.test <- corr %>%
    dunn_test(biologic_production ~ hac, p.adjust.method = "bonferroni") %>%
    add_significance() %>% 
    add_xy_position(x = "hac")
  
  # Boxplots with added p-values
  bio_bxp <- ggplot() +
    geom_boxplot(aes(x = factor(hac), 
                     y = biologic_production, 
                     fill = factor(hac)), 
                 data = corr) +
    ggpubr::stat_pvalue_manual(stat.test, 
                               label = "p.adj.signif", size = 8) +
    scale_fill_manual(name = "HAC cluster", 
                      values = RColorBrewer::brewer.pal(n = length(unique(coord_sp$hac)), 
                                                        name = "Set2")) +
    theme_bw(base_size = 24) +
    theme(axis.title.x = element_blank()) +
    labs(y = "Mean semi-partial correlation coefficient")
  
  # Saving graph 
  ggsave(filename = here::here("figures", 
                               "post_mod_figures",
                               "biologic_production_mean_boxplot.png"),
         plot = bio_bxp,
         width = 4000,
         height = 3660,
         units = "px",
         dpi = "retina")
  
  # Combining plots into one 
  library(patchwork)
  
  combined_bxp <- water_bxp + 
    (dynamics_bxp + theme(axis.title = element_blank())) + 
    (bio_bxp + theme(axis.title = element_blank())) +
    plot_annotation(tag_levels = "a")
  
  # Saving combined graph
  ggsave(filename = here::here("figures",
                               "post_mod_figures",
                               "euphau_bxp_mean_corr.png"),
         plot = combined_bxp,
         height = 3660,
         width = 10000,
         units = "px",
         dpi = "retina")
  
  return(coord_sp)

}