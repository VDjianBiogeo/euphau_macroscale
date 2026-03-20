library(targets)
library(tarchetypes)
library(crew)

# Configuration ----
# tar_option_set(controller = crew_controller_local(workers = 12, 
#                                                   seconds_idle = 60))
tar_source()
set.seed(4951)

# Defining pipeline ----
list(
  
  # Semi-correlations pipeline ----
 
  # Importing raw occurrences 
  tar_target(name = grouped_occ,
             command = read.csv(here::here("outputs", 
                                           "pre_mod_outputs", 
                                           "euphau_longer_occ.csv")),
             priority = 1),
  
  # Creating species list
  tar_target(name = sp_list,
             command = unique(grouped_occ$species[order(grouped_occ$species)]), 
             priority = 0.99),
  
  # Creating directories
  tar_target(name = dir_creation,
             command = post_mod_dir(sp_list), priority = 0.98),
  
  # Importing parsimonious raster layers
  tar_target(name = parsim_stack,
             command = read_parsim(sp_list, grouped_occ), priority = 0.97),
  
  # Importing environmental data
  tar_target(name = env_data,
             command = read.csv(file = here::here("data",
                                                  "treated",
                                                  "env",
                                                  "env_euphau_dataset.csv")),
             priority = 0.96),
  
  # Joining to environmental data
  tar_target(name = parsim_and_env,
             command = join_rasters(parsim_stack, env_data),
             priority = 0.95),
  
  # Plotting model metrics according to species raw occurrences
  tar_target(name = metric_plots,
             command = plot_model_metrics(sp_list)),
  
  # Joining binary distribution with environment
  tar_target(name = binary_and_env,
             command = join_rasters(binary_sp, env_data)),
  
  # Computing species environmental ranges and optimums
  tar_target(name = sp_ranges_optimums,
             command = compute_range_and_optimum(binary_and_env)),
  
  # Computing semi-correlations for each species
  tar_target(name = corr_table,
             command = semi_corr(sp_list),
             priority = 0.94),

  # Making a CA on semi-correlations data
  tar_target(name = corr_pca,
             command = make_pca_on_corr(corr_table, 
                                        sp_ranges_optimums),
             priority = 0.93),
  
  # Plotting environmental density ridges for all species
  tar_target(name = range_plots,
             command = plotting_env_ranges(binary_and_env,
                                           corr_pca)),
  
  # Computing diversity indices
  tar_target(name = div_indices,
             command = compute_div(binary_sp, parsim_stack),
             priority = 0.91),

  # Ecoregionalisation pipeline ----
  
  # Turning probabilities into binary data
  tar_target(name = binary_sp,
             command = turn_binary(parsim_stack, grouped_occ, sp_list),
             priority = 0.9),
  
  # Computing ecoregions
  tar_target(name = ecoregions,
             command = ecoregionalisation(binary_sp)),
  
  # Mapping ecoregions 
  tar_target(name = ecoreg_map,
             command = mapping_ecoregions(ecoregions)),
  
  # Computing biogeographic metrics for resulting ecoregions ----
  tar_target(name = ecoreg_metrics,
             command = biogeographic_metrics(ecoregions = ecoregions,
                                             binary_sp = binary_sp))

)
