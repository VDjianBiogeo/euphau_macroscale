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
  
  # Importing species list
  tar_target(name = sp_list,
             command = readRDS(here::here("outputs",
                                          "explor_outputs",
                                          "species_lists.rds"))),
  
  # Importing occurrence data
  tar_target(name = bio_occ,
             command = read.csv(here::here("data",
                                           "treated",
                                           "bio",
                                           "euphau_0.25_stack.csv"))),
  
  # Importing environmental data
  tar_target(name = env_data,
             command = read.csv(file = here::here("data",
                                                  "treated",
                                                  "env",
                                                  "env_euphau_dataset.csv"))),
  
  # Filtering species occurrence according to prevalence criteria
  tar_target(name = filtered_occ,
             command = species_filtering(dataset = bio_occ, 
                                         sp_list = sp_list, 
                                         xy = T, 
                                         xy_name = c("x","y"))),
  
  # Branching targets by species (parallelisation)
  tar_group_by(name = grouped_occ,
                command = na.omit(tidyr::pivot_longer(data = filtered_occ, 
                                           cols = 3:ncol(filtered_occ), 
                                           names_to = "species",
                                           values_to = "occurrence")) |>
                 as.data.frame(),
               species),
  
  # Writing longer occurrence dataframe as a csv file
  tar_target(name = write_grouped_occ,
             command = write.csv(grouped_occ, 
                                 here::here("outputs",
                                            "pre_mod_outputs",
                                            "euphau_longer_occ.csv"),
                                 row.names = F)),
  
  # Selecting environmental data for each species based on collinearity
  tar_target(name = bio_and_env_data,
           command = env_select_branch(filtered_occ = grouped_occ, 
                                       env_data = env_data,
                                       coord_names = c("x","y")), 
           pattern = map(grouped_occ), 
           iteration = "list"),
  
  # Optimising modelling options to use
  tar_target(name = mod_prep,
             command = model_tuning(bio_and_env_data),
             pattern = map(bio_and_env_data)),

  # Optimising modelling options to use for parsimonious modelling
  tar_target(name = parsim_prep,
             command = parsim_tuning(bio_and_env_data),
             pattern = map(bio_and_env_data)),

  # Optimising modelling options to use for environmental partitioning
  tar_target(name = env_part_prep,
             command = partition_tuning(bio_and_env_data),
             pattern = map(bio_and_env_data))
  
)
