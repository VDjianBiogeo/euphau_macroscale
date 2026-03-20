library(targets)
library(tarchetypes)

# Configuration ----
tar_source()
set.seed(4951)

# Defining pipeline ----
list(
  
  tar_target(name = create_output_dir,
             command = create_dir(here::here("outputs",
                                             "explor_outputs")), 
             priority = 1),
  
  tar_target(name = create_figure_dir,
             command = create_dir(here::here("figures",
                                             "explor_figures")),
             priority = 0.99),
  
  ## Importing bio data
  tar_target(name = occ_list,
             command = readRDS(here::here("data",
                                          "treated",
                                          "bio",
                                          "euphau_stacks_several_res.rds"))),
  
  tar_target(name = list.sp,
             command = species_freq(occ_list)),
  
  tar_target(name = ca_analysis,
             command = do_ca(occ_list, list.sp)),
  
  tar_target(name = sp_rich,
             command = species_richness(occ_list, list.sp))
  
)
