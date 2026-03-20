# Setup ----
set.seed(4951)
library(tidyverse)

# Raw data treatment ----
source(here::here("analyses", "raw_euphau_data_treatment.R"))

# Exploratory analyses ----
## Load package and function ----
targets::tar_config_set(store = here::here("outputs", 
                                           "expl_pipeline"),
                        script = here::here("analyses", 
                                            "exploratory_analyses.R"))

## Pre visualisation ----
targets::tar_visnetwork(targets_only = T)

## Launch pipeline ----
targets::tar_make()

## Post visualisation ----
targets::tar_visnetwork(targets_only = T)

# Pre-modelling pipeline ----
targets::tar_config_set(store = here::here("outputs",
                                           "pre_mod_pipeline"),
                        script = here::here("analyses",
                                            "pre_mod_analyses.R"))

## Pre visualisation ----
targets::tar_visnetwork(targets_only = T)

## Launch pipeline ----
targets::tar_make()

## Post visualisation ----
targets::tar_visnetwork(targets_only = T)

# Post modelling pipeline ----
targets::tar_config_set(store = here::here("outputs",
                                           "post_mod_pipeline"),
                        script = here::here("analyses",
                                            "post_mod_analyses.R"))

## Pre visualisation ----
targets::tar_visnetwork(targets_only = T)

## Launch pipeline ----
targets::tar_make()

## Post visualisation ----
targets::tar_visnetwork(targets_only = T)
