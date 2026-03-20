set.seed(4951)

# Import libraries ####
library(tidyverse)

# Importing all datasets ####

## OBIS ####
source(here::here("R", "obis_read.R"))

obis_db <- obis_read(here::here("data",
                                "raw",
                                "bio",
                                "euphau_raw_obis.csv"))

## GBIF ####
source(here::here("R", "gbif_read.R"))

gbif_db <- gbif_read(here::here("data",
                                'raw',
                                "bio",
                                "euphau_raw_gbif.txt"))

## Biodiversity.aq ####
biodiv_db <- gbif_read(here::here("data",
                                  'raw',
                                  'bio',
                                  'euphau_raw_biodiversity.txt'))

## Historic data ####
source(here::here("R", "historic_data_treatment.R"))
historic_db <- read_historic(here::here("data",
                                        "raw",
                                        "bio",
                                        "euphau_raw_historic.csv"))

## REPCCOAI survey data ####
source(here::here("R", "read_repccoai.R"))
repccoai_db <- read_repccoai(here::here("data",
                                        "raw",
                                        "bio",
                                        "euphau_repccoai.csv"))

# Combining all raw data ####
raw_db <- rbind(gbif_db,
                obis_db,
                biodiv_db,
                historic_db,
                repccoai_db) 

# Correcting some issues with coordinates
raw_db$lon <- as.numeric(raw_db$lon)
raw_db$lat <- as.numeric(raw_db$lat)

raw_db <- raw_db[raw_db$lon != 115000.0000,]

write.csv(x = raw_db, 
          file = here::here("data",
                            "raw",
                            "bio",
                            "raw_euphausiid_database.csv"),
          row.names = F)

# Mapping raw euphausiid database ####
source(here::here("R", "raw_mapping.R"))
raw_map <- raw_mapping(db_to_map = raw_db)

# Cleaning dataset ####
source(here::here("R", "data_cleaning.R"))
cleaned_db <- data_cleaning(raw_db)

# Mapping cleaned occurrences ####
source(here::here("R", "cleaned_mapping.R"))
cleaned_map <- cleaned_mapping(cleaned_db)

# Finding frequent species in our area ####


# Gridding biological data at different resolutions ####
source(here::here("R", "data_gridding.R"))
sp_stack <- data_gridding(cleaned_db, resolution = c(0.25, 0.5, 1))