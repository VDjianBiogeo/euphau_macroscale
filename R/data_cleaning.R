#'@description
#'This function cleans a biodiversity database resulting from the previous
#'steps.
#'
#'@return
#'This function returns a dataframe.
#'
#'@export

data_cleaning <- function(raw_db){
  
  # Eliminate duplicates
  cleaned_db <- CoordinateCleaner::cc_dupl(x = raw_db, 
                                                 lon = "lon",
                                                 lat = "lat",
                                                 species = "Species")
  
  # Eliminate points with same lon/lat
  cleaned_db <- CoordinateCleaner::cc_equ(x = cleaned_db, 
                                                lon = "lon",
                                                lat = "lat")
  
  # Eliminate points with 0 at lon/lat
  cleaned_db <- CoordinateCleaner::cc_zero(x = cleaned_db,
                                                 lon = "lon",
                                                 lat = "lat")
  
  # Eliminate terrestrial points
  sea_occ <- CoordinateCleaner::cc_sea(x = cleaned_db,
                                       lon = "lon",
                                       lat = "lat", 
                                       value = "flagged")
  
  cleaned_db <- cleaned_db[sea_occ == F,]
  
  # Eliminate invalid coordinates
  cleaned_db <- CoordinateCleaner::cc_val(cleaned_db,
                                                lon = "lon",
                                                lat = "lat")
  
  # Eliminate outliers in species distribution
  cleaned_db <- CoordinateCleaner::cc_outl(cleaned_db,
                                           lon = "lon",
                                           lat = "lat",
                                           species = "Species",
                                           method = "quantile")
  
  # Eliminate occurrences before 1950
  cleaned_db <- cleaned_db[cleaned_db$year > 1950,]
  
  # Eliminate genus-only occurrences
  cleaned_db_genus <- cleaned_db[cleaned_db$Species %in% 
                                               str_split(pattern = " ", 
                                                         string = cleaned_db$Species, 
                                                         simplify = T)[,1],]
  cleaned_db <- cleaned_db[!cleaned_db$Species %in% 
                                         str_split(pattern = " ", 
                                                   string = cleaned_db$Species, 
                                                   simplify = T)[,1],]
  
  # Eliminate NA in species
  cleaned_db <- cleaned_db[!is.na(cleaned_db$Species) &
                                         !is.na(cleaned_db$year),]
  
  # Eliminate DNA based occurrences
  cleaned_db <- cleaned_db[-grep(pattern = "DNA|Sequence", 
                                             cleaned_db$organismquantitytype),]
  
  # Add ID to occurrences
  cleaned_db$id <- 1:nrow(cleaned_db)
  
  write.csv(x = cleaned_db, 
            file = here::here("data",
                              "treated",
                              "bio",
                              "cleaned_euphausiid_database.csv"),
            row.names = F)
  return(cleaned_db)
  
}