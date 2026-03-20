#'@description
#'This function reads a given GBIF txt file and gives in the right format for
#'later combining with other biological data from different sources.
#'@return
#'This function returns a dataframe with biological data and associated metadata
#'from an GBIF txt file.
#'@export

gbif_read <- function(file_path){
  
  # Importing all datasets ####
  
  ## GBIF ####
  gbif_db <- read.delim(file_path, sep = "\t")
  
  # Filtering according to occurrence status to keep only presences
  gbif_db <- gbif_db[gbif_db$occurrenceStatus == "PRESENT",]
  
  # Adding a source flag
  gbif_db$source <- "GBIF"
  
  # Filtering which columns to keep
  gbif_db <- gbif_db[, c("decimalLongitude",
                         "decimalLatitude",
                         "species",
                         "eventDate",
                         "year",
                         "month",
                         "day",
                         "basisOfRecord",
                         "datasetName",
                         "bibliographicCitation",
                         "associatedReferences",
                         "organismQuantityType",
                         "samplingProtocol",
                         "source")]
  
  # Renaming kept columns for later merging
  colnames(gbif_db) <- c("lon",
                         "lat",
                         "Species",
                         "date",
                         "year",
                         "month",
                         "day",
                         "BasisOfRecord",
                         "dataset",
                         "citation",
                         "associatedreferences",
                         "organismquantitytype",
                         "samplingprotocol",
                         "source")
  
  return(gbif_db)
}