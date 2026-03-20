#'@description
#'This function reads a given OBIS csv file and gives in the right format for
#'later combining with other biological data from different sources.
#'@return
#'This function returns a dataframe with biological data and associated metadata
#'from an OBIS csv file.
#'@export

obis_read <- function(file_path){
  # Importing all datasets ####
  
  ## OBIS ####
  obis_db <- read.csv(file_path)
  
  # Adding a source flag
  obis_db$source <- "OBIS"
  
  # Resuming all year information into one column
  obis_db[["year"]][is.na(obis_db[["year"]])] <- obis_db[["date_year"]][is.na(obis_db[["year"]])]
  
  # Filtering flags and columns kept
  obis_db <- obis_db[obis_db$flags == "{}",
                     c("decimallongitude",
                       "decimallatitude",
                       "scientificname",
                       "eventdate",
                       "year",
                       "month",
                       "day",
                       "basisofrecord",
                       "datasetname",
                       "bibliographiccitation",
                       "associatedreferences",
                       "organismquantitytype",
                       "samplingprotocol",
                       "source")]
  
  # Changing colnames for later homogeneisation
  colnames(obis_db) <- c("lon",
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
  
  return(obis_db)
}