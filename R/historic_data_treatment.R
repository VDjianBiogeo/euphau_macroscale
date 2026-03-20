#'@description
#'This function reads historical data file and gives in the right format for
#'later combining with other biological data from different sources.
#'@return
#'This function returns a dataframe with biological data and associated metadata
#'from historical surveys file.
#'@export

read_historic <- function(file_path){
  
  # Importing data
  historic <- read.csv(file_path, 
                       sep = ";")
  
  # Renaming colnames to the same species names from other files
  colnames(historic)[8:18] <- c("Euphausia_cristallorophias", 
                                "Euphausia_frigida", 
                                "Euphausia_longirostris", 
                                "Euphausia_similis", 
                                "Euphausia_superba", 
                                "Euphausia_triacantha", 
                                "Euphausia_vallentini", 
                                "Euphausia_sp", 
                                "Thysanoessa_gregaria", 
                                "Thysanoessa_macrura", 
                                "Euphausiacea")
  
  # Pivoting dataframe to have all species under one column
  historic_pivot <- pivot_longer(data = historic, cols = c(8:18), 
                                 names_to = "Species", 
                                 values_to = "Occurrence")
  
  # Adding a data source flag
  historic_pivot$source <- "Historic surveys"
  
  # Simplifying and decomposing temporal information
  historic_pivot$year <- str_split(string = historic_pivot$Date.prelevements, 
                                   pattern = "/", 
                                   simplify = T)[,3]
  historic_pivot$month <- str_split(string = historic_pivot$Date.prelevements, 
                                    pattern = "/", 
                                    simplify = T)[,2]
  historic_pivot$day <- str_split(string = historic_pivot$Date.prelevements, 
                                  pattern = "/", 
                                  simplify = T)[,1]
  
  # Completing missing columns
  historic_pivot$BasisOfRecords <- "Sampled_specimens"
  historic_pivot$bibliographicCitation <- ""
  historic_pivot$associatedReferences <- ""
  historic_pivot$OrganismQuantityType <- "Individuals"
  historic_pivot$samplingProtocol <- ""
  
  # Filtering columns to keep
  historic_pivot <- historic_pivot[historic_pivot$Occurrence == 1,
                                   c("longitude",
                                     "latitude",
                                     "Species",
                                     "Date.prelevements",
                                     "year",
                                     "month",
                                     "day",
                                     "BasisOfRecords",
                                     "campagne",
                                     "bibliographicCitation",
                                     "associatedReferences",
                                     "OrganismQuantityType",
                                     "samplingProtocol",
                                     "source")]
  
  # Changing species names
  historic_pivot$Species <- str_replace_all(string = historic_pivot$Species, 
                                            pattern = "_", 
                                            replacement = " ")
  
  # Changing column names for future merging
  colnames(historic_pivot) <- c("lon",
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
  
  return(historic_pivot)
  
}
