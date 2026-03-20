#'@description
#'This function reads REPCCOAI data file and gives in the right format for
#'later combining with other biological data from different sources.
#'@return
#'This function returns a dataframe with biological data and associated metadata
#'from REPCCOAI surveys file.
#'@export

read_repccoai <- function(file_path){
  # Importing REPCCOAI database
  repccoai_db <- read.csv(file_path,
                          sep = ";", dec = ".", fileEncoding = "latin1")
  
  # Changing var armata into E.similis
  repccoai_db$Espèce_abond[repccoai_db$Espèce_abond == "Euphausia similis v. armata"] <- "Euphausia similis"
  
  # Building dataframe with information to keep
  repccoai_db <- data.frame(lon = repccoai_db$lon,
                            lat = repccoai_db$lat,
                            Species = repccoai_db$Espèce_abond,
                            date = repccoai_db$Year,
                            BasisOfRecord = "HUMAN_OBSERVATION",
                            year = repccoai_db$Year,
                            month = NA,
                            day = NA,
                            dataset = "REPCCOAI program",
                            citation = NA,
                            associatedReferences = NA,
                            organismquantitytype = "Individuals",
                            samplingProtocol = "IKMT nets",
                            source = "REPCCOAI")
  
  # Changing "Nematoscelis" to "Hansarsia"
  taxon_names <- as.data.frame(str_split(string = repccoai_db$Species,
                                         pattern = " ",
                                         simplify = T))
  taxon_names$V1[taxon_names$V1 == "Nematoscelis"] <- "Hansarsia"
  taxon_names$V2[taxon_names$V2 == "boÃ¶pis"] <- "boopis"
  
  repccoai_db$Species <- paste(taxon_names$V1,
                               taxon_names$V2,
                               sep = " ")
  
  # Changing colnames for future merging
  colnames(repccoai_db) <- c("lon",
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
  
  return(repccoai_db)
}