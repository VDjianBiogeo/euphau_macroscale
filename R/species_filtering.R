#' @description
#' This function takes an occurrence dataset and filter it based on a list of 
#' species to keep.
#' @param sp_list A list of species to select in the given dataset.
#' @param dataset A dataset containing species occurrences as columns
#' @param xy Logical. Indicates if spatial coordinates are in the dataset.
#' @param xy_name Names of the columns corresponding to spatial coordinates (ex: c("x","y")).
#' 
#' @export
#' 

species_filtering <- function(sp_list, dataset, xy, xy_name){
  
  if(xy == T){
    df <- dataset[, c(xy_name, sp_list[[1]])]
  } else{
    df <- dataset[, sp_list[[1]]]
  }
  colnames(df)[3:ncol(df)] <- str_split(colnames(subset(x = df, 
                                            select = -c(x,y))), 
                            pattern = "_", 
                            simplify = T)[,1]
  return(df)
}