#'@description
#'This function tests if the given directory exists or not, and if not,
#'creates it.
#'
#'@param path Path to the directory.
#'
#'@export

create_dir <- function(path){
  
  if(dir.exists(path) == F){
    
    dir.create(path)
  }
}