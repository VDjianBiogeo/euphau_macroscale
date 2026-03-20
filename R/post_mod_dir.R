post_mod_dir <- function(sp_list){
  
  dir.create(here::here("outputs",
                        "post_mod_outputs"))
  
  dir.create(here::here("figures",
                        "post_mod_figures"))
  
  for(sp in sp_list){
    dir.create(here::here("figures",
                          "post_mod_figures",
                          sp))
    
    dir.create(here::here("outputs",
                          "post_mod_outputs",
                          sp))
  }
  return("done")
}