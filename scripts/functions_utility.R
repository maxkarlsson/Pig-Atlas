
savepath <- 
  function(savename) { 
    wd <- getwd()
    result_folder <- paste0(str_extract(wd, ".*Pig Atlas/"), "results/", Sys.Date())
    
    dir.create(result_folder, showWarnings = FALSE)
    
    paste0(result_folder, "/", savename)
  }
