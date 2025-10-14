read_npz <- function(npz_path){
  # https://rstudio.github.io/reticulate/
  np <- import("numpy", convert=TRUE )
  py_npz <- np$load(npz_path, allow_pickle = TRUE)
  
  results <- list()
  for (key in py_npz$files){
    results[[key]] = py_npz$f[[key]]
  }
  
  return(results)
}

