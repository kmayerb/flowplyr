# override.R
# I wrote this for cases when we might want 
# to provide a params file but also have the ability to over-ride 
# some parameters all over-rides must come as full flags --flag

json_params_with_overrides <- function(){
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  # Initialize an empty list to store the parameters
  raw_params <- list()
  # Loop through the arguments
  for (i in seq(1, length(args), by = 2)) {
    if (i < length(args)) {
      flag <- gsub("^--", "", args[i])
      raw_params[[flag]] <- args[i + 1]
    }
  }
  if ('params' %in% names(raw_params)){
    # Get the JSON params file path from the command line
    params_file <- raw_params$params
    # Check if the file exists
    if (!file.exists(params_file)) {
      cat("Error: The specified JSON params file does not exist.\n")
      quit(status = 1)
    }
    # load params
    params <- jsonlite::fromJSON(params_file)
    # override and suppliment with raw_params
    for (arg in names(raw_params)){
      print(arg)
      if (arg %in% names(params)){
        value = raw_params[[arg]]
        target_class = class(params[[arg]])
        #print(value)
        value <- switch(target_class,
                        "numeric" = as.numeric(value ),
                        "character" = as.character(value ),
                        "logical" = as.logical(value ),
                        stop("Unsupported class type"))
        #print(value)
        params[[arg]] = value
      } else {
        params[[arg]] = raw_params[[arg]]
      }
    }
  } else {
    # here is case where all params come as raw. I don't want to allow this
    params <- raw_params
    stop('params file not suppiled')
  }
  return(params)
}
