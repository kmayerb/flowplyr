validate_args_or_params_file<-function(req_args){
  # VALIDATE INPUT ARGUMENTS

  # If no params json file is provide check for required command line arguments
  if (is.na(args$params)){
    # Check that either the required arguments or a params file was provided
    for (req in req_args){
      if (is.na(args[[req]])){
        message(paste0("Error: Either --params or --", req, " is required.\n"))
        quit(status = 1)
      }
    }
    params = args
  }

  # Alternative, if a params .json is provided, check that it exists, and that
  # the required argument were provided
  if (!is.na(args$params)){
    # otherwise read from param file
    params_file <- args$params
    # Check if the file exists
    if (!file.exists(params_file)) {
      message("Error: The specified JSON params file does not exist.\n")
      quit(status = 1)
    }
    # Read and parse the JSON params file
    params <- jsonlite::fromJSON(params_file)
    for (req in req_args){
      if (is.na(params[[req]])){
        message(paste0("Error: Either", req, " is required in --params json input.\n"))
        quit(status = 1)
      }
    }
  }

  return(params)
}
