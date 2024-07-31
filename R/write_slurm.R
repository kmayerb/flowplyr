# Write slurm

generate_slurm_script <- function(filename, 
                                  job_name, 
                                  wd = 'cd /fh/fast/gilbert_p/fg_data/flowplyr/',
                                  ml = 'ml fhR/4.2.0-foss-2021b',
                                  command = "Rscript extract_flow_events",
                                  params = "--params /fh/fast/gilbert_p/fg_data/flowplyr/2023-11-15_5C717641_TEST_RUN/extract_3264-B-HVTN137.json",
                                  ntasks = 1, 
                                  time = "60:00", 
                                  mem = '16G') {
  # Create the script content
  script_content <- paste(
    "#!/bin/bash",
    paste("#SBATCH --job-name='", job_name,"'",sep=""),
    paste("#SBATCH --ntasks=", ntasks, sep=""),
    paste("#SBATCH --time=", time, sep=""),
    paste("#SBATCH --mem=", mem, sep=""),
    "", # empty line
    paste(wd),
    paste(ml),
    paste(command, params),
    sep="\n"
  )
  
  # Write the script to a file
  writeLines(script_content, con = filename)
  
  # Make the script executable
  system(paste("chmod +x", filename))
  
  return(paste("SLURM script", filename, "generated successfully."))
}

