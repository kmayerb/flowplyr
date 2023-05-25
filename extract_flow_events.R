# May 24, 2023
# commandline version of extrator
source("R/extract.R")

parser <- argparser::arg_parser(
  description='extract_flow_events - extract events from FlowJo to .hdf5 format')
parser <- argparser::add_argument(parser, arg="--params", type="character",default = NA, help = "path to json.param")
args   <- argparser::parse_args(parser)
print(args)
# Get the JSON params file path from the command line
params_file <- args$params
# Check if the file exists
if (!file.exists(params_file)) {
  cat("Error: The specified JSON params file does not exist.\n")
  quit(status = 1)
}
# Read and parse the JSON params file
params <- jsonlite::fromJSON(params_file)
print(params)



# First, assemble the gating set using xml file, fcs folder and xml_keywords
# NOTE: this step take up to a few minutes depending on batch and filesize
message(sprintf("%s Loading a Gating Set from %s", format(Sys.time(), "%X"), params$fcs_folder_path))
my_gs = get_gs(xml_path        = params$xml_path,
               fcs_folder_path = params$fcs_folder_path,
               xml_keywords    = params$xml_keywords)
message(sprintf("%s Loaded a Gating Set with %s Samples", format(Sys.time(), "%X"), length(my_gs)))

# Next we can analyze the associated marker tree
message(sprintf("%s Gating Paths:", format(Sys.time(), "%X")))
message(paste0("\n",flowWorkspace::gs_get_pop_paths(my_gs)))
#
message(sprintf("%s Makers:", format(Sys.time(), "%X")))
print(params$markers)
message(sprintf("%s Key Makers:", format(Sys.time(), "%X")))
print(params$key_markers)


# Loop through all of the samples in the gated_set:
# <store> will hold each samples results
message(sprintf("%s Loaded a Gating Set with %s Samples", format(Sys.time(), "%X"), length(my_gs)))

store = list()
# Loop through each gated entry in the gated set
for (i in 1:length(my_gs)){
  pd = flowWorkspace::pData(my_gs[[i]])
  exp_name <- pd$`EXPERIMENT NAME`
  fcs_name <- paste(pd$`EXPERIMENT NAME`, pd$"Sample Order", pd$Replicate, pd$Stim, pd$name, sep = "|")

  my_result = extract_events(g = my_gs[[i]],
                             parent_gate = params$parent_gate,
                             markers = params$markers,
                             key_markers = params$key_markers)
  message(sprintf("(%s of %s) %s\t%s\t%s", i, length(my_gs), format(Sys.time(), "%X"), params$batch_name, fcs_name))
  # each result is a list() with 3 data entries ('pos', fi', and 'fcs_index')
  store[[i]] = my_result
}

# Concatenate the 'pos', 'fi' and 'fcs_index' data.frames from stored list
# single data.frames.
pos_ = do.call(rbind, purrr::map(store, ~.x$pos))*1
fi_  = do.call(rbind, purrr::map(store, ~.x$fi))
fcs_ = do.call(rbind, purrr::map(store, ~.x$fcs_index))
cols_fi = colnames(fi_)
cols_pos = colnames(pos_) 
# Section on output

# Optionally write .csv files
if (params$write_csv){
  message(sprintf("Writing output to .csv files"))
  write.csv(file = file.path(params$output_path, paste0(params$batch_name, "_pos.csv")),
            x = pos_,
            row.names = FALSE)
  write.csv(file = file.path(params$output_path, paste0(params$batch_name, "_fi.csv")),
            x = fi_,
            row.names = FALSE)
  write.csv(file = file.path(params$output_path, paste0(params$batch_name, "_fcs_index.csv")),
            x = fcs_,
            row.names = FALSE)
}

# Optionally write, hdf5 files
if (params$write_h5){
  message(sprintf("Writing output to .h5 files"))
  h5_name =  file.path(params$output_path, paste0(params$batch_name, ".h5"))
  if (!dir.exists(params$output_path)){dir.create(params$output_path)}
  if (file.exists(h5_name)){file.remove(h5_name)}
  rhdf5::h5createFile(h5_name)
  rhdf5::h5createGroup(h5_name ,"data")
  rhdf5::h5write(pos_, h5_name, "data/pos")
  rhdf5::h5write(fi_,  h5_name, "data/fi")
  rhdf5::h5write(fcs_, h5_name, "data/fcs_index")
  rhdf5::h5write(cols_pos, h5_name, "data/cols_pos")
  rhdf5::h5write(cols_fi, h5_name, "data/cols_fi")
  
  rhdf5::h5closeAll()
}

