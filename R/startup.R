# startup.R
# START HERE TO CONFIGURE THE PIPELINE
# The startup module contains code for setting up the pipeline.
# 1. It will create the /output_folder/ to hold all the outputs of the pipeline
# 2. It will create the critical /output_folder/metadata.csv file that is needed to 
  # link sample_order to ptid, visit information
# 3. It will create the /output_folder/pipeline.json that will be used to compile
# 4. It will attempt to make extraction .json files that will be used extract 
  # relevant events from batches of .fcs files
require(dplyr)
source('R/start.R')
source('R/extract.R')
# USER PROVIDES CRITICAL MINIMAL SET OF INPUTS
# USER MUST IDENTIFY fcm08 file, .fcs folder paths, and FJ.xml files describing 
# each .fcs folder paths
# SPEC 0: All files and directories must be subfolders of the base dir
# SPEC 1: all .fcs file for a single batch must be in a single folder. 
# SPEC 2: No other .fcs files may be in that folder
# SPEC 3: One and only one xml_file ending "FJ.xml" must be in the same directory as the xml files
# SPEC 4: The parent_gate must exist in fcm08 file 
# SPEC 5: Specify a set of functional makers based on careful inspection of the fcm08 raw file.

params = list()
params[["run_name"]]        = "TEST"
params[["output_path"]]     =  NULL # Will be created if none if specified
params[["base_dir"]]        = '/fh/fast/gilbert_p/fg_data/VTN137/ics/facs_raw/'
params[["fcm08_path"]]      = '/fh/fast/gilbert_p/fg_data/VTN137/ics/fcm08_raw/vtn137_AP51_fcm08_FH_AP51_20220907.txt'
params[["fcs_folder_path"]] = "/fh/fast/gilbert_p/fg_data/VTN137/ics/facs_raw/3264-B-HVTN137/"
params[["xml_file_path"]]   = file.path(params[["fcs_folder_path"]], list.files(params[["fcs_folder_path"]], pattern = "FJ.xml")[[1]])
params[["parent_gate"]]     = 'Time/S/Lv/K1/K2/K3/K4/K5/K6/K7/K8/14-/S/L/19-/3+/3+excl 16br/56-16-/4+'
params[["stim_exclusion_terms"]] = c( "phactrl",
                                  "sebctrl",
                                  "posctrl")
params[["functional_markers"]] = c("IFNg+","IL2+","TNFa+","154+","IL4+","IL5_OR_IL13","IL17a+")
params[['xml_keywords']]       = c( "$FIL","Stim","Sample Order","EXPERIMENT NAME","Replicate")

# <marker_io> 
  # This package relies on marker io to reconcile multiple possible names 
  # for each marker. A standard marker is provided in the output, with 
  # variable input. The positive/negative status of the marker is ignored
marker_io = readxl::read_xlsx("marker_mapping.xlsx")$output
names(marker_io) = readxl::read_xlsx("marker_mapping.xlsx")$input
marker_io

# Step 1: Create a folder for the run
# If no output path is created on will be created with a random date/string/run_name
if (is.null(params$output_path)){
  # Get the current date
  current_date <- Sys.Date()
  # Generate a random -character hexadecimal string
  num_chars   <- 8
  random_hex  <- paste(sample(c(0:9, LETTERS[1:6]), num_chars, replace = TRUE), collapse = "")
  # Combine the date and random hex string
  output_path <- file.path( getwd(), paste0(current_date, "_", random_hex, "_", params$run_name))
  run_name    <- params$run_name
  dir.create(output_path)
}else{
  dir.create(params$output_path)
  output_path = params$output_path
  run_name = params$run_name
}

# Step 2: Get metadata, write metadata
metadata = generate_metadata(fcm08_filepath       = params$fcm08_path,
                             parent_gate          = params$parent_gate,
                             stim_exclusion_terms = params$stim_exclusion_terms)
configs$metadata_path = file.path(configs$output_path, "metadata.csv")
metadata %>% readr::write_csv(configs$metadata_path)
batches = metadata %>% pull( experiment_name) %>% unique() %>% sort()

# Step 3: get all the markers used in parent, these will likely be excluded 
# from downstream clustering
parent_markers = get_within_parent_gate(parent_gate = params$parent_gate) %>% 
  purrr::map_chr(., ~get_or_default(dict = marker_io, key = .x))
# Step 4: get all child elements
all_child_elements = get_child_paths(parent_gate = params$parent_gate) %>% 
  purrr::map_chr(., ~extract_last_element(.x)) 
# Step 5:
# Check that all functional markers are in elements
stopifnot(params$functional_markers %in% all_child_elements)

# get only the single items markers not complex combinations
all_child_markers_raw = all_child_elements[ purrr::map_lgl(all_child_elements, ~is_single_marker(.x)) ]
all_child_markers = purrr::map_chr(all_child_markers_raw, ~get_or_default(dict = marker_io, key = .x))

# OUT pipeline.json is based on configurations <configs>

configs = list()
configs$run_name                 = params$run_name
configs$parent_gate              = params$parent_gate
configs$output_path              = output_path
configs$batches                  = batches
configs$parent_markers           = parent_markers 
configs$functional_markers       = params$functional_markers
configs$continuous_markers_fcm08 = all_child_markers


message(sprintf("%s Loading a Gating Set from %s", format(Sys.time(), "%X"), params$fcs_folder_path))
messsage("This step will likely take 10-20 minutes, got get coffee!")
my_gs = get_gs(xml_path        = params$xml_path,
               fcs_folder_path = params$fcs_folder_path,
               xml_keywords    = params$xml_keywords)
fi_channels <- flowWorkspace::pData(flowCore::parameters(flowWorkspace::gh_pop_get_data(my_gs[[1]]))[ ,c("name", "desc")])
fi_channel_short = fi_channels$desc
fi_channel_long = paste(seq(1,length(fi_channels$name)),fi_channels$name, fi_channels$desc, sep = "|")
ix = (fi_channel_short %in% configs$continuous_markers_fcm08) & (!fi_channel_short %in% configs$parent_markers)

configs$all_channels_desc     = fi_channel_short
configs$all_channels          = fi_channel_long
configs$cluster_channels      = fi_channel_long[ix] 
configs

json_string = jsonlite::toJSON(configs, pretty = TRUE, auto_unbox = TRUE)
print(json_string)
write(json_string, file = file.path(configs$output_path, 'pipeline.json'))




for (batch in batches){
  base_dir = params[["base_dir"]] 
  batch_dir  = output_path
  batch_json = paste0("extract_",gsub(batch, pattern = " ", replacement = "_"), ".json")
  batch_path = file.path(output_path, paste0("extract_",gsub(batch, pattern = " ", replacement = "_"), ".json"))
  
  bconfigs[["json_dir"]]  = batch_dir
  bconfigs[["json_filename"]] = paste0("extract_",gsub(batch, pattern = " ", replacement = "_"), ".json")
  bconfigs[["batch_name"]] = batch
  bconfigs[["fcs_folder_path"]] = file.path(base_dir, batch)
  bconfigs[["xml_path"]] = file.path(base_dir, batch, paste0(batch, "FJ.xml"))
  bconfigs[["output_path"]] = output_path
  bconfigs[["parent_gate"]] = params$parent_gate
  bconfigs[["write_csv"]] = FALSE
  bconfigs[["write_h5"]] = TRUE
  bconfigs[["xml_keywords"]] = params$xml_keywords 
  bconfigs[["functional_markers"]] = params$functional_markers
  bconfigs[["stim_exclusion_terms"]] = params$stim_exclusion_terms
}

  




"IFNg+",
"IL2+",
"TNFa+",
"154+",
"IL4+",
"IL5_OR_IL13",
"IL17a+"
]
"stim_exclusion_terms"]] [
  "phactrl",
  "sebctrl",
  "posctrl"
]
}


















# my_result = extract_events(g                  = my_gs[[1]],
#                            sample_name        = sample_name,
#                            parent_gate        = params$parent_gate,
#                            markers            = params$functional_markers,
#                            functional_markers = params$functional_markers,
#                            experiment_name    = params$experiment_name, 
#                            sample_order       = params$sample_order, 
#                            replicate          = params$replicate, 
#                            stim               = params$stim,
#                            name               = params$stim)
# 
# 
# 
# 
# 
# # 1. generate metadata
# metadata = generate_metadata() 
# print(jsonlite::toJSON(metadata))
# 
# 
# data$fcs$sample_order <- as.character(data$fcs$sample_order)
# metadata$sample_order <- as.character(metadata$sample_order)
# 
# data$fcs %>% select(experiment_name, sample_order, stim) %>% 
#   group_by(experiment_name, sample_order, stim) %>% 
#   slice(1) %>% 
#   dim()
# 
# data$fcs %>% left_join(metadata, by = c("experiment_name", "sample_order", "stim")) %>% 
#   filter(is.na(ptid)) %>% 
#   group_by(experiment_name, stim, sample_order, sample_name) %>% 
#   slice(1)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# generate_pos_map <-function(){
#   
# }
# 
# generaate_batch_list <- function(){
#   
# }
# 
# generate_markers <- function(){
#   
# }