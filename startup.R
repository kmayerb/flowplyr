# startup.R -- start here
cat('..d8888b...888.....................888.....................\n')
cat('d88P..Y88b.888.....................888.....................\n')
cat('Y88b.......888.....................888.....................\n')
cat('."Y888b....888888..8888b...888d888.888888.888..888.88888b..\n')
cat('...."Y88b..888........"88b.888P"...888....888..888.888."88b\n')
cat('......"888.888.....d888888.888.....888....888..888.888..888\n')
cat('Y88b..d88P.Y88b...888..888.888.....Y88b...Y88b.888.888.d88P\n')
cat('."Y8888P"..."Y888."Y888888.888......"Y888.."Y88888.88888P".\n')
cat('...................................................888.....\n')
cat('.....FLOWPLYR V.1..................................888.....\n')
cat('...................................................888.....\n')
testing_only = FALSE
# THIS IS HOW YOU RUN 
"
cd /fh/fast/gilbert_p/fg_data/flowplyr
ml fhR/4.2.0-foss-2021b
Rscript startup.R \
--run_name 'VTN137_Part_A_CD4_TESTRUN' \
--base_dir '/fh/fast/gilbert_p/fg_data/VTN137/ics/facs_raw/' \
--fcm08_path '/fh/fast/gilbert_p/fg_data/VTN137/ics/fcm08_raw/vtn137_AP51_fcm08_FH_AP51_20220907.txt' \
--fcs_folder_path '/fh/fast/gilbert_p/fg_data/VTN137/ics/facs_raw/3264-B-HVTN137/' \
--xml_file_path '/fh/fast/gilbert_p/fg_data/VTN137/ics/facs_raw/3264-B-HVTN137/3264-B-HVTN137 FJ.xml' \
--parent_gate 'Time/S/Lv/K1/K2/K3/K4/K5/K6/K7/K8/14-/S/L/19-/3+/3+excl 16br/56-16-/4+' \
--stim_exclusion_terms 'phactrl,sebctrl,posctrl' \
--functional_markers 'IFNg+,IL2+,TNFa+,154+,IL4+,IL5_OR_IL13,IL17a+' \
--xml_keywords '$FIL,Stim,Sample Order,EXPERIMENT NAME,Replicate'
"

parser <- argparser::arg_parser(
  description='startup - Create a project output folder and .json control files for the pipeline')
parser <- argparser::add_argument(parser, arg="--run_name", 
                                  type="character",
                                  default = "TEST_RUN",
                                  help = "run name -- user supplied informative naem for the run")
parser <- argparser::add_argument(parser, arg="--base_dir", 
                                  type="character",
                                  help = "Base directory -- all input batch folders must be sub directories")
parser <- argparser::add_argument(parser, arg="--output_path", 
                                  type="character",
                                  default = NA,
                                  help = "Full file path where all the run's parameters and output will be saved")
parser <- argparser::add_argument(parser, arg="--fcm08_path", 
                                  type="character",
                                  help = "fcm08 file path")
parser <- argparser::add_argument(parser, arg="--fcs_folder_path", 
                                  type="character",
                                  help = "fcs folder path (one representative for learning channels)")
parser <- argparser::add_argument(parser, arg="--xml_file_path", 
                                  type="character",
                                  help = "xml file path (one representative for learnign channels)")
parser <- argparser::add_argument(parser, arg="--parent_gate", 
                                  type="character",
                                  default = 'Time/S/Lv/K1/K2/K3/K4/K5/K6/K7/K8/14-/S/L/19-/3+/3+excl 16br/56-16-/4+',
                                  help = "String specifying the parent gate for this run")
parser <- argparser::add_argument(parser, arg="--stim_exclusion_terms", 
                                  type="character",
                                  default = "phactrl,sebctrl,posctrl",
                                  help = "Comma separated list of stims we want to ignore")
parser <- argparser::add_argument(parser, arg="--functional_markers", 
                                  type="character",
                                  default = "IFNg+,IL2+,TNFa+,154+,IL4+,IL5_OR_IL13,IL17a+",
                                  help = "")
parser <- argparser::add_argument(parser, arg="--xml_keywords", 
                                  type="character",
                                  default = "$FIL,Stim,Sample Order,EXPERIMENT NAME,Replicate",
                                  help = "Comma separated list of xml keywords")
params   <- argparser::parse_args(parser)
# Convert to "," sep strings to lists
params$stim_exclusion_terms = stringr::str_split(params$stim_exclusion_terms,",")[[1]]
params$xml_keywords       = stringr::str_split(params$xml_keywords, ",")[[1]]
params$functional_markers = stringr::str_split(params$functional_markers, ",")[[1]]
params$flowplyr_dir = 'fh/fast/gilbert_p/fg_data/flowplyr'
print(params)

require(dplyr)
require(stringr)
source('R/start.R')
source('R/extract.R')
source('R/write_slurm.R')

if (testing_only){
  params = list()
  params[["run_name"]]        = "TEST"
  params[["output_path"]]     =  NA # Will be created if none if specified
  params[["base_dir"]]        = '/fh/fast/gilbert_p/fg_data/VTN137/ics/facs_raw/'
  # These can be one representative of the full batch
  params[["fcm08_path"]]      = '/fh/fast/gilbert_p/fg_data/VTN137/ics/fcm08_raw/vtn137_AP51_fcm08_FH_AP51_20220907.txt'
  params[["fcs_folder_path"]] = "/fh/fast/gilbert_p/fg_data/VTN137/ics/facs_raw/3264-B-HVTN137/"
  params[["xml_file_path"]]   = file.path(params[["fcs_folder_path"]], list.files(params[["fcs_folder_path"]], pattern = "FJ.xml")[[1]])
  params[["parent_gate"]]     = 'Time/S/Lv/K1/K2/K3/K4/K5/K6/K7/K8/14-/S/L/19-/3+/3+excl 16br/56-16-/4+'
  params[["stim_exclusion_terms"]] = c( "phactrl",
                                        "sebctrl",
                                        "posctrl")
  # positivivity for one or more of these markers is needed to retain the event
  params[["functional_markers"]] = c("IFNg+","IL2+","TNFa+","154+","IL4+","IL5_OR_IL13","IL17a+")
  params[['xml_keywords']]       = c( "$FIL","Stim","Sample Order","EXPERIMENT NAME","Replicate")
  params$flowplyr_dir = 'fh/fast/gilbert_p/fg_data/flowplyr'
}

# <marker_io> 
  # This package relies on marker io to reconcile multiple possible names 
  # for each marker. A standard marker is provided in the output, with 
  # variable input. The positive/negative status of the marker is ignored
marker_io = readxl::read_xlsx("marker_mapping.xlsx")$output
names(marker_io) = readxl::read_xlsx("marker_mapping.xlsx")$input

# Step 1.

# If no output path is created on will be created with a random date/string/run_name
if (is.na(params$output_path)){
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

cat(paste0("\nStep 1:  Step 1: Create a folder for the run: ", params$run_name, "\n"))
cat(paste0("\t", output_path,"\n"))

# Step 2: Get metadata, write metadata
metadata = generate_metadata(fcm08_filepath       = params$fcm08_path,
                             parent_gate          = params$parent_gate,
                             stim_exclusion_terms = params$stim_exclusion_terms)
metadata_path = file.path(output_path, "metadata.csv")
metadata %>% readr::write_csv(metadata_path)
cat(paste0("\nStep 2: Generated metadata file from: ", 
           params$fcm08_path, "\n"))
cat(paste0("\n Wrote metadata file to: ", 
           metadata_path, "\n"))


# Step 3: 
cat(paste0("\nStep 3: Identify the number of batches in metadata file\n")) 
batches = metadata %>% pull( experiment_name) %>% unique() %>% sort()
for (i in 1:length(batches)){
  cat(paste0("\tBatch ", i, " -- ", batches[i], "\n"))
}

cat(paste0("\nStep 4: Identify the markers in the parent gate\n")) 
# Step 3: get all the markers used in parent, these will likely be excluded 
# from downstream clustering
parent_markers = get_within_parent_gate(parent_gate = params$parent_gate) %>% 
  purrr::map_chr(., ~get_or_default(dict = marker_io, key = .x))

cat(paste0("\nStep 5: Identify child markers in parent gate\n")) 
# Step 4: get all child elements
all_child_elements = get_child_paths(parent_gate = params$parent_gate) %>% 
  purrr::map_chr(., ~extract_last_element(.x)) 
# Step 5:
# Check that all functional markers are in elements
stopifnot(params$functional_markers %in% all_child_elements)
# get only the single items markers not complex combinations
all_child_markers_raw = all_child_elements[ purrr::map_lgl(all_child_elements, ~is_single_marker(.x)) ]
continuous_markers_fcm08 = purrr::map_chr(all_child_markers_raw, ~get_or_default(dict = marker_io, key = .x))

cat(paste0("\nStep 6: writing batch specific .json and .sh scripts\n")) 
# OUT pipeline.json is based on configurations <configs>
cat(sprintf("%s Loading a repressentaive gating set from %s\n", format(Sys.time(), "%X"), params$fcs_folder_path))
cat("\nThis step will likely take 10-20 minutes\n")
cat("Maybe, go get coffee\n")
my_gs = get_gs(xml_path        = params$xml_file_path,
               fcs_folder_path = params$fcs_folder_path,
               xml_keywords    = params$xml_keywords)
cat(sprintf("%s Finished loading", format(Sys.time(), "%X") ))

# CONDITIONALS 
g = my_gs[[1]]
# get all conditional gates that follow the parent_gate
all_conditionals = flowWorkspace::gh_get_leaf_nodes(g, params$parent_gate)
# extract the final part of the gate
all_conditionals = purrr::map_chr(all_conditionals, ~extract_last_element(.x))
# keep on the simple conditions (i.e. no more than one +/-/AND/OR)
all_simple_conditionals = all_conditionals[ purrr::map_lgl(all_conditionals, ~is_single_marker(.x)) ]
# Check that the conditionals in the fcm08 match those derived from .xml workspace
stopifnot(all(all_simple_conditionals %in% all_child_markers_raw))
# Check that all the functional markers are in the simple conditionals.
# If they are not the pipeline will fail at the extraction step.
stopifnot(all(params$functional_markers %in% all_simple_conditionals))


# From the representative workspace identify the 
# channels.
  # <fi_channel_short>
  # <fi_channel_long>
  # <ix> - index those that will be used for clustering, based on (i) being in
  # the child nodes of teh parent gate in the fcm08 (e.g.,. FSC and SSC
  # are channels but not relevant for clustering), and (ii) exclude those
  # markers aht
fi_channels <- flowWorkspace::pData(flowCore::parameters(flowWorkspace::gh_pop_get_data(my_gs[[1]]))[ ,c("name", "desc")])
fi_channel_short = fi_channels$desc
fi_channel_long = paste(seq(1,length(fi_channels$name)),fi_channels$name, fi_channels$desc, sep = "|")
# channel that is in fcm08 file but not a marker used in the parent gate
ix = (fi_channel_short %in% continuous_markers_fcm08) & (!fi_channel_short %in% parent_markers)

# At this 

for (batch in batches){
  base_dir = params[["base_dir"]] 
  batch_dir  = output_path
  batch_json = paste0("extract_",gsub(batch, pattern = " ", replacement = "_"), ".json")
  batch_path = file.path(output_path, paste0("extract_",gsub(batch, pattern = " ", replacement = "_"), ".json"))
  
  bconfigs = list()
  bconfigs[["json_dir"]]  = batch_dir
  bconfigs[["json_filename"]] = paste0("extract_",gsub(batch, pattern = " ", replacement = "_"), ".json")
  bconfigs[["batch_name"]] = gsub(batch, pattern = " ", replacement = "_")
  bconfigs[["fcs_folder_path"]] = file.path(base_dir, batch)
  bconfigs[["xml_path"]] = file.path(base_dir, batch, paste0(batch, " FJ.xml"))
  bconfigs[["output_path"]] = output_path
  bconfigs[["parent_gate"]] = params$parent_gate
  bconfigs[["write_csv"]] = FALSE
  bconfigs[["write_h5"]] = TRUE
  bconfigs[["xml_keywords"]] = params$xml_keywords 
  bconfigs[["functional_markers"]] = params$functional_markers
  bconfigs[["markers"]] = all_simple_conditionals
  bconfigs[["stim_exclusion_terms"]] = params$stim_exclusion_terms
  bjson_string = jsonlite::toJSON(bconfigs, pretty = TRUE, auto_unbox = TRUE)
  write(bjson_string, file =  batch_path)
  
  # Example usage of the function
  json_filepath = file.path(bconfigs[["json_dir"]], bconfigs[["json_filename"]])
  bash_filename = gsub(json_filepath, pattern = ".json", replacement = ".sh")
  
  # see R/write_slurm.R
  generate_slurm_script(filename =  bash_filename, 
                        job_name = paste0('extract_', batch), 
                        wd = paste('cd', params$flowplyr_dir),
                        ml = 'ml fhR/4.2.0-foss-2021b',
                        command = "Rscript extract_flow_events.R",
                        params = paste0("--params ", json_filepath),
                        ntasks = 1, 
                        time = "60:00", 
                        mem = '16G')
  
}

batches_no_spaces = purrr::map_chr(batches, ~gsub(.x, pattern = " ", replacement = "_"))

cat(paste0("\nStep 7: writing full pipeline.json\n")) 
# WRITING GLOBAL CONFIGS FOR OVERALL PIPELINE
configs = list()
configs$run_name                 = params$run_name
configs$metadata_file            = metadata_path
configs$parent_gate              = params$parent_gate
configs$output_path              = output_path
configs$batches                  = batches
configs$batches_no_spaces        = batches_no_spaces
configs$parent_markers           = parent_markers 
configs$functional_markers       = params$functional_markers # USER MUST SUPPLY
configs$markers                  = all_simple_conditionals
configs$continuous_markers_fcm08 = continuous_markers_fcm08
configs$all_channels_desc        = fi_channel_short
configs$all_channels             = fi_channel_long
configs$cluster_channels         = fi_channel_long[ix] 

# compilations based elements
configs$h5_output_path_all_events                     = file.path(configs$output_path, "batches_compiled_all_events.h5")
configs$h5_output_path_subset_events                  = file.path(configs$output_path, "batches_compiled_subset_events.h5")
configs$h5_output_path_subset_umap                    = file.path(configs$output_path, "batches_compiled_subset_umap.h5")
configs$h5_output_path_subset_events_clustering       = file.path(configs$output_path, "batches_compiled_subset_clustering.h5")
configs$h5_output_path_subset_events_annoy_clustering = file.path(configs$output_path, "batches_compiled_subset_annoy_clustering.h5")
configs$h5_output_path_subset_events_summary          = file.path(configs$output_path, "batches_compiled_subset_summary.h5")
configs$csv_output_path_subset_events_summary         = file.path(configs$output_path, "batches_compiled_subset_summary.csv")
configs$subset_frac          = 1
configs$subset_seed          = 1
configs$annoy_threads        = 4
configs$umap_n_neighors      = 31
configs$umap_learning_rate   = 0.25
configs$min_cell_per_cluster = 100
configs$cluster_name_string  = "Leiden "
configs$batch_list           = purrr::map_chr(configs$batches_no_spaces, ~file.path(configs$output_path, paste0(.x, ".h5")))
configs$ptid_visit_exclude   = c()
configs$group_vars           = c("batch","ptid","visit","stim","parent_ct")
configs$join_vars            = c("batch", "ptid", "visit","key")
configs$negative_stim_string = "negctrl"
# WRITE IT OUT
json_string = jsonlite::toJSON(configs, pretty = TRUE, auto_unbox = TRUE)
print(json_string)
write(json_string, file = file.path(configs$output_path, 'pipeline.json'))
