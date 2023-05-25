# Help identify the parent path
# May 24, 2023
# Use this commandline tool to inspect potential marker paths

source('R/validate_input.R')
source('R/extract.R')

# Examples:
# 1. Example with direct commandline inputs
# Rscript extract_marker_paths.R \
#--xml_path /fh/fast/gilbert_p/fg_data/ics_test/xml_files/ICS_PROFICIENTCY_DATA_FJ.xml \
#--fcs_folder_path /fh/fast/gilbert_p/fg_data/ics_test/fcs_folders/3639-L-ICS_PROFICIENCY_3/

# 2. Example with .json inputs
# Rscript extract_marker_paths.R --params tests/test_params.json

# Get commandline arguments
parser <- argparser::arg_parser(description='inspect marker paths')
parser <- argparser::arg_parser(
  description='extract_flow_events - extract events from FlowJo to .hdf5 format')

parser <- argparser::add_argument(parser,
                                  arg="--params",
                                  type="character",
                                  help = "path to params as .json file")
parser <- argparser::add_argument(parser,
                                  arg="--xml_path",
                                  type="character",
                                  help = "path to xml file, specifying FlowJo workspace")
parser <- argparser::add_argument(parser,
                                  arg="--fcs_folder_path",
                                  type="character",
                                  help = "path to folder with fcs files referenced in xml")

args  <- argparser::parse_args(parser)

# validate commandline arguments
# what are the required arguments
req_args <- c("xml_path",
              "fcs_folder_path")
# check that the user successfully passed these arguments
# in a params file or did so in .json params file
params = validate_args_or_params_file(req_args = req_args)

print(params)

# HARD CODE xml_keywords for now, until we know if these will ever vary.
params$xml_keywords = c("$FIL", "Stim", "Sample Order", "EXPERIMENT NAME", "Replicate")

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
# <markers>
all_marker_paths = extract_marker_paths(my_gs)
all_markers      = extract_markers(my_gs)
key_markers      = c(extract_cd4_markers(my_gs), extract_cd8_markers(my_gs))

message(sprintf("%s All Markers", format(Sys.time(), "%X")))
print(dput(all_markers))
message(sprintf("%s key CD4 Markers", format(Sys.time(), "%X")))
print(dput(extract_cd4_markers(my_gs)))
message(sprintf("%s key CD8 Markers", format(Sys.time(), "%X")))
print(dput(extract_cd8_markers(my_gs)))

