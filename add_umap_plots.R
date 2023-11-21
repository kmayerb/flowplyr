# add_umap_plots.R
# add umap plots to output_directory 
cat('888.....888.888b.....d888........d8888.8888888b..\n')
cat('888.....888.8888b...d8888.......d88888.888...Y88b\n')
cat('888.....888.88888b.d88888......d88P888.888....888\n')
cat('888.....888.888Y88888P888.....d88P.888.888...d88P\n')
cat('888.....888.888.Y888P.888....d88P..888.8888888P".\n')
cat('888.....888.888..Y8P..888...d88P...888.888.......\n')
cat('Y88b...d88P.888..."...888..d8888888888.888.......\n')
cat('."Y88888P"..888.......888.d88P.....888.888.......\n')
parser <- argparser::arg_parser(
  description='compile_flow_events - extract events from mulitple batches
  from .hdf5 formatted data')
parser <- argparser::add_argument(parser, arg="--params", 
                                  type="character",
                                  default = NA, help = "path to parameters as a .json file")
args   <- argparser::parse_args(parser)
# Get the JSON params file path from the command line
params_file <- args$params
# Check if the file exists
if (!file.exists(params_file)) {
  cat("Error: The specified JSON params file does not exist.\n")
  quit(status = 1)
}

source("R/open_hdf5.R")
source("R/visualize.R")

# Use the standard pipeline.json file
params <- jsonlite::fromJSON('2023-11-20_EEA9FFEE_VTN137_Part_A_CD4_TESTRUN/pipeline.json')

# prepare_umap_data
# SEE R/visualize.R for function showing how to assemble the various parts of .h5
# file into one file
umap_data = prepare_visualize_umap_data(params = params)


cat("PLOT 1: UMAP by experimental name\n")
p = visualize_umap(data = umap_data, by = "experiment_name")
# NOTE ONE CAN FACET BY ANY VALID VARIABLE IN <UMAP_DATA>
p = p + facet_wrap(~experiment_name) + 
  guides(col = guide_legend(override.aes = list(size = 2)))
pdf_name = file.path(params$output_path, "experimental_name.umap.pdf")
cat(paste("Writing to pdf: " ,pdf_name, "\n"))
pdf(pdf_name, 
    width = 8, 
    height = 8)
p
dev.off()

cat("PLOT 2: UMAP by antigen stimulation\n")
p = visualize_umap(data = umap_data, by = "stim")
# NOTE ONE CAN FACET BY ANY VALID VARIABLE IN <UMAP_DATA>
p = p + facet_wrap(~stim) +
  guides(col = guide_legend(override.aes = list(size = 2)))
pdf_name = file.path(params$output_path, "stimulation.umap.pdf")
cat(paste("Writing to pdf: " ,pdf_name, "\n"))
pdf(pdf_name, 
    width = 8, 
    height = 8)
p
dev.off()

cat("PLOT 3: UMAP by phenographic cluster\n")
p = visualize_umap(data = umap_data, by = "cluster")
# NOTE ONE CAN FACET BY ANY VALID VARIABLE IN <UMAP_DATA>
p = p + facet_wrap(~stim) +
  guides(col = guide_legend(override.aes = list(size = 2)))
pdf_name = file.path(params$output_path, "cluster.umap.pdf")
cat(paste("Writing to pdf: " ,pdf_name, "\n"))
pdf(pdf_name, 
    width = 8, 
    height = 8)
p
dev.off()


cat("PLOT 4: UMAP by channel intensity\n")

extract_last_element <- function(x){
  split_string <- stringr::str_split(string = x, pattern = stringr::fixed("|"))[[1]]
  # Extract the last element
  final_element <- tail(split_string, 1)
  return(final_element)
}

for (channel in params$cluster_channels){
  cat("PLOT 4: UMAP for: ", channel,"\n")
  p = visualize_umap(data = umap_data, 
                     by = "channel", 
                     channel = channel)
  channel_short = extract_last_element(channel)
  pdf_name = file.path(params$output_path,  paste0(channel_short, ".umap.pdf"))
  cat(paste("Writing to pdf: " ,pdf_name, "\n"))
  pdf(pdf_name, 
      width = 4, 
      height =4)
  print(p)
  dev.off()
}





