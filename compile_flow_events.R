# Compile flow events
# This script is used after extracting flow events
# parameters arguments are provided in a .json file

# Rscript compile_flow_events.R --params tests/test_summarize.json --annoy_cluster TRUE --umap TRUE
#params = jsonlite::fromJSON('/fh/fast/gilbert_p/fg_data/flowplyr/tests/test_summarize.json')
# read a .json file with run paramters
parser <- argparser::arg_parser(
  description='compile_flow_events - extract events from mulitple batches
  from .hdf5 formatted data')

parser <- argparser::add_argument(parser, arg="--params", 
                                  type="character",
                                  default = NA, help = "path to parameters as a .json file")
parser <- argparser::add_argument(parser, arg="--cluster", 
                                  type="boolean",
                                  default = FALSE, help = "Run Leiden Phenograph")

parser <- argparser::add_argument(parser, arg="--annoy_cluster", 
                                  type="boolean",
                                  default = FALSE, help = "Run Leiden Phenograph with Annoy")

parser <- argparser::add_argument(parser, arg="--umap", 
                                  type="boolean",
                                  default = FALSE, help = "Run Umap")
parser <- argparser::add_argument(parser, arg="--verbose", 
                                  type="boolean",
                                  default = FALSE, help = "extra verbosity")
args   <- argparser::parse_args(parser)

# Get the JSON params file path from the command line
params_file <- args$params
# Check if the file exists
if (!file.exists(params_file)) {
  cat("Error: The specified JSON params file does not exist.\n")
  quit(status = 1)
}

# LOAD 
source('R/open_hdf5.R')
suppressMessages(require(magrittr))
suppressMessages(require(dplyr))

# Start a Log file
log_file <- paste0(".complile_flow_events.log")
sink(log_file, append = FALSE, split = TRUE)
current_time <- Sys.time()
cat("Running compile_flow_event.R\n")
cat(paste0("See log file: ",log_file,"\n"))
cat(paste0("Run started :", format(current_time, format = "%Y-%m-%d %H:%M:%S")))

# Read and parse the JSON params file
params <- jsonlite::fromJSON(params_file)
# A function for messaging the contents of the params.json
verb <- function(params){
  cat('\n"h5_output_path_all_events:\n\t')
  cat(params$h5_output_path_all_events)
  cat('\n"h5_output_path_subset_events:\n\t')
  cat(params$h5_output_path_subset_events)
  cat('\nbatch_list:\n\t')
  cat(params$batch_list)
  cat('\nmetadata_file:\n\t')
  cat(params$metadata_file)
  cat('\nsubset_frac:\n\t')
  cat(params$subset_frac)
  cat('\nsubset_seed:\n\t')
  cat(params$subset_seed)
  cat('\nptid_visit_exclude:\n\t')
  cat(params$ptid_visit_exclude)
  cat('\ncluster_sample_exclusion:\n\t')
  cat(params$cluster_sample_exclusion)
  cat('\numap_markers:\n\t')
  cat(params$umap_markers )
  cat('\ncluster_markers:\n\t')
  cat(params$cluster_markers)
  cat('\ngroup_vars:\n\t')
  cat(params$group_vars)
  cat('\njoin_vars:\n\t')
  cat(params$join_vars)
  cat('\nboolean_strings:\n\t')
  cat(params$boolean_strings)
  cat('\nany_of_strings:\n\t')
  cat(params$any_of_strings)
  cat('\nnegative_stim_string:\n\t')
  cat(params$negative_stim_string)
}
if (args$verbose ){verb(params)}


is_list_not_null <- function(input_list) {
  if (!is.null(input_list) && length(input_list) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

cat("\n")
cat('.........................................................\n')
cat('..####....####...##...##..#####...######..##......######.\n')
cat('.##..##..##..##..###.###..##..##....##....##......##.....\n')
cat('.##......##..##..##.#.##..#####.....##....##......####...\n')
cat('.##..##..##..##..##...##..##........##....##......##.....\n')
cat('..####....####...##...##..##......######..######..######.\n')
cat('.........................................................\n')
# Step 1: Load metadata across batches from <metadata_file> param
cat("\nStep 1: Load metadata across batches from <metadata_file> param\n")
cat(paste0("\t", params$metadata_file,"\n"))

metadata = read.csv(params$metadata_file)
stopifnot('batch' %in% names(metadata))
stopifnot('sample_order' %in% names(metadata))
stopifnot('ptid' %in% names(metadata))
stopifnot('visit' %in% names(metadata))
metadata = metadata %>% mutate(ptid_visit = paste0(ptid, "_", visit))

# Step 2: Compile all data across batches based on <batch_list> param
cat("Step 2: Compile events across batches\n")
cat(paste0("\tBatch : ", params$batch_list,"\n"))
data = unpack_hdf5s(paths = params$batch_list)

cat("Step 3: Fix the names in the data$pos matrix using <pos_map>\n")
pos_map = params$pos_map
if (!is.null(params$pos_map)){
  colnames(data$pos) <- as.vector(pos_map[colnames(data$pos)])
}
if (!is.null(params$fi_map)){
  colnames(data$raw) <- as.vector(fi_map[colnames(data$raw)])
  colnames(data$fi)  <- as.vector(fi_map[colnames(data$fi)])
}
# Step 3: Add metadata to the data$fcs dataframe, by batch and sample order
cat("Step 4: Add metadata to the data$fcs dataframe, by batch and sample order\n")
data$fcs$sample_order <- as.character(data$fcs$sample_order)
metadata$sample_order <- as.character(metadata$sample_order)
data$fcs_ptid = data$fcs %>% 
  left_join(metadata, 
            by = c("batch","sample_order"))
# defensive, make sure that left join was 1 to 1
if (!( dim(data$fcs_ptid)[1] == dim(data$fcs)[1] )){
  stop("Your metadata file did not perfectly match the batch, sample_order info")
}
# defensive, make sure no NA in ptid column, issue warning if ther is
if(any(is.na(data$fcs_ptid$ptid))){
  warning("You have some NAs in the ptid column, suggesting metadata file was incomplete")
}

# !!!! # 
# We need to store unique sample count
unique_samples = data$fcs_ptid %>% 
  group_by(dummy, 
           batch, 
           sample_order,
           sample_name, 
           ptid,
           visit,
           parent_ct,
           total_ct) %>% 
  slice(1)

# Step 5: subset data , potentially exclude ptid_visits
cat("Step 5: Subset data, potentially excluding ptid_visits\n")
cat(paste0("\tSubset fraction: ",params$subset_frac, "\n"))
cat(paste0("\tSubset seed: ",params$subset_seed, "\n"))    
data_subset = get_data_subset(data, 
                              ptid_visit_exclude = params$ptid_visit_exclude, 
                              frac = params$subset_frac, 
                              seed = params$subset_seed)

unique_samples_subset <- data_subset$fcs_ptid %>% 
  group_by(dummy, 
           batch, 
           sample_order,
           sample_name, 
           ptid,
           visit,
           parent_ct,
           total_ct) %>% 
  slice(1)

batches_to_h5 <- function(data, h5_name, all_samples, subset_samples){
  if (!is.null(h5_name)){
    message(sprintf("Writing output to .h5 files"))
    message(h5_name)
    if (file.exists(h5_name)){file.remove(h5_name)}
    rhdf5::h5createFile(h5_name)
    rhdf5::h5createGroup(h5_name, "data")
    rhdf5::h5write(data$pos,             h5_name, "data/pos")
    rhdf5::h5write(data$fi,              h5_name, "data/fi")
    rhdf5::h5write(data$raw,             h5_name, "data/raw")
    rhdf5::h5write(data$fcs,             h5_name, "data/fcs")
    rhdf5::h5write(data$fcs_ptid,        h5_name, "data/fcs_ptid")
    rhdf5::h5write(colnames(data$pos),   h5_name, "data/cols_pos")
    rhdf5::h5write(colnames(data$fi),    h5_name, "data/cols_fi")
    rhdf5::h5write(all_samples,          h5_name, "data/all_samples")
    rhdf5::h5write(subset_samples,       h5_name, "data/subset_samples")
    rhdf5::h5closeAll()
  }
}

cat("Step 6a: Write all events to h5\n")
cat(paste0("\tEvents: ", dim(data$fcs_ptid)[1],"\n"))
ptid_visit_n = data$fcs_ptid %>% 
  pull(ptid_visit) %>% 
  unique() %>% 
  length()
cat(paste0("\tPtid_Visit (n) : ", ptid_visit_n,"\n"))
cat(paste0("\tSamples : ", dim(unique_samples)[1],"\n"))
cat(paste0("\tFile: ", params$h5_output_path_all_events,"\n"))
# write to H5
suppressMessages(batches_to_h5(data = data, 
                               params$h5_output_path_all_events,
                               all_samples = unique_samples, 
                               subset_sample = unique_samples))

cat("Step 6b: Write subset of events to h5\n")
cat(paste0("\tEvents: ", dim(data_subset$fcs_ptid)[1],"\n"))
ptid_visit_n_subset = data_subset$fcs_ptid %>% 
  pull(ptid_visit) %>% 
  unique() %>% 
  length()
cat(paste0("\tPtid_Visit (n) : ", ptid_visit_n_subset,"\n"))
cat(paste0("\tDropped (ptid_visit) : ", params$ptid_visit_exclude,"\n"))
cat(paste0("\tDropped (samples) : ", dim(unique_samples)[1]-dim(unique_samples_subset)[1],"\n"))
cat(paste0("\tFile: ", params$h5_output_path_subset_events,"\n"))
# write to H5
suppressMessages(batches_to_h5(data = data_subset, 
                               params$h5_output_path_subset_events,
                               all_samples = unique_samples, 
                               subset_sample = unique_samples_subset))


if (args$umap){
  cat("Step 7a: Finding UMAP coordinates for the subset data\n")
  cat(paste0("\tFinding UMAP on ", dim(data_subset$fi)[1], " cells"))
  set.seed(1)
  start_time <- Sys.time()
  data_umap = uwot::umap(data_subset$fi[,params$umap_markers], 
                         n_neighbors   = params$umap_n_neighors,
                         learning_rate = params$umap_learning_rate,
                         init = "random") 
  data_umap %>% write.csv("umap_test.csv")
  
  data_umap2 = data_umap[] %>% 
      as.data.frame()%>%
      dplyr::transmute(UMAP1 = V1,
                       UMAP2 = V2) %>% 
      as.matrix()

  
  end_time <- Sys.time()
  elap = end_time - start_time
  cat("\n\tUMAP: ")
  cat(elap)
  cat(" seconds \n")

  umap_to_h5 <- function(x,
                         h5_path, 
                         fcs_ptid = NULL){
    cat("\n\tWriting umap output to .h5 files:\n")
    cat(h5_path)
    if ( file.exists(h5_path) ){
      cat("\n\tcleaning up old file\n")
      file.remove(h5_path)
    }
    rhdf5::h5createFile(h5_path)
    rhdf5::h5createGroup(h5_path, "umap")
    rhdf5::h5write(x, h5_path, "umap/coordinates")
    if (!is.null(fcs_ptid) ){
      rhdf5::h5write(fcs_ptid,  h5_path, "umap/fcs_ptid")
    }
    rhdf5::h5closeAll()
  }
  
  umap_to_h5(x = data_umap2, 
             h5_path  = params$h5_output_path_subset_umap,
             fcs_ptid = data_subset$fcs_ptid)
}

## LEIDEN PHENOGRAPH CLUSTERING USING APPROXIMATE NN BASED CLUSTERING USING ANNOY (Spotify)
if (args$annoy_cluster){
  suppressPackageStartupMessages(loadNamespace("Rcpp"))
  suppressPackageStartupMessages(require(Rphenograph))
  jaccard_coeff <- function(idx) {
    .Call('Rphenograph_jaccard_coeff', PACKAGE = 'Rphenograph', idx)
  }
  cat("Step 7b: Cluster data subset with Annoy\n")
  k = 30
  data_mat = data_subset$fi[,params$cluster_markers]
  cat(paste0("\tFinding nearest neighbors on ", dim(data_mat)[1], " cells"))
  start_time <- Sys.time()
  snn <- uwot:::annoy_nn(data_mat, 
                   k = k + 1, 
                   n_threads = params$annoy_threads)$idx[, -1]
  end_time <- Sys.time()
  elap = end_time - start_time
  cat(paste0("\n\tuwot:::annoy_nn with ", params$annoy_threads, " threads "))
  cat(elap)
  cat("seconds\n")
  # Compute jaccard coefficient between nearest-neighbor set
  # Undirected graph from the weighted links
  #snn = readRDS('/fh/fast/gilbert_p/fg_data/flowplyr-cluster/flowplyrcluster/snn_copy.Rds')
  start_time <- Sys.time()
  links <- jaccard_coeff(snn)
  end_time <- Sys.time()
  elap = end_time - start_time
  cat("\tJaccard: ")
  cat(elap)
  cat(" seconds \n")
  start_time <- Sys.time()
  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")
  g <- igraph::graph.data.frame(relations, directed = FALSE)
  r <- quantile(igraph::strength(g))[2] / (igraph::gorder(g) - 1)
  set.seed(1)
  ldc <- igraph::cluster_leiden(g, resolution_parameter=r)
  end_time <- Sys.time()
  cat("\tLeiden: ")
  elap = end_time - start_time
  cat(elap)
  cat(" seconds \n")
  stopifnot(dim(data_subset$fcs_ptid)[1] == length(ldc$membership))
  clustering_to_h5 <- function(x,
                               h5_path, 
                               fcs_ptid = NULL){
    
    #x = ldc 
    #h5_path = '/fh/fast/gilbert_p/fg_data/POI_COR/data/flowplyr_ics_gs/all_batches_clustered.h5'
    #fcs_ptid = data_subset$fcs_ptid
    message("Writing clustering output to .h5 files")
    message(h5_path)
    if ( file.exists(h5_path) ){
      file.remove(h5_path)
    }
    rhdf5::h5createFile(h5_path)
    rhdf5::h5createGroup(h5_path, "clustering")
    rhdf5::h5write(x$membership,   h5_path, "clustering/membership")
    rhdf5::h5write(x$nb_clusters,  h5_path, "clustering/nb_clusters")
    rhdf5::h5write(x$quality,      h5_path, "clustering/quality")
    rhdf5::h5write(x$algorithm,    h5_path, "clustering/algorithm")
    rhdf5::h5write(x$vcount,       h5_path, "clustering/vcount")
    rhdf5::h5write(x$names,        h5_path, "clustering/names")
    if (!is.null(fcs_ptid) ){
      rhdf5::h5write(fcs_ptid,        h5_path, "clustering/fcs_ptid")
    }
    rhdf5::h5closeAll()
  }
  clustering_to_h5(x = ldc, 
                   h5_path = params$h5_output_path_subset_events_annoy_clustering,
                   fcs_ptid = data_subset$fcs_ptid)

  }

## LEIDEN PHENOGRAPH CLUSTERING USING EXACT NN (RANN) BASED CLUSTERING
if (args$cluster){
  loadNamespace("Rcpp")
  require(Rphenograph)
  jaccard_coeff <- function(idx) {
    .Call('Rphenograph_jaccard_coeff', PACKAGE = 'Rphenograph', idx)
  }
  
  cat("Step 7: Cluster data subset5\n")
  k = 30
  data_mat = data_subset$fi[,params$cluster_markers]
  cat(paste0("Finding nearest neighbors on ", dim(data_mat)[1], " cells\n"))
  cat(paste0("Expect 20-30 minutes for 100,000 cells\n"))
  
  start_time <- Sys.time()
  snn <- RANN::nn2(data_mat, searchtype = "standard", 
                   k = k + 1, 
                   eps = params$cluster_eps)$nn.idx[, -1]
  end_time <- Sys.time()
  elap = end_time - start_time
  cat("RANN\n")
  cat(elap)
  cat("\n")
  # Compute jaccard coefficient between nearest-neighbor set
  # Undirected graph from the weighted links
  #snn = readRDS('/fh/fast/gilbert_p/fg_data/flowplyr-cluster/flowplyrcluster/snn_copy.Rds')
  start_time <- Sys.time()
  links <- jaccard_coeff(snn)
  end_time <- Sys.time()
  elap = end_time - start_time
  cat("Jaccard\n")
  cat(elap)
  cat("\n")
  start_time <- Sys.time()
  links <- links[links[, 1] > 0, ]
  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")
  g <- igraph::graph.data.frame(relations, directed = FALSE)
  r <- quantile(igraph::strength(g))[2] / (igraph::gorder(g) - 1)
  set.seed(1)
  ldc <- igraph::cluster_leiden(g, resolution_parameter=r)
  end_time <- Sys.time()
  cat("Leiden\n")
  elap = end_time - start_time
  cat(elap)
  cat("\n")
  stopifnot(dim(data_subset$fcs_ptid)[1] == length(ldc$membership))
  
  clustering_to_h5 <- function(x,
                               h5_path, 
                               fcs_ptid = NULL){
    

    cat("Step 7c: Writing clustering output to .h5 files:\n")
    cat(paste0("File: ", h5_path))
    cat("\n")
    if ( file.exists(h5_path) ){
      file.remove(h5_path)
    }
    rhdf5::h5createFile(h5_path)
    rhdf5::h5createGroup(h5_path, "clustering")
    rhdf5::h5write(x$membership,   h5_path, "clustering/membership")
    rhdf5::h5write(x$nb_clusters,  h5_path, "clustering/nb_clusters")
    rhdf5::h5write(x$quality,      h5_path, "clustering/quality")
    rhdf5::h5write(x$algorithm,    h5_path, "clustering/algorithm")
    rhdf5::h5write(x$vcount,       h5_path, "clustering/vcount")
    rhdf5::h5write(x$names,        h5_path, "clustering/names")
    if (!is.null(fcs_ptid) ){
      rhdf5::h5write(fcs_ptid,     h5_path, "clustering/fcs_ptid")
    }
    rhdf5::h5closeAll()
  }
  clustering_to_h5(x = ldc, 
                   h5_path = params$h5_output_path_subset_events_clustering,
                   fcs_ptid = data_subset$fcs_ptid)

}
sink()

