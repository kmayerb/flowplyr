# start.R 
# Nov 14, 2023

# Replicate dict.get(x,x) functionality in R 
get_or_default <- function(dict, key, default = key) {
  if (key %in% names(dict)) {
    return(dict[key])
  } else {
    warning(stringr::str_glue("\n!!!\tThe marker: ", key, " was unchanged"))
    return(default)
  }
}

get_within_parent_gate <- function(parent_gate = 'Time/S/Lv/K1/K2/K3/K4/K5/K6/K7/K8/14-/S/L/19-/3+/3+excl 16br/56-16-/4+'){
  parts = stringr::str_split(parent_gate , pattern = stringr::fixed("/"))
  return(parts[[1]])
}

get_child_paths <- function(fcm08_filepath        = params$fcm08_path,
                            parent_gate           = params$parent_gate){
  fcm08a = readr::read_tsv(fcm08_filepath)
  all_paths = fcm08a %>% pull(SUBSET) %>% unique()
  child_paths = all_paths[all_paths %>% stringr::str_detect(., pattern = stringr::fixed(parent_gate))]
  return(child_paths)
  
}

extract_last_element <- function(x){
  split_string <- stringr::str_split(string = x, pattern = stringr::fixed("/"))[[1]]
  # Extract the last element
  final_element <- tail(split_string, 1)
  return(final_element)
}

is_single_marker<-function(x){
  n = stringr::str_count(string = x, pattern = "\\+|\\-|AND|OR|NOT|and|or|not")
  if (n > 1){
    return(FALSE)
  }else
    return(TRUE)
}

generate_metadata <- function(fcm08_filepath       = params$fcm08_path,
                              parent_gate          = params$parent_gate,
                              stim_exclusion_terms = params$stim_exclusion_terms){
  
  fcm08a = readr::read_tsv(fcm08_filepath)
  metadata = fcm08a %>% 
    filter(SUBSET == parent_gate) %>%
    select(experiment_name = ASSAYID, 
           sample_order = SAMP_ORD, 
           ptid   = PTID, 
           visit  = VISITNO, 
           guspec = GUSPEC, 
           stim   = ANTIGEN) %>% 
    group_by(experiment_name, sample_order, ptid, visit, stim, guspec) %>% 
    filter(!stim %in% stim_exclusion_terms)
  return(metadata)
}
