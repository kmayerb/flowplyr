######
# May 23, 2023
# Dependencies contained within ml fhR/4.2.0-foss-2021b

#' get_gated_set
#'
#' get_gated_set from inputs
#'
#' @param xml_path - string to file path flowJo .xml output
#' @param fcs_folder_path - file path to folder with .fcs files
#' @param xml_keywords - vector of keywords used to extract work space elements
#'
#' @return gatingset object
#' @export
#'
#' @examples
get_gs <- function(xml_path,
                   fcs_folder_path,
                   xml_keywords){

  # xml_keywords <- c("$FIL",
  #                   "Stim",
  #                   "Sample Order",
  #                   "EXPERIMENT NAME",
  #                   "Replicate")

  # Make a wksp from the xlm_path
  wksp <- CytoML::open_flowjo_xml(xml_path,
                                  options = 1)

  # Make a set from the wksp + folder of fcs files
  gs   <- CytoML::flowjo_to_gatingset(wksp,
                                      name="Samples",
                                      path=fcs_folder_path,
                                      channel.ignore.case=TRUE,
                                      keywords=xml_keywords,
                                      keywords.source = "XML")
  return(gs)
}

#' extract_events
#'
#' @param g - a sample object within a FlowWorkspace gated set object
#' @param parent_gate - string specifying the parent gate of all markers
#' @param markers - vector of all marker to retain in output
#' @param functional_markers - vector of key marker names, used for subseting
#' to only events gated within one or more key markers
#'
#' @return list with three data.frame elements 'pos' , 'fi' , and 'fcs_index'
#' @export
#'
#' @examples
extract_events <- function(g,
                           parent_gate,
                           markers,
                           functional_markers){


  # #xtract meta-data, on the experiment name <exp_name> and <fcs_name>
  pd = flowWorkspace::pData(g)
  exp_name <- pd$`EXPERIMENT NAME`
  fcs_name <- paste(pd$`EXPERIMENT NAME`, pd$"Sample Order", pd$Replicate, pd$Stim, pd$name, sep = "|")

  # Get <total_events> identify the number of total_events - integer sepcifying
  # total number of events recorded
  total_events <- length(flowWorkspace::gh_pop_get_indices(g, '/'))

  # Get <pos> is a matrix of booleans whether an event falls in a specific gate
  #   it has row dimensions events, column dimensions markers
  pos <- vapply(X         = markers,
                FUN       = function(x) flowWorkspace::gh_pop_get_indices(g, file.path(parent_gate, x)),
                FUN.VALUE = logical(total_events))

  # Get <ix> is a logical index, True if the event coincided with any of the
    # functional_markers' gates
  ix = rowSums(pos[,functional_markers]) >= 1
  # Get <fi> matrix of florescent intensity
  fi <- flowCore::exprs(flowWorkspace::gh_pop_get_data(g))

  # Subset to only those rows where event falls in gate(s) of a key marker
  # <pos> subset
  pos = pos[ix,]
  # <fi> subset
  fi   = fi[ix,]
  # Names of fi channels
  fi_channels <- flowWorkspace::pData(flowCore::parameters(flowWorkspace::gh_pop_get_data(g))[ ,c("name", "desc")])
  # Combine the names e.g., (<G780-A>) with desc  e.g., IL-17a
  colnames(fi) <- paste(seq_len(ncol(fi)),fi_channels$name, fi_channels$desc, sep = "|")
  # Assign rownams based on fcs_name and seq counter
  rownames(fi) <- paste(fcs_name, seq_len(nrow(fi)), sep = "|")
  rownames(pos) <- paste(fcs_name, seq_len(nrow(pos)), sep = "|")

  fcs_cells = data.frame(list("cell"=seq_len(nrow(fi)),
                              "experiment_name"=pd$`EXPERIMENT NAME`,
                              "sample_order"=pd$`Sample Order`,
                              "stim"=pd$Stim,
                              "name"=pd$name,
                              "parent_ct"=total_events))
  result = list('pos' = pos,
                'fi' = fi,
                'fcs_index' = fcs_cells)

  return(result)
}



#' extract_marker_paths
#'
#' @param gs
#'
#' @return list of full paths (e.g., /Time/S/14-/Lv/L/3+)
#' @export
#'
#' @examples
extract_marker_paths <-function(gs){
  paths = flowWorkspace::gs_get_pop_paths(gs)
  return(paths)
}


#' extract_markers
#'
#' @param gs
#'
#' @return list of short markers names using path = "auto" (e.g, "154+" )
#' @export
#'
#' @examples
extract_markers <-function(gs){
  paths = flowWorkspace::gs_get_pop_paths(gs, path = "auto")
  return(paths)
}


# extract all cd4 markers
#' Title
#'
#' @param gs
#'
#' @return a list of markers with parent 4+ e.g., c("4+/154+", "4+/GzB+")
#' @export
#'
#' @examples
extract_cd4_markers <-function(gs){
  paths = flowWorkspace::gs_get_pop_paths(gs)
  cd4_paths = paths[grepl(paths , pattern = '/4\\+/', perl = F)]
  cd4_markers = stringr::str_match(cd4_paths, pattern = '/(4\\+/.*)')[,2]
  return(cd4_markers)
}


# extract all cd8 markers
#' Title
#'
#' @param gs
#'
#' @return  a list of markers with parent 4+ e.g., c("8+/154+", "8+/GzB+")
#' @export
#'
#' @examples
extract_cd8_markers <-function(gs){
  paths = flowWorkspace::gs_get_pop_paths(gs)
  cd8_paths = paths[grepl(paths , pattern = '/8\\+/', perl = F)]
  cd8_markers = stringr::str_match(cd8_paths, pattern = '/(8\\+/.*)')[,2]
  return(cd8_markers)
}
