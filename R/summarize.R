# summarize.R

# Main functions:
#   get_standard_pctpos_pctposadj - gets values for all markers
#   get_boolean_pctpos_pctposadj  - gets values based on boolean combinations
#   get_x_of_n_pctpos_pctposadj   - gets values for x_of_n

# Helper functions:
#   validate_grammar
#   parse_boolean_string
#   return_x_of_n 
#   all_equal_one
#   all_equal_zero
#   get_all_index_key_combinations


#' get_standard_pctpos_pctposadj
#'
#' @param data 
#' @param group_vars 
#' @param join_vars 
#' @param negative_stim_string 
#'
#' @return
#' @export
#'
#' @examples
#' require(ggplot2)
#' data_pos_freq_summary <- get_standard_pctpos_pctposadj(data)
#' data_pos_freq_summary %>% 
#'   ggplot(aes(x = ptid, y = key, col = log10(pctposadj))) + 
#'   geom_point() + 
#'   facet_wrap(~stim.pos, ncol= 1) + 
#'   theme(axis.text.x = element_text(angle = 90, size = 6, hjust = 1)) +
#'   theme(axis.text.y = element_text(angle = 0, size = 6, hjust = 1)) +
#'   scale_color_viridis_c(limits = c(-5,-2), option = "D", oob = scales::squish)
get_standard_pctpos_pctposadj <- function(data,
                                          group_vars = c('batch', 
                                                         'ptid',
                                                         'visit', 
                                                         'stim',
                                                         'parent_ct'),
                                          join_vars = c('batch', 'ptid', 'visit','key'),
                                          negative_stim_string = 'negctrl'){
  
  # All variables that we have yes or now values for
  variables =  colnames(data$pos)
  # Fuse the fcs_ptid dataframe  with in-gate positive call matrix, 
  # so that we can groupby index columns
  data_pos  =  cbind(data$fcs_ptid, data$pos)
  # we will summarise over all of variables in the data$pos matrix
  sum_vars <- colnames(data$pos)
  # For each sample, and for each variable we take the sum of positive 
  # divided by the parent_gate count to get the 'freq', then, after we've 
  # pivoted to long form we can regenerated the count by taking the product
  # of the frequency and parent_ct
  data_pos_freq <- data_pos %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>% 
    dplyr::summarize(dplyr::across(dplyr::all_of(sum_vars), sum), .groups = "keep") %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of(sum_vars), ~ . / parent_ct)) %>% 
    tidyr::pivot_longer(cols = !all_of(group_vars), 
                        names_to = "key", values_to = "freq") %>% 
    mutate(count = freq*parent_ct)
  # Next,  we need to extract the negative control data
  data_pos_freq_negctrl = 
    data_pos_freq %>% 
    dplyr::filter(stim == negative_stim_string)
  # Left join negative control frequency and count to positive frequency and count
  # Then compute a (i) pctpos and (ii) pctposadj
  data_pos_freq_summary = data_pos_freq %>% 
    dplyr::left_join(data_pos_freq_negctrl, 
                     by = join_vars, 
                     suffix = c(".pos",".neg")) %>%
    dplyr::mutate(count.neg = tidyr::replace_na(count.neg, 0))%>%
    dplyr::mutate(freq.neg = tidyr::replace_na(freq.neg, 0)) %>%
    dplyr::mutate(pctpos = freq.pos, 
                  pctposadj = freq.pos-freq.neg)
  
  return(data_pos_freq_summary)
}



get_cluster_pctpos_pctposadj <- function(data_pos,
                                         group_vars = c('batch', 
                                                         'ptid',
                                                         'visit', 
                                                         'stim',
                                                         'parent_ct',
                                                         'cluster'),
                                          min_cell_per_cluster = 100,
                                          populate_zeros = TRUE,
                                          cluster_string = "Leiden ",
                                          join_vars = c('batch', 'ptid', 'visit','key'),
                                          negative_stim_string = 'negctrl'){
  

  data_pos_freq <- data_pos %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>% 
    dplyr::summarize(count = n(), .groups = "keep") %>% 
    dplyr::mutate(freq = count/parent_ct) %>% 
    dplyr::rename(key = cluster)

  # Next,  we need to extract the negative control data
  data_pos_freq_negctrl = 
    data_pos_freq %>% 
    dplyr::filter(stim == negative_stim_string)
  # Left join negative control frequency and count to positive frequency and count
  # Then compute a (i) pctpos and (ii) pctposadj
  data_pos_freq_summary = data_pos_freq %>% 
    dplyr::left_join(data_pos_freq_negctrl, 
                     by = join_vars, 
                     suffix = c(".pos",".neg")) %>%
    dplyr::mutate(count.neg = tidyr::replace_na(count.neg, 0))%>%
    dplyr::mutate(freq.neg = tidyr::replace_na(freq.neg, 0)) %>%
    dplyr::mutate(freq.neg = ifelse(is.na(stim.neg,
                                          params$negative_stim_string,
                                          stim.neg))) %>%
    dplyr::mutate(pctpos = freq.pos, 
                  pctposadj = freq.pos-freq.neg)
  

  return(data_pos_freq_summary)
}















parse_boolean_string<-function(string){
  x = stringr::str_split(string, "/",)
  x = x[[1]]
  posneg = purrr::map_chr(x, ~stringr::str_extract(.x, ".$"))
  stopifnot(all(posneg %in% c("+","-")))
  x      = purrr::map_chr(x, ~stringr::str_remove(.x, ".$"))
  #print(x)
  #print(posneg)
  stopifnot(all(posneg %in% c("+","-")))
  pos_markers = x[posneg == "+"]
  neg_markers = x[posneg == "-"]
  #print(length(pos_markers))
  #print(length(neg_markers))
  if (length(pos_markers) < 1){
    pos_markers = NA
  }
  if (length(neg_markers)< 1){
    neg_markers = NA
  }
  
  return(list(pos_markers = pos_markers, 
              neg_markers = neg_markers))
  
}
# Tests
identical(parse_boolean_string("IL2+/IFNg-"),list('pos_markers'=c("IL2"), "neg_markers"=c("IFNg")))
identical(parse_boolean_string("IL2+/IFNg+"),list('pos_markers'=c("IL2","IFNg"), "neg_markers"=NA))
identical(parse_boolean_string("IL2-/IFNg-"),list('pos_markers'=NA, "neg_markers"=c("IL2","IFNg")))
identical(parse_boolean_string("IL2+/IFNg+/CD45RA-"),list('pos_markers'=c("IL2","IFNg"), "neg_markers"=c("CD45RA")))

all_equal_one <- function(row) {
  all(row == 1)
}
all_equal_zero <- function(row) {
  all(row == 0)
}

get_boolean_value <- function(data, 
                              pos_markers=NA, 
                              neg_markers=NA){
  # Can handle only positives
  if (all(is.na(pos_markers))){
    neg = data$pos[,neg_markers, drop=F]
    all_neg = apply(neg, 1, all_equal_zero)
    all_pos_and_neg = all_neg
    # Can handle only negatives
  }else if (all(is.na(neg_markers))){
    pos = data$pos[,pos_markers, drop=F]
    all_pos = apply(pos, 1, all_equal_one)
    all_pos_and_neg = all_pos
    # Will return NAs if neither is provided
  }else if ( all(is.na(pos_markers)) & all(is.na(neg_markers))){
    all_pos_and_neg = rep(NA, diim(data)[1])
    # Can handle mixture of positivs and negatives 
  }else{
    pos = data$pos[,pos_markers, drop=F]
    neg = data$pos[,neg_markers, drop=F]
    all_pos = apply(pos, 1, all_equal_one)
    all_neg = apply(neg, 1, all_equal_zero)
    all_pos_and_neg = all_pos & all_neg
  }
  return(as.numeric(all_pos_and_neg))
}

get_boolean_pctpos_pctposadj <- function(data,
                                         boolean_strings = c("IFNg+",
                                                             "IL2+",
                                                             "IFNg+/IL2+",
                                                             "IFNg+/IL2-",
                                                             "IFNg-/IL2+",
                                                             "IFNg+/EM+"),
                                         group_vars = c('batch', 
                                                        'ptid',
                                                        'visit', 
                                                        'stim',
                                                        'parent_ct'),
                                         join_vars = c('batch', 'ptid', 'visit','key'),
                                         negative_stim_string = 'negctrl'){
  #browser()
  result = list()
  for (i in boolean_strings){
    ps = parse_boolean_string(i)
    #print(ps)
    result[[i]] = get_boolean_value(data=data,
                                    pos_markers = ps$pos_markers,
                                    neg_markers = ps$neg_markers)
  }
  
  # so that we can group by index columns
  data_pos  =  cbind(data$fcs_ptid, tibble::as_tibble(result))
  # we will summaries over all of variables in the data$pos matrix
  sum_vars <- boolean_strings
  # For each sample, and for each variable we take the sum of positive 
  # divided by the parent_gate count to get the 'freq', then, after we've 
  # pivoted to long form we can regenerated the count by taking the product
  # of the frequency and parent_ct
  data_pos_freq <- data_pos %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>% 
    dplyr::summarize(dplyr::across(dplyr::all_of(sum_vars), sum), .groups = "keep") %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of(sum_vars), ~ . / parent_ct)) %>% 
    tidyr::pivot_longer(cols = !all_of(group_vars), 
                        names_to = "key", values_to = "freq") %>% 
    mutate(count = freq*parent_ct)
  # Next,  we need to extract the negative control data
  data_pos_freq_negctrl = 
    data_pos_freq %>% 
    dplyr::filter(stim == negative_stim_string)
  # Left join negative control frequency and count to positive frequency and count
  # Then compute a (i) pctpos and (ii) pctposadj
  data_pos_freq_summary = data_pos_freq %>% 
    dplyr::left_join(data_pos_freq_negctrl, 
                     by = join_vars, 
                     suffix = c(".pos",".neg")) %>%
    dplyr::mutate(count.neg = tidyr::replace_na(count.neg, 0))%>%
    dplyr::mutate(freq.neg = tidyr::replace_na(freq.neg, 0)) %>%
    dplyr::mutate(pctpos = freq.pos, 
                  pctposadj = freq.pos-freq.neg)
  
  return(data_pos_freq_summary)
}


validate_any_of_grammar <- function(x,check_names){
  # must have commas
  if(!grepl(x, pattern='/')) 
    stop(paste0("Error: ",x," does not contain '/'"))
  if(!grepl(x, pattern="^[0-9]/")) 
    stop(paste0("Error: ",x," does not start with an integer"))
  x_of_n_list= strsplit(x , '/')
  x = as.numeric(x_of_n_list[[1]][1])
  n = x_of_n_list[[1]][2:length(x_of_n_list[[1]])]
  #print(n)
  #print(check_names)
  if(!all(n %in% c(check_names))){
    missing = n[!n %in% c(check_names)]
    stop(paste0("\nError:",missing," is not in input data columns"))
  }
  return(x)
}

return_x_of_n <- function(data_pos,
                          x_of_n, 
                          split = "/"){
  validate_any_of_grammar(x = x_of_n, check_names=colnames(data_pos))
  x_of_n_list= strsplit(x_of_n , split)
  #print(x_of_n)
  x = as.numeric(x_of_n_list[[1]][1])
  n = x_of_n_list[[1]][2:length(x_of_n_list[[1]])]
  if (length(n) == 1){
    return(as.numeric(data$pos[,n] >= x))
  }else{
    sum_over_n = apply(data$pos[, n, drop = F], 1, sum)
    return(as.numeric(sum_over_n >=x))
  }
}


#' get_x_of_n_pctpos_pctposadj
#'
#' @param data must be a data object with data$pos, data$fcs_ptid
#' @param combinations list of in comma delim format e.g., c("1,IL2+,IFNg+")
#' indicating at least 1 of following list
#' @param group_vars variables that define a flow sample
#' @param join_vars variables that would join experimental and neg. ctrl stim
#' @param negative_stim_string a string matching negative control stimulation
#'
#' @return data_pos_freq_summary - long form dataframe with pct positive
#' with new columns for each combination
#' @export
#'
#' @examples
#' 
#' 
#' get_x_of_n_pctpos_pctposadj(
#'     combinations = c("1,IFNg+",
#'                      "1,IL2+",
#'                      "1,TNFa+",
#'                      "1,IFNg+,IL2+,TNFa+",
#'                      "1,IFNg+,IL2+",
#'                      "2,IFNg+,IL2+"),
#'     group_vars = c('batch', 
#'                    'ptid',
#'                    'visit', 
#'                    'stim',
#'                    'parent_ct'),
#'     join_vars = c('batch', 'ptid', 'visit','key'),
#'     negative_stim_string = 'negctrl')
get_x_of_n_pctpos_pctposadj <- function(data,
                                        combinations,
                                        group_vars = c('batch', 
                                                       'ptid',
                                                       'visit', 
                                                       'stim',
                                                       'parent_ct'),
                                        join_vars = c('batch', 'ptid', 'visit','key'),
                                        negative_stim_string = 'negctrl'){
  names(combinations) = combinations
  #print(combinations)
  k = purrr::map_df(combinations,
                    ~return_x_of_n(data_pos = data$pos, x_of_n = .x))
  variables = combinations
  # so that we can group by index columns
  data_pos  =  cbind(data$fcs_ptid, k)
  # we will summaries over all of variables in the data$pos matrix
  sum_vars <- combinations
  # For each sample, and for each variable we take the sum of positive 
  # divided by the parent_gate count to get the 'freq', then, after we've 
  # pivoted to long form we can regenerated the count by taking the product
  # of the frequency and parent_ct
  data_pos_freq <- data_pos %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>% 
    dplyr::summarize(dplyr::across(dplyr::all_of(sum_vars), sum), .groups = "keep") %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of(sum_vars), ~ . / parent_ct)) %>% 
    tidyr::pivot_longer(cols = !all_of(group_vars), 
                        names_to = "key", values_to = "freq") %>% 
    mutate(count = freq*parent_ct)
  # Next,  we need to extract the negative control data
  data_pos_freq_negctrl = 
    data_pos_freq %>% 
    dplyr::filter(stim == negative_stim_string)
  # Left join negative control frequency and count to positive frequency and count
  # Then compute a (i) pctpos and (ii) pctposadj
  data_pos_freq_summary = data_pos_freq %>% 
    dplyr::left_join(data_pos_freq_negctrl, 
                     by = join_vars, 
                     suffix = c(".pos",".neg")) %>%
    dplyr::mutate(count.neg = tidyr::replace_na(count.neg, 0)) %>%
    dplyr::mutate(freq.neg = tidyr::replace_na(freq.neg, 0)) %>%
    dplyr::mutate(pctpos = freq.pos, 
                  pctposadj = freq.pos-freq.neg)
  
  return(data_pos_freq_summary)
}

#' get_all_index_key_combinations
#' 
#' deal with cases where a cluster is not present
#'
#' @param d data frame 
#' @param index_cols index cols specifying combinations
#'
#' @return
#' @export
#'
#' @examples
get_all_index_key_combinations <- function(d, 
                                           index_cols, 
                                           unique_keys,
                                           key_name){
  unique_indices = d %>% 
    dplyr::select(dplyr::all_of(index_cols)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(index_cols))) %>%
    dplyr::tally() %>% 
    dplyr::select(dplyr::all_of(index_cols))
  
  foo<-function(i, ui){
    ui[[key_name]] = i
    return(ui)
  }
  comb_ = purrr::map_df(unique_keys, ~foo(.x, ui= unique_indices))
  return(comb_)
}
