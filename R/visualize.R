require(ggplot2)


prepare_visualize_umap_data <-function(params){
  d = unpack_data_hdf5(params$h5_output_path_subset_events)
  u = unpack_umap_hdf5(params$h5_output_path_subset_umap)
  cl = unpack_cluster_hdf5(params$h5_output_path_subset_events_clustering) 
  stopifnot(dim(cl$membership)[1] ==  dim(u$coordinates)[1])
  cl$membership = data.frame(cluster = cl$membership)
  colnames(u$coordinates) <- c("UMAP1","UMAP2")
  umap_data = cbind(u$fcs_ptid, u$coordinates,cl$membership,d$fi)
  return(umap_data)
}

visualize_umap <- function(data, 
                           by = "experiment_name", 
                           channel = NULL){
  if (by == "experiment_name"){
    cat("Visualizing by experiment_name (a.k.a. batch)\n")
    p= ggplot(data, 
           aes(x = UMAP1,
               y = UMAP2,
               col = factor(experiment_name))) + 
    geom_point(size = .01)+ 
    theme_classic()+ 
    theme(legend.position = "top")
  }
  if (by == "stim"){
    cat("Visualizing by stimm\n")
    p= ggplot(data, 
              aes(x = UMAP1,
                  y = UMAP2,
                  col = factor(stim))) + 
      geom_point(size = .01)+ 
      theme_classic()+ 
      theme(legend.position = "top")
  }
  if (by == "cluster"){
    cat("Visualizing by phenographic cluster\n")
    p= ggplot(data, 
              aes(x = UMAP1,
                  y = UMAP2,
                  col = factor(cluster))) + 
      geom_point(size = .01)+ 
      theme_classic()+ 
      theme(legend.position = "top")
  }
  if (by == "channel"){
    if (is.null(channel)){
      stop("No channel specified")
    }
    cat("Visualizing by channel\n")
    min_col = quantile(data[[channel]], .01)
    max_col = quantile(data[[channel]], .99)
    print(min_col)
    print(max_col)
    col_var = paste0("`", channel, "`")

    p= ggplot(data, 
              aes_string(
                  x = 'UMAP1',
                  y = 'UMAP2',
                  col = col_var)) + 
      geom_point(size = .01)+ 
      theme_classic()+ 
      theme(legend.position = "top") + 
      scale_color_viridis_c(limits = c(min_col, max_col))
  }
  
  print("Visualize UMAP")
  return(p)
  
}