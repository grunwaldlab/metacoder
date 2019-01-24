#' Reformat heat_tree arguments to taxmap
#'
#' Reformat heat_tree arguments to taxmap
#' 
#' @param obj The input taxmap object
#' @param arguments The validated arguments for heat_tree in a list.
#' 
#' @return a taxmap
#' 
#' @keywords internal
heat_tree_init_taxmap <- function(obj, arguments)
{
  # Transfer the taxonomy to a new taxmap object
  output <- taxmap()
  output$taxa <- obj$taxa
  output$edge_list <- obj$edge_list
  
  # Put all plot aesthetics in a list
  is_aes_arg <- grepl(names(arguments), pattern = "(node|edge|tree).+(size|color)$")
  aes_args <- arguments[is_aes_arg]
  # aes_table <- tibble::tibble(taxon_id = unlist(purrr::map(aes_args[! purrr::map_lgl(aes_args, is.null)], names)),
  #                             setting = rep(names(aes_args), purrr::map_int(aes_args, length)),
  #                             user = unname(do.call(c, aes_args)))
  
  # Put all labels in a list
  is_label_arg <- names(arguments) %in% c("node_label", "edge_label", "tree_label")
  label_args <- arguments[is_label_arg]
  # label_table <- tibble::tibble(taxon_id = unlist(purrr::map(label_args[! purrr::map_lgl(label_args, is.null)], names)),
  #                               setting = rep(names(label_args), purrr::map_int(label_args, length)),
  #                               user = unlist(label_args))
  
  # Put all transforamtion arguments in a list
  is_trans_arg <- grepl(names(arguments), pattern = "trans$")
  trans_list <-  arguments[is_trans_arg]
  names(trans_list) <- sub(names(trans_list), pattern = "_trans$", replacement = "")
  
  # Put all range arguments in a list
  is_range_arg <- grepl(names(arguments), pattern = "range$")
  range_list <-  arguments[is_range_arg]
  names(range_list) <- sub(names(range_list), pattern = "_range$", replacement = "")
  
  # Put all interval arguments in a list
  is_interval_arg <- grepl(names(arguments), pattern = "interval$")
  interval_list <-  arguments[is_interval_arg]
  names(interval_list) <- sub(names(interval_list), pattern = "_interval$", replacement = "")
  
  # Add all data to output taxmap
  output$data <- list(
    aes = aes_args,
    labels = label_args,
    trans = trans_list,
    ranges = range_list,
    intervals = interval_list
  )
  
  return(output)
}


#' Determine sizes used in heat_tree
#' 
#' Determine the sizes used in a heat_tree
#' 
#' @param obj a taxmap
#' 
#' @return a taxmap
#' 
#' @keywords internal
heat_tree_sizes <- function(obj) {
  size_params <- grep("size", names(obj$data$aes_trans), value = TRUE)
  obj$data$aes_table <- tibble::tibble(.rows = length(obj$taxa))
  obj$data$aes_table[size_params] <- purrr::map(obj$data$aes_trans[size_params], function(x) group_list_by_taxa(obj, x, as.numeric))
}


#' Determine colors used in heat_tree
#' 
#' Determine the colors used in a heat_tree
#' 
#' @param obj a taxmap
#' 
#' @return a taxmap
#' 
#' @keywords internal
heat_tree_colors <- function(obj) {
  print(1:3)
}


#' Group list items by taxon id
#' 
#' Group list items by taxon id and fill empty ids with NULL
#' 
#' @param obj a taxmap
#' @param a_list a list to convert
#' @param def_type_func An as.* type function to corece values to a certain type. 
#' 
#' @return A list corresponding to taxa in obj
#' 
#' @keywords internal
group_list_by_taxa <- function(obj, a_list, def_type_func) {
  output <- purrr::map(obj$taxon_ids(), function(tid) {
    if (tid %in% names(a_list)) {
      tid_items <- def_type_func(unlist(a_list[tid == names(a_list)]))
      # unique_types <- unique(purrr::map_chr(obj$data, function(x) class(x)[1]))
      # if (length(unique_types) != 1) {
      #   stop(call. = FALSE, '')
      # }
    } else {
      tid_items <- def_type_func(character(0))
    }
  })
  names(output) <- obj$taxon_ids()
  return(output)
}

