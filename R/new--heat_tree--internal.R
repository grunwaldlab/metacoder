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
  
  # Put all plot aesthetics in a table
  is_aes_arg <- grepl(names(arguments), pattern = "(node|edge|tree).+(size|color)$")
  aes_args <- arguments[is_aes_arg]
  aes_table <- tibble::tibble(taxon_id = unlist(purrr::map(aes_args[! purrr::map_lgl(aes_args, is.null)], names)),
                              setting = rep(names(aes_args), purrr::map_int(aes_args, length)),
                              user = unlist(aes_args))
  
  # Put all labels in a table
  is_label_arg <- names(arguments) %in% c("node_label", "edge_label", "tree_label")
  label_args <- arguments[is_label_arg]
  label_table <- tibble::tibble(taxon_id = unlist(purrr::map(label_args[! purrr::map_lgl(label_args, is.null)], names)),
                                setting = rep(names(label_args), purrr::map_int(label_args, length)),
                                user = unlist(label_args))
  
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
    aes = aes_table,
    labels = label_table,
    trans = trans_list,
    ranges = range_list,
    intervals = interval_list
  )
  
  return(output)
}

