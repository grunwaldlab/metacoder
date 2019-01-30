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
  is_aes_arg <- grepl(names(arguments), pattern = "(node|edge|tree).+(size|color|shape)$")
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
  
  # make a table with standardized versions of size parameters
  size_params <- grep("size", names(obj$data$aes_trans), value = TRUE)
  obj$data$aes_table <- tibble::as.tibble(purrr::map(obj$data$aes_trans[size_params],
                                                     function(x) group_list_by_taxa(obj, x, as.numeric)))
  
  # transform to proportion of plot size and apply intervals
  obj$data$aes_table[size_params] <- purrr::map(size_params, function(param) {
    if (param %in% names(obj$data$intervals)) {
      interval <- obj$data$intervals[[param]]
    } else {
      interval <- range(get_numerics(obj$data$aes_table[[param]]), na.rm = TRUE, finite = TRUE)
    }
    purrr::map(obj$data$aes_table[[param]], function(values) {
      rescale(values, to = c(0, 1), from = interval)
    })
  })
  
  return(obj)
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
  
  # Identify color parameters
  color_params <- grep("color", names(obj$data$aes_trans), value = TRUE)

  # Apply intervals and rescale to 0-1 for numeric values
  color_data <- purrr::map(color_params, function(param) {
    purrr::map(obj$data$aes_trans[[param]], function(values) {
      if (is.numeric(values)) {
        return(apply_color_scale(values,
                                 color_series = obj$data$ranges[[param]],
                                 interval = obj$data$intervals[[param]])) 
      } else {
        return(values)
      }
    })
  })
  names(color_data) <- color_params
  
  # Reformat to per-taxon list 
  color_data <- purrr::map(color_data, function(x) group_list_by_taxa(obj, x, as.character))
  
  # Add results to aesthetics table
  obj$data$aes_table[color_params] <- color_data
  return(obj)
}


#' Determine labels used in heat_tree
#' 
#' Determine the labels used in a heat_tree
#' 
#' @param obj a taxmap
#' 
#' @return a taxmap
#' 
#' @keywords internal
heat_tree_labels <- function(obj) {
  
  # Get label data
  label_data <- obj$data$labels
  
  # Reformat to per-taxon list 
  label_data <- purrr::map(label_data, function(x) group_list_by_taxa(obj, x, as.character))
  
  # Add to aesthetics table
  obj$data$aes_table[names(label_data)] <- label_data
  return(obj)
}


#' Determine shpaes used in heat_tree
#' 
#' Determine the shapes used in a heat_tree
#' 
#' @param obj a taxmap
#' 
#' @return a taxmap
#' 
#' @keywords internal
heat_tree_shapes <- function(obj) {
  
  node_shapes <- obj$data$aes$node_shape
  edge_shapes <- obj$data$aes$edge_shape
  
  # Convert factors to shapes
  is_factor <- purrr::map_lgl(node_shapes, is.factor)
  node_shapes[is_factor] <- purrr::map_chr(node_shapes[is_factor], function(x) {
    names(obj$data$ranges$node_shape)[as.integer(x)]
  })
  is_factor <- purrr::map_lgl(edge_shapes, is.factor)
  edge_shapes[is_factor] <- purrr::map_chr(edge_shapes[is_factor], function(x) {
    names(obj$data$ranges$edge_shap)[as.integer(x)]
  })
  
  # Reformat to per-taxon list 
  node_shapes <- group_list_by_taxa(obj, node_shapes, as.character)
  edge_shapes <- group_list_by_taxa(obj, edge_shapes, as.character)
  
  # Add results to aesthetics table
  obj$data$aes_table[c("node_shape", "edge_shape")] <- list(node_shapes, edge_shapes)
  return(obj)
  
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


#' Return numeric values from a list or vector
#' 
#' Return numeric values from a list or vector
#' 
#' @param x A list or vector
#' 
#' @return a numeric vector
#' 
#' @keywords internal
get_numerics <- function(x) {
  if (is.vector(x) && is.numeric(x)) {
    return(x)
  }
  if (is.list(x)) {
    is_num <- purrr::map_lgl(x, is.numeric)
    if (sum(is_num) > 0) {
      return(unlist(x[is_num]))
    }
  }
  return(numeric(0))
}


#' Return character values from a list or vector
#' 
#' Return character values from a list or vector
#' 
#' @param x A list or vector
#' 
#' @return a character vector
#' 
#' @keywords internal
get_characters <- function(x) {
  if (is.vector(x) && is.character(x)) {
    return(x)
  }
  if (is.list(x)) {
    is_char <- purrr::map_lgl(x, is.character)
    if (sum(is_char) > 0) {
      return(unlist(x[is_char]))
    }
  }
  return(character(0))
}


#' Return factors from a list or vector
#' 
#' Return factors from a list or vector
#' 
#' @param x A list or vector
#' 
#' @return a factor vector
#' 
#' @keywords internal
get_factors <- function(x) {
  if (is.vector(x) && is.factor(x)) {
    return(x)
  }
  if (is.list(x)) {
    is_num <- purrr::map_lgl(x, is.factor)
    if (sum(is_num) > 0) {
      return(unlist(x[is_num]))
    }
  }
  return(factor())
}


#' Check if a list/vector is categorical
#' 
#' Check if a list/vector that is the value for a color aesthetic is composed of factors and perhaps characters, but not numbers.
is_categorical <- function(x) {
  return(length(get_factors(x)) > 1 && length(get_numerics(x)) == 0)
}