#' Calculate proportions from observation counts
#'
#' For a given table in a taxmap object, convert one or more columns containing
#' counts to proportions. This is meant to be used with counts associated with
#' observations (e.g. OTUs), as opposed to counts that have already been summed
#' per taxon. If counts have already been summed per taxon, use
#' [calc_tax_props()] instead.
#'
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{dataset} to use. By
#'   default, all numeric columns are used. Takes one of the following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All/No columns will used.}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Use the columns
#'   corresponding to \code{TRUE} values.}
#'     \item{Character vector:}{The names of columns to use}
#'     \item{Numeric vector:}{The indexes of columns to use}
#'   }
#' @param other_cols Preserve in the output non-target columns present in the
#'   input data. The "taxon_id" column will always be preserved. Takes one of
#'   the following inputs:
#'   \describe{
#'     \item{TRUE/FALSE:}{All non-target columns will be preserved or not.}
#'     \item{Vector of TRUE/FALSE of length equal to the number of columns:}{Preserve the columns
#'   corresponding to \code{TRUE} values.}
#'     \item{Character vector:}{The names of columns to preserve}
#'     \item{Numeric vector:}{The indexes of columns to preserve}
#'   }
#' @param new_names If supplied, rename the output proportion columns. Must be
#'   the same length as \code{cold}.
#'
#' @return A tibble
#'
#' @family calculations
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Parse dataset for examples
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#'                    
#' # Calculate proportions for all numeric columns
#' calc_obs_props(x, "tax_data")
#' 
#' # Calculate proportions for a subset of columns
#' calc_obs_props(x, "tax_data", cols = c("700035653", "700097433", "700100489"))
#' calc_obs_props(x, "tax_data", cols = 4:6)
#' calc_obs_props(x, "tax_data", cols = startsWith(colnames(x$data$tax_data), "7"))
#' 
#' # Including all other columns in ouput
#' calc_obs_props(x, "tax_data", other_cols = TRUE)
#' 
#' # Inlcuding specific columns in output
#' calc_obs_props(x, "tax_data", cols = c("700035653", "700097433", "700100489"),
#'                other_cols = 2:3)
#' 
#' 
#' }
calc_obs_props <- function(obj, dataset, cols = NULL, other_cols = FALSE, new_names = NULL) {
  
  # Get count table
  count_table <- get_taxmap_table(obj, dataset)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  }
  
  # Check that count columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  # Check that new column names are the same length as calculation columns
  if (is.null(new_names)) {
    new_names <- cols
  } else if (length(new_names) != length(cols)) {
    stop(paste0('The `new_names` option (length = ', length(new_names), 
                ') must be the same length as the `cols` used (length = ',
                length(cols), ').'))
  }
  
  # Find other columns
  #   These might be added back to the output later
  cols_to_keep <- get_taxmap_other_cols(obj, dataset, cols, other_cols)
  
  # Calculate proportions
  prop_data <- do.call(cbind, lapply(count_table[cols], function(x) x / sum(x)))
  colnames(prop_data) <- new_names
  
  # Add back other columns if specified
  prop_data <- cbind(count_table[, cols_to_keep], prop_data)
  
  # Convert to tibble
  prop_data <- dplyr::as_tibble(prop_data)
  
  return(prop_data)
}


#' Calculate proportions from observation counts
#' 
#' For a given table in a taxmap object, rarefy counts to a constant total. This
#' is a wrapper around \code{\link[vegan]{rrarefy}} that automatically detects
#' which columns are numeric and handles the reformatting needed to use tibbles.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} that have counts. By 
#'   Default, all numeric columns in \code{data} are used.
#' @param sample_size The sample size counts will be rarefied to. This can be 
#'   either a single integer or a vector of integers of equal length to the 
#'   number of columns.
#' @param other_cols If \code{TRUE}, keep non-count cols in the input data.
#'   The "taxon_id" column will always be preserved.
#'   
#' @return A tibble
#' 
#' @family calculations
#'   
#' @export
rarefy_obs <- function(obj, dataset, cols = NULL, sample_size = NULL, other_cols = FALSE) {
  
  # Get count table
  count_table <- get_taxmap_table(obj, dataset, cols)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  } else { # remove any columns that do not exist
    cols <- cols[! cols %in% get_invalid_cols(count_table, cols)]
  }
  
  # Check that columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  # Calculate minimum count if no sample size is given
  if (is.null(sample_size)) {
    sample_size <- min(colSums(count_table[cols]))
  }
  
  # Rarefy
  count_table[cols] <- dplyr::as.tbl(as.data.frame(t(vegan::rrarefy(t(count_table[cols]), sample = sample_size))))

  # Remove other columns if specified
  if (! other_cols) {
    cols_to_keep <- c(colnames(count_table[cols]), "taxon_id")
    count_table <- count_table[colnames(count_table) %in% cols_to_keep]
  }
  
  return(count_table)
}


#' Compare treatments
#' 
#' Apply a function to compare data, usually abundance, from pairs of 
#' treatments. By default, every pairwise combination of treatments are 
#' compared. A custom function can be supplied to perform the comparison.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains data for each 
#'   sample in columns.
#' @param sample_ids The names of sample columns in \code{dataset}
#' @param treatments The treatment associated with each sample. Must be the same
#'   order and length as \code{sample_ids}.
#' @param func The function to apply for each comparison. For each row in 
#'   \code{dataset}, for each combination of treatments, this function will 
#'   recieve the data for each treatment, passed a two character vecotors.
#'   Therefor the function must take at least 2 arguments corresponding to the
#'   two treatments compared. The function should return a vector or list or
#'   results of a fixed length. If named, the names will be used in the output.
#'   The names should be consistent as well. A simple example is
#'   \code{function(x, y) mean(x) - mean(y)}. By default, the following function
#'   is used:
#'   \preformatted{
#'   function(abund_1, abund_2) {
#'     log_ratio <- log2(median(abund_1) / median(abund_2))
#'     if (is.nan(log_ratio)) {
#'       log_ratio <- 0
#'     }
#'     list(log2_median_ratio = log_ratio,
#'          median_diff = median(abund_1) - median(abund_2),
#'          mean_diff = mean(abund_1) - mean(abund_2),
#'          wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)
#'   }
#'   }
#' @param combinations Which combinations of treatments to use. Must be a list 
#'   of vectors, each containing the names of 2 treatments to compare. By 
#'   default, all pairwise combinations of treatments are compared.
#' @param other_cols If \code{TRUE}, preserve all columns not in 
#'   \code{sample_ids} in the output. If \code{FALSE}, dont keep other columns. 
#'   If a column names or indexes are supplied, only preserve those columns.
#'   
#' @return A tibble
#'   
#' @family calculations
#' 
#' @examples
#' \dontrun{
#' # Parse dataset for plotting
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#' 
#' # Convert counts to proportions
#' x$data$otu_table <- calc_obs_props(x, dataset = "tax_data", cols = hmp_samples$sample_id)
#' 
#' # Get per-taxon counts
#' x$data$tax_table <- calc_taxon_abund(x, dataset = "otu_table", cols = hmp_samples$sample_id)
#' 
#' # Calculate difference between treatments
#' x$data$diff_table <- compare_treatments(x, dataset = "tax_table",
#'                                         sample_ids = hmp_samples$sample_id,
#'                                         treatments = hmp_samples$body_site)
#'
#' # Plot results (might take a few minutes)
#' heat_tree_matrix(x,
#'                  dataset = "diff_table",
#'                  node_size = n_obs,
#'                  node_label = taxon_names,
#'                  node_color = log2_median_ratio,
#'                  node_color_range = diverging_palette(),
#'                  node_color_trans = "linear",
#'                  node_color_interval = c(-3, 3),
#'                  edge_color_interval = c(-3, 3),
#'                  node_size_axis_label = "Number of OTUs",
#'                  node_color_axis_label = "Log2 ratio median proportions")
#' 
#' }
#' 
#' @export
compare_treatments <- function(obj, dataset, sample_ids, treatments,
                               func = NULL, combinations = NULL,
                               other_cols = FALSE) {
  # Get abundance by sample data
  abund_data <- get_taxmap_table(obj, dataset, sample_ids)
  
  # Define defualt function
  if (is.null(func)) {
    func <- function(abund_1, abund_2) {
      log_ratio <- log2(stats::median(abund_1) / stats::median(abund_2))
      if (is.nan(log_ratio)) {
        log_ratio <- 0
      }
      list(log2_median_ratio = log_ratio,
           median_diff = stats::median(abund_1) - stats::median(abund_2),
           mean_diff = mean(abund_1) - mean(abund_2),
           wilcox_p_value = stats::wilcox.test(abund_1, abund_2)$p.value)
    }
  }
  
  # Account for columns that do not exist
  invalid_cols <- get_invalid_cols(abund_data, sample_ids)
  treatments <- treatments[! sample_ids %in% invalid_cols]
  sample_ids <- sample_ids[! sample_ids %in% invalid_cols]
  
  # Parse "other_cols" option
  if (is.logical(other_cols)) {
    if (length(other_cols) != 1) {
      stop('The "other_cols" option must either be TRUE/FALSE or a vector of valid column names/indexes.', call. = FALSE)
    } else if (other_cols) {
      other_cols <- colnames(abund_data)[! colnames(abund_data) %in% sample_ids]
    } else {
      other_cols <- c()
    }
  } else {
    not_used <- get_taxmap_table(obj, dataset, other_cols) # Checks that columns specified are valid
  }
  if (! "taxon_id" %in% other_cols) {
    other_cols <- c("taxon_id", other_cols)
  }
  
  # Get every combination of treatments to compare
  if (is.null(combinations)) {
    combinations <- t(utils::combn(unique(treatments), 2))
    combinations <- lapply(seq_len(nrow(combinations)), function(i) combinations[i, ])
  }
  
  # Make function to compare one pair of treatments
  one_comparison <- function(treat_1, treat_2) {
    output <- lapply(seq_len(nrow(abund_data)), function(i) {
      # get samples ids for each treatment
      samples_1 <- sample_ids[treatments == treat_1]
      samples_2 <- sample_ids[treatments == treat_2]
      
      # get abundance data for each treatment
      abund_1 <- unlist(abund_data[i, samples_1])
      abund_2 <- unlist(abund_data[i, samples_2])
      
      # Run user-supplied function on abundance info
      result <- as.list(func(abund_1, abund_2))
      
      # Complain if nothing is returned
      if (length(result) == 0) {
        stop(paste0("The function supplied returned nothing when given the following data:\n",
                    '  Treatments compared: "', treat_1, '" vs "', treat_2, '"\n',
                    '  Row index: ', i, '\n',
                    limited_print(abund_1, prefix = paste0('  "', treat_1, '" data: '), type = "silent"),
                    limited_print(abund_2, prefix = paste0('  "', treat_2, '" data: '), type = "silent")))
      }
      
      # If the output does not have names, add defaults
      if (is.null(names(result))) {
        if (length(result) == 1) {
          names(result) <- "value"
        } else {
          names(result) <- paste("value_", seq_along(result)) 
        }
      }
      
      do.call(data.frame, result)
    })
    
    # Combine the results for each row into a data frame
    output <- as.data.frame(do.call(rbind, output))
    
    # Add treatments compared
    output <- cbind(data.frame(stringsAsFactors = FALSE,
                               treatment_1 = treat_1,
                               treatment_2 = treat_2),
                    output)
    
    # Add in other columns if specified
    output <- cbind(abund_data[, other_cols], output)
    
    return(output)
  }
  
  # Compare all pairs of treatments and combine
  output <- lapply(seq_len(length(combinations)), function(i) {
    one_comparison(combinations[[i]][1], combinations[[i]][2])
  })
  output <- do.call(rbind, output)
  
  # Convert to tibble and return
  dplyr::as.tbl(output)
}


#' Sum observation values for each taxon
#' 
#' For a given table in a taxmap object, sum the values in each column for each 
#' taxon. This is useful to convert per-observation counts (e.g. OTU counts) to 
#' per-taxon counts.
#' 
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} that have counts. By 
#'   Default, all numeric columns in \code{data} are used.
#' @param groups Group counts of multiple columns per treatment/group. This
#'   should be a vector of group IDs (e.g. character, integer) the same
#'   length as \code{cols} that defines which samples go in which group. When
#'   used, there will be one column in the output for each unique value in
#'   \code{groups}.
#' @param out_names The names of count columns in the output. Must be the same
#'   length as \code{cols} (or \code{unique(groups)}, if used).
#'   
#' @return A tibble
#'   
#' @family calculations
#'   
#' @export
calc_taxon_abund <- function(obj, dataset, cols = NULL, groups = NULL,
                             out_names = NULL) {
  
  # Get count table
  count_table <- get_taxmap_table(obj, dataset, cols)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  }  else { # remove any columns that do not exist
    cols <- cols[! cols %in% get_invalid_cols(count_table, cols)]
  }

  # Check that count columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  # Check that groups and output names make sense
  if (is.null(groups)) {
    if (is.null(out_names)) { # groups and out_names are NULL
      out_names <- colnames(count_table[cols])
    } else { # groups is NULL, but out_names set
      if (length(out_names) != length(cols)) {
        stop("The length of 'cols' and 'out_names' are not equal")
      }
    }
    groups <- seq_along(cols)
  } else {
    if (length(groups) != length(cols)) {
      stop("`groups` must be the same length as cols.")
    }
    if (is.null(out_names)) { # groups is set, but out_names is NULL
      if (is.numeric(groups)) {
        warning("Numeric groups used without supplying out_name. This will result in numeric column names.")
      }
      out_names <- unique(groups)
    } else { # groups and out_names are both set
      if (length(out_names) != length(unique(groups))) {
        stop("The length of 'unique(groups)' and 'out_names' are not equal")       
      }
    }
  }
  
  # Check that out_names is a character
  if (! is.null(out_names) && is.numeric(out_names)){
    warning("`out_names` is numeric. This will result in numeric column names.")
  }
  
  # Sum counts per taxon for each sample
  obs_indexes <- obj$obs(dataset)
  output <- lapply(split(cols, groups), function(col_index) {
    vapply(obs_indexes, function(i) sum(unlist(count_table[i, col_index])), numeric(1))
  })
  output <- as.data.frame(output, stringsAsFactors = FALSE)
  
  # Add taxon_id column
  output <- cbind(data.frame(taxon_id = obj$taxon_ids(), stringsAsFactors = FALSE),
                  output)
  
  # Rename cols
  colnames(output)[-1] <- out_names
  
  # Convert to tibble and return
  dplyr::as.tbl(output)
}



#' Count the number of samples with reads
#'
#' For a given table in a taxmap object, count the number of samples with
#' greater than zero reads.
#'
#' @param obj A taxmap object
#' @param dataset The name of a table in \code{obj} that contains counts.
#' @param cols The names/indexes of columns in \code{data} that have counts. By
#'   Default, all numeric columns in \code{data} are used.
#' @param groups Group counts of multiple columns per treatment/group. This
#'   should be a vector of group IDs (e.g. character, integer) the same length
#'   as \code{cols} that defines which samples go in which group. When used,
#'   there will be one column in the output for each unique value in
#'   \code{groups}.
#' @param out_names The names of count columns in the output. Must be the length
#'   1 or same length as \code{unique(groups)}, if used.
#' @param drop If \code{groups} is not used, return a vector of the results instead
#'   of a table with one column.
#' @param append If \code{TRUE}, append results to input table and return.
#'
#' @return A tibble
#'
#' @family calculations
#'
#' @examples
#' \dontrun{
#' # Parse dataset for example
#' x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#'                    class_key = c(tax_rank = "info", tax_name = "taxon_name"),
#'                    class_regex = "^(.+)__(.+)$")
#'                    
#' # Count samples with reads
#' calc_n_samples(x, dataset = "tax_data")
#' 
#' # Return a vector instead of a table
#' calc_n_samples(x, dataset = "tax_data", drop = TRUE)
#' 
#' # Only use some columns
#' calc_n_samples(x, dataset = "tax_data", cols =  hmp_samples$sample_id[1:5])
#' 
#' # Return a count for each treatment
#' calc_n_samples(x, dataset = "tax_data", groups = hmp_samples$body_site)
#' 
#' # Rename output columns 
#' calc_n_samples(x, dataset = "tax_data", groups = hmp_samples$body_site,
#'                out_names = c("A", "B", "C", "D", "E"))
#' 
#' # Add results to input table
#' calc_n_samples(x, dataset = "tax_data", append = TRUE)
#' }
#' 
#' @export
calc_n_samples <- function(obj, dataset, cols = NULL, groups = NULL,
                           out_names = NULL, drop = FALSE, append = FALSE) {
  # Get count table
  count_table <- get_taxmap_table(obj, dataset, cols)
  
  # Find default columns if needed
  if (is.null(cols)) {
    cols <- which(vapply(count_table, is.numeric, logical(1)))
  }  else { # remove any columns that do not exist
    cols <- cols[! cols %in% get_invalid_cols(count_table, cols)]
  }
  
  # Check that count columns are numeric
  col_is_num <- vapply(count_table[cols], is.numeric, logical(1))
  if (! all(col_is_num)) {
    stop(paste0("All columns must be numeric. The following columns are not numeric:  ",
                limited_print(cols[!col_is_num], type = "silent")))
  }
  
  # Check that groups and output names make sense
  if (is.null(groups)) {
    if (is.null(out_names)) { # groups and out_names are NULL
      out_names <- "n_samples"
    } else { # groups is NULL, but out_names set
      if (length(out_names) != length(cols)) {
        stop("The length of 'cols' and 'out_names' are not equal")
      }
    }
    groups <- rep("placeholder", length(cols))
  } else {
    if (length(groups) != length(cols)) {
      stop("`groups` must be the same length as cols.")
    }
    if (length(unique(groups)) != 1 && drop) {
      stop("Cannot drop second dimension since there are multiple groups specified.")
    }
    if (is.null(out_names)) { # groups is set, but out_names is NULL
      if (is.numeric(groups)) {
        warning("Numeric groups used without supplying out_name. This will result in numeric column names.")
      }
      out_names <- unique(groups)
    } else { # groups and out_names are both set
      if (length(out_names) != length(unique(groups))) {
        stop("The length of 'unique(groups)' and 'out_names' are not equal")       
      }
    }
  }
  
  # Check that out_names is a character
  if (! is.null(out_names) && is.numeric(out_names)){
    warning("`out_names` is numeric. This will result in numeric column names.")
  }
  
  # Count number of samples
  output <- lapply(split(cols, groups), function(col_index) {
    vapply(seq_len(nrow(count_table)), function(i) sum(count_table[i, col_index] > 0), integer(1))
  })
  output <- as.data.frame(output, stringsAsFactors = FALSE)
  
  # Drop second dimension
  if (drop) {
    output <- output[[1]]
    names(output) <- count_table$taxon_id
    return(output)
  } else {
    # Rename cols
    colnames(output) <- out_names
    
    # Add taxon_id column
    if (append) {
      output <- cbind(get_taxmap_table(obj, dataset),
                      output)
    } else {
      output <- cbind(data.frame(taxon_id = count_table$taxon_id, stringsAsFactors = FALSE),
                      output)
    }
    
    # Convert to tibble and return
    return(dplyr::as.tbl(output))
  }
  
}



