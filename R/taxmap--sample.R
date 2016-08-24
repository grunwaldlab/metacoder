#' Sample n observations from \code{\link{taxmap}}
#' 
#' Randomly sample some number of observations from a \code{\link{taxmap}} object. Weights can be 
#' specified for observations or the taxa they are taxmap by.
#' See \link[dplyr]{sample_n} for the inspiration for this function.
#' 
#' @param .data (\code{\link{taxmap}}) The object to sample from.
#' @param size (\code{numeric} of length 1) The number of observations to sample.
#' @param replace (\code{logical} of length 1) If \code{TRUE}, sample with replacement.
#' @param taxon_weight (\code{numeric}) Non-negative sampling weights of each taxon. If 
#'   \code{use_supertaxa} is \code{TRUE}, the weights for each taxon in an observation's classification are
#'   supplied to \code{collapse_func} to get the observation weight. The expression given is evaluated in 
#'   the context of \code{\link{taxon_data}}. In other words, any column name that appears in 
#'   \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. If 
#'   \code{obs_weight} is also specified, the two weights are multiplied (after \code{taxon_weight}
#'   for each observation is calculated).
#' @param obs_weight (\code{numeric}) Sampling weights of each observation. The expression given is 
#'   evaluated in the context of \code{\link{obs_data}}. In other words, any column name that 
#'   appears in \code{\link{obs_data}(.data)} can be used as if it was a vector on its own. If 
#'   \code{taxon_weight} is also specified, the two weights are multiplied (after 
#'   \code{taxon_weight} for each observation is calculated).
#' @param use_supertaxa (\code{logical} of length 1) Affects how the \code{taxon_weight} is used. If
#'   \code{TRUE}, the weights for each taxon in an observation's classification are multiplied to get the 
#'   observation weight. Otherwise, just the taxonomic level the observation is assign to it considered.
#' @param collapse_func (\code{function} of length 1) If \code{taxon_weight} option is used and 
#'   \code{supertaxa} is \code{TRUE}, the weights for each taxon in an observation's classification are 
#'   supplied to \code{collapse_func} to get the observation weight. This function should take  numeric 
#'   vector and return a single number.
#' @param ... Additional options are passed to \code{\link{filter_obs}}.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#' 
#' @examples
#' # Subsample without replacement, keeping all taxa
#' sample_n_obs(unite_ex_data_3, 100)
#' # Subsample without replacement and remove unsampled taxa
#' sample_n_obs(unite_ex_data_3, 100, unobserved = FALSE)
#' # Subsample with taxon weight 
#' sample_n_obs(unite_ex_data_3, 100, unobserved = FALSE, taxon_weight = 1 / n_obs)
#' # Sample with replacement
#' sample_n_obs(unite_ex_data_3, 10000, replace = TRUE)
#'   
#' @export
sample_n_obs <- function(.data, size, replace = FALSE, taxon_weight = NULL, obs_weight = NULL,
                           use_supertaxa = TRUE, collapse_func = mean, ...) {
  # Calculate taxon component of taxon weights -----------------------------------------------------
  my_taxon_data <- taxon_data(.data)
  taxon_weight <- lazyeval::lazy_eval(lazyeval::lazy(taxon_weight), data = my_taxon_data)
  if (is.null(taxon_weight)) {
    obs_taxon_weight <- rep(1, nrow(.data$obs_data))
  } else {
    obs_index <- match(.data$obs_data$obs_taxon_ids, .data$taxon_data$taxon_ids)
    my_supertaxa <- supertaxa(.data, recursive = use_supertaxa, simplify = FALSE,
                              include_input = TRUE, index = TRUE, na = FALSE)
    taxon_weight_product <- vapply(my_supertaxa, function(x) collapse_func(taxon_weight[x]), numeric(1))
    obs_taxon_weight <- taxon_weight_product[obs_index]
  }
  obs_taxon_weight <- obs_taxon_weight / sum(obs_taxon_weight)
  
  # Calculate observation component of observation weights -------------------------------------------------------
  my_obs_data <- obs_data(.data)
  obs_weight <- lazyeval::lazy_eval(lazyeval::lazy(obs_weight), data = my_obs_data)
  if (is.null(obs_weight)) {
    obs_weight <- rep(1, nrow(.data$obs_data)) 
  }
  obs_weight <- obs_weight / sum(obs_weight)
  
  # Combine observation and taxon weight components  ------------------------------------------------------
  combine_func <- prod
  weight <- mapply(obs_taxon_weight, obs_weight, FUN = function(x, y) combine_func(c(x,y)))
  weight <- weight / sum(weight)
  
  # Sample observations -----------------------------------------------------------------------------------
  sampled_rows <- sample.int(nrow(my_obs_data), size = size, replace = replace, prob = weight)
  filter_obs(.data, sampled_rows, ...)
}

#' Sample a proportion of observations from \code{\link{taxmap}}
#' 
#' Randomly sample some propoortion of observations from a \code{\link{taxmap}} object. Weights can be 
#' specified for observations or the taxa they are taxmap by.
#' See \link[dplyr]{sample_frac} for the inspiration for this function.
#' 
#' @inheritParams sample_n_obs
#' @param size (\code{numeric} of length 1) The proportion of observations to sample.
#' 
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#'   
#' @examples
#' # Subsample without replacement, keeping all taxa
#' sample_frac_obs(unite_ex_data_3, 0.1)
#' # Subsample without replacement and remove unsampled taxa
#' sample_frac_obs(unite_ex_data_3, 0.1, unobserved = FALSE)
#' # Subsample with taxon weight 
#' sample_frac_obs(unite_ex_data_3, 0.1, unobserved = FALSE, taxon_weight = 1 / n_obs)
#' # Sample with replacement
#' sample_frac_obs(unite_ex_data_3, 10, replace = TRUE)
#' 
#' @export
sample_frac_obs <- function(.data, size = 1, replace = FALSE, taxon_weight = NULL, obs_weight = NULL,
                           use_supertaxa = TRUE, collapse_func = mean, ...) {
  sample_n_obs(.data = .data, size = size * nrow(.data$obs_data), replace = replace,
                 taxon_weight = taxon_weight, obs_weight = obs_weight,
                 use_supertaxa = use_supertaxa, collapse_func = collapse_func, ...)
}


#' Sample n taxa from \code{\link{taxmap}}
#' 
#' Randomly sample some number of taxa from a \code{\link{taxmap}} object. Weights can be 
#' specified for taxa or the observations assigned to them.
#' See \link[dplyr]{sample_n} for the inspiration for this function.
#' 
#' @param .data (\code{\link{taxmap}}) The object to sample from.
#' @param size (\code{numeric} of length 1) The number of taxa to sample.
#' @param taxon_weight (\code{numeric}) Non-negative sampling weights of each taxon. The expression 
#'   given is evaluated in the context of \code{\link{taxon_data}}. In other words, any column name 
#'   that appears in \code{\link{taxon_data}(.data)} can be used as if it was a vector on its own. 
#'   If \code{obs_weight} is also specified, the two weights are multiplied (after 
#'   \code{obs_weight} for each taxon is calculated).
#' @param obs_weight (\code{numeric}) Sampling weights of each observation. The weights for each observation 
#'   assigned to a given taxon are supplied to \code{collapse_func} to get the taxon weight. If 
#'   \code{use_subtaxa} is \code{TRUE} then the observations assigned to every subtaxa are also used. The 
#'   expression given is evaluated in the context of \code{\link{obs_data}}. In other words, any 
#'   column name that appears in \code{\link{obs_data}(.data)} can be used as if it was a vector on
#'   its own. If \code{taxon_weight} is also specified, the two weights are multiplied (after 
#'   \code{obs_weight} for each observation is calculated).
#' @param use_subtaxa (\code{logical} of length 1) Affects how the \code{obs_weight} option is
#'   used. If \code{TRUE}, the weights for each taxon in an observation's classification are multiplied to
#'   get the observation weight. Otherwise, just the taxonomic level the observation is assign to it considered.
#' @param collapse_func (\code{function} of length 1) If \code{taxon_weight} is used and 
#'   \code{supertaxa} is \code{TRUE}, the weights for each taxon in an observation's classification are 
#'   supplied to \code{collapse_func} to get the observation weight. This function should take  numeric 
#'   vector and return a single number.
#' @param ... Additional options are passed to \code{\link{filter_taxa}}.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#' 
#' @examples 
#' # subsample taxa, preserving shared supertaxa
#' sample_n_taxa(unite_ex_data_3, 100, supertaxa = TRUE)
#' # subsample taxa using weights, preserving subtaxa
#' sample_n_taxa(unite_ex_data_3, 10, subtaxa = TRUE,
#'               taxon_weight = ifelse(unite_rank == "g" & n_subtaxa > 3, 1, 0))
#'
#' @export
sample_n_taxa <- function(.data, size, taxon_weight = NULL, obs_weight = NULL,
                          use_subtaxa = TRUE, collapse_func = mean, ...) {
  # Calculate observation component of taxon weights ------------------------------------------------------
  my_obs_data <- obs_data(.data)
  obs_weight <- lazyeval::lazy_eval(lazyeval::lazy(obs_weight), data = my_obs_data)
  if (is.null(obs_weight)) {
    taxon_obs_weight <- rep(1, nrow(.data$taxon_data))
  } else {
    my_obs <- obs(.data, recursive = use_subtaxa, simplify = FALSE)
    taxon_obs_weight <- vapply(my_obs, function(x) collapse_func(obs_weight[x]), numeric(1))
  }
  taxon_obs_weight <- taxon_obs_weight / sum(taxon_obs_weight)
  
  # Calculate taxon component of taxon weights -------------------------------------------------------
  my_taxon_data <- taxon_data(.data)
  taxon_weight <- lazyeval::lazy_eval(lazyeval::lazy(taxon_weight), data = my_taxon_data)
  if (is.null(taxon_weight)) {
    taxon_weight <- rep(1, nrow(.data$taxon_data)) 
  }
  taxon_weight <- taxon_weight / sum(taxon_weight)
  
  # Combine observation and taxon weight components  ------------------------------------------------------
  combine_func <- prod
  weight <- mapply(taxon_weight, taxon_obs_weight, FUN = function(x, y) combine_func(c(x,y)))
  weight <- weight / sum(weight)
  
  # Sample observations -----------------------------------------------------------------------------------
  sampled_rows <- sample.int(nrow(.data$taxon_data), size = size, replace = FALSE, prob = weight)
  filter_taxa(.data, sampled_rows, ...)
}


#' Sample a proportion of taxa from \code{\link{taxmap}}
#' 
#' Randomly sample some propoortion of taxa from a \code{\link{taxmap}} object. Weights can be 
#' specified for taxa or the observations assigned to them. See \link[dplyr]{sample_frac} for the
#' inspiration for this function.
#' 
#' @inheritParams sample_n_taxa
#' @param size (\code{numeric} of length 1) The proportion of taxa to sample.
#'   
#' @return An object of type \code{\link{taxmap}}
#'   
#' @family dplyr-like functions
#'   
#' @examples 
#' # subsample taxa, preserving shared supertaxa
#' sample_frac_taxa(unite_ex_data_3, 0.1, supertaxa = TRUE)
#' # subsample taxa using weights, preserving subtaxa
#' sample_frac_taxa(unite_ex_data_3, 0.01, subtaxa = TRUE,
#'                  taxon_weight = ifelse(unite_rank == "g" & n_subtaxa > 3, 1, 0))
#'               
#' @export
sample_frac_taxa <- function(.data, size = 1, taxon_weight = NULL, obs_weight = NULL,
                              use_subtaxa = TRUE, collapse_func = mean, ...) {
  sample_n_taxa(.data = .data, size = size * nrow(.data$obs_data),
                 taxon_weight = taxon_weight, obs_weight = obs_weight,
                use_subtaxa = use_subtaxa, collapse_func = collapse_func, ...)
}

