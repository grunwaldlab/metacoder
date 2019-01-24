#' Apply transformations
#'
#' Apply transformations to user-supplied data contained in a taxmap object.
#' Used in the heat_tree function.
#' 
#' @param obj The input taxmap object
#' 
#' @return a taxmap
#' 
#' @keywords internal
heat_tree_transform_data <- function(obj) {
  
  obj$data$aes_trans <- purrr::map(names(obj$data$aes), function(n) {
    if (n %in% names(obj$data$trans)) {
      if (is.vector(obj$data$aes[[n]]) && is.numeric(obj$data$aes[[n]])) {
        return(obj$data$trans[[n]](obj$data$aes[[n]]))
      }
      if (is.list(obj$data$aes[[n]])) {
        return(purrr::map(obj$data$aes[[n]], obj$data$trans[[n]]))
      }
      warning(call. = FALSE, 'Could not transform parameter "', n, '"')
    }
    return(obj$data$aes[[n]])
  })
  names(obj$data$aes_trans) <- names(obj$data$aes)
  
  
  # obj$data$aes$transformed <- purrr::map2_dbl(obj$data$aes$setting, obj$data$aes$user,
  #                                             function(setting, value) {
  #                                               if (is.numeric(value)) {
  #                                                 obj$data$trans[[setting]](value)
  #                                               } else {
  #                                                 as.numeric(NA)
  #                                               }
  #                                             })
  return(obj)
}
