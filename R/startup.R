#' Return startup message
#' 
#' Return startup message
#' 
#' @keywords internal
startup_msg <- function() {
  my_version <- utils::packageVersion("metacoder")
  is_devel <- stringr::str_count(as.character(my_version), "\\.") == 3
  paste0('This is metacoder version ', my_version, ' ',
         ifelse(is_devel, crayon::bold("(development version)"), "(stable)"))#,
         # '. If you use metacoder for published research, please cite our paper:\n\n',
         # 'Foster Z, Sharpton T and Grunwald N (2017). "Metacoder: An R package for',
         # ' visualization and manipulation of community taxonomic diversity data." ',
         # crayon::italic('PLOS Computational Biology'), ', ', crayon::bold('13'),
         # '(2), pp. 1-15. doi: 10.1371/journal.pcbi.1005404\n\n',
         # 'Enter `citation("metacoder")` for a BibTeX entry for this citation.')
}


#' Run when package loads
#' 
#' Run when package loads
#' 
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(startup_msg())
}
