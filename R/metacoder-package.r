#===================================================================================================
#' Metacoder
#' 
#' A package for planning and analysis of amplicon metagenomics research projects.
#' 
#' The goal of the \code{metacoder} package is to provide a set of tools for:
#' 
#' \itemize{
#'   \item Standardized parsing of taxonomic information from diverse resources.
#'   \item Visualization of statistics distributed over taxonomic classifications and phylogenies.
#'   \item Evaluating potential metabarcoding primers for taxonomic specificity.
#'   \item Evaluating potential metabarcoding loci for taxonomic discrimiation (in development).
#' }
#' 
#' To accomplish these goals, `metacoder` leverages resources from other R packages, interfaces with
#' external programs, and provides novel functions where needed to allow for entire analyses within R.
#'
#' To learn how to use the package and what it can do view the package vignettes by typing: 
#' \code{browseVignettes("metacoder")}
#'
#' @section Most important functions:
#' 
#' \strong{Parseing taxonomy information:}
#' 
#' \itemize{
#'   \item \code{\link{extract_taxonomy}}
#' }
#' 
#' \strong{Dplyr-style manipulations of taxonomic data:}
#' 
#' \itemize{
#'   \item \code{\link{arrange_items}}
#'   \item \code{\link{arrange_taxa}}
#'   \item \code{\link{filter_items}}
#'   \item \code{\link{filter_taxa}}
#'   \item \code{\link{mutate_items}}
#'   \item \code{\link{mutate_taxa}}
#'   \item \code{\link{transmute_items}}
#'   \item \code{\link{transmute_taxa}}
#'   \item \code{\link{sample_n_items}}
#'   \item \code{\link{sample_n_taxa}}
#'   \item \code{\link{sample_frac_items}}
#'   \item \code{\link{sample_frac_taxa}}
#'   \item \code{\link{select_items}}
#'   \item \code{\link{select_taxa}}
#' }
#' 
#' \strong{Taxonomically balanced sub-sampling:}
#' 
#' \itemize{
#'   \item \code{\link{taxonomic_sample}}
#' }
#' 
#' \strong{Plotting taxonomic distribution of data:}
#' 
#' \itemize{
#'   \item \code{\link{plot.classified}}
#' }
#'
#' \strong{In silico PCR:}
#' 
#' \itemize{
#'   \item \code{\link{primersearch}}
#' }
#'
#' @author Zachary Foster
#' @name metacoder
#' @docType package
NULL
