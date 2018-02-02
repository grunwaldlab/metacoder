#' A HMP subset
#'
#' A subset of the Human Microbiome Project abundance matrix produced by QIIME.
#' It contains OTU ids, taxonomic lineages, and the read counts for 50 samples.
#' See \code{\link{hmp_samples}} for the matching dataset of sample information.
#' 
#' The 50 samples were randomly selected such that there were 10 in each of 5
#' treatments: "Saliva", "Throat", "Stool", "Right_Antecubital_fossa",
#' "Anterior_nares". For each treatment, there were 5 samples from men and 5
#' from women. 
#'
#' @name hmp_otus
#' @format A 1,000 x 52 tibble.
#' @source Subset from data available at https://www.hmpdacc.org/hmp/HMQCP/
#' @family hmp_data
#' @keywords data
NULL


#' Sample information for HMP subset
#'
#' The sample information for a subset of the Human Microbiome Project data. It
#' contains the sample ID, sex, and body site for each sample in the abundance
#' matrix stored in \code{\link{hmp_otus}}. The "sample_id" column corresponds
#' to the column names of \code{\link{hmp_otus}}.
#'
#' The 50 samples were randomly selected such that there were 10 in each of 5
#' treatments: "Saliva", "Throat", "Stool", "Right_Antecubital_fossa",
#' "Anterior_nares". For each treatment, there were 5 samples from men and 5
#' from women. "Right_Antecubital_fossa" was renamed to "Skin" and
#' "Anterior_nares" to "Nose".
#'
#' @name hmp_samples
#' @format A 50 x 3 tibble.
#' @source Subset from data available at https://www.hmpdacc.org/hmp/HMQCP/
#' @family hmp_data
#' @keywords data
NULL
