
#' @keywords internal
parse_ft <- function(text) {
  text <- gsub(text, pattern = "\n\t\t\t", replacement = "\t", fixed = TRUE)
  parts <- strsplit(text, "\n\n", fixed = TRUE)[[1]]
  part_data <- lapply(parts, function(x) {
    first_line <- sub(x, pattern = "\\n.+$", replacement = "")
    acc <- stringr::str_match(first_line, pattern = "\\|(.+)\\|")[,2]
    the_rest <- sub(x, pattern = "^>.+?\\n", replacement = "")
    # replace extra \t with , when there are multiple features
    lines <- strsplit(the_rest, "\n")[[1]]
    lines <- purrr::map2(stringr::str_locate_all(lines, "\t"), lines, function(matches, line) {
      if (nrow(matches) > 4) {
        for (i in matches[3:(nrow(matches) - 2), 1]) {
          substr(line, i, i) <- ","
        }
      }
      return(line)
    })
    the_rest <- paste0(lines, collapse = "\n")
    output <- readr::read_tsv(paste0(the_rest, "\n"), col_names = c("start", "end", "feature", "type", "name"), col_types = "ccccc")
    output <- tibble::as_tibble(cbind(list(acc = acc), output, stringsAsFactors = FALSE))
  })
  output <- dplyr::bind_rows(part_data)
  output$complete <- ifelse(startsWith(as.character(output$start), "<") | startsWith(as.character(output$end), ">"), FALSE, TRUE)
  output$start <- as.integer(gsub(output$start, pattern = "<", replacement = ""))
  output$end <- as.integer(gsub(output$end, pattern = ">", replacement = ""))
  return(output)
}

#' @keywords internal
parse_seqs <- function(text) {
  xml <- xml2::read_xml(text)
  tibble::tibble(acc = xml2::xml_text(xml2::xml_find_all(xml, "//TSeq_accver")),
                 seq = xml2::xml_text(xml2::xml_find_all(xml, "//TSeq_sequence")),
                 header = xml2::xml_text(xml2::xml_find_all(xml, "//TSeq_defline")),
                 length = xml2::xml_text(xml2::xml_find_all(xml, "//TSeq_length")))
}

#' Lookup gene sequences from NCBI
#'
#' Look for sequences of a particular gene for a list of species/isolates/genes from the Genbank
#' nucleotide database.
#'
#' @param species The names of species to look up.
#' @param genes The names of the genes to look up.
#' @param isolates The names of isolates to look up. Must be the same length as \code{species} if
#'   used.
#' @param extract_features If TRUE, return the sequence for each feature in the sequence annotation
#'   instead of the whole sequence.
#' @param gene_name_in_feature If TRUE, only return features that have one of the gene names
#'   somewhere in their description. Only has an effect if extract_features is TRUE.
#' @param flanking A vector of length 2. The number of base pairs before and after the target gene
#'   to include in the sequence returned. If the flanking sequence is not available, the sequence
#'   will be considered incomplete.
#' @param db The name of the NCBI database to query. Only tested with "nucleotide", but a few others
#'   might work.
#' @param pause The number of seconds to pause between each query. This avoids annoying NCBI and
#'   having them block you IP address. Should be at least 0.35 seconds if you dont have an NCBI API
#'   key and at least 0.1 seconds if you do.
#' @param ... Additional terms to add to the search request for each species/isolate, see NCBI
#'   documentation for a complete list:
#'   http://www.ncbi.nlm.nih.gov/books/NBK25499/#_chapter4_ESearch_
#'
#' @return
#'
#' A tibble (a type of data.frame) with the following columns, depending on settings:
#'
#' * species: The species search term associated with the result
#' * isolate: The isolate search term associated with the result
#' * query: The query used to search for sequences on NCBI
#' * acc: The Genbank accession number
#' * start: The index of the first base pair of the target gene in the original sequence
#' * end: The index of the last base pair of the target gene in the original sequence
#' * feature: The type of locus
#' * type: The type of annotation for the sequence
#' * name: The name of the locus
#' * complete: If the locus was complete in the original sequence, according to the annotation
#' * seq: The whole sequence the gene was found in
#' * header: The name of the overall sequence the gene was found in
#' * length: The length of the whole sequence
#' * flank_start: The index of the first base pair of the target gene plus the flanking region in the original sequence
#' * flank_end: The index of the last base pair of the target gene plus the flanking region in the original sequence
#' * flank_complete: If the locus plus flanking region was complete in the original sequence
#' * gene_seq: The target gene sequence plus the flanking region
#' * gene_length: The length of the target gene plus the flanking region
#'
#' @examples \dontrun{
#'
#' # Search for the whole seqeunces for with P. infestans Cox I
#' get_isolate_seqs(species = c("Phytophthora infestans"),
#'                  genes = c("cox I", "cox 1", "cox1", "coxI", "cytochrome oxidase I", "cytochrome oxidase 1"),
#'                  retmax = 100)
#'
#' # Search for the just the gene sequence for  P. infestans Cox I
#' get_isolate_seqs(species = c("Phytophthora infestans"),
#'                  genes = c("cox I", "cox 1", "cox1", "coxI", "cytochrome oxidase I", "cytochrome oxidase 1"),
#'                  retmax = 100,
#'                  extract_features = TRUE)
#'
#' # Search for all the gene sequences in whole sequences that contain P. infestans Cox I
#' get_isolate_seqs(species = c("Phytophthora infestans"),
#'                  genes = c("cox I", "cox 1", "cox1", "coxI", "cytochrome oxidase I", "cytochrome oxidase 1"),
#'                  retmax = 100,
#'                  extract_features = TRUE,
#'                  gene_name_in_feature = FALSE)
#'
#' # Search for whole sequences for P. infestans Cox I for just some isolates
#' get_isolate_seqs(species = c("Phytophthora infestans", "Phytophthora infestans", "Phytophthora infestans"),
#'                  isolates = c("44", "580", "180"),
#'                  genes = c("cox I", "cox 1", "cox1", "coxI", "cytochrome oxidase I", "cytochrome oxidase 1"))
#'
#' # Search for just the gene sequences for P. infestans Cox I for just some isolates
#' get_isolate_seqs(species = c("Phytophthora infestans", "Phytophthora infestans", "Phytophthora infestans"),
#'                  isolates = c("44", "580", "180"),
#'                  genes = c("cox I", "cox 1", "cox1", "coxI", "cytochrome oxidase I", "cytochrome oxidase 1"),
#'                  extract_features = TRUE)
#'
#' # Search for P infestans ITS with flanking regions and subset for complete results
#' result <- get_isolate_seqs(species = c("Phytophthora"),
#'                            genes = c("internal transcribed spacer"),
#'                            retmax = 300,
#'                            extract_features = TRUE,
#'                            flanking = c(300, 300))
#' result[result$complete, ]
#' result[result$flank_complete, ]
#'
#' }
#'
#' @export
get_isolate_seqs <- function(species, genes, isolates = NULL, extract_features = FALSE,
                             gene_name_in_feature = TRUE, flanking = c(0, 0),
                             db = "nucleotide", pause = 0.5, ...) {

  if (! is.numeric(flanking) || length(flanking) != 2) {
    stop('The "flanking" option must be a numeric vector of length 2')
  }

  get_one <- function(name, isolate = NULL) {

    # Wait a bit so NCBI doesnt get unhappy
    Sys.sleep(pause)

    # Search for sequences
    if (is.null(isolate)) {
      query <- paste0('"', name, '"[Organism] AND (', paste0('"', genes, '"[All Fields]', collapse = " OR "), ')')
    } else {
      query <- paste0('"', name, '"[Organism] AND ("', isolate, '"[Isolate] OR "', isolate, '"[Strain]) AND (', paste0('"', genes, '"[All Fields]', collapse = " OR "), ')')
    }
    search <- rentrez::entrez_search(db, term = query, ...)
    if (length(search$ids) == 0) {
      return(NULL)
    }

    if (extract_features) {
      # Parse features
      features <- parse_ft(rentrez::entrez_fetch(db, id = search$ids, rettype = "ft", retmode = "text"))
      if (gene_name_in_feature) {
        gene_in_feature <- purrr:::map_lgl(features$name, function(text) {
          purrr:::reduce(lapply(genes, function(gene) grepl(tolower(text), pattern = tolower(gene), fixed = TRUE)), `|`)
        })
        features <- features[gene_in_feature, ]
      }
      if (nrow(features) == 0) {
        return(NULL)
      }

      # Parse sequences
      sequences <- parse_seqs(rentrez::entrez_fetch(db, id = search$ids, rettype = "fasta", retmode = "xml"))

      # Join feature and sequence data
      output <- dplyr::left_join(features, sequences, by = "acc")

      # Get positions to return
      output$flank_start <- purrr::map2_dbl(output$start, output$end, function(x, y) {
        min(c(x, y)) - flanking[1]
      })
      output$flank_end <- purrr::map2_dbl(output$start, output$end, function(x, y) {
        max(c(x, y)) + flanking[2]
      })
      output$flank_complete <- output$complete & output$flank_start >= 1 & output$flank_end <= nchar(output$seq)
      output$flank_start[output$flank_start < 1] <- 1
      output$flank_end[output$flank_end > nchar(output$seq)] <- nchar(output$seq)[output$flank_end > nchar(output$seq)]

      # Subset sequences to fetures
      output$gene_seq <- substr(output$seq, output$flank_start, output$flank_end)
      output$gene_length <- nchar(output$gene_seq)

    } else {
      output <- parse_seqs(rentrez::entrez_fetch(db, id = search$ids, rettype = "fasta", retmode = "xml"))
    }

    # Add query info
    if (is.null(isolate)) {
      output <- tibble::as_tibble(cbind(list(species = name, query = query), output, stringsAsFactors = FALSE))
    } else {
      output <- tibble::as_tibble(cbind(list(species = name, isolate = isolate, query = query), output, stringsAsFactors = FALSE))
    }

    return(output)
  }

  if (is.null(isolates)) {
    return(dplyr::bind_rows(purrr::pmap(list(species), get_one)))
  } else {
    return(dplyr::bind_rows(purrr::pmap(list(species, isolates), get_one)))
  }

}
