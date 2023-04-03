#' Print a item
#'
#' Used to print each item in the `taxmap` print method.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data The item to be printed
#' @param max_rows (`numeric` of length 1) The maximum number of rows in
#'   tables to print.
#' @param max_items (`numeric` of length 1) The maximum number of list
#'   items to print.
#' @param max_width (`numeric` of length 1) The maximum number of
#'   characters to print.
#' @param prefix (`numeric` of length 1) What to print in front of each
#'   line.
#'
#' @keywords internal
print_item <- function(obj, data, name = NULL, max_rows = 3, max_items = 3,
                       max_width = getOption("width") - 10, prefix = "") {

  # Find best print method
  print_methdods <- c(paste0("print__", class(data)),
                      "print__default_")
  applicable_methods <- print_methdods[vapply(print_methdods, function(x) exists(x), logical(1))]
  best_method <- applicable_methods[1]

  # Call print method
  get(best_method)(obj, data, name = name, prefix = prefix, max_width = max_width,
                   max_rows = max_rows)

  invisible(data)
}


#' Print a object with a prefix
#'
#' Print a object with a prefix. Uses the standard print method of the object.
#'
#' @param x What to print.
#'
#' @keywords internal
prefixed_print <- function(x, prefix, ...) {
  output <- paste0(prefix, utils::capture.output(print(x, ...)))
  cat(paste0(paste0(output, collapse = "\n"), "\n"))
}

#' Highlight taxon ID column
#'
#' Changes the font of a taxon ID column in a table print out.
#'
#' @param table_text The print out of the table in a character vector, one element per line.
#' @param header_index The row index that contains the table column names
#' @param row_indexes The indexes of the rows to be formatted.
highlight_taxon_ids <- function(table_text, header_index, row_indexes) {
  tax_id_bounds <- stringr::str_locate(table_text[header_index], "taxon_id[[:space:]]*")[1,]
  tax_id_part <- substr(table_text[row_indexes],
                        start = tax_id_bounds[1], stop = tax_id_bounds[2])
  tax_id_part <- sub(tax_id_part,
                     pattern = "^(\\S+)",
                     replacement = tid_font("\\1"),
                     perl = TRUE)
  table_text[row_indexes] <- paste0(substr(table_text[row_indexes],
                                           start = 1,
                                           stop = tax_id_bounds[1] - 1),
                                    tax_id_part,
                                    substr(table_text[row_indexes],
                                           start = tax_id_bounds[2] + 1,
                                           stop = nchar(table_text[row_indexes])))
  return(table_text)
}


#' Print a tibble
#'
#' Print a table for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__tbl_df <- function(obj, data, name, prefix, max_width, max_rows) {

  # used for tibble print methods
  loadNamespace("dplyr")

  # Make data set name
  if (length(name) > 0 && ! is.na(name)) {
    cat(paste0(prefix, name_font(name), ":\n"))
  }

  # Increment prefix
  prefix <- paste0(prefix, "  ")

  # Capture tibble print output
  output <- paste0(prefix, utils::capture.output(print(data, n = max_rows,
                                                       width = max_width - nchar(prefix),
                                                       max_extra_cols = 10)))
  # Remove any existing fonts
  output <- crayon::strip_style(output)

  # Modify font of tibble dims
  output[1] <- desc_font(output[1])

  # Highlight taxon IDs if they exist
  if (! is.null(obj$get_data_taxon_ids(name))) {
    if (nrow(data) > 0) {
      output <- highlight_taxon_ids(output,
                                    header_index = 2,
                                    row_indexes = 4:(3 + min(c(max_rows, nrow(data)))))
    }
    output[2] <- sub(output[2], pattern = "(^|\\W)taxon_id($|\\W)", replacement = tid_font("\\1taxon_id\\2"))
  } else {
    output[2] <- sub(output[2], pattern = "(^|\\W)taxon_id($|\\W)", replacement = error_font("\\1taxon_id\\2"))
  }

  # Change font of column classes
  if (nrow(data) > 0) {
    output[3] <- desc_font(output[3])
  }

  # Print output
  cat(paste0(paste0(output, collapse = "\n"), "\n"))
}


#' Print a data.frame
#'
#' Print a data.frame for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__data.frame <- function(obj, data, name, prefix, max_width, max_rows) {
  # Make data set name
  cat(paste0(prefix, name_font(name), ":\n"))

  # Increment prefix
  prefix <- paste0(prefix, "  ")

  # Make description of data.frame
  if (nrow(data) > max_rows) {
    cat(desc_font(paste0(prefix, "# A data.frame: ", nrow(data), " x ", ncol(data),
                         " (first ", max_rows, " rows shown)\n")))
    data <- data[1:max_rows, , drop = FALSE]
  } else {
    cat(desc_font(paste0(prefix, "# A data.frame: ", nrow(data), " x ", ncol(data), "\n")))
  }

  # Highlight taxon IDs if they exist
  output <- paste0(prefix, utils::capture.output(print(data)))
  if (! is.null(obj$get_data_taxon_ids(name))) {
    output <- highlight_taxon_ids(output,
                                  header_index = 1,
                                  row_indexes = 2:(1 + min(c(max_rows, nrow(data)))))
    output[1] <- sub(output[1], pattern = "(^|\\W)taxon_id($|\\W)", replacement = tid_font("\\1taxon_id\\2"))
  } else {
    output[1] <- sub(output[1], pattern = "(^|\\W)taxon_id($|\\W)", replacement = error_font("\\1taxon_id\\2"))

  }

  # Print output
  cat(paste0(paste0(output, collapse = "\n"), "\n"))
}


#' Print a list
#'
#' Print a list for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__list <- function(obj, data, name, prefix, max_width, max_rows) {
  if (length(data) < 1) {
    cat(paste0(prefix, name_font(name), ": empty list\n"))
  } else {
    cat(paste0(prefix, name_font(name), ": a list of ", length(data),
               ifelse(length(data) == 1, " item", " items")))
    if (is.null(names(data))) {
      cat("\n")
    } else {
      if (is.null(obj$get_data_taxon_ids(name))) {
        cat(paste0(" with names:\n  ",
                   limited_print(names(data),
                                 prefix = prefix,
                                 sep = punc_font(", "),
                                 mid = punc_font(" ... "),
                                 trunc_char = punc_font("[truncated]"),
                                 type = "silent")))

      } else {
        cat(paste0(" named by taxa:\n  ",
                   limited_print(tid_font(names(data)),
                                 prefix = prefix,
                                 sep = punc_font(", "),
                                 mid = punc_font(" ... "),
                                 trunc_char = punc_font("[truncated]"),
                                 type = "silent")))
      }
    }
  }
}


#' Generic vector printer
#'
#' Print a vector for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#' @param type The name of the type of vector to print (e.g. numeric).
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__vector <- function(obj, data, name, prefix, max_width, max_rows, type = class(data)[1]) {
  cat(paste0(prefix, name_font(name), ": ", ifelse(is.null(names(data)), "a vector of '", "a named vector of '"), type, "' with ", length(data),
             " item", ifelse(length(data) == 1, "", "s"), "\n"))
  if (length(data) > 0) {
    if (is.null(names(data))) { # Not named
      limited_print(data, prefix = paste0(prefix, "  "), max_chars = max_width, sep = punc_font(", "),
                    mid = punc_font(" ... "), trunc_char = punc_font("[truncated]"),
                    type = "cat")
    } else { # Is named
      if (is.null(obj$get_data_taxon_ids(name))) { # no taxon ids
        limited_print(paste0(names(data), punc_font(". "), data),
                      prefix = paste0(prefix, "  "), max_chars = max_width, sep = punc_font(", "),
                      mid = punc_font(" ... "), trunc_char = punc_font("[truncated]"),
                      type = "cat")
      }
      else { # has taxon ids
        limited_print(paste0(tid_font(names(data)), punc_font(". "), data),
                      prefix = paste0(prefix, "  "), max_chars = max_width, sep = punc_font(", "),
                      mid = punc_font(" ... "), trunc_char = punc_font("[truncated]"),
                      type = "cat")
      }
    }
  }
}


#' Print an integer
#'
#' Print an integer for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__integer <- function(obj, data, name, prefix, max_width, max_rows) {
  print__vector(obj, data, name, prefix, max_width, max_rows)
}



#' Print a numeric
#'
#' Print a numeric vector for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__numeric <- function(obj, data, name, prefix, max_width, max_rows) {
  print__vector(obj, data, name, prefix, max_width, max_rows)
}


#' Print a character
#'
#' Print a character for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__character <- function(obj, data, name, prefix, max_width, max_rows) {
  print__vector(obj, data, name, prefix, max_width, max_rows)
}


#' Print a logical
#'
#' Print a logical for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__logical <- function(obj, data, name, prefix, max_width, max_rows) {
  print__vector(obj, data, name, prefix, max_width, max_rows)
}


#' Print a factor
#'
#' Print a factor for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__factor <- function(obj, data, name, prefix, max_width, max_rows) {
  print__vector(obj, data, name, prefix, max_width, max_rows)
}


#' Print a ordered factor
#'
#' Print a ordered factor for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__ordered <- function(obj, data, name, prefix, max_width, max_rows) {
  print__vector(obj, data, name, prefix, max_width, max_rows, type = "ordered factor")
}


#' Print a matrix
#'
#' Print a matrix for the print method of taxmap objects.
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__matrix <- function(obj, data, name, prefix, max_width, max_rows) {
  cat(paste0(prefix, name_font(name), ":\n"))
  if (nrow(data) > max_rows) {
    cat(paste0(prefix, "  A ", nrow(data), " by ", ncol(data), " matrix (first ",
               max_rows, " rows shown)\n"))
    data <- data[1:max_rows, , drop = FALSE]
  } else {
    cat(paste0(prefix, "  A ", nrow(data), " by ", ncol(data), " matrix\n"))
  }
  prefixed_print(data, prefix = paste0(prefix, "  "))
}


#' Print method for unsupported
#'
#' Print method for unsupported classes for taxmap objects
#'
#' Which print method is called is determined by its name, so changing the name
#' of this function will change when it is called.
#'
#' @param obj The taxmap object containing the thing to print
#' @param data Something to print
#' @param name The name of the thing to print
#' @param prefix What to put before the thing printed. Typically a space.
#' @param max_width Maximum width in number of characters to print
#' @param max_rows Maximum number of rows to print
#'
#' @family taxmap print methods
#'
#' @keywords internal
print__default_ <- function(obj, data, name, prefix, max_width, max_rows) {
  cat(paste0(prefix, name_font(name), ":\n"))
  prefixed_print(data, prefix = paste0(prefix, "  "))
}


#' Taxon id formatting in print methods
#'
#' A simple wrapper to make changing the formatting of text printed easier.
#'
#' @param text What to print
#'
#' @family printer fonts
#'
#' @keywords internal
tid_font <- function(text) {
  crayon::green(text)
}


#' Punctuation formatting in print methods
#'
#' A simple wrapper to make changing the formatting of text printed easier.
#' This is used for non-data, formatting characters
#'
#' @param text What to print
#'
#' @family printer fonts
#'
#' @keywords internal
punc_font <- function(text) {
  crayon::silver(text)
}


#' Description formatting in print methods
#'
#' A simple wrapper to make changing the formatting of text printed easier.
#' This is used for non-data, formatting characters
#'
#' @param text What to print
#'
#' @family printer fonts
#'
#' @keywords internal
desc_font <- function(text) {
  crayon::italic(text)
}


#' Variable name formatting in print methods
#'
#' A simple wrapper to make changing the formatting of text printed easier.
#' This is used for non-data, formatting characters
#'
#' @param text What to print
#'
#' @family printer fonts
#'
#' @keywords internal
name_font <- function(text) {
  crayon::bold(text)
}


#' Font to indicate an error
#'
#' A simple wrapper to make changing the formatting of text printed easier.
#' This is used for non-data, formatting characters
#'
#' @param text What to print
#'
#' @family printer fonts
#'
#' @keywords internal
error_font <- function(text) {
  crayon::bgRed(text)
}
