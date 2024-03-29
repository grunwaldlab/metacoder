% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_taxa--hierarchy.R
\name{hierarchy}
\alias{hierarchy}
\title{The Hierarchy class}
\usage{
hierarchy(..., .list = NULL)
}
\arguments{
\item{...}{Any number of object of class `Taxon` or taxonomic names as
character strings}

\item{.list}{An alternate to the `...` input. Any number of object of class
[taxon()] or character vectors in a list. Cannot be used with `...`.}
}
\value{
An `R6Class` object of class `Hierarchy`
}
\description{
A class containing an ordered list of [taxon()] objects that represent a
hierarchical classification.
}
\details{
On initialization, taxa are sorted if they have ranks with a known
  order.

**Methods**
  \describe{
    \item{`pop(rank_names)`}{
      Remove `Taxon` elements by rank name, taxon name or taxon ID. The
      change happens in place, so you don't need to assign output to a new
      object. returns self - rank_names (character) a vector of rank names
    }
    \item{`pick(rank_names)`}{
      Select `Taxon` elements by rank name, taxon name or taxon ID. The
      change happens in place, so you don't need to assign output to a new
      object. returns self - rank_names (character) a vector of rank names
    }
  }
}
\examples{
(x <- taxon(
  name = taxon_name("Poaceae"),
  rank = taxon_rank("family"),
  id = taxon_id(4479)
))

(y <- taxon(
  name = taxon_name("Poa"),
  rank = taxon_rank("genus"),
  id = taxon_id(4544)
))

(z <- taxon(
  name = taxon_name("Poa annua"),
  rank = taxon_rank("species"),
  id = taxon_id(93036)
))

(res <- hierarchy(z, y, x))

res$taxa
res$ranklist

# null taxa
x <- taxon(NULL)
(res <- hierarchy(x, x, x))
## similar to hierarchy(), but `taxa` slot is not empty
}
\seealso{
Other classes: 
\code{\link{hierarchies}()},
\code{\link{taxa}()},
\code{\link{taxmap}()},
\code{\link{taxon}()},
\code{\link{taxon_database}()},
\code{\link{taxon_id}()},
\code{\link{taxon_name}()},
\code{\link{taxon_rank}()},
\code{\link{taxonomy}()}
}
\concept{classes}
