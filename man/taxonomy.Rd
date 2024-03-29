% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/old_taxa--taxonomy--class.R
\name{taxonomy}
\alias{taxonomy}
\title{Taxonomy class}
\usage{
taxonomy(..., .list = NULL, named_by_rank = FALSE)
}
\arguments{
\item{...}{Any number of object of class [hierarchy()] or character
vectors.}

\item{.list}{An alternate to the `...` input. Any number of object of class
[hierarchy()] or character vectors in a list. Cannot be used with `...`.}

\item{named_by_rank}{(`TRUE`/`FALSE`) If  `TRUE` and the input is a list of
vectors with each vector named by ranks, include that rank info in the
output object, so it can be accessed by `out$taxon_ranks()`. If `TRUE`,
taxa with different ranks, but the same name and location in the taxonomy,
will be considered different taxa.}
}
\value{
An `R6Class` object of class `Taxonomy`
}
\description{
Stores a taxonomy composed of [taxon()] objects organized in a tree
structure. This differs from the [hierarchies()] class in how the [taxon()]
objects are stored. Unlike [hierarchies()], each taxon is only stored once
and the relationships between taxa are stored in an [edge
list](https://en.wikipedia.org/wiki/Adjacency_list).
}
\examples{

# Making a taxonomy object with vectors
taxonomy(c("mammalia", "felidae", "panthera", "tigris"),
         c("mammalia", "felidae", "panthera", "leo"),
         c("mammalia", "felidae", "felis", "catus"))

# Making a taxonomy object from scratch
#   Note: This information would usually come from a parsing function.
#         This is just for demonstration.
x <- taxon(
  name = taxon_name("Notoryctidae"),
  rank = taxon_rank("family"),
  id = taxon_id(4479)
)
y <- taxon(
  name = taxon_name("Notoryctes"),
  rank = taxon_rank("genus"),
  id = taxon_id(4544)
)
z <- taxon(
  name = taxon_name("Notoryctes typhlops"),
  rank = taxon_rank("species"),
  id = taxon_id(93036)
)

a <- taxon(
  name = taxon_name("Mammalia"),
  rank = taxon_rank("class"),
  id = taxon_id(9681)
)
b <- taxon(
  name = taxon_name("Felidae"),
  rank = taxon_rank("family"),
  id = taxon_id(9681)
)

cc <- taxon(
  name = taxon_name("Puma"),
  rank = taxon_rank("genus"),
  id = taxon_id(146712)
)
d <- taxon(
  name = taxon_name("Puma concolor"),
  rank = taxon_rank("species"),
  id = taxon_id(9696)
)

m <- taxon(
  name = taxon_name("Panthera"),
  rank = taxon_rank("genus"),
  id = taxon_id(146712)
)
n <- taxon(
  name = taxon_name("Panthera tigris"),
  rank = taxon_rank("species"),
  id = taxon_id(9696)
)

(hier1 <- hierarchy(z, y, x, a))
(hier2 <- hierarchy(cc, b, a, d))
(hier3 <- hierarchy(n, m, b, a))

(hrs <- hierarchies(hier1, hier2, hier3))

ex_taxonomy <- taxonomy(hier1, hier2, hier3)
}
\seealso{
Other classes: 
\code{\link{hierarchies}()},
\code{\link{hierarchy}()},
\code{\link{taxa}()},
\code{\link{taxmap}()},
\code{\link{taxon}()},
\code{\link{taxon_database}()},
\code{\link{taxon_id}()},
\code{\link{taxon_name}()},
\code{\link{taxon_rank}()}
}
\concept{classes}
