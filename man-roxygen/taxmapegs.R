#' @examples
#' # The code below shows how to contruct a taxmap object from scratch.
#' # Typically, taxmap objects would be the output of a parsing function,
#' #  not created from scratch, but this is for demostration purposes.
#'
#' notoryctidae <- taxon(
#' name = taxon_name("Notoryctidae"),
#' rank = taxon_rank("family"),
#' id = taxon_id(4479)
#' )
#' notoryctes <- taxon(
#'   name = taxon_name("Notoryctes"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(4544)
#' )
#' typhlops <- taxon(
#'   name = taxon_name("typhlops"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(93036)
#' )
#' mammalia <- taxon(
#'   name = taxon_name("Mammalia"),
#'   rank = taxon_rank("class"),
#'   id = taxon_id(9681)
#' )
#' felidae <- taxon(
#'   name = taxon_name("Felidae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(9681)
#' )
#' felis <- taxon(
#'   name = taxon_name("Felis"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(9682)
#' )
#' catus <- taxon(
#'   name = taxon_name("catus"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(9685)
#' )
#' panthera <- taxon(
#'   name = taxon_name("Panthera"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(146712)
#' )
#' tigris <- taxon(
#'   name = taxon_name("tigris"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(9696)
#' )
#' plantae <- taxon(
#'   name = taxon_name("Plantae"),
#'   rank = taxon_rank("kingdom"),
#'   id = taxon_id(33090)
#' )
#' solanaceae <- taxon(
#'   name = taxon_name("Solanaceae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(4070)
#' )
#' solanum <- taxon(
#'   name = taxon_name("Solanum"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(4107)
#' )
#' lycopersicum <- taxon(
#'   name = taxon_name("lycopersicum"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(49274)
#' )
#' tuberosum <- taxon(
#'   name = taxon_name("tuberosum"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(4113)
#' )
#' homo <- taxon(
#'   name = taxon_name("homo"),
#'   rank = taxon_rank("genus"),
#'   id = taxon_id(9605)
#' )
#' sapiens <- taxon(
#'   name = taxon_name("sapiens"),
#'   rank = taxon_rank("species"),
#'   id = taxon_id(9606)
#' )
#' hominidae <- taxon(
#'   name = taxon_name("Hominidae"),
#'   rank = taxon_rank("family"),
#'   id = taxon_id(9604)
#' )
#' unidentified <- taxon(
#'   name = taxon_name("unidentified")
#' )
#'
#' tiger <- hierarchy(mammalia, felidae, panthera, tigris)
#' cat <- hierarchy(mammalia, felidae, felis, catus)
#' human <- hierarchy(mammalia, hominidae, homo, sapiens)
#' mole <- hierarchy(mammalia, notoryctidae, notoryctes, typhlops)
#' tomato <- hierarchy(plantae, solanaceae, solanum, lycopersicum)
#' potato <- hierarchy(plantae, solanaceae, solanum, tuberosum)
#' potato_partial <- hierarchy(solanaceae, solanum, tuberosum)
#' unidentified_animal <- hierarchy(mammalia, unidentified)
#' unidentified_plant <- hierarchy(plantae, unidentified)
#'
#' info <- data.frame(stringsAsFactors = FALSE,
#'                    name = c("tiger", "cat", "mole", "human", "tomato", "potato"),
#'                    n_legs = c(4, 4, 4, 2, 0, 0),
#'                    dangerous = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
#'
#' abund <- data.frame(code = rep(c("T", "C", "M", "H"), 2),
#'                     sample_id = rep(c("A", "B"), each = 2),
#'                     count = c(1,2,5,2,6,2,4,0),
#'                     taxon_index = rep(1:4, 2))
#'
#' phylopic_ids <- c("e148eabb-f138-43c6-b1e4-5cda2180485a",
#'                   "12899ba0-9923-4feb-a7f9-758c3c7d5e13",
#'                   "11b783d5-af1c-4f4e-8ab5-a51470652b47",
#'                   "9fae30cd-fb59-4a81-a39c-e1826a35f612",
#'                   "b6400f39-345a-4711-ab4f-92fd4e22cb1a",
#'                   "63604565-0406-460b-8cb8-1abe954b3f3a")
#'
#' foods <- list(c("mammals", "birds"),
#'               c("cat food", "mice"),
#'               c("insects"),
#'               c("Most things, but especially anything rare or expensive"),
#'               c("light", "dirt"),
#'               c("light", "dirt"))
#'
#' reaction <- function(x) {
#'   ifelse(x$data$info$dangerous,
#'          paste0("Watch out! That ", x$data$info$name, " might attack!"),
#'          paste0("No worries; its just a ", x$data$info$name, "."))
#' }
#'
#' ex_taxmap <- taxmap(tiger, cat, mole, human, tomato, potato,
#'                     data = list(info = info,
#'                                 phylopic_ids = phylopic_ids,
#'                                 foods = foods,
#'                                 abund = abund),
#'                     funcs = list(reaction = reaction))
