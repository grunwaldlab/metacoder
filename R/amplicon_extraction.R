#' Extract expected amplicons
#'
#' Extract expected amplicons from sequences using primer sequences.
#'
#' @param seqs The DNA sequences of the target to be amplified. Letter case is ignored but
#'   preserved.
#' @param forward The DNA sequence of the forward primer in the 5' to 3' orientation. IUPAC codes
#'   can be present.
#' @param reverse The DNA sequence of the forward primer in the 5' to 3' orientation. IUPAC codes
#'   can be present.
#' @param max_mismatch The maximum number of mismatching bases allowed. See
#'   [Biostrings::matchPattern] for more information.
#' @param min_mismatch The minimum number of mismatching bases allowed. See
#'   [Biostrings::matchPattern] for more information.
#' @param allow_indels If TRUE then indels are allowed. See [Biostrings::matchPattern] for more
#'   information.
#' @param algorithm One of: "auto", "naive-exact", "naive-inexact", "boyer-moore" or "shift-or". See
#'   [Biostrings::matchPattern] for more information.
#' @param infer_amps If TRUE, performe an additional step to find amplicons in sequences that might
#'   not contain the primer binding sites due to being produced with similar primers. This is done
#'   by aligning amplicons found using primers to sequences not found to contain the primers. For
#'   each sequence, the region aligned to the best-matching amplicon will be included as an amplion
#'   if:
#'   *  The aligned region extends to the boundry of the sequence, suggesting the missing primer site is not in the sequence. 
#'   *  If the aligned region covers the minimum proportion of the best matching amplicon defined by `infer_thresh_cov`.
#'   *  If the aligned region has the minimum percent identity to the best matching amplicon defined by `infer_thresh_pid`.
#' @param infer_thresh_cov See the `infer_amps` option description.
#' @param infer_thresh_pid See the `infer_amps` option description.
#' 
#' @export
calc_amplicons <- function(seqs, forward, reverse, min_coverage = 0.9, ...) {
  
  get_amplicon_chr <- function(seq) {
    map2_chr(seq@ranges@start, seq@ranges@width, function(s, w) {
      as.character(seq@subject[seq(from = s, length.out = w)])
    })
  }
  
  
  # Get simulated amplicons using primers
  full_amps <- future_map(seqs, function(s) {
    get_amplicon_chr(matchProbePair(DNAString(s), Fprobe = forward, Rprobe = reverse))
  })
  
  # Check for multiple possible amplicons per input
  if (any(map_int(full_amps, length) > 1)) {
    stop('Some inputs have more than one amplicons.')
  }
  full_amps <- unlist(full_amps[map_int(full_amps, length) == 1])
  
  # Remove primers from amplicons
  full_amps <- substr(full_amps, start = nchar(forward) + 1, stop = nchar(full_amps) - nchar(reverse))
  
  # Align unamplified sequences with best matching amplicons
  unamped <- seqs[! names(seqs) %in% names(full_amps)]
  unamped_aligned <- future_map_chr(unamped, function(s) {
    aligned <- pairwiseAlignment(pattern = full_amps, subject = s, 
                                 type = 'global-local',
                                 # type = 'overlap',
                                 gapOpening = 10, gapExtension = 4)
    best_score_i <- which.max(aligned@score)
    best_align <- aligned[best_score_i]
    aligned_subject_char <- gsub(as.character(best_align@subject), pattern = '-', replacement = '')
    aligned_amp_char <- gsub(as.character(best_align@pattern), pattern = '-', replacement = '')
    if (best_align@pattern@range@width < best_align@pattern@unaligned@ranges@width * min_coverage || nchar(aligned_subject_char) == 0) {
      aligned_subject_char <- NA
    }
    aligned_subject_char
  })
  
  # Add full amplicons to inferred amplicons from alignment
  all_amps <- c(full_amps, unamped_aligned[! is.na(unamped_aligned)])
  all_amps <- all_amps[is.na(all_amps) | nchar(all_amps) > 0]
  setNames(all_amps[names(seqs)], names(seqs))
}