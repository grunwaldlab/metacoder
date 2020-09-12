library(metacoder)
library(testthat)
context("Simulated PCR")

# Make test data set
primer_1_site <- "AAGTACCTTAACGGAATTATAG"
primer_2_site <- "ATTCGTTTCGTAGGTGGAGC"
amplicon <- "NNNAGTGGATAGATAGGGGTTCTGTGGCGTTTGGGAATTAAAGATTAGAGANNN"
seq_1 <- paste0("AA", primer_1_site, amplicon, primer_2_site, "AAAA")
seq_2 <- rev_comp(seq_1)

f_primer <- "ACGTACCTTAACGGAATTATAG" # Note the "C" mismatch at position 2
r_primer <- rev_comp(primer_2_site)
seqs <- c(a = seq_1, b = seq_2)



test_that("primersearch works", {
  skip_on_cran()
  skip_if_not(metacoder:::primersearch_is_installed(must_be_installed = FALSE))
  
  result <- primersearch_raw(seqs, 
                             forward = c("p1" = f_primer),
                             reverse = c("p2" = r_primer))
  
  # Primer start indexes
  expect_equal(result$f_start[1], nchar("AA") + 1)
  expect_equivalent(result$r_start[1], nchar(paste0("AA", primer_1_site, amplicon)) + 1)
  expect_equal(result$f_start[2], nchar("AAAA") + 1)
  expect_equivalent(result$r_start[2], nchar(paste0("AAAA", primer_2_site, amplicon)) + 1)
  
  # Primer end indexes
  expect_true(all(result$f_end == result$f_start + nchar(result$f_primer) - 1))
  expect_true(all(result$r_end == result$r_start + nchar(result$r_primer) - 1))
  
  # PCR product that would be produced
  expect_equal(result$product[1], paste0(f_primer, amplicon, primer_2_site))
  expect_equal(result$product[2], rev_comp(paste0(f_primer, amplicon, primer_2_site)))
  
  # Sequence between primers
  expect_equal(result$amplicon[1], amplicon)
  expect_equal(result$amplicon[2], rev_comp(amplicon))
  
  # Primer binding sites
  expect_equivalent(result$f_match[1], primer_1_site)
  expect_equivalent(result$f_match[2], primer_2_site)
  expect_equivalent(result$r_match[1], primer_2_site)
  expect_equivalent(result$r_match[2], primer_1_site)
  
  # Mismatches
  expect_equal(result$f_mismatch[1], 1)
  expect_equal(result$f_mismatch[2], 0)
  expect_equal(result$r_mismatch[1], 0)
  expect_equal(result$r_mismatch[2], 1)
  
  # Can read from a file
  in_path <- system.file("extdata/silva_subset.fa", package = "metacoder")
  file_result <- primersearch_raw(file = in_path,
                                  forward = c("U519F" = "CAGYMGCCRCGGKAAHACC"),
                                  reverse = c("Arch806R" = "GGACTACNSGGGTMTCTAAT"),
                                  mismatch = 10)
  
  # Can read ape format (Current example file has U instead of T and this breaks ape)
  # ape_result <- primersearch_raw(ape::read.FASTA(in_path),
  #                                forward = c("U519F" = "CAGYMGCCRCGGKAAHACC"),
  #                                reverse = c("Arch806R" = "GGACTACNSGGGTMTCTAAT"),
  #                                mismatch = 10)

  # Can read seqinr format
  seqinr_result <- primersearch_raw(seqinr::read.fasta(in_path, forceDNAtolower = FALSE),
                                    forward = c("U519F" = "CAGYMGCCRCGGKAAHACC"),
                                    reverse = c("Arch806R" = "GGACTACNSGGGTMTCTAAT"),
                                    mismatch = 10)
  
  # Check that all input types return the same object
  expect_equivalent(lapply(file_result, unname), lapply(seqinr_result, unname))
  # expect_equal(ape_result, seqinr_result) # (Current example file has U instead of T and this breaks ape)
  
  })
