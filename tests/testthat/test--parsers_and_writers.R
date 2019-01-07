library(metacoder)
library(testthat)
context("Input parsing")

test_that("Mothur classify.seqs *.taxonomy parsing", {
  raw_data <-
"AY457915	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;Johnsonella_et_rel.;Johnsonella_et_rel.;Eubacterium_eligens_et_rel.;Lachnospira_pectinoschiza;
AY457914	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;Johnsonella_et_rel.;Johnsonella_et_rel.;Eubacterium_eligens_et_rel.;Eubacterium_eligens;Eubacterium_eligens;
AY457913	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;Johnsonella_et_rel.;Roseoburia_et_rel.;Roseoburia_et_rel.;Eubacterium_ramulus_et_rel.;uncultured;
AY457912	Bacteria;Firmicutes;Clostridiales;Johnsonella_et_rel.;Johnsonella_et_rel.;
AY457911	Bacteria;Firmicutes;Clostridiales;Ruminococcus_et_rel.;Anaerofilum-Faecalibacterium;Faecalibacterium;Faecalibacterium_prausnitzii;
"

  result <- parse_mothur_taxonomy(text = raw_data)
  expect_equal(length(result$taxa), 18)
  expect_equal(length(roots(result)), 1)
  expect_true(all(c("Bacteria", "Firmicutes") %in% result$taxon_names()))
  
  # Check that the input can be replicated
  out_path <- "test_mothur_tax_output.txt"
  write_mothur_taxonomy(result, file = out_path)
  expect_equal(readLines(out_path), strsplit(raw_data, split = "\n")[[1]])
  expect_error(write_mothur_taxonomy(result))
  
  # Delete files used for tests
  file.remove(out_path)
})


test_that("Mothur classify.seqs *.taxonomy parsing w/ scores", {
  raw_data <-
    "AY457915\tBacteria(100);Firmicutes(99);Clostridiales(99);Johnsonella_et_rel.(99);Johnsonella_et_rel.(99);Johnsonella_et_rel.(91);Eubacterium_eligens_et_rel.(89);Lachnospira_pectinoschiza(80);
AY457914\tBacteria(100);Firmicutes(100);Clostridiales(100);Johnsonella_et_rel.(100);Johnsonella_et_rel.(100);Johnsonella_et_rel.(95);Eubacterium_eligens_et_rel.(92);Eubacterium_eligens(84);Eubacterium_eligens(81);
AY457913\tBacteria(100);Firmicutes(100);Clostridiales(100);Johnsonella_et_rel.(100);Johnsonella_et_rel.(100);Roseoburia_et_rel.(97);Roseoburia_et_rel.(97);Eubacterium_ramulus_et_rel.(90);uncultured(90);
AY457912\tBacteria(100);Firmicutes(99);Clostridiales(99);Johnsonella_et_rel.(99);Johnsonella_et_rel.(99);
AY457911\tBacteria(100);Firmicutes(99);Clostridiales(98);Ruminococcus_et_rel.(96);Anaerofilum-Faecalibacterium(92);Faecalibacterium(92);Faecalibacterium_prausnitzii(90);
"
  
  result <- parse_mothur_taxonomy(text = raw_data)
  expect_equal(length(result$taxa), 18)
  expect_equal(length(roots(result)), 1)
  expect_true(all(c("Bacteria", "Firmicutes") %in% result$taxon_names()))
  expect_equal(nrow(result$data$class_data), stringr::str_count(raw_data, ";"))
  expect_true("score" %in% colnames(result$data$class_data))
  
  
  # Check that the input can be replicated
  out_path <- "test_mothur_tax_output.txt"
  write_mothur_taxonomy(result, file = out_path)
  expect_equal(readLines(out_path), strsplit(raw_data, split = "\n")[[1]])
  expect_error(write_mothur_taxonomy(result))
  
  # Delete files used for tests
  file.remove(out_path)
})


test_that("Mothur classify.seqs *.tax.summary  detailed parsing", {
  raw_data <-
"taxlevel	 rankID	 taxon	 daughterlevels	 total	A	B	C	
0	0	Root	2	242	84	84	74	
1	0.1	Bacteria	50	242	84	84	74	
2	0.1.2	Actinobacteria	38	13	0	13	0	
3	0.1.2.3	Actinomycetaceae-Bifidobacteriaceae	10	13	0	13	0	
4	0.1.2.3.7	Bifidobacteriaceae	6	13	0	13	0	
5	0.1.2.3.7.2	Bifidobacterium_choerinum_et_rel.	8	13	0	13	0	
6	0.1.2.3.7.2.1	Bifidobacterium_angulatum_et_rel.	1	11	0	11	0	
7	0.1.2.3.7.2.1.1	unclassified	1	11	0	11	0	
8	0.1.2.3.7.2.1.1.1	unclassified	1	11	0	11	0	
9	0.1.2.3.7.2.1.1.1.1	unclassified	1	11	0	11	0	
10	0.1.2.3.7.2.1.1.1.1.1	unclassified	1	11	0	11	0	
11	0.1.2.3.7.2.1.1.1.1.1.1	unclassified	1	11	0	11	0	
12	0.1.2.3.7.2.1.1.1.1.1.1.1	unclassified	1	11	0	11	0	
6	0.1.2.3.7.2.5	Bifidobacterium_longum_et_rel.	1	2	0	2	0	
7	0.1.2.3.7.2.5.1	unclassified	1	2	0	2	0	
8	0.1.2.3.7.2.5.1.1	unclassified	1	2	0	2	0	
9	0.1.2.3.7.2.5.1.1.1	unclassified	1	2	0	2	0"  
  result <- parse_mothur_tax_summary(text = raw_data)
  result_from_file <- parse_mothur_tax_summary(file = "example_data/mothur_summary.txt")
  
  expect_equal(result, result_from_file)
  expect_equal(length(result$taxa), 17)
  expect_equal(length(roots(result)), 1)
  expect_true(all(c("Bacteria", "Actinobacteria") %in% result$taxon_names()))
})


test_that("Mothur classify.seqs *.tax.summary simple parsing", {
  raw_data <- 
'taxon	total	A	B	C
"k__Bacteria";"p__Actinobacteria";"c__Actinobacteria";"o__Bifidobacteriales";"f__Bifidobacteriaceae";"g__Bifidobacterium";"s__";	1	0	1	0
"k__Bacteria";"p__Actinobacteria";"c__Actinobacteria";"o__Bifidobacteriales";"f__Bifidobacteriaceae";"g__Bifidobacterium";"s__adolescentis";	1	0	1	0
"k__Bacteria";"p__Actinobacteria";"c__Actinobacteria";"o__Bifidobacteriales";"f__Bifidobacteriaceae";"g__Bifidobacterium";"s__longum";	1	0	1	0
'
  result <- parse_mothur_tax_summary(text = raw_data)
  expect_equal(length(result$taxa), 9)
  expect_equal(length(roots(result)), 1)
  expect_true(all(c("k__Bacteria", "p__Actinobacteria") %in% result$taxon_names()))
})



test_that("Newick parsing", {
  expect_warning(result <- parse_newick("example_data/newick_example_1.txt"))
  expect_equal(length(result$taxa), 21)
  expect_equal(length(roots(result)), 2)
  expect_true(all(c("node_1", "node_2") %in% result$taxon_names()))
})


test_that("Parsing the UNITE general release fasta", {
  # Reading
  seq_in_path <- "example_data/unite_general.fa"
  result <- parse_unite_general(file = seq_in_path)
  expect_equal(length(result$taxa), 183)
  expect_equal(length(roots(result)), 1)
  expect_equivalent(result$taxon_names()[result$data$tax_data$taxon_id[5]], "Orbilia_sp")
  expect_equal(result$data$tax_data$organism[5], "Orbilia_sp")
  expect_equal(result$data$tax_data$unite_seq[5], "CCAAATCATGTCTCCCGGCCGCAAGGCAGGTGCAGGCGTTTAACCCTTTGTGAACCAAAAAACCTTTCGCTTCGGCAGCAGCTCGGTTGGAGACAGCCTCTGTGTCAGCCTGCCGCTAGCACCAATTATCAAAACTTGCGGTTAGCAACATTGTCTGATTACCAAATTTTCGAATGAAAATCAAAACTTTCAACAACGGATCTCTTGGTTCCCGCATCGATGAAGAACGCAGCGAAACGCGATAGTTAATGTGAATTGCAGAATTCAGTGAATCATCGAGTCTTTGAACGCACATTGCGCCCATTGGTATTCCATTGGGCATGTCTGTTTGAGCGTCATTACAACCCTCGGTCACCACCGGTTTTGAGCGAGCAGGGTCTTCGGATCCAGCTGGCTTTAAAGTTGTAAGCTCTGCTGGCTGCTCGGCCCAACCAGAACATAGTAAAATCATGCTTGTTCAAGGTTCGCGGTCGAAGCGGTACGGCCTGAACAATACCTACCACCTCTTAGG")
  
  # Check that the input can be replicated
  seq_out_path <- "test_unite_output.fa"
  write_unite_general(result, file = seq_out_path)
  expect_equal(readLines(seq_out_path), readLines(seq_in_path))
  expect_error(write_unite_general(result))
  
  # Delete files used for tests
  file.remove(seq_out_path)
})


test_that("Parsing the RDP fasta release", {
  # Reading
  seq_in_path <- "example_data/rdp_example.fa"
  result <- parse_rdp(file = seq_in_path)
  expect_equal(length(result$taxa), 26)
  expect_equal(length(roots(result)), 1)
  expect_equivalent(result$taxon_names()[result$data$tax_data$taxon_id[3]], "Saccharomyces")
  expect_equal(result$data$tax_data$rdp_id[3], "S004468774")
  expect_true(startsWith(result$data$tax_data$rdp_seq[3], "gtttgacctcaaatcaggtaggagtacccgctgaacttaagcatatcaataagcggaggaaaagaaaccaaccgggattg"))
  
  # Check that the input can be replicated
  seq_out_path <- "test_rdp_output.fa"
  write_rdp(result, file = seq_out_path)
  expect_equal(readLines(seq_out_path), readLines(seq_in_path))
  expect_error(write_greengenes(result))
  
  # Delete files used for tests
  file.remove(seq_out_path)
  })


test_that("Parsing the SILVA fasta release", {
  # Reading
  seq_in_path <- "example_data/silva_example.fa"
  result <- parse_silva_fasta(file = seq_in_path)
  expect_equal(length(result$taxa), 164)
  expect_equal(length(roots(result)), 2)
  expect_equivalent(result$taxon_names()[result$data$tax_data$taxon_id[5]], "Physalis peruviana")
  expect_equal(result$data$tax_data$ncbi_id[5], "GEET01005309")
  expect_true(startsWith(result$data$tax_data$silva_seq[5], "GAUGGAUGCCUUGGCUUCAUCAGGCGAAGAAGGACGCAGCAAGCUGCGAUAAGCUUCGGGGAGCGGCACGCACGCUUUGA"))

   # Check that the input can be replicated
  seq_out_path <- "test_rdp_output.fa"
  write_silva_fasta(result, file = seq_out_path)
  # expect_equal(readLines(seq_out_path)[c(-89, -2580)],
  #              readLines(seq_in_path)[c(-89, -2580)])
  expect_error(write_greengenes(result))
  
  # Delete files used for tests
  file.remove(seq_out_path)
})


test_that("Parsing/writing the greengenes database", {
  # Reading
  tax_in_path <- "example_data/gg_tax_example.txt"
  seq_in_path <- "example_data/gg_seq_example.fa"
  result <- parse_greengenes(tax_file = tax_in_path, seq_file = seq_in_path)
  expect_equal(length(result$taxa), 119)
  expect_equal(length(roots(result)), 1)
  expect_equivalent(result$taxon_names()[result$data$tax_data$taxon_id[5]], "Rhodobacteraceae")
  expect_equal(result$data$tax_data$gg_id[5], "1111758")
  expect_true(startsWith(result$data$tax_data$gg_seq[5], "TTAGAGTTTGATCCTGGCTCAGAACGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGCGCCCTTCGGGGTGAGCGGCGGACGGGTGAGTAACGCGTGGGAACGTGCCCTCTTCTGCGGGATAGCC"))
  
  # Check that the input can be replicated
  tax_out_path <-  "test_gg_output.txt"
  seq_out_path <- "test_gg_output.fa"
  write_greengenes(result, tax_file = tax_out_path, seq_file = seq_out_path)
  expect_equal(readLines(tax_out_path), readLines(tax_in_path))
  expect_equal(readLines(seq_out_path), readLines(seq_in_path))
  expect_error(write_greengenes(result))
  
  # Delete files used for tests
  file.remove(tax_out_path)
  file.remove(seq_out_path)
})


test_that("Converting to phyloseq", {
  # test round-trip
  library(phyloseq)
  data(enterotype)
  x <- parse_phyloseq(enterotype)
  y <- as_phyloseq(x)
  expect_equivalent(enterotype, y)
})


test_that("Parsing/writing dada2 output", {
  # test round-trip
  load("example_data/dada2.RData")
  obj <- parse_dada2(seq_table = seqtab.nochim, tax_table = taxa)
  seq_table <- make_dada2_asv_table(obj)
  tax_table <- make_dada2_tax_table(obj)
  expect_equal(seqtab.nochim, seq_table)
  expect_equal(taxa, tax_table)
})