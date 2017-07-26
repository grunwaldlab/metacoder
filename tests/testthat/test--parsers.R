library(metacoder)
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
  
})


test_that("Mothur classify.seqs *.tax.summary parsing", {
  raw_data <-
"taxlevel	 rankID	 taxon	 daughterlevels	 total	
0	0	Root	2	242	
1	0.1	Bacteria	50	242	
2	0.1.2	Actinobacteria	38	13	
3	0.1.2.3	Actinomycetaceae-Bifidobacteriaceae	10	13	
4	0.1.2.3.7	Bifidobacteriaceae	6	13	
5	0.1.2.3.7.2	Bifidobacterium_choerinum_et_rel.	8	13	
6	0.1.2.3.7.2.1	Bifidobacterium_angulatum_et_rel.	1	11	
7	0.1.2.3.7.2.1.1	unclassified	1	11	
8	0.1.2.3.7.2.1.1.1	unclassified	1	11	
9	0.1.2.3.7.2.1.1.1.1	unclassified	1	11	
10	0.1.2.3.7.2.1.1.1.1.1	unclassified	1	11
11	0.1.2.3.7.2.1.1.1.1.1.1	unclassified	1	11	
12	0.1.2.3.7.2.1.1.1.1.1.1.1	unclassified	1	11	
6	0.1.2.3.7.2.5	Bifidobacterium_longum_et_rel.	1	2		
"  
  result <- parse_mothur_taxonomy(text = raw_data)
})
