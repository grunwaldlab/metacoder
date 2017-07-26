library(metacoder)
context("Tree plotting")

test_that("Mothur .taxonomy parsing works", {
  raw_data <-
"AY457915\tBacteria(100);Firmicutes(99);Clostridiales(99);Johnsonella_et_rel.(99);Johnsonella_et_rel.(99);Johnsonella_et_rel.(91);Eubacterium_eligens_et_rel.(89);Lachnospira_pectinoschiza(80);
AY457914\tBacteria(100);Firmicutes(100);Clostridiales(100);Johnsonella_et_rel.(100);Johnsonella_et_rel.(100);Johnsonella_et_rel.(95);Eubacterium_eligens_et_rel.(92);Eubacterium_eligens(84);Eubacterium_eligens(81);
AY457913\tBacteria(100);Firmicutes(100);Clostridiales(100);Johnsonella_et_rel.(100);Johnsonella_et_rel.(100);Roseoburia_et_rel.(97);Roseoburia_et_rel.(97);Eubacterium_ramulus_et_rel.(90);uncultured(90);
AY457912\tBacteria(100);Firmicutes(99);Clostridiales(99);Johnsonella_et_rel.(99);Johnsonella_et_rel.(99);
AY457911\tBacteria(100);Firmicutes(99);Clostridiales(98);Ruminococcus_et_rel.(96);Anaerofilum-Faecalibacterium(92);Faecalibacterium(92);Faecalibacterium_prausnitzii(90);
"
  
  result <- parse_mothur_taxonomy(text = raw_data)
  
})