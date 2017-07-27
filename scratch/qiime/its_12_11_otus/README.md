UNITE/QIIME ITS reference OTUs
------------------------------

These reference sequence sets represent clustered versions (at 99% and 97% sequence similarity) of all fungal rDNA ITS sequences in the current UNITE+INSD (International Nucleotide Sequence Databases: NCBI, EMBL, DDBJ) release of circa 300,000 sequences (http://unite.ut.ee/repository.php ; filename: UNITE_public_24.09.12.fasta). The taxonomy mapping files provided here were created from the Index Fungorum ranked classification schema provided by UNITE in association with each sequence in their database. These files have been designed to facilitate use with [QIIME](http://www.qiime.org).

Efforts have been made by UNITE to improve the taxonomic information associated with some of the sequences in their database.  The QIIME reference sequence sets linked here have not been subject to any other form of curation (manual or automated) and certainly include incorrectly identified sequences, chimeras, and other problematic sequences.

**These fungal ITS reference sets are alpha versions.** Responsibility lies with the user to verify the accuracy of associated taxonomic information as well as quality of sequences in these datasets.  Improved fungal rDNA ITS reference sets based on the semi-curated centroids for sequence clusters in the UNITE Global Key Annotations module will soon be available here.

Acknowledgements
----------------

The creation of these reference sets is the direct result of a Fungal ITS Workshop held on 19–20 October, 2012 in Boulder, Colorado, USA.  The workshop would not have been possible without the intellectual and financial support of Paula Olsiewski and the [Alfred P. Sloan Foundation](http://www.sloan.org/).  We thank Jason Stajich, for organizing the meeting, as well as all the workshop participants for their valuable contributions.  We specially acknowledge Scott Bates, Greg Caporaso, Noah Fierer, Rob Knight, Jonathan Leff, Daniel McDonald, Henrik Nilsson, Jai Ram Rideout, and Lee Taylor for their expertise, support, and/or advice in constructing the reference sets.

Questions about these data should be directed to the [Qiime Forum](http://forum.qiime.org).

Find information on [UNITE here](http://unite.ut.ee/index.php) and on [QIIME here](http://www.qiime.org).

Technical notes
---------------

Technical notes were compiled by Jon Leff, Jai Ram Rideout and Greg Caporaso.

All commands were run using [QIIME](https://github.com/qiime/qiime) 1.5.0-dev (git commit hash  ``1ed0a4f3cea5e6da04896ff9e5ac337d40f3358``) and the [nested_reference_otu](https://github.com/qiime/nested_reference_otus) workflow (git commit hash ``821a98df6773ea4e4d209af20b9a8cf34d00324e``).

All input files were downloaded November of 2012, and were identified as the ``UNITE_public_24.09.12`` release on the UNITE website. The data in this directory was compiled and published on 27 Nov 2012.

Generating the current reference OTUs
-------------------------------------

```
export WORKING_DIR=/Users/leffj/ITS_ref/UNITE/UNITE_public_240912/

sort_seqs.py -i $WORKING_DIR/UNITE_public_24.09.12.fasta -t $WORKING_DIR/UNITE_public_taxonomy_mapping_24.09.12.txt -o $WORKING_DIR/UNITE_public_24.09.12_sorted.fasta

nested_reference_workflow.py -i $WORKING_DIR/UNITE_public_24.09.12_sorted.fasta -o /Users/leffj/ITS_ref/UNITE/NRWout/ -r 20121105 -s 97

nested_reference_workflow.py -i $WORKING_DIR/UNITE_public_24.09.12_sorted.fasta -o /Users/leffj/ITS_ref/UNITE/NRWout99/ -r 20121105 -s 99
```

The representative sequence files were subsequently renamed to ``97_otus.fasta`` and ``99_otus.fasta``, respectively.


Filter the UNITE taxonomy to files containing only the 97% and 99% OTUs
-----------------------------------------------------------------------

Steps performed in IPython 0.13

```
from cogent.parse.fasta import MinimalFastaParser

g = open('./99_otu_taxonomy.txt','w')
otu_ids = []
for seq_id,_ in MinimalFastaParser(open('../rep_set/99_otus.fasta','U')):
    otu_ids.append(seq_id)
otu_ids = set(otu_ids)

for line in open('otu_taxonomy.txt','U'):
    if line.strip().split()[0] in otu_ids:
        g.write(line)
g.close()

g = open('./97_otu_taxonomy.txt','w')
otu_ids = []
for seq_id,_ in MinimalFastaParser(open('../rep_set/97_otus.fasta','U')):
    otu_ids.append(seq_id)
otu_ids = set(otu_ids)

for line in open('otu_taxonomy.txt','U'):
    if line.strip().split()[0] in otu_ids:
        g.write(line)
g.close()
```

Using these reference OTUs to process ITS Data
-------------------------------------------------

Define a [parameters file](http://qiime.org/documentation/file_formats.html#qiime-parameters) (``params.txt``) containing these values:

```
pick_otus:otu_picking_method    uclust_ref
pick_otus:enable_rev_strand_match       True
```

If you want to include taxonomy assignment in this process, you should have the following settings in your QIIME config (note: [in the near future](https://github.com/qiime/qiime/issues/468) you will not have to do this - these options will be hooked up as command line options).

```
assign_taxonomy_reference_seqs_fp	<path_to_its-reference-otus>/rep_set/97_otus.fasta
assign_taxonomy_id_to_taxonomy_fp	<path_to_its-reference-otus>/taxonomy/97_taxonomy.txt
```

Demultiplex sequences and apply [QIIME's subsampled open-reference OTU picking protocol](http://qiime.org/tutorials/open_reference_illumina_processing.html#subsampled-otu-picking-workflow-evaluation). The demultiplexing step shown here illustrates how to work with barcoded Illumina reads, but this protocol illustrates [how to demultiplex 454 data](http://qiime.org/tutorials/tutorial.html).

```
export REFSEQS_FP=<path_to_its-reference-otus>/rep_set/97_otus.fasta

split_libraries_fastq.py -i $PWD/sequence_reads.fastq -b $PWD/barcode_reads.fastq -m $PWD/mapping_file.txt -o $PWD/slout/ --rev_comp_mapping_barcodes

pick_subsampled_reference_otus_through_otu_table.py -i $PWD/slout/seqs.fna -o $PWD/ucrss/ -p $PWD/params.txt -aO <number_of_cores> --prefilter_percent_id=0.75 -r $REFSEQS_FP
```

That's it! After this step you will have an OTU table.

Citing this data
----------------

If you use this dataset, please cite the following to credit the QIIME and UNITE development groups:

The UNITE database for molecular identification of fungi - recent updates and future perspectives.
Abarenkov, Kessy; Nilsson, R. Henrik; Larsson, Karl-Henrik; Alexander, Ian J.; Eberhardt, Ursula; Erland, Susanne; Høiland, Klaus; Kjøller, Rasmus; Larsson, Ellen; Pennanen, Taina; Sen, Robin; Taylor, Andy F. S.; Tedersoo, Leho; Ursing, Björn M.; Vrålstad, Trude; Liimatainen, Kare; Peintner, Ursula; Kõljalg, Urmas (2010). New Phytologist, 186(2), 281-285.

QIIME allows integration and analysis of high-throughput community sequencing data.
J Gregory Caporaso; Justin Kuczynski; Jesse Stombaugh; Kyle Bittinger; Frederic D Bushman; Elizabeth K Costello; Noah Fierer; Antonio Gonzalez Pena; Julia K Goodrich; Jeffrey I Gordon; Gavin A Huttley; Scott T Kelley; Dan Knights; Jeremy E Koenig; Ruth E Ley; Catherine A Lozupone; Daniel McDonald; Brian D Muegge; Meg Pirrung; Jens Reeder; Joel R Sevinsky; Peter J Turnbaugh; William A Walters; Jeremy Widmann; Tanya Yatsunenko; Jesse Zaneveld; Rob Knight. (2010). Nature Methods 7: 335-336.
