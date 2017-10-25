Pausing detection in Ribosome Profiling data
============================================

Detection of ribosome pausing was performed in the following way.

  1. The P-sites were obtained from the reads by shifting the reads by a fixed amount for each read length (shift by 11, 12, 13, 12 nucleotides for RPF reads of length 27, 28, 29, 30, respectively).  The amount of offset was determined by observing the position of the peak at the start of ORF for the given nucleotide length.
  2. Counts were combined for each codon in the reading frame, thus allowing for the reduction of systematic difference between the nucleotides in the codon, as well as the observed three-nucleotide periodicity, and making the algorithm less sensitive to errors in position of 1 nucleotide.
  3. For each gene individually the distribution of the counts for all the codons was modeled by a negative binomial distribution.
  4. For each codon the ratio of its count to the average over the gene and the P-value according to the fitted negative binomial distributions were obtained. These two numbers allow to quantify the level of pausing and the statistical reliability of the pausing.

Installation
------------

### Requirements

- R > 3.3.1
  + GenomicAlignments
  + GenomicFeatures
  + RiboProfiling
  + rtracklayer

Usage
-----

Example of the use is given in example.R

Note, that initial analysis/alignment of the data can be made using the pipeline https://github.com/fedxa/RPF_pipeline
