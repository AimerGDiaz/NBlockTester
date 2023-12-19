NBlockTester: A bioinformatic tool to detect small ncRNA-derived
fragments
================
Aimer G. Diaz

## Description of the tool

NBlockTester is a perl code (in proccess to be a cpan module) written by
Aimer G. Diaz, with the supervision and comments of Clara Bermudez and
Steve Hoffmann. It’s an adaptation of the previously developed tool
[blockbuster
code](https://www.bioinf.uni-leipzig.de/~david/LIFE/LIFE/blockbuster.html)
\[[1](#ref-langenberger2009evidence)\], which detect blocks of
overlapping reads using a gaussian-distribution approach. NBlockTester
it’s specially adapted for the detection of small ncRNAs fragments
derived from longer ncRNAs like tRNAs, rRNAs, snoRNAs, or even mRNAs
\[[2](#ref-tuck2011rna)\], using mapped small RNA-seq data (bam files)
and a set of genomic coordinates which would be used as a search space.
The main purpose of the code is to discriminate small fragments (like
tRFs -tRNA derived fragments-) who has a miRNA-like expression pattern
from background-source expression (the tRNA source), recover in the
small RNAseq protocol of sequencing which might be degradation or
technical byproducts \[[3](#ref-gutierrez2018deteccion)\].

## Installation

By now just download this repository.

## Running with test data

To Run the code with a small test data set, please just run the shell
file

    bash run.sh

Or alternatively directly launch the perl script:

    perl nBlock_tester_npv.pl test_data/Transcribed-non-protein-genes_regions/ncRNA-or-intergenic_regions.bed  test_data/Mapped_tag_reduced_data/12d4_S7.mapped.ncRNA.bam

## Dependencies

The code it’s also written to be run in parallel (fork), however it does
requires the cpan module “Parallel::ForkManager”. To install this module
just run:

    cpan upgrade Test::More

    cpan Parallel::ForkManager

Test installation by running

    perl -e 'use Parallel::ForkManager;'

To run parallel version of the code just assing the number of threads to
be use:

    perl nBlock_tester.pl  2  test_data/Transcribed-non-protein-genes_regions/ncRNA-or-intergenic_regions.bed  test_data/Mapped_tag_reduced_data/12d4_S7.mapped.ncRNA.bam

## License

The code is freely available to download and run, but it’s protected and
licensed under a [Creative Commons Attribution-ShareAlike 4.0
International License](https://creativecommons.org/licenses/by-nc/4.0/),
meaning you can use it but citing it’s source.

[![License: CC BY-NC
4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

## Viewers

[![Image of
Viewers](https://github.com/AimerGDiaz/Viewers/blob/master/svg/409164432/badge.svg)](https://github.com/AimerGDiaz/Viewers/blob/master/readme/409164432/week.md)

## References

<div id="refs" class="references csl-bib-body">

<div id="ref-langenberger2009evidence" class="csl-entry">

1\. Langenberger D, Bermudez-Santana C, Hertel J, Hoffmann S, Khaitovich
P, Stadler PF. Evidence for human microRNA-offset RNAs in small RNA
sequencing data. Bioinformatics. 2009;25:2298–301.

</div>

<div id="ref-tuck2011rna" class="csl-entry">

2\. Tuck AC, Tollervey D. RNA in pieces. Trends in genetics.
2011;27:422–32.

</div>

<div id="ref-gutierrez2018deteccion" class="csl-entry">

3\. Gutiérrez Dı́az AA. Detección automatizada de pequeños fragmentos
derivados de RNAs no-codificantes expresados diferencialmente frente a
la infección del virus dengue. PhD thesis. 2018.

</div>

</div>
