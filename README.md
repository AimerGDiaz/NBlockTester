NBlockTester: A bioinformatic tool to detect small ncRNA-derived
fragments
================
Aimer G. Diaz

## Description of the tool

NBlockTester is a perl code (in proccess to be a cpan module) written by
Aimer G. Diaz, with the supervision and comments of Clara Bermudez and
Steve Hoffmann. It’s an adaptation of the previously developed tool
[blockbuster
code](http://hoffmann.bioinf.uni-leipzig.de/LIFE/blockbuster.html),
which detect blocks of overlapping reads using a gaussian-distribution
approach. NBlockTester it’s specially adapted for the detection of small
ncRNAs fragments derived from longer ncRNAs like tRNAs, rRNAs, snoRNAs,
or even mRNAs [1](https://pubmed.ncbi.nlm.nih.gov/21741109/), using
mapped small RNA-seq data (bam files) and a set of genomic coordinates
which would be used as a search space. The main purpose of the code is
to discriminate small fragments (like tRFs -tRNA derived fragments-) who
has a miRNA-like expression pattern from background-source expression
(the tRNA source), recover in the small RNAseq protocol of sequencing
which might be degradation or technical byproducts.

## Installation

By now just download this repository.

## Running with test data

To Run the code with a small test data set, please just run the shell
file

    bash run.sh

Or alternativelly directly launch the perl script:

    perl nBlock_tester_npv.pl test_data/Transcribed-non-protein-genes_regions/ncRNA-or-intergenic_regions.bed  test_data/Mapped_tag_reduced_data/12d4_S7.mapped.ncRNA.bam

## Dependencies

The code it’s also written to be run in parallel (fork), however it does
requieres the cpan module “Parallel::ForkManager”. To install this
module just run:

    cpan upgrade Test::More

    cpan Parallel::ForkManager

Test installation by running

    perl -e 'use Parallel::ForkManager;'

## References

By now the tool only have been used as part of my master thesis,
[avaliable in spanish
here](https://repositorio.unal.edu.co/handle/unal/63688)

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
