#!/bin/bash

perl nBlock_tester_npv.pl Transcribed-non-protein-genes_regions/ncRNA-or-intergenic_regions.bed  Mapped_tag_reduced_data/12d4_S7.mapped.ncRNA.bam


echo -ne "\n################################### Sampling of Results ###################################\n"

head -n 2 SmallRNA_patterns_unnanotated/ncRNA-or-intergenic_regions.sm-blocks.bed 


echo -ne "\n################################### Interpretation ###################################\n"

echo " 
#### Columns output bed file  #
##### Subcolums or colums split by ;  ## 
# chr
# start whole expression block 
# end whole expresion block 
# name intergenic expressed region 
# Total de-duplicate read counds eg 2301
# Strand
# Info Column 
## sl= coordinates intergenic expressed region 
## dl= block of reads following a normal distribution pattern 
## sc= repetitivity of start count  (how many reads support the same start coordinate?) 
## ec= repetitivity of end count 
## soc = count of total redundant reads or extract total read count 
## tc = tag count or non redundant reads counts
## rm = raw block middle position calculation  
## bm = block mean calculation using total reads population
#----------- > coincidence of these two measures is required to ensure quality of prediction
## sd = standart deviation
# length small predicted fragment
" 


# Data sample 7 and 8 column" 

##sl=564875:564997;dl=564921:564998;sc=1095;ec=2262;soc=2713;tc=211;rm=564937;bm=564936.5;sd=14.00        29 
