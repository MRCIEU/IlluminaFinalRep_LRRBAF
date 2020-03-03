# IlluminaFinalRep_LRRBAF

This script was used to estimate log R Ratio (LRR) and B allele frequency (BAF) from
Illumina Final Report files with the headings:

> SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,X,Y,X Raw,Y Raw

The purpose of this script was to generate estimated LRR and BAF in the absence of these
values in the Final Report (and in the absence of original IDAT files and GenomeStudio files).

The estimated LRR and BAF values are intended as input for PennCNV.

The approach to estimating LRR and BAF is based on [Peiffer et al. "High-resolution genomic profiling of chromosomal aberrations using Infinium whole-genome genotyping". Genome Res. 2006. 16: 1136-1148](https://genome.cshlp.org/content/16/9/1136.full)
