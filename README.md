# SoDpipe

## Description
  We developed SoDpipe, an automated bioinformatic pipeline for the identification of redundant genes and translation regulatory elements. The basic protocols consist of redundant gene and translation initiation signal annotation, which in conjunction with virulence factors and resistance genes annotation facilitate large-scale analysis of evolutionary events dominated by redundant genes for prokaryotes, especially for pathogens.
  The pipeline can handle different types of genomic data for various implementation requirements, including but not limited to genome assembly data and whole genome sequencing data. For genome assembly data from GenBank/RefSeq, the pipeline mainly includes the identification of redundant genes, translation initiation site correction, translation initiation signal annotation, and VF/AMR annotation. Quality control, genome assembly, and genome annotation were also implemented when processing whole genome sequencing data. The final report generated contains detailed information on the redundant clusters, types of translation initiation signals (SD-like, TA-like, or no signal), signal motifs, and the start site of the signal.


### For installation



1.Docker 
```
docker pull peihw/sodpipe:1.00
```

2.Download from the website

