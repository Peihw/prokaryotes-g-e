# SoDpipe

## Description
  We developed SoDpipe, an automated bioinformatic pipeline for the identification of redundant genes and translation regulatory elements. The basic protocols consist of redundant gene and translation initiation signal annotation, which in conjunction with virulence factors and resistance genes annotation facilitate large-scale analysis of evolutionary events dominated by redundant genes for prokaryotes, especially for pathogens. The pipeline can handle different types of genomic data for various implementation requirements, including but not limited to genome assembly data and whole genome sequencing data. For genome assembly data from GenBank/RefSeq, the pipeline mainly includes the identification of redundant genes, translation initiation site correction, translation initiation signal annotation, and VF/AMR annotation. Quality control, genome assembly, and genome annotation were also implemented when processing whole genome sequencing data. The final report generated contains detailed information on the redundant clusters, types of translation initiation signals (SD-like, TA-like, or no signal), signal motifs, and the start site of the signal.

## Workflow
<img src="./workflow.png" alt="Alt Text" width="500" height="700">


## Installation

### Download from the website
1. Download and Install the Miniconda (Optional, depends on usage).
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b
```
2. Download SoDpipe.
```
cd SoDpipe
conda env create --name SoDpipe --file Install/environment.yaml
conda activate SoDpipe
```
3. Install the dependencies listed in Install/Requirements and setup dbs.
```
prokka –setupdb
diamond makedb --in bin/db/CARD.faa -d bin/db/CARD
diamond makedb --in bin/db/VFDB.faa -d bin/db/VFDB
```

### Docker 
```
docker pull peihw/sodpipe:1.20
```
