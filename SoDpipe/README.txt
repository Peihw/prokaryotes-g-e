 ==================================
	Introduction
 ===================================
SoDpipe is an integration pipeline for automatically analyzing redundant genes and their regulations in prokaryotes. SoDpipe provides a framework to support genomic surveillance of the occurrence, gene expression and adaptive evolution from the perspective of gene redundancy for prokaryotes, especially for pathogen. SoDpipe takes both genome assembly and whole genome sequencing data as input and automatically performs the analysis of the data all through a single command-line instruction. It generates a detailed report of the duplicated gene clusters, function annotation, virulence factors, antimicrobe resistance, types of translation initiation mechanisms (SD-like, TA-like, Atypical, or no signal), translation initiation signal motifs, and possible promoter sequences.
 ===================================
	Step-by-Step tutorial
 ===================================
# Installation

## Download from the website
   Download and Install the Miniconda (Optional, depends on usage).
	- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	- bash miniconda.sh -b
   Download SoDpipe.
	- cd SoDpipe
   Install conda environment from the prebuilt environment YAML file.
	- conda env create --name SoDpipe --file Install/environment.yaml
	- conda activate SoDpipe
   Install the dependencies listed in Install/Requirements and setup dbs.
	- prokka –setupdb
	- diamond makedb --in bin/db/CARD.faa -d bin/db/CARD
	- diamond makedb --in bin/db/VFDB.faa -d bin/db/VFDB

## Docker
	- docker pull peihw/sodpipe:1.20 

# Invoking SoDpipe

## For genome assembly from GenBank/RefSeq

   1. Input files
	SoDpipe takes FASTA format of the genomic sequences (*_genomic.fna), tab-delimited text file reporting annotated features (*_feature_table.txt), FASTA format of the sequences corresponding to all CDS (*_translated_cds.faa or other translated amino acid sequence) as input.

   2. Configuration
	The config_a.yaml file allows you to define the parameters and their values.
	----------------------------------------------------------------------------
	#Input genome assembly files
	Genomefna: path to your genomic sequence file, e.g., testdata/GCF_000005845.2_ASM584v2_genomic.fna.
	Annotation: path to your annotation file, e.g., testdata/GCF_000005845.2_ASM584v2_feature_table.txt.
	cds: path to your CDS sequence file, e.g., testdata/GCF_000005845.2_ASM584v2_translated_cds.faa.
	List: path to your file with a list of accession numbers in the genomic file, e.g., testdata/GCF_000005845.2_ASM584v2.list.

	#Parameters
	Reference: bin/parameters/Bacteria.ref.dat	#Bacteria.ref.dat for bacteria and Archaea.ref.dat for archaea.
	Length: 20	#length of untranslated regions upstream translation initiation site. 20 was recommended for bacteria and 30 for archaea.
	identity: 80	#identity threshold of duplicated genes.
	coverage: 0.8		#coverage threshold of duplicated genes.
	query_id: 60	#minimum identity% for VF/AMR annotation.
	query_cover: 50	#minimum query cover% for VF/AMR annotation.
        ----------------------------------------------------------------------------

   3. Software usage
	#Getting translation initiation signal annotation of all genes in, e.g., GCF_000005845.2_ASM584v2 with 20 cores.
	- snakemake -s Snakefile.a result/records/GCF_000005845.2_ASM584v2.tis.rec.dat --jobs 20
	#Getting possible promoter sequence of all operon in, e.g., GCF_000005845.2_ASM584v2.
	- snakemake -s Snakefile.a result/promoter/GCF_000005845.2_ASM584v2.mp.proseq --jobs 20
	#Getting duplicated gene clusters
	- snakemake -s Snakefile.a result/duplication/GCF_000005845.2_ASM584v2.dup --jobs 20
	#Getting translation initiation signal annotation of duplicated genes.
	- snakemake -s Snakefile.a result/duplication/GCF_000005845.2_ASM584v2.dup.tis --jobs 20
	#Getting detailed report of duplicated genes.
	- snakemake -s Snakefile.a result/duplication/GCF_000005845.2_ASM584v2.info --jobs 20

	#Translation initiation site reannotation was enabled by default, users can enable or depress the rule of tritisa or without_tritisa by adding or removing the "#" at the beginning of the sentence in Snakefile.a.

   4. Output files
	result/Tritisa/GCF_000005845.2_ASM584v2.tritisa.rec.dat	#reannotation feature file. 
	result/TISseq/GCF_000005845.2_ASM584v2.tis.fa	#FASTA format of untranslated sequence upstream translation initiation site.
	result/records/GCF_000005845.2_ASM584v2.tis.rec.dat	#translation initiation signal annotation of all genes within the genome.
	result/duplication/GCF_000005845.2_ASM584v2.dup	#duplicated gene clustered at setting threshold.
	result/duplication/GCF_000005845.2_ASM584v2.dup.tis	#translation initiation signal annotation of duplicated genes within the genome.
	result/TSseq/GCF_000005845.2_ASM584v2.tis.fa	#untranslated sequence for promoter detection
	result/promoter/GCF_000005845.2_ASM584v2.mp.proseq	#promoter sequence of the maximum posibility.
	result/duplication/GCF_000005845.2_ASM584v2.info	#a detailed report of the duplicated gene clusters, function annotation, virulence factors, antimicrobe resistance, types of translation initiation mechanisms, translation initiation signal motifs, and possible promoter sequences.


## For whole genome sequencing data

   1. Input files
	SoDpipe also takes paired-end whole genome sequencing data (*_1.fastq, *_2.fastq) as input.

   2. Configuration
	The config_s.yaml file allows you to define the parameters and their values.
        ----------------------------------------------------------------------------
	#Input Whole Genome Sequencing files
	Fastq1: path to your FASTQ format file, e.g., testdata/SRR19707997_1.fastq.
	Fastq2: path to your FASTQ format file, e.g., testdata/SRR19707997_2.fastq.
	Fasta: path to your FASTA format file, e.g., result/spades_out/SRR19707997/contigs.fasta. It can be the output file of Spades or your FASTA format file.

	#Parameters
	Genus: Acinetobacter	#specify a genus for prokka.
	Species: baumannii	#specify a species for prokka.
	Reference: bin/parameters/Bacteria.ref.dat	#Bacteria.ref.dat for bacteria and Archaea.ref.dat for archaea.
	Length: 20	#length of untranslated regions upstream translation initiation site. 20 was recommended for bacteria and 30 for archaea.
	identity: 80	#identity threshold of duplicated genes.
	coverage: 0.8		#coverage threshold of duplicated genes.
	query_id: 60    #minimum identity% for VF/AMR annotation.
        query_cover: 50 #minimum query cover% for VF/AMR annotation.
        ----------------------------------------------------------------------------

   3. Software usage
	#Getting genome annotation of sequencing data with 20 cores.
	- snakemake -s Snakefile.s result/prokka_out/SRR19707997/SRR19707997.fna --jobs 20
	#Getting translation initiation signal annotation of all genes with 20 cores.
	- snakemake -s Snakefile.s result/records/SRR19707997.tis.rec.dat --jobs 20
	#Getting translation initiation signal annotation of duplicated genes with 20 cores.
	- snakemake -s Snakefile.s result/duplication/SRR19707997.dup.tis --jobs 20
	#Getting detailed report of duplicated genes.
	- snakemake -s Snakefile.s result/duplication/SRR19707997.info --job 20
	#Translation initiation site reannotation was disabled by default when analyzing contigs data, users can enable tritisa by using the “Snakefile.st”file. e.g.,
	- snakemake -s Snakefile.st result/records/SRR19707997.tis.rec.dat --jobs 20
	Users can also take scaffolds as input by editing the config_s.yaml file, e.g.,
		Fasta: result/spades_out/SRR19707997/scaffolds.fasta

   4. Output files
	result/spades_out/SRR19707997/	#genome assembly files.
	result/prokka_out/SRR19707997/	#genome annotation files.
	result/Tritisa/SRR19707997.tritisa.rec.dat	#reannotation feature file. 
	result/TISseq/SRR19707997.tis.fa	#FASTA format of untranslated sequence upstream translation initiation site.
	result/records/SRR19707997.tis.rec.dat	#translation initiation signal annotation of all genes within the genome.
	result/duplication/SRR19707997.dup	#duplicated gene clustered at setting thresholds.
	result/duplication/SRR19707997.dup.tis	#translation initiation signal annotation of duplicated genes within the genome.
	result/duplication/SRR19707997.info	#a detailed report of the duplicated gene clusters, function annotation, virulence factors, antimicrobe resistance, types of translation initiation mechanisms, translation initiation signal motifs, and possible promoter sequences.

   
## For genomic data in other formats
	We recommend reannotating the genome with Prokka and continuing the following analyses with "Snakefile.s" or "Snakefile.st". Users can edit the "Fasta" option in "config_s.yaml" with the information of your FASTA file.

