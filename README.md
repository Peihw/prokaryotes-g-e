# Prokaryotic genome evolution: a moderate profusion of genetic elements

## Description
Ferther details of the codes and results related to the paper can be found herewithin. We provide scripts that allow you to reproduce the results or manage your own data.


### Analysis of translation initiation signal 
*Genome1_genomic.fna and Genome1_feature_table.txt are files of genome sequences and annotations, respectively, which can be download from the Refseq database. (ftp://ftp.ncbi.nlm.nih.gov/genomes/)*


1.Splitting the fasta file of genome assembly into fasta file containing single chromosome 
```
Prokaryotes -o fna/ -f Genome1_genomic.fna -N
```

2.Generating med file from genome annotation provided by Refseq database
```
./Prokaryotes -o med/ -nc Genome1 -a Genome1_feature_table.txt -M
```

3.Translation initiation site correction
```
./TriTISA Genome1.fna Genome1.med Tritisa/Genome1
```

4.Extracting TIS upstream sequences
```
./extract -o TISseq/ -g Genome1 -p tis -r result/Genome1.tritisa.rec.dat -n Genome1_genomic.fna -l 20 -R 
```

5.Probabilistic modeling for translation initiation signal
```
./TIS -o model -g Genome1 -f result/Genome1.tis.fa -pb Genome1_genomic.fna -mj 0.52 -mn 0.16 -ep 1.0e-4 -tr 4 -k 10 -r 20 -m 2 -w 8 -M
```

6.Signal classification 
```
./TIS -o signal/Bacteria.sig -g Genome1 -r Bacteria.ref.dat -q model/ -sd 6 -ta 0.75  -l 5 -C
```

7.Signal scanning
```
./TIS -o signal/ -g Genome1  -q model/ -f TISseq/ -d signal/Bacteria.sig -R
```
