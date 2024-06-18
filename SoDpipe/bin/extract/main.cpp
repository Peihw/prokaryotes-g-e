#include "extract.h"
using namespace std;

/*
0.	Information preparation for all prokaryotes.Done!
0.1	Sequence and annotation file association for all prokaryotes.Done!
0.2	Taxonomic annotation for all prokaryotes(both automatically and manually).Done!
0.3	Genome ID association for all prokaryotes(manually).Done!

1.	TriTISA correction for all translation initiation sites(TIS).Done!
1.0	TriTISA input preparation.Done!
1.0.1	Ptt format input preparation.Done!
1.0.2	Fasta format input for each chromosome and plasmid.Done!
1.1	TriTISA correction for all TISs.Done!
1.2	TriTISA output inclusion in records.Done!

2.	TIS model learning.
2.0	TIS sequence extraction for model learning.
2.1	EM learning on all TIS sequences.
2.2	Type classification for all signals.
2.2.1	Calculation of signal distances.
2.2.2	Type classification by distances.
2.3	Signal scanning for all genes.
2.4	TIS signal information inclusion in records.

3.	Ribosome binding site(RBS) translation rate prediction
3.0	RBS	input sequence preparation.
3.0.1	5' upstream untranslational region sequences.
3.0.2	5' region of mRNA.
3.1	mRNA folding free energy calculated by ViennaRNA.
3.2	Folding results integration and calculate full energy.
3.3	RBS translation rate inclusion in records.

4.	Experimental Leaderless annotation inclusion.
4.0	RNA - seq data inclusion and standarization.
4.1	Leaderless TIS identification.

5.	COG annotation inclusion.

6.	Phylogenetic tree construction.
6.0	16S rRNA sequence extraction for all prokaryotes.
6.1	Phylogenetic distance calculation by Clustal Omega.
*/

int main(int argc, char* argv[]) {
	if (argc > 1) {
		Str Out, Geno, Opt;
		Ma_Str_Str fa, rfa, fna, codon, label;
		Ma_Str_Ve_Str key;
		Ma_Str_T tax;
		rec red;
		ftt feat;
		int len = 0;
		int num = 0;
		bool END = false;
		for (int i = 0; i < argc; i++) {
			Str Temp = argv[i];
			Str In;
			if (i != argc - 1)
				In = argv[i + 1];
			//Parameter Input
			if (Temp == "-o")
				Out = In;
			else if (Temp == "-g")
				Geno = In;
			else if (Temp == "-p")
				Opt = In;
			else if (Temp == "-a")
				feat = read_in_ftt(In);
			else if (Temp == "-n")
				fna = read_in_fna(In);
			else if (Temp == "-f")
				fa = read_in_fa(In, num);
			else if (Temp == "-rf")
				rfa = read_in_rfa(In, Geno);
			else if (Temp == "-r")
				red = read_in_record(In, Geno);
			else if (Temp == "-k")
				key = read_in_keywords(In);
			else if (Temp == "-c")
				codon = read_in_codon(In);
			else if (Temp == "-b")
				label = read_in_label(In);
			else if (Temp == "-l")
				convertFromString(len, In);
			else if (Temp == "-m")
				convertFromString(num, In);
			else if (Temp == "-t")
				tax = read_in_taxonomy(In);
			else if (Temp == "-e") {
				if (In == "5\'")
					END = true;
			}
			//Functional Options
			else if (Temp == "-G")
				cout << gene(Out, Geno, key, feat, fna, codon, label, len) << endl;
			else if (Temp == "-R")
				cout << regions(Out, Opt, Geno, red, fna, codon, len, num) << endl;
			else if (Temp == "-I")
				cout << ips(Out, Geno, fa, fna) << endl;
			else if (Temp == "-N")
				cout << nonredundant(Out, rfa) << endl;
			else if (Temp == "-16")
				cout << rRNA(Out, fa, tax, len) << endl;
		}
		return 0;
	}
	else {
		/*Str Ftt = "E:\\ProTISA2.0\\ftt\\GCF_002025205.1_ASM202520v1_feature_table.txt";
		ftt feat = read_in_ftt(Ftt);
		Str Fna = "E:\\ProTISA2.0\\ftt\\GCF_002025205.1_ASM202520v1_genomic.fna";
		Ma_Str_Str fna = read_in_fna(Fna);
		Str Key = "E:\\ProTISA2.0\\ftt\\16s.dat";
		Ma_Str_Ve_Str key = read_in_keywords(Key);
		Str Geno = "GCF_002025205.1_ASM202520v1";
		Str Out = "E:\\ProTISA2.0\\ftt\\";
		Str Cod = "E:\\ProTISA2.0\\ftt\\codon.dat";
		Ma_Str_Str codon = read_in_codon(Cod);
		Str Label = "E:\\ProTISA2.0\\ftt\\trans.dat";
		Ma_Str_Str label = read_in_label(Label);
		int len = 0;
		cout << gene(Out, Geno, key, feat, fna, codon, label, len) << endl;
		Str Out = "E:\\ProTISA2.0\\records\\";
		Str Opt = "tis";
		Str Geno = "GCF_000734895.2_ASM73489v2";
		Str Rec = "E:\\ProTISA2.0\\records\\GCF_000734895.2_ASM73489v2.tritisa.rec.dat";
		rec red = read_in_record(Rec, Geno);
		Str Fna = "E:\\ProTISA2.0\\records\\GCF_000734895.2_ASM73489v2_genomic.fna";
		Ma_Str_Str fna = read_in_fna(Fna);
		Str Cod = "E:\\ProTISA2.0\\records\\codon.dat";
		Ma_Str_Str codon = read_in_codon(Cod);
		int len = 70;
		int num = 35;
		len = 20;
		cout<< regions(Out, Opt, Geno, red, fna, codon, len, num) << endl;*/
		cout << "#################################################################" << endl
			<< "# Welcome using this sequence extraction tool!                  #" << endl
			<< "# Developed by Longshu Yang.                                    #" << endl
			<< "# Copyright: Center for Quantitative Biology, Peking University.#" << endl
			<< "# Version 1.0 released at 2017/12/18.                           #" << endl
			<< "# Latest modification released at 2017/12/18.                   #" << endl
			<< "#################################################################" << endl
			<< "This is a brief introduction!" << endl
			<< "1. Functional options:" << endl
			<< "1)	-G:	extraction of selected gene sequence;" << endl
			<< "		Necessary inputs: -o, -n, -k, -c, -b, -l and -a" << endl
			//<< "2)	-U:	Translation initiation rate calculation;" << endl
			//<< "		Necessary inputs: -o, -n, -k, -e, -l, -b and -a" << endl
			<< "2)	-B:	extraction of RBS upstream sequence;" << endl
			<< "		Necessary inputs: -o, -r, -n, -g, -b and -l" << endl
			<< "3)	-M:	extraction of mRNA sequence;" << endl
			<< "		Necessary inputs: -o, -r, -n, -g, -b and -l" << endl
			<< "2. Parameters for inputs:" << endl
			<< "1)	Numerical inputs:" << endl
			<< "	-l: <the length of RBS upstream untranslated sequences>" << endl
			<< "	-m: <the number of column as the label of input sequences>" << endl
			<< "2)	String inputs:" << endl
			<< "	-g: <the assembly accession ID of a genome>" << endl
			<< "	-e: <the end option, 5' or 3'>" << endl
			<< "3)	Input full directory without the expansion name of a file: " << endl
			<< "	-o:	<output directory>" << endl
			<< "	-n: <fna file directory>" << endl
			<< "	-a: <ftt file directory>" << endl
			<< "	-f: <fasta file directory>" << endl
			<< "	-r: <record file directory>" << endl
			<< "	-k: <file directory of keywords>" << endl
			<< "	-c: <file directory of codon list>" << endl
			<< "	-b: <file directory of transferring label>" << endl;
		return 1;
	}
}