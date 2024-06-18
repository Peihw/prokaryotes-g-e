#include "prokaryotes.h"
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
		Str Pa, In, Out, Fna, Ftt, NC;
		Ma_Str_Str fa, trans, cog, csv, exftt;
		ftt ano;
		Ma_Str_K keg;
		Ve_Str idf, idt;
		Ma_I_M med;
		Ma_Str_L bla;
		//Ma_I_E exp;
		Ve_E exp;
		Ma_Str_E trs;
		rec red;
		Ma_Str_T tu;
		Ve_Str genes;
		Ve_Pa_Str_I opr;
		for (int i = 1; i < argc; i++) {
			Pa = argv[i];
			if (i != argc - 1)
				In = argv[i + 1];
			//Parameters
			if (Pa == "-o")
				Out = In;
			else if (Pa == "-nc")
				NC = In;
			else if (Pa == "-n")
				fa = read_in_fna(In);
			else if (Pa == "-f")
				fa = read_in_fa(In);
			else if (Pa == "-a")
				ano = read_in_ftt(In);
			else if (Pa == "-k")
				keg = read_in_kegg(In);
			else if (Pa == "-e") {
				In += Fna + ".exp.dat";
				exp = read_in_exp(In);
			}
			else if (Pa == "-if")
				idf = read_in_column(In);
			else if (Pa == "-it")
				idt = read_in_column(In);
			else if (Pa == "-fn")
				Fna = In;
			else if (Pa == "-ft")
				Ftt = In;
			else if (Pa == "-m")
				med = read_in_med(In);
			else if (Pa == "-r") {
				In += Fna + ".tritisa.rec.dat";
				red = read_in_record(In, Fna);
			}
			else if (Pa == "-re") {
				In += Fna + ".exp.rec.dat";
				red = read_in_record(In, Fna);
			}
			else if (Pa == "-tf")
				trans = read_in_trans(In);
			else if (Pa == "-v")
				csv = read_in_cog_csv(In);
			else if (Pa == "-b")
				bla = read_in_blast(In);
			//else if (Pa == "-tr")
			//	trs = read_in_trans_record(In);
			else if (Pa == "-tu") {
				In += Fna + ".ope.dat";
				tu = read_in_tu(In);
			}
			else if (Pa == "-ef")
				exftt = read_in_exftt(In);
			else if (Pa == "-tg")
				genes = read_in_tu_genes(In);
			else if (Pa == "-op")
				opr = read_in_opr(In);
			//Functional options
			//0.1	Sequence and annotation file association for all prokaryotes.
			else if (Pa == "-R")
				cout << correlate(Out, idf, idt) << endl;
			//0.2	Prokaryotes taxonomic annotation.
			else if (Pa == "-P")
				cout << prokaryotes(Fna, Ftt, Out, fa, keg, ano) << endl;
			//1.0.1	Ptt format input preparation.
			else if (Pa == "-M")
				cout << ptt_output(Out, NC, ano) << endl;
			//1.0.2	Fasta format input for each chromosome and plasmid.
			else if (Pa == "-N")
				cout << nc_fasta(Out, fa) << endl;
			//1.2	TriTISA output postprocession as records.
			else if (Pa == "-T")
				cout << tritisa_record(Out, NC, med, ano, trans) << endl;
			//4.1	Leaderless TIS identification.
			else if (Pa == "-E")
				cout << exp_map(Out, Fna, red, fa, exp) << endl;
			//5.2	BLAST results classified into different COG ID.
			else if (Pa == "-B2C")
				cout << blast2cog(Out, NC, ano, bla, csv) << endl;
			else if (Pa == "-RA")
				cout << reannotate(Out, Fna, fa, ano) << endl;
			else if (Pa == "-TU")
				cout << transcript_map(Out, Fna, tu, red) << endl;
			else if (Pa == "-OP")
				cout << operon(Out, Fna, genes, opr, exftt) << endl;
		}
		return 0;
	}
	else {
		/*Str Out = "E:\\ProTISA2.0\\conditional and internal\\ope\\";
		Str Geno = "GCF_000008445.1_ASM844v1";
		//Str Genes = "E:\\ProTISA2.0\\conditional and internal\\ope std\\GCF_000008425.1_ASM842v1.pre.dat";
		//Ve_Str genes = read_in_tu_genes(Genes);
		Ve_Str genes;
		Str Ftt = "E:\\ProTISA2.0\\conditional and internal\\ftt\\GCF_000008445.1_ASM844v1_feature_table.txt";
		Ma_Str_Str exftt = read_in_exftt(Ftt);
		Str Opr = "E:\\ProTISA2.0\\conditional and internal\\ope std\\GCF_000008445.1_ASM844v1.pre.dat";
		Ve_Pa_Str_I opr = read_in_opr(Opr);
		//Ve_Pa_Str_I opr;
		cout << operon(Out, Geno, genes, opr, exftt) << endl;*/
		Str Fna = "GCF_000005845.2_ASM584v2";
		Str Out = "E:\\ProTISA2.0\\conditional and internal\\exp\\";
		//Str Rec = "E:\\ProTISA2.0\\conditional and internal\\tritisa\\" + Fna + ".tritisa.rec.dat";
		Str Rec = "E:\\ProTISA2.0\\Transcriptome\\experiment\\Bacteria\\" + Fna + ".exp.rec.dat";
		rec red = read_in_record(Rec, Fna);
		//Str Fa = "E:\\ProTISA2.0\\conditional and internal\\exp\\" + Fna + "_genomic.fna";
		//Ma_Str_Str fna = read_in_fna(Fa);
		Str Exp = "E:\\ProTISA2.0\\conditional and internal\\exp\\Bacteria\\" + Fna + ".exp.dat";
		Ve_E exp = read_in_exp(Exp);
		Str TU = "E:\\ProTISA2.0\\conditional and internal\\exp\\" + Fna + ".ope.dat";
		Ma_Str_T tu = read_in_tu(TU);
		cout << transcript_map(Out, Fna, tu, red) << endl;
		//cout<< exp_map(Out, Fna, red, fna, exp) << endl;
		/*Str Ano = "E:\\ProTISA2.0\\records\\debug\\GCF_000011005.1_ASM1100v1_feature_table.txt";
		ftt ano = read_in_ftt(Ano);
		Str Rec = "E:\\ProTISA2.0\\records\\GCF_000007985.2_ASM798v2.tritisa.rec.dat";
		rec red = read_in_tritisa_record(Rec, "GCF_000007985.2_ASM798v2");
		Str Exp = "E:\\ProTISA2.0\\experiments\\Bacteria transcriptome\\Data set\\TSS\\GCF_000007985.2_ASM798v2.exp.tss.dat";
		Ma_I_E exp = read_in_exp_tss(Exp);
		Str Out = "E:\\ProTISA2.0\\experiments\\";
		Str Fa = "E:\\ProTISA2.0\\experiments\\GCF_000007985.2_ASM798v2_genomic.fna";
		Ma_Str_Str fna = read_in_fna(Fa);
		cout << tss_map(Out, "GCF_000007985.2_ASM798v2", red, fna, exp)<<endl;*/
		cout << "###########################################################################" << endl
			<< "# Welcome using this Translation Initiation Site Analysis Program!        #" << endl
			<< "# Developed by Longshu Yang.                                              #" << endl
			<< "# Copyright: Center for Quantitative Biology, Peking University.          #" << endl
			<< "# Version 2.0 released at 2018/01/01.                                     #" << endl
			<< "# Latest modification released at 2018/01/03.                             #" << endl
			<< "###########################################################################" << endl
			<< "This is a brief introduction!" << endl
			<< "1.	Functional options:" << endl
			<< "	1)	-R:	correlate fna and ftt files for the same assembly;" << endl
			<< "		Necessary inputs: -if, -it and -o;" << endl
			<< "	2)	-P:	list all prokaryotes with description;" << endl
			<< "		Necessary inputs: -fn, -ft, -f, -t, -k and -o" << endl
			<< "2.	Parameters for inputs:" << endl
			<< "	1)	File directory inputs:" << endl
			<< "		-o:	<directory of output file>" << endl
			<< "		-f: <directory of fna file>" << endl
			<< "		-t: <directory of ftt file>" << endl
			<< "		-k: <directory of kegg file>" << endl
			<< "		-if: <directory of fna file name list>" << endl
			<< "		-it: <directory of ftt file name list>" << endl
			<< "	2)	String inputs:" << endl
			<< "		-fn: <fna file name>" << endl
			<< "		-ft: <ftt file name>" << endl
			<< "Please input necessary parameters!" << endl;
		return 1;
	}
}