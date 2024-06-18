#include "EM.h"
#include "rps.h"
#include "record.h"
#include "estimate.h"
#include "classify.h"
using namespace std;


int main(int argc, char *argv[]) {
	Str Out, GenoID, Fa;
	int wid, cand, mod, step, len;
	wid = cand = mod = step = len = 0;
	double eps, maj, min, sd, ta;
	eps = maj = min = sd = ta = 0;
	Ma_Str_X tax;
	Ma_Str_Str fa;
	Ma_C_D pb;
	Ma_Str_M ref, que;
	Ma_Str_Str cog;
	Ma_Str_Pa_Str_Str rps1;
	Ma_Str_Ve_Str dis;
	rec red;
	if (argc > 2) {
		for (int i = 0; i < argc; i++) {
			//Input arguments
			Str Temp = argv[i];
			if (Temp == "-o")
				Out = argv[i + 1];
			else if (Temp == "-g")
				GenoID = argv[i + 1];
			else if (Temp == "-b") {
				Temp = argv[i + 1];
				Temp += GenoID + "_genomic.fna";
				pb = background_probability(Temp);
			}
			else if (Temp == "-f") {
				Fa = argv[i + 1];
				//Temp += GenoID + ".tis.fa";
				//fa = read_in_fa(Temp);
			}
			else if (Temp == "-t") {
				Temp = argv[i + 1];
				tax = read_in_taxonomy(Temp);
			}
			else if (Temp == "-r") {
				Temp = argv[i + 1];
				ref = read_in_motif(Temp);
			}
			else if (Temp == "-q") {
				Temp = argv[i + 1];
				Temp += GenoID + ".sig.dat";
				que = read_in_motif(Temp);
			}
			else if (Temp == "-p") {
				Temp = argv[i + 1];
				rps1 = read_in_rps1(Temp);
			}
			else if (Temp == "-cog") {
				Temp = argv[i + 1];
				Temp += GenoID + ".cog";
				cog = read_in_cog(Temp);
			}
			else if (Temp == "-m") {
				Temp = argv[i + 1];
				convertFromString(mod, Temp);
			}
			else if (Temp == "-d") {
				Temp = argv[i + 1];
				dis = read_in_distance(Temp);
			}
			else if (Temp == "-w") {
				Temp = argv[i + 1];
				convertFromString(wid, Temp);
			}
			else if (Temp == "-s") {
				Temp = argv[i + 1];
				convertFromString(step, Temp);
			}
			else if (Temp == "-c") {
				Temp = argv[i + 1];
				convertFromString(cand, Temp);
			}
			else if (Temp == "-l") {
				Temp = argv[i + 1];
				convertFromString(len, Temp);
			}
			else if (Temp == "-e") {
				Temp = argv[i + 1];
				convertFromString(eps, Temp);
			}
			else if (Temp == "-j") {
				Temp = argv[i + 1];
				convertFromString(maj, Temp);
			}
			else if (Temp == "-n") {
				Temp = argv[i + 1];
				convertFromString(min, Temp);
			}
			else if (Temp == "-rd") {
				Temp = argv[i + 1];
				Temp += GenoID + ".tritisa.rec.dat";
				red = read_in_record(Temp, GenoID);
			}
			else if (Temp == "-re") {
				Temp = argv[i + 1];
				Temp += GenoID + ".exp.rec.dat";
				red = read_in_record(Temp, GenoID);
			}
			else if (Temp == "-rt") {
				Temp = argv[i + 1];
				Temp += GenoID + ".tis.rec.dat";
				red = read_in_record(Temp, GenoID);
			}
			else if (Temp == "-ts") {
				Temp = argv[i + 1];
				Temp += GenoID + ".trans.rec.dat";
				red = read_in_record(Temp, GenoID);
			}
			else if (Temp == "-sd") {
				Temp = argv[i + 1];
				convertFromString(sd, Temp);
			}
			else if (Temp == "-ta") {
				Temp = argv[i + 1];
				convertFromString(ta, Temp);
			}
			//Functional parameters
			else if (Temp == "-M") {
				Fa += GenoID + ".tis.fa";
				fa = read_in_fa(Fa);
				cout << TIS_motif(Out, GenoID, fa, cand, mod, wid, maj, min, eps, step, pb) << endl;
			}
			else if (Temp == "-ME") {
				fa = read_in_fa(Fa);
				cout << TIS_motif(Out, GenoID, fa, cand, mod, wid, maj, min, eps, step, pb) << endl;
			}
			else if (Temp == "-C")
				cout << dist(Out, ref, que, GenoID, tax, ta, sd, rps1) << endl;
			else if (Temp == "-R") {
				Fa += GenoID + ".tis.fa";
				fa = read_in_fa(Fa);
				cout << tis_record(Out, GenoID, fa, red, que[GenoID], dis, len) << endl;
			}
			else if (Temp == "-E") {
				cout << estimate(Out, GenoID, red) << endl;
			}
			else if (Temp == "-P")
				cout << rps(Out, GenoID, red, tax, cog) << endl;
		}
		return 0;
	}
	else {
		cout << "########################################################################" << endl
			<< "# An EM based translation initiation site (TIS) prediction program.    #" << endl
			<< "# Developed by Longshu Yang.                                           #" << endl
			<< "# Copyright: Center for Quantitative Biology, Peking University.       #" << endl
			<< "# Version 1.0 released at 2018/02/02.                                  #" << endl
			<< "########################################################################" << endl
			<< "This is a brief introduction!" << endl
			<< "1. Functional options:" << endl
			<< "1)   -M: Motif training for TIS;" << endl
			<< "         Necessary inputs:   -o, -f, -k, -m, -w and -r;" << endl
			<< "2. Parameters for input:" << endl
			<< " 1) Numerical inputs:" << endl
			<< "     -k: <rank of possible initiation kmer>" << endl
			<< "     -m: <models of TIS>" << endl
			<< "     -r: <possible region of a TIS>" << endl
			<< "     -w: <width of signal>" << endl
			<< " 2) Input full directory and expansion name of a file:" << endl
			<< "     -f: <readable fasta file>" << endl
			<< "     -o: <output directory and file name>" << endl;
		return 1;
	}
}
