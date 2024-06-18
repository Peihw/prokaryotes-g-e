#include "EM.h"
#include "scoring.h"
#include "tune.h"
using namespace std;


int main(int argc, char *argv[]) {
	Str Out, Fa, GenoID, Ori, Que, Ref;
	int mod, wid, reg, trial, thr;
	reg = wid = mod = trial = thr = 0;
	double tthr, sthr, pthr, maj, min, eps;
	tthr = sthr = pthr = maj = min = eps = 0;
	Ma_Str_M que, ref;
	ftt ft;
	Ma_Str_Str fa, cog;
	Ma_Str_Pa_Str_Str check;
	Ma_Str_Ve_Str dis;
	Ma_Str_Ve_L o2q, q2o;
	Ma_Str_T tss;
	//Ma_Str_R orr, qur;
	Ma_Str_Ma_Str_O orth;
	Ma_Str_X tax;
	Ma_C_D pb;
	rec red;
	if (argc > 2) {
		for (int i = 0; i < argc; i++) {
			//Input arguments
			Str Temp = argv[i];
			if (Temp == "-o")
				Out = argv[i + 1];
			else if (Temp == "-f")
				Fa = argv[i + 1];
			else if (Temp == "-a") {
				Temp = argv[i + 1];
				ft = read_in_ftt(Temp);
			}
			else if (Temp == "-w") {
				Str Wid = argv[i + 1];
				convertFromString(wid, Wid);
			}
			else if (Temp == "-r") {
				Str Reg = argv[i + 1];
				convertFromString(reg, Reg);
			}
			else if (Temp == "-k") {
				Str Rank = argv[i + 1];
				convertFromString(thr, Rank);
			}
			else if (Temp == "-m") {
				Str Mod = argv[i + 1];
				convertFromString(mod, Mod);
			}
			else if (Temp == "-t") {
				Str TA = argv[i + 1];
				convertFromString(tthr, TA);
			}
			else if (Temp == "-s") {
				Str SD = argv[i + 1];
				convertFromString(sthr, SD);
			}
			else if (Temp == "-p") {
				Str SP = argv[i + 1];
				convertFromString(pthr, SP);
			}
			else if (Temp == "-g")
				GenoID = argv[i + 1];
			else if (Temp == "-rf")
				Ref = argv[i + 1];
			else if (Temp == "-ts") {
				Temp = argv[i + 1];
				Temp += GenoID + ".trans.rec.dat";
				red = read_in_record(Temp, GenoID);
			}
			else if (Temp == "-ex") {
				Temp = argv[i + 1];
				Temp += GenoID + ".exp.rec.dat";
				red = read_in_record(Temp, GenoID);
			}
			else if (Temp == "-q") {
				Temp = argv[i + 1];
				que = read_in_motif(Temp);
				//cout << que.rwm.size() << endl;
			}
			else if (Temp == "-e") {
				Temp = argv[i + 1];
				ref = read_in_motif(Temp);
			}
			else if (Temp == "-d") {
				Temp = argv[i + 1];
				dis = read_in_distance(Temp);
			}
			else if (Temp == "-rd") {
				Temp = argv[i + 1];
				red = read_in_record(Temp, GenoID);
			}
			else if (Temp == "-c") {
				Temp = argv[i + 1];
				Temp += GenoID + ".cog";
				cog = read_in_cog(Temp);
			}
			/*else if (Temp == "-or") {
				Temp = argv[i + 1];
				orr = read_in_record(Temp);
			}
			else if (Temp == "-qr") {
				Temp = argv[i + 1];
				qur = read_in_record(Temp);
			}*/
			else if (Temp == "-ot") {
				Temp = argv[i + 1];
				orth = read_in_ortholog(Temp);
			}
			else if (Temp == "-og")
				Ori = argv[i + 1];
			else if (Temp == "-qg")
				Que = argv[i + 1];
			/*else if (Temp == "-oro") {
				Temp = argv[i + 1];
				orr = read_in_record_v2(Temp);
			}
			else if (Temp == "-qro") {
				Temp = argv[i + 1];
				qur = read_in_record_v2(Temp);
			}*/
			else if (Temp == "-pb") {
				Temp = argv[i + 1];
				pb = background_probability(Temp);
			}
			else if (Temp == "-mj") {
				Temp = argv[i + 1];
				convertFromString(maj, Temp);
			}
			else if (Temp == "-mn") {
				Temp = argv[i + 1];
				convertFromString(min, Temp);
			}
			else if (Temp == "-ep") {
				Temp = argv[i + 1];
				convertFromString(eps, Temp);
			}
			else if (Temp == "-tr") {
				Temp = argv[i + 1];
				convertFromString(trial, Temp);
			}
			else if (Temp == "-rp") {
				Temp = argv[i + 1];
				check = read_in_rps1(Temp);
			}
			else if (Temp == "-x") {
				Temp = argv[i + 1];
				tax = read_in_taxonomy(Temp);
			}
			//Functional parameters
			else if (Temp == "-M") {
				fa = read_in_fa(Fa);
				cout << TIS_motif(Out, GenoID, fa, thr, mod, wid, maj, min, eps, trial, pb) << endl;
			}
			else if (Temp == "-C")
				cout << dist(Out, ref, que, GenoID, reg, tthr, sthr, pthr, check) << endl;
			else if (Temp == "-R"){
				fa = read_in_fas(Fa);
				cout << tis_record(Out, GenoID, fa, red, que[GenoID], dis, reg) << endl;
			}
			else if (Temp == "-T") {
				fa = read_in_fa(Fa);
				cout << tune(Out, fa, thr, mod, wid, reg, maj, min, eps, pb) << endl;
			}
			else if (Temp == "-M2F") 
				cout << motif2fasta(Out, GenoID, que[GenoID], reg) << endl;
			//else if(Temp == "-SP")
			//	cout << specialize(Out, ref, que, GenoID, reg, tthr, sthr) << endl;
			//else if(Temp == "-E")
			//	cout << max_entropy(Out, ref, GenoID) << endl;
			else if(Temp == "-P")
				cout << RPS1(Out, GenoID, red, tax) << endl;
			else if (Temp == "-ST") {
				fa = read_in_fas(Fa);
				cout << tune_record(Out, GenoID, fa, red, que[GenoID], dis) << endl;
			}
		}
		return 0;
	}
	else {
		cout<<"########################################################################"<<endl
            <<"# An EM based translation initiation site (TIS) prediction program.    #"<<endl
            <<"# Developed by Longshu Yang.                                           #"<<endl
            <<"# Copyright: Center for Quantitative Biology, Peking University.       #"<<endl
            <<"# Version 1.0 released at 2018/02/02.                                  #"<<endl
            <<"########################################################################"<<endl
            <<"This is a brief introduction!"<<endl
            <<"1. Functional options:"<<endl
            <<"1)   -M: Motif training for TIS;"<<endl
            <<"         Necessary inputs:   -o, -f, -k, -m, -w and -r;"<<endl
            <<"2. Parameters for input:"<<endl
            <<" 1) Numerical inputs:"<<endl
            <<"     -k: <rank of possible initiation kmer>"<<endl
            <<"     -m: <models of TIS>"<<endl
            <<"     -r: <possible region of a TIS>"<<endl
            <<"     -w: <width of signal>"<<endl
            <<" 2) Input full directory and expansion name of a file:"<<endl
            <<"     -f: <readable fasta file>"<<endl
            <<"     -o: <output directory and file name>"<<endl;
		return 1;
	}
}
