#ifndef SCORING_H
#define SCORING_H

#include "common.h"

struct scan {
	int start;
	double pb, sc;
	Str Type, Seq, Name;
};

typedef std::map <Str, scan> Ma_Str_S;

Str dist(Str Out, Ma_Str_M ref, Ma_Str_M que, Str GenoID, int len, double tthr, double sthr, double pthr, Ma_Str_Pa_Str_Str check);

double SD_distance(Str SD, Ve_Ma_C_D que, Ma_C_D pb, int len, bool EXC);

double TA_distance(Ve_Ma_C_D ref, Ma_C_D rpb, Ve_Ma_C_D que, Ma_C_D qpb, int len, bool EXC);

Str tis_record(Str Out, Str Que, Ma_Str_Str fa, rec red, motif mot, Ma_Str_Ve_Str dis, int len);

Ma_Str_S prob(Ma_Str_Str fa, motif mot, Ve_Str dis);

Str ortholog_statistics(Str Out, Str Ori, Str Que, Ma_Str_R ori, Ma_Str_R que, Ma_Str_Ma_Str_O orth);

Str signal_statistics(Str Out, Str Ori, Ma_Str_R ori);

Str motif2fasta(Str Out, Str GenoID, motif mot, int depth);

//Str specialize(Str Out, motif ref, motif que, Str GenoID, int len, double tthr, double sthr);

Str max_entropy(Str Out, motif ref, Str GenoID);

Str RPS1(Str Out, Str GenoID, rec red, Ma_Str_X tax);

#endif
