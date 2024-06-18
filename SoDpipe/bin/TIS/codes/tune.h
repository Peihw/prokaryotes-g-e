#pragma once
#ifndef  TUNE_H

#define TUNE_H

#include "common.h"

struct Cmp {
	int operator()(const Pa_Str_I& lhs, const Pa_Str_I& rhs) {
		return lhs.second > rhs.second;
	}
};
Str tune_nonredundant(Str Kmer, bool AG, bool CG, int hit, Ma_Str_I kmer);

//bool tune_similarity(Str Kmer, Ve_Str kmer);

Ve_Str kmer_tune(Ma_Str_Str fa, int wid, int thr, Str Kmer);

Ve_Ma_C_D pwm_initialize(Str Seq, Ma_C_D pb, double maj, double min);

Ma_Str_D background_likelihood(Ma_Str_Str fa, Ma_C_D pb);

motif EM_tune(motif test, Ma_Str_Str fa, int mod, int reg, double eps, int step = 5000);

Str tune(Str Out, Ma_Str_Str fa, int rank, int mod, int wid, int reg, double maj, double min, double eps, Ma_C_D pb);

Str tune_record(Str Out, Str Que, Ma_Str_Str fa, rec red, motif mot, Ma_Str_Ve_Str dis);

#endif // ! TUNE_H
