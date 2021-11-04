#ifndef EM_H
#define EM_H

#include "common.h"

struct dp {
	Pa_I_I trace;
	double score;
	int continuous;
};

typedef std::vector<dp> Ve_P;

typedef std::vector<std::vector<dp> > Ve_Ve_P;

struct align {
	double score;
	Ve_Pa_I_I trace;
};

typedef std::vector<align> Ve_A;

typedef std::vector<std::vector<align> > Ve_Ve_A;

typedef std::map<Str, std::set<Str> > Ma_Str_Se_Str;

Ve_Ma_C_D pwm_assign(Str Seq, Ma_C_D pb, double maj, double min);

align dynamic_program(Str Ref, Str Que);

Ve_Str alignment(Ma_Str_Str fa, int wid, int thr);

Ma_C_D background(Ma_Str_Str fa);

typedef std::vector< motif > Ve_M;

motif EM(motif test, Ma_Str_Str fa, int mod, int reg, double eps, int step = 500);

Str TIS_motif(Str Out, Str GenoID, Ma_Str_Str fa, double thr, int mod, int wid, double maj, double min, double eps, int trial, Ma_C_D pb);

struct CmpByValue {
	int operator()(const Pa_Str_I& lhs, const Pa_Str_I& rhs) {
		return lhs.second > rhs.second;
	}
};

Ma_Str_D background_score(Ma_Str_Str fa, Ma_C_D pb);


#endif // !EM_H


#pragma once
