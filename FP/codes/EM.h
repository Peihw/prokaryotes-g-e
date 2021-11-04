#pragma once
#pragma once
#ifndef KMER_H
#define KMER_H

#include "common.h"

struct candidate {
	Str Consensus;
	Ve_Ma_C_D pwm;
};

struct comparison {
	align alg;
	Str Best;
	double sim;
	int slide;
	int maxlen;
	int start;
	int end;
};

typedef std::map <Str, std::map<Str, comparison> > Ma_Str_Ma_Str_P;

typedef std::vector<candidate> Ve_E;

typedef std::map <Str, std::map<Str, align> > Ma_Str_Ma_Str_A;

typedef std::map<Str, align> Ma_Str_A;

typedef std::vector<std::set<Str> > Ve_Se_Str;

typedef std::map<Str, std::vector<std::map<char, double> > > Ma_Str_Ve_Ma_C_D;

typedef std::map<Str, std::map<Str, double> > Ma_Str_Ma_Str_D;

typedef std::map<Str, std::map<char, double> > Ma_Str_Ma_C_D;

typedef std::vector<std::map<Str, int> > Ve_Ma_Str_I;

typedef std::vector< motif > Ve_M;

typedef std::map<Str, std::vector<double> > Ma_Str_Ve_D;

struct category {
	Se_Str cate;
	Pa_D_D best;
};

struct kmer {
	Ma_Str_A matr;
	Ve_D pos;
	int count;
	Ma_Str_D dis;
	Ma_Str_D corr;
};

typedef std::map<Str, kmer> Ma_Str_K;

struct difference {
	Str Cent;
	int pos;
	int size;
	double dis;
	double cont;
	double corr;
	Se_Str mem;
};

typedef std::map<Str, difference> Ma_Str_F;

bool check_g(Str In, Str Type);

bool check_promoter(Str In, Str Type);

align dynamic_program(Str Ref, Str Que);

double composition_correlation(Str Ref, Str Que);

Ve_Str alignment(Ma_Str_Str fa, int wid, int lim, int mod, Ma_C_D pb);

Ma_C_D background(Ma_Str_Str fa);

Ma_Str_D background_score(Ma_Str_Str fa, Ma_C_D pb);

Ve_Ma_C_D pwm_assign(Str Seq, Ma_C_D pb, double maj, double min);

motif EM(motif test, Ma_Str_Str fa, int mod, int reg, double eps, int step = 500);

Str TIS_motif(Str Out, Str GenoID, Ma_Str_Str fa, double thr, int mod, int wid, double maj, double min, double eps, int trial, Ma_C_D pb);

#endif