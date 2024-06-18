#pragma once
#ifndef RECORD_H
#define RECORD_H

#include "common.h"

struct scan {
	int start;
	double pb, sc;
	Str Type, Seq, Name;
};

typedef std::map <Str, scan> Ma_Str_S;

Ma_Str_S prob(Ma_Str_Str fa, motif mot, Ve_Str dis);

Str tis_record(Str Out, Str GenoID, Ma_Str_Str fa, rec red, motif mot, Ma_Str_Ve_Str dis, int reg);

#endif