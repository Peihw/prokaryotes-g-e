#pragma once
#ifndef EXTRACT_H

#define EXTRACT_H

#include "common.h"

Str gene(Str Out, Str Geno, Ma_Str_Ve_Str keyword, ftt feat, Ma_Str_Str fa, Ma_Str_Str cod, Ma_Str_Str label, int reg);

Str regions(Str Out, Str Geno, Str Opt, rec red, Ma_Str_Str fna, Ma_Str_Str cod, int up, int down);

Str ips(Str Out, Str Geno, Ma_Str_Str fa, Ma_Str_Str fna);

Str nonredundant(Str Out, Ma_Str_Str rfa);

Str rRNA(Str Out, Ma_Str_Str fa, Ma_Str_T tax, int len);

#endif