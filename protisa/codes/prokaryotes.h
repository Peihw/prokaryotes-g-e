#pragma once
#ifndef PROKARYOTES_H
#define PROKARYOTES_H

#include "common.h"

Ma_Str_Str identify(Ve_Str id, Str Key);

Str correlate(Str Out, Ve_Str fna, Ve_Str ftt);

Str prokaryotes(Str Fna, Str Ftt, Str Out, Ma_Str_Str fna, Ma_Str_K keg, ftt ano);

Str nc_fasta(Str Out, Ma_Str_Str fna);

Str ptt_output(Str Out, Str NC, ftt ano);

Str tritisa_record(Str Out, Str Chr, Ma_I_M tisa, ftt ano, Ma_Str_Str trans);

Str blast2cog(Str Out, Str Geno, ftt ano, Ma_Str_L bla, Ma_Str_Str csv);

//Str tss_map(Str Out, Str Geno, rec red, Ma_Str_Str fna, Ma_I_E tss);

Str reannotate(Str Out, Str Geno, Ma_Str_Str fna, ftt ano);

Str exp_map(Str Out, Str Geno, rec red, Ma_Str_Str fna, Ve_E ep);

//Str tu_map(Str Out, Str Geno, rec red, Ma_Str_E trs, Ma_Str_T tu);

Str transcript_map(Str Out, Str Geno, Ma_Str_T tu, rec red);

Str operon(Str Out, Str Geno, Ve_Str genes, Ve_Pa_Str_I opr, Ma_Str_Str exp);

#endif