#pragma once
#ifndef COMMON_H

#define COMMON_H

#include "TypeDefBase.h"


template <typename T>
void convertFromNumber(std::string &s, T &value) {
	std::ostringstream ss(s);
	ss << value;
	s = ss.str();
}

template <typename S>
void convertFromString(S &value, const std::string &s) {
	std::stringstream ss(s);
	ss >> value;
}

//Data Structure

typedef std::map<Str, std::vector<Str> > Ma_Str_Ve_Str;

struct site {
	int utr;
	Se_Str ref, phase;
};

typedef std::map<int, site> Ma_I_S;

struct record {
	Str Replicon, Strand, Name, Product, TISRef, Signal, Seq, Type, TSS, TTS, Feature;
	int start, end, shift, sig, opr, tis;
	Ma_I_S tss, tts;
	Ma_I_I ldl;
	double prob, rate;
};

typedef std::map<Str, record> Ma_Str_R;

struct rec {
	Ma_Str_R record;
	Ve_Str order;
	Str Genome;
};

struct feature {
	Str Feature, Class, Unit, Seq_type, Chromosome, GenoID, Strand, ProdID, Nr, RelID, Product, Symbol, GeneID, Synonym, Attributes;
	int start, end, interval, product;
};

typedef std::map<Str, feature> Ma_Str_F;

struct ftt {
	Str AssemID;
	Se_Str Chr;
	Ma_Str_F S2F;
	Ve_Str order;
};

/*struct med {
	Str Strand;
	int sta, end;
};

typedef std::map<int, med> Ma_I_M;*/

struct taxonomy {
	Str Domain, Phylum, Genus, Species, AssemID, Sub;
};

typedef std::map<Str, taxonomy> Ma_Str_T;

//Commonly Used Functions

Ve_Str string_parse(con_Str Line, con_Str Key);

De_I subscript(long long int i, Ve_I limit);

char complement(char S, bool NUC);

double gc_content(Ma_Str_Str fna);

Str antisense(Str Forward, bool NUC);

Str transcribe(Str DNA);

Str translate(Str Cds, Ma_Str_Str cod);

Str transverse(Str Forward);

Str cleave(Str Seq, Str Strand, int start, int end, bool NUC);

//Data Loading Functions

Ma_Str_Str read_in_fna(Str In);

Ma_Str_Str read_in_fa(Str In, int num);

rec read_in_record(Str In , Str Geno);

ftt read_in_ftt(Str In);

Ma_Str_Ve_Str read_in_keywords(Str In);

Ma_Str_Str read_in_codon(Str In);

Ma_Str_Str read_in_label(Str In);

Ma_Str_Str read_in_rfa(Str In, Str Geno);

Se_Str read_in_list(Str In);

Ma_Str_T read_in_taxonomy(Str In);

#endif // !COMMON_H
