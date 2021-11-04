#pragma once
#ifndef COMMON_H
#define COMMON_H

#include "TypeDefBase.h"

//Constants

con_Str DNA = "ACGT";

con_Str RNA = "ACGU";

//Data structure

typedef std::map<Str, std::vector<Str> > Ma_Str_Ve_Str;

struct feature {
	Str Feature, Class, Unit, Seq_type, Chromosome, GenoID, Strand, ProdID, Nr, RelID, Product, Symbol, GeneID, Synonym, Attributes;
	int start, end, interval, product;
};

typedef std::map<Str, feature> Ma_Str_F;

struct ftt {
	Str AssemID;
	Ma_Str_F S2F;
	Ve_Str order;
	Se_Str Chr;
	//Ma_Str_Ve_Str forward;
	//Ma_Str_Ve_Str reverse;
};

struct kegg {
	Str Domain, Kingdom, Phylum, Genus;
};

typedef std::map<Str, kegg> Ma_Str_K;

struct med {
	Str Strand;
	int p5, p3;
};

typedef std::map<int, med> Ma_I_M;

struct blast{
	Str Que, Ref;
	int len, mis, gap, qst, qed, rst, red;
	double id, eval, bit;
};

typedef std::map<Str, blast> Ma_Str_L;

struct site {
	int utr;
	Se_Str ref, phase;
};

typedef std::map<int, site> Ma_I_S;

struct record {
	Str Replicon, Strand, Name, Product, TISRef, Feature;
	int start, end, shift;
	Ma_I_S tss, tts;
};

typedef std::map<Str, record> Ma_Str_R;

struct rec {
	Ma_Str_R record;
	Ve_Str order;
	Str Genome;
	Ma_Str_Ve_Str forward;
	Ma_Str_Ve_Str reverse;
};

struct ex {
	Str Strand, Replicon;
	Se_I tss;
	int tts;
	Ma_I_S ts, tt;
	Se_Str ref, phase;
};

typedef std::map<int, ex> Ma_I_E;

typedef std::map<Str, ex> Ma_Str_E;

typedef std::vector<ex> Ve_E;

struct transcript {
	Str TSS, TTS, ID;
};

typedef std::map<Str, transcript> Ma_Str_T;

typedef std::map<Str, std::pair<Str, Str> > Ma_Str_Pa_Str_Str;

typedef std::map<Str, std::set<Str> > Ma_Str_Se_Str;

//Commonly used functions

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

Ve_Str string_parse(con_Str Line, con_Str Key);

char complement(char S, bool NUC);

Str antisense(Str Forward, bool NUC);

Str transcribe(Str DNA);

Pa_I_I locate(Ve_Str order, int loc, Ma_Str_R red, Str Strand, int len, bool UTR);//UTR is true for 5'UTR;

Str cleave(Str Seq, Str Strand, int start, int end, bool NUC);

//Data loading functions

Ma_Str_Str read_in_codon(Str In);

Ma_Str_Str read_in_fna(Str In);

Ma_Str_Str read_in_fa(Str In);

ftt read_in_ftt(Str In);

Ma_Str_K read_in_kegg(Str In);

Ve_Str read_in_column(Str In);

Ma_I_M read_in_med(Str In);

Ma_Str_Str read_in_trans(Str In);

Ma_Str_Str read_in_cog(Str In);

Ma_Str_L read_in_blast(Str In);

Ma_Str_Str read_in_cog_csv(Str In);

rec read_in_record(Str In, Str Geno);

//Ma_I_E read_in_exp_tss(Str In);

Ve_E read_in_exp(Str In);

Ma_Str_T read_in_tu(Str In);

//Ma_Str_E read_in_trans_record(Str In);

Ma_Str_Str read_in_exftt(Str In);

Ve_Str read_in_tu_genes(Str In);

Ve_Pa_Str_I read_in_opr(Str In);

#endif
