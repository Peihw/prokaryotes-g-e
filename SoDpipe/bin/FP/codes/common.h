#ifndef COMMON_H
#define COMMON_H

#include "TypeDefBase.h"

//Data Structures
struct taxonomy {
	Str AssemID, Domain, Super, Gram, Phylum, Class, Order, Family, Genus, Species;
};

typedef std::map<Str, taxonomy> Ma_Str_X;

typedef std::map<Str, std::vector<Str> > Ma_Str_Ve_Str;

typedef std::vector<std::vector<std::vector<double> > > Ve_Ve_Ve_D;

typedef std::vector<std::vector<std::map<char, double> > > Ve_Ve_Ma_C_D;

typedef std::map<Str, std::pair<Str, Str> > Ma_Str_Pa_Str_Str;

struct motif {
	Ve_Ve_Ma_C_D rwm;//m weight matrices for w width of bases;
	Ve_Ve_D pj;//probability of signal initiates at jth base upstream of start codon for m signals;
	Ma_C_D pb;//probability of background base usage;
	Ma_C_D pc;//Counts of background bases;
	double loglike, last, p0, pd;
	Ve_D entropy;
	Ve_D trace;
	Ve_Str cons;
};

typedef std::map<Str, motif> Ma_Str_M;

struct site {
	int utr;
	Se_Str ref, phase;
};

typedef std::map<int, site> Ma_I_S;

struct record {
	Str Replicon, Strand, Name, Product, TISRef, Signal, Type, TSS, TTS, Feature, Seq;
	int start, end, shift, sig, opr;
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


//Commonly Used Functions

const Str Base = "ACGT";

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

struct align {
	int cont;
	double dis;
	double pcorr;
	double ccorr;
	double score;
	Ve_Pa_I_I trace;
};

typedef std::vector<align> Ve_A;

typedef std::vector<std::vector<align> > Ve_Ve_A;

struct dp {
	Pa_I_I trace;
	double score;
	int continuous;
	int best;
};

typedef std::vector<dp> Ve_P;

typedef std::vector<std::vector<dp> > Ve_Ve_P;

typedef std::map<Str, std::set<Str> > Ma_Str_Se_Str;

struct CmpByValue {
	int operator()(const Pa_Str_I& lhs, const Pa_Str_I& rhs) {
		return lhs.second > rhs.second;
	}
};

Ve_Str string_parse(con_Str Line, con_Str Key);

De_I subscript(long long int i, Ve_I limit);

//Data Loading Functions

Ma_Str_Str read_in_fa(Str In);

Ma_Str_Str read_in_cog(Str In);

Ma_Str_X read_in_taxonomy(Str In);

Ma_Str_Pa_Str_Str read_in_rps1(Str In);

Ma_C_D background_probability(Str In);

Ma_Str_M read_in_motif(Str In);

Ma_Str_Ve_Str read_in_distance(Str In);

rec read_in_record(Str In, Str Geno);

#endif

