#ifndef COMMON_H
#define COMMON_H

#include "TypeDefBase.h"

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

typedef std::map<Str, std::vector<Str> > Ma_Str_Ve_Str;

typedef std::vector<std::vector<std::vector<double> > > Ve_Ve_Ve_D;

typedef std::vector<std::vector<std::map<char, double> > > Ve_Ve_Ma_C_D;

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

struct ortholog {
	Str OID, QID;
	double o2qi, q2oi, o2qe, q2oe;
	int olen, qlen;
};

typedef std::map<Str, ortholog> Ma_Str_O;

typedef std::map<Str, std::map <Str, ortholog> > Ma_Str_Ma_Str_O;

struct blast {
	Str SID;
	Str QID;
	double id;
	int len;
	int mis;
	int gap;
	int qstart;
	int qend;
	int sstart;
	int send;
	double evalue;
	double bit;
};

typedef std::map<Str, std::vector<blast> > Ma_Str_Ve_L;

struct site {
	int utr;
	Se_Str ref, phase;
};

typedef std::map<int, site> Ma_I_S;

struct record {
	Str Replicon, Strand, Name, Product, TISRef, Signal, Type, TSS, TTS, Feature;
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

struct tss_exp {
	Str Replicon, Strand, Synonym;
	int start;
	Se_Str tssref;
	Ve_Pa_I_I tss;//Absolute, Relative.
};

typedef std::map<Str, tss_exp> Ma_Str_T;

struct taxonomy {
	Str Domain, Phylum, Genus, Species, AssemID, Sub;
};

typedef std::map<Str, taxonomy> Ma_Str_X;

typedef std::map<Str, std::pair<Str, Str> > Ma_Str_Pa_Str_Str;

ftt read_in_ftt(Str In);

Ve_Str string_parse(con_Str Line, con_Str Key);

De_I subscript(long long int i, Ve_I limit);

//Data Loading Functions

Ma_Str_Str read_in_fa(Str In);

Ma_Str_Str read_in_fas(Str In);

Ma_Str_M read_in_motif(Str In);

Ma_Str_Ve_Str read_in_distance(Str In);

Ma_Str_Ma_Str_O read_in_ortholog(Str In);

Ma_C_D background_probability(Str In);

rec read_in_record(Str In, Str Geno);

Ma_Str_Str read_in_cog(Str In);

Ma_Str_T read_in_exp_tss(Str In);

Ma_Str_Pa_Str_Str read_in_rps1(Str In);

Ma_Str_X read_in_taxonomy(Str In);

#endif

