#include "rps.h"
using namespace std;


Str rps(Str Out, Str GenoID, rec red, Ma_Str_X tax, Ma_Str_Str cog) {
	Ofs out(Out.data(), ios_base::app);
	bool FOUND = false;
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		Ve_Str tmp = string_parse(red.record[Synonym].Product, " ");
		bool R, P, S, N, B, E, D, A, NB;
		R = P = S = N = B = E = D = A = NB = false;
		for (int j = 0; j < tmp.size(); j++) {
			if (tmp[j] == "ribosomal" || tmp[j] == "Ribosomal" || tmp[j] == "Ribosome" || tmp[j] == "ribosome")
				R = true;
			if (tmp[j] == "protein" || tmp[j] == "Protein")
				P = true;
			if (tmp[j] == "S1" || tmp[j] == "S1p" || tmp[j] == "S1P")
				S = true;
			if (tmp[j] == "RNA" || tmp[j] == "rna")
				N = true;
			if (tmp[j] == "binding" || tmp[j] == "Binding")
				B = true;
			if (tmp[j] == "RNA-binding")
				NB = true;
			if (tmp[j] == "4-hydroxy-3-methylbut-2-enyl")
				E = true;
			if (tmp[j] == "diphosphate")
				D = true;
			if (tmp[j] == "reductase")
				A = true;
		}
		if ((R & P & S) || (((N && B) || NB) && S) || (cog[Synonym] == "COG0539")) {
			out << GenoID << "\t" << Synonym << "\t" << red.record[Synonym].Product << "\t" << cog[Synonym] << "\t" << tax[GenoID].Phylum << endl;
			FOUND = true;
			break;
		}
		else if (E && D && A && (tax[GenoID].Phylum == "Fusobacteria" || tax[GenoID].Phylum == "Thermotogae" || tax[GenoID].Phylum == "Firmicutes")) {
			out << GenoID << "\t" << Synonym << "\t" << red.record[Synonym].Product << "\t" << cog[Synonym] << "\t" << tax[GenoID].Phylum << endl;
			FOUND = true;
			break;
		}
	}
	if (!FOUND)
		out << GenoID << "\tNA\tNA\tNA\t" << tax[GenoID].Phylum << endl;
	out.close();
	return "RPS1 of " + GenoID + " have been searched!";
}