#include "record.h"
using namespace std;

Ma_Str_S prob(Ma_Str_Str fa, motif mot, Ve_Str dis) {
	Ma_Str_S result;
	for (Ma_Str_Str::iterator iter = fa.begin(); iter != fa.end(); iter++) {
		scan temp;
		temp.pb = temp.sc = 0;
		double sum, lsc;
		int lst = 0;
		sum = lsc = 0;
		for (int i = 0; i < mot.rwm.size(); i++) {
			double pb = 0;
			for (int j = 0; j < iter->second.length() - mot.rwm[i].size() +1; j++) {
				int pos = mot.pj[i].size() - iter->second.length() + mot.rwm[i].size() - 1;
				//cout << iter->first << endl << pos << endl;
                long double sc = 1;
				sc *= mot.pj[i][pos + j];
				//cout << sc << "\t" << i << "\t" << j << endl;
				for (int k = 0; k < mot.rwm[i].size(); k++)
					sc *= mot.rwm[i][k][iter->second[j + k]] / mot.pb[iter->second[j + k]];
				//cout << sc << "\t" << i << "\t" << j << endl;
				if (sc > lsc) {
					lsc = sc;
					lst = j;
					//cout << lsc << "\t" << lst << endl;
				}
				if (sc > 1.0e-10)
					pb += sc;
			}
			sum += pb;
			//cout << sum << "\t" << pb << endl;
			if (pb > temp.pb) {
				temp.pb = pb;
				temp.sc = lsc;
				int w = mot.rwm[i].size();
				int d = iter->second.length();
				temp.start = lst - d;
				//temp.start = lst - w - d + 1;
				temp.Seq = iter->second.substr(lst, w);
				temp.Type = dis[i];
				int a = i + 1;
				Str N;
				convertFromNumber(N, a);
				temp.Name = "Signal" + N;
			}
		}
		//cout << temp.pb << "\t" << temp.sc << "\t" << temp.Name << endl;
		sum += mot.p0;
		if (temp.pb > mot.p0)
			temp.pb /= sum;
		else {
			temp.pb = mot.p0 / sum;
			temp.sc = 0;
			temp.Seq = "NA";
			temp.start = 0;
			temp.Name = "None";
			temp.Type = "NA";
		}
		Ve_Str sub = string_parse(iter->first, "|");
		result[sub[2]] = temp;
		//cout << iter->first << endl;
	}
	return result;
}

Str tis_record(Str Out, Str GenoID, Ma_Str_Str fa, rec red, motif mot, Ma_Str_Ve_Str dis, int reg) {
	Out += GenoID + ".tis.rec.dat";
	Ofs out(Out.data());
	int utr = 0;
	Ma_Str_S sinfo = prob(fa, mot, dis[GenoID]);
	out << "Genomic Accession\tStrand\tLocus\tGene Name\tProduct\tGene Start\tGene End\tTIS Annotation Reference\tTIS shift (5\'-,3\'+)\tFeature"
		<< "\tExperimentally Verified TSS (Position|5\'UTR Length|Phase|Reference)\tExperimentally Verified TTS (Position|Phase|Reference)"
		<< "\tExperimentally Verified Leaderless TIS (Position|5\'End Length)\tOperon ID\tRole of TSSs\tRole of TTSs"
		<< "\tTranslation Initiation Signal\tTranslation Signal Type\tTranslation Signal Probability"
		<< "\tTranslation Signal Start\tTranslation Signal Sequence" << endl;
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		//cout << Synonym << "\t";
		//if (!fa[Synonym].empty() && fa[Synonym] != "Error!")
		out << red.record[Synonym].Replicon << "\t" << red.record[Synonym].Strand << "\t" << Synonym
			<< "\t" << red.record[Synonym].Name << "\t" << red.record[Synonym].Product
			<< "\t" << red.record[Synonym].start << "\t" << red.record[Synonym].end
			<< "\t" << red.record[Synonym].TISRef;
		if (red.record[Synonym].shift)
			out << "\t" << red.record[Synonym].shift;
		else
			out << "\tNA";
		out << "\t" << red.record[Synonym].Feature;
		if (!red.record[Synonym].tss.empty()) {
			for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
				if (iter == red.record[Synonym].tss.begin())
					out << "\t" << iter->first << "|" << iter->second.utr;
				else
					out << ";" << iter->first << "|" << iter->second.utr;
				for (Se_Str::const_iterator ite = iter->second.phase.begin(); ite != iter->second.phase.end(); ite++) {
					if (ite == iter->second.phase.begin())
						out << "|" << *ite;
					else
						out << "&" << *ite;
				}
				for (Se_Str::const_iterator ite = iter->second.ref.begin(); ite != iter->second.ref.end(); ite++) {
					if (ite == iter->second.ref.begin())
						out << "|" << *ite;
					else
						out << "&" << *ite;
				}
			}
		}
		else
			out << "\tNA";
		if (!red.record[Synonym].tts.empty()) {
			for (Ma_I_S::const_iterator iter = red.record[Synonym].tts.begin(); iter != red.record[Synonym].tts.end(); iter++) {
				if (iter == red.record[Synonym].tts.begin())
					out << "\t" << iter->first;
				else
					out << ";" << iter->first;
				for (Se_Str::const_iterator ite = iter->second.phase.begin(); ite != iter->second.phase.end(); ite++) {
					if (ite == iter->second.phase.begin())
						out << "|" << *ite;
					else
						out << "&" << *ite;
				}
				for (Se_Str::const_iterator ite = iter->second.ref.begin(); ite != iter->second.ref.end(); ite++) {
					if (ite == iter->second.ref.begin())
						out << "|" << *ite;
					else
						out << "&" << *ite;
				}
			}
		}
		else
			out << "\tNA";
		if (!red.record[Synonym].ldl.empty()) {
			for (Ma_I_I::const_iterator ite = red.record[Synonym].ldl.begin(); ite != red.record[Synonym].ldl.end(); ite++) {
				if (ite == red.record[Synonym].ldl.begin())
					out << "\t" << ite->first << "|" << ite->second;
				else
					out << ";" << ite->first << "|" << ite->second;
			}
		}
		else
			out << "\tNA";
		if (red.record[Synonym].opr)
			out << "\t" << red.record[Synonym].opr;
		else
			out << "\tNA";
		if (!red.record[Synonym].TSS.empty())
			out << "\t" << red.record[Synonym].TSS;
		else if (!red.record[Synonym].tss.empty())
			out << "\tVerified TSS";
		else
			out << "\tNA";
		if (!red.record[Synonym].TTS.empty())
			out << "\t" << red.record[Synonym].TTS;
		else
			out << "\tNA";
		if (!sinfo[Synonym].Name.empty()) {
			out << "\t" << sinfo[Synonym].Name << "\t" << sinfo[Synonym].Type << "\t" << sinfo[Synonym].pb
				<< "\t" << sinfo[Synonym].start << "\t" << sinfo[Synonym].Seq << endl;
		}
		else
			out << "\tNA\tNA\tNA\tNA\tNA" << endl;
	}
	out.close();
	return "All information about " + GenoID + " have been recorded!";
}
