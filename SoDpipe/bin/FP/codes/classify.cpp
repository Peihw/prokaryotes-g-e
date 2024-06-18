#include "classify.h"
using namespace std;


Str dist(Str Out, Ma_Str_M ref, Ma_Str_M que, Str GenoID, Ma_Str_X taxo, double tthr, double sthr, Ma_Str_Pa_Str_Str check) {
	bool IN = true;
	Ifs in(Out.data());
	if (!in)
		IN = false;
	in.close();
	Ofs out(Out.data(), ios_base::app);
	bool RBS = false;
	if (!IN) {
		out << "Genome Accession\tGram Test\tPhylum\tClass\tOrder\tFamily\tGenus\tGenome\tRPS1 Check\tMotif\tConsensus\tEntropy\tType"; // \tDistance to Canonical RBS";
		for (Ma_Str_M::const_iterator it = ref.begin(); it != ref.end(); it++) {
			out << "\t Distance to " << it->first;
			//if (it->first == "Canonical Leaderless")
			//	out << "\tDistance to Canonical RBS";
		}
		out << endl;
	}
	for (int i = 0; i < que[GenoID].rwm.size(); i++) {
		out << GenoID << "\t" << taxo[GenoID].Gram << "\t" << taxo[GenoID].Phylum << "\t" << taxo[GenoID].Class << "\t" << taxo[GenoID].Order
			<< "\t" << taxo[GenoID].Family << "\t" << taxo[GenoID].Genus << "\t" << taxo[GenoID].Species;
		if (check[GenoID].first.empty())
			out << "\tNA";
		else
			out << "\t" << check[GenoID].first;
		out << "\tSignal" << i + 1 << "\t" << que[GenoID].cons[i] << "\t" << que[GenoID].entropy[i];
		Ma_Str_D score;
		double dis = 0;
		Str Type;
		double tmp = 1000;
		for (Ma_Str_M::const_iterator iter = ref.begin(); iter != ref.end(); iter++) {
			for (int j = 0; j < iter->second.rwm.size(); j++) {
				int len = 6;
				if (iter->first == "TANNTTTG-like")
					tmp = TA_distance(iter->second.rwm[j], iter->second.pb, que[GenoID].rwm[i], que[GenoID].pb, 8);
				else if (iter->first == "SD-like")// || iter->first == "Non-canonical RBS recognized by RPS1")
					tmp = TA_distance(iter->second.rwm[j], iter->second.pb, que[GenoID].rwm[i], que[GenoID].pb, 5);
				else
					tmp = TA_distance(iter->second.rwm[j], iter->second.pb, que[GenoID].rwm[i], que[GenoID].pb, len);
				score[iter->first] = tmp;
				if (iter->first == "RPS1 RBS") {
					if (!check[GenoID].first.empty() && check[GenoID].first != "NA") {
						if (tmp > dis && tmp > tthr) {
							dis = tmp;
							Type = iter->first;
						}
					}
				}
				else if (iter->first == "TANNTTTG-like") {
					if (tmp > dis && tmp > 0.8 && taxo[GenoID].Domain == "Bacteria") {
						dis = tmp;
						Type = iter->first;
					}
				}
				/*else if (iter->first == "Non-canonical RBS in Halophilic Archaea") {
					if (tmp > dis && tmp > 0.8) {
						dis = tmp;
						Type = iter->first;
					}
				}*/
				else {
					if (tmp > dis && tmp > tthr) {
						dis = tmp;
						Type = iter->first;
					}
				}
			}
		}
		if (score["TANNTTTG-like"] > 0.8 && taxo[GenoID].Domain == "Bacteria")
			Type = "TANNTTTG-like";
		if (!Type.empty())
			out << "\t" << Type;
		else
			out << "\tAtypical";
		//out << "\t" << score["Canonical RBS"];
		for (Ma_Str_D::const_iterator ite = score.begin(); ite != score.end(); ite++)
			out << "\t" << ite->second;
		out << endl;
	}
	out.close();
	return "The distances to all six standard signals have been calculated!";
}

double SD_distance(Str SD, Ve_Ma_C_D que, Ma_C_D pb, int len) {
	double result = 1000;
	for (int i = 0; i < que.size() - len + 1; i++) {
		for (int j = 0; j < SD.length() - len + 1; j++) {
			double score = 0;
			for (int k = 0; k < len; k++)
				score -= log2(que[i + k][SD[j + k]]);
			if (score < result)
				result = score;
			//cout << i << "\t" << SD.substr(j, len) << "\t" << score << "\t" << result << endl;
		}
	}
	return result;
}

double TA_distance(Ve_Ma_C_D ref, Ma_C_D rpb, Ve_Ma_C_D que, Ma_C_D qpb, int len) {
	double result = -1;
	bool PURE = false;
	for (int m = 0; m < qpb.size(); m++) {
		if (qpb[m] > 0.4)
			PURE = true;
	}
	for (int i = 0; i < ref.size() - len + 1; i++) {
		for (int j = 0; j < que.size() - len + 1; j++) {
			double score = 0;
			double scr = 0;
			for (int k = 0; k < len; k++) {
				double dis, moda, modb, disr, modra, modrb;
				dis = moda = modb = disr = modra = modrb = 0;
				for (int l = 0; l < Base.length(); l++) {
					double a, b;
					//if (que[j + k][Base[l]] / qpb[Base[l]] > 1.0e-10)
						//a = log2(que[j + k][Base[l]] / qpb[Base[l]]) * que[j + k][Base[l]];
					//else
						//a = 0;
					//if (ref[i + k][Base[l]] / rpb[Base[l]] > 1.0e-10)
						//b = log2(ref[i + k][Base[l]] / rpb[Base[l]]) * ref[i + k][Base[l]];
					//else
						//b = 0;
                    a = que[j + k][Base[l]] / qpb[Base[l]];
                    b = ref[i + k][Base[l]] / rpb[Base[l]];
					if (PURE) {
						disr += que[j + k][Base[l]] * ref[i + k][Base[l]];
						modra += que[j + k][Base[l]] * que[j + k][Base[l]];
						modrb += ref[i + k][Base[l]] * ref[i + k][Base[l]];
					}
					dis += a * b;
					moda += a * a;
					modb += b * b;
				}
				//dis = sqrt(dis);
				moda = sqrt(moda);
				modb = sqrt(modb);
				dis = dis / moda / modb;
				score += dis;
				if (PURE) {
					modra = sqrt(modra);
					modrb = sqrt(modrb);
					disr = disr / modra / modrb;
					scr += disr;
					if (scr > score)
						score = scr;
				}
			}
			if (score > result)
				result = score;
		}
	}
	result /= len;
	return result;
}
