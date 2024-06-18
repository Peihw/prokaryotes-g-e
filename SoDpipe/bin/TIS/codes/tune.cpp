#include "tune.h"
using namespace std;

/*Str tune_nonredundant(Str Kmer, int hit, Ma_Str_I kmer) {
	Str Res = "Candidate";
	Ma_C_D kvec;
	double n = Kmer.length();
	for (int i = 0; i < Kmer.length(); i++)
		kvec[Kmer[i]] += 1 / n;
	for (Ma_Str_I::iterator iter = kmer.begin(); iter != kmer.end(); iter++) {
		Ma_C_D cvec;
		for (int i = 0; i < iter->first.length(); i++)
			cvec[iter->first[i]] += 1 / n;
		double prod, normk, normc;
		prod = normk = normc = 0;
		for (int j = 0; j < Base.length(); j++) {
			prod += kvec[Base[j]] * cvec[Base[j]];
			normk += kvec[Base[j]] * kvec[Base[j]];
			normc += cvec[Base[j]] * cvec[Base[j]];
		}
		normk = sqrt(normk);
		normc = sqrt(normc);
		double corr = prod / normk / normc;
		if (corr >= 0.75 && iter->second > hit) {
			Res = "Abandoned";
			break;
		}
	}
	return Res;
}*/

Str tune_nonredundant(Str Kmer, bool AG, bool CG, int hit, Ma_Str_I kmer) {
	Str Res = "Candidate";
	int ag, cg, g;
	ag = cg = g = 0;
	//AG or CG rich check!
	for (int i = 0; i < Kmer.length(); i++) {
		if (Kmer[i] == 'A') {
			ag++;
			if (cg >= 5 && g > 0) {
				if (CG)
					Res = "Abandoned";
				else
					Res = "CG rich";
			}
			cg = 0;
		}
		else if (Kmer[i] == 'C') {
			cg++;
			if (ag >= 5 && g > 0) {
				if (AG)
					Res = "Abandoned";
				else
					Res = "AG rich";
			}
			ag = 0;
		}
		else if (Kmer[i] == 'G') {
			ag++;
			g++;
			cg++;
		}
		else {
			if (cg >= 5 && g > 0) {
				if (CG)
					Res = "Abandoned";
				else
					Res = "CG rich";
			}
			if (ag >= 5 && g > 0) {
				if (AG)
					Res = "Abandoned";
				else
					Res = "AG rich";
			}
			ag = g = cg = 0;
		}
	}
	if (cg >= 5 && g > 0) {
		if (CG)
			Res = "Abandoned";
		else
			Res = "CG rich";
	}
	if (ag >= 5 && g > 0) {
		if (AG)
			Res = "Abandoned";
		else
			Res = "AG rich";
	}
	//Checking for SD covered or covered by other candidates!
	for (int j = 0; j < Kmer.length() - 4; j++) {
		Str SD = "AAGGAGGTGA";
		Str Sub = Kmer.substr(j, 5);
		for (Ma_Str_I::const_iterator iter = kmer.begin(); iter != kmer.end(); iter++) {
			if (iter->first.find(Sub) != Str::npos && iter->second > hit)
				Res = "Abandoned";
		}
		if (SD.find(Sub) != Str::npos && AG)
			Res = "Abandoned";
        if (Sub == "TTATT" && j < 2)
            Res = "Abandoned";
        if (Sub == "TTTTG" && j < 2)
            Res = "Abandoned";
	}
	//Checking for TATA box shifting!
	if (Kmer[0] == 'A' && Kmer[1] == 'T' && Kmer[4] == 'T')
		Res = "Abandoned";
	int tat = Kmer.find("TAT");
	while (tat != Str::npos) {
		if (tat >= Kmer.length() - 5 && tat < Kmer.length() - 3)
			Res = "Abandoned";
		tat = Kmer.find("TAT", tat + 1);
	}
    Str TTTG = Kmer.substr(0, 4);
    if(TTTG == "TTTG")
        Res = "Abandoned";
	return Res;
}

Ve_Str kmer_tune(Ma_Str_Str fa, int wid, int thr, Str Kmer) {
	Ve_Str res;
	Ma_Str_I count, cand;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		//cout << iter->first << "\t" << iter->second << endl;
		for (int i = 0; i < iter->second.length() - wid + 1; i++) {
			//cout << iter->second.length() << "\t" << i;
			Str Kmer = iter->second.substr(i, wid);
			//cout << "\t" << Kmer << endl;
			count[Kmer] ++;
		}
	}
	cout << "All K-mer types: " << count.size() << endl;
	Ve_Pa_Str_I rank(count.begin(), count.end());
	sort(rank.begin(), rank.end(), Cmp());
	int chosen = 0;
	bool AG = false;
	bool CG = false;
	Ofs out(Kmer.data());
	for (int i = 0; i < rank.size(); i++) {
		Str Tmp = tune_nonredundant(rank[i].first, AG, CG, rank[i].second, cand);
		out << rank[i].first << "\t" << rank[i].second << "\t" << Tmp << endl;
		if (!i) {
			cand[rank[i].first] = rank[i].second;
			res.push_back(rank[i].first);
			chosen++;
		}
		else if (Tmp != "Abandoned") {
			cand[rank[i].first] = rank[i].second;
			res.push_back(rank[i].first);
			chosen++;
		}
		if (Tmp == "AG rich")
			AG = true;
		if (Tmp == "CG rich")
			CG = true;
		if (chosen == thr)
			break;
	}
	out.close();
	return res;
}

Ve_Ma_C_D pwm_initialize(Str Seq, Ma_C_D pb, double maj, double min) {
	Ve_Ma_C_D res;
	for (int i = 0; i < Seq.size(); i++) {
		Ma_C_D temp;
		for (int j = 0; j < Base.size(); j++) {
			if (Base[j] == Seq[i])
				temp[Seq[i]] = maj / pb[Seq[i]];
			else
				temp[Base[j]] = min / pb[Seq[i]];
		}
		res.push_back(temp);
	}
	return res;
}

Ma_Str_D background_likelihood(Ma_Str_Str fa, Ma_C_D pb) {
	Ma_Str_D res;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		double b = 0;
		for (int i = 0; i < iter->second.size(); i++)
			b += log2(pb[iter->second[i]]);
		res[iter->first] = b;
	}
	return res;
}

Ma_C_D base_distribution(Ma_Str_Str fa){
	Ma_C_D res;
	int sum = 0;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		for (int i = 0; i < iter->second.length(); i++) {
			if (Base.find(iter->second[i] != Str::npos)) {
				res[iter->second[i]] ++;
				sum++;
			}
		}
	}
	return res;
}

motif EM_tune(motif test, Ma_Str_Str fa, int mod, int reg, double eps, int step) {
	motif res = test;
	for (int n = 0; n < step; n++) {
		Ma_C_D B;
		for (int v = 0; v < Base.length(); v++)
			B[Base[v]] = 0;
		Ve_Ve_Ma_C_D D;
		for (int m = 0; m < mod; m++) {
			Ve_Ma_C_D tem;
			for (int i = 0; i < res.rwm[m].size(); i++) {
				tem.push_back(B);
			}
			D.push_back(tem);
		}
		Ve_Ve_D C(mod, Ve_D(reg - res.rwm[0].size() + 1, 0));
		Ve_Ve_Ve_D lam(mod, Ve_Ve_D(fa.size(), Ve_D(reg - res.rwm[0].size() + 1, 0)));
		Ma_Str_D bs = background_likelihood(fa, res.pb);
		int k = 0;
		double p0 = 0;
		for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
			double sumk = res.p0;
			for (int j = 0; j < reg - res.rwm[0].size() + 1; j++) {
				for (int m = 0; m < mod; m++) {
					double prod = 1;
					for (int i = 0; i < res.rwm[m].size(); i++) {
						prod *= res.rwm[m][i][iter->second[i + j]];
					}
					lam[m][k][j] = prod * res.pj[m][j];
					sumk += lam[m][k][j];
				}
			}
			res.loglike += log2(sumk) + bs[iter->first];
			p0 += res.p0 / sumk;
			for (int j = 0; j < reg - res.rwm[0].size() + 1; j++) {
				for (int m = 0; m < mod; m++) {
					lam[m][k][j] /= sumk;
					C[m][j] += lam[m][k][j];
					for (int i = 0; i < res.rwm[m].size(); i++) {
						D[m][i][iter->second[i + j]] += lam[m][k][j];
						B[iter->second[i + j]] -= lam[m][k][j];
					}
				}
			}
			k++;
			//cout << "Background log likelihood:\t" << bs[iter->first] << "\tLog likelihood:\t" << res.loglike << endl;
		}
		if (abs(res.loglike - res.last) > eps * abs(res.loglike)) {
			res.trace.push_back(res.loglike);
			res.last = res.loglike;
			res.loglike = 0;
			double sumb = 0;
			if (n != 0 && res.loglike < res.last)
				cout << n << " step backward!" << endl;
			for (int v = 0; v < Base.length(); v++) {
				B[Base[v]] += res.pc[Base[v]];
				sumb += B[Base[v]];
			}
			//pb estimation
			for (int v = 0; v < Base.length(); v++) 
				res.pb[Base[v]] = B[Base[v]] / sumb;
			double sumc = 0;
			for (int m = 0; m < mod; m++) {
				for (int i = 0; i < res.rwm[m].size(); i++) {
					double sumd = 0;
					for (int v = 0; v < Base.length(); v++)
						sumd += D[m][i][Base[v]];
					//RWM estimation
					for (int v = 0; v < Base.length(); v++) {
						D[m][i][Base[v]] /= sumd;
						res.rwm[m][i][Base[v]] = D[m][i][Base[v]] / res.pb[Base[v]];
						D[m][i][Base[v]] = 0;
					}
				}
				for (int j = 0; j < C[m].size(); j++)
					sumc += C[m][j];
			}
			for (int m = 0; m < mod; m++) {
				//pj estimation
				for (int j = 0; j < C[m].size(); j++)
					res.pj[m][j] = C[m][j] / (sumc + p0);
			}
			//p0 estimation
			res.p0 = p0 / (sumc + p0);
		}
		else {
			res.trace.push_back(res.loglike);
			res.entropy.clear();
			for (int m = 0; m < mod; m++) {
				double entropy = 0;
				for (int i = 0; i < res.rwm[m].size(); i++) {
					for (int v = 0; v < Base.length(); v++)
						entropy += res.pb[Base[v]] * res.rwm[m][i][Base[v]] * log2(res.rwm[m][i][Base[v]]);
				}
				res.entropy.push_back(entropy);
			}
			break;
		}
	}
	return res;
}

Str tune(Str Out, Ma_Str_Str fa, int rank, int mod, int wid, int reg, double maj, double min, double eps, Ma_C_D pb) {
	Str Log = Out + ".log";
	Str Kmer = Out + ".kmer.log";
	Str Mod;
	convertFromNumber(Mod, mod);
	Ve_Str km = kmer_tune(fa, wid, rank, Kmer);
	cout << "Selected K-mer types: " << km.size() << endl;
	long long int max = 1;
	Ve_I limit;
	for (int m = 0; m < mod; m++) {
		limit.push_back(km.size() - mod + 1);
		max *= km.size() - mod + 1;
	}
	for (int k = 0; k < max; k++) {
		De_I sub = subscript(k, limit);
		bool Check = true;
		for (int m = 1; m < mod; m++) {
			cout << sub[m - 1] + m - 1 << "\t" << sub[m] + m << endl;
			if (sub[m - 1] + m - 1 >= sub[m] + m)
				Check = false;
		}
		if (Check) {
			Ofs log(Log.data(), ios_base::app);
			Ve_Str tem;
			motif res;
			res.last = -1.0e8;
			res.pb = pb;
			res.pc = base_distribution(fa);
			res.p0 = 1.0 / ((reg - wid + 1) * mod + 1);
			for (int m = 0; m < mod; m++) {
				log << "\t" << m << ":" << sub[m] << ":" << km[sub[m] + m];
				res.rwm.push_back(pwm_initialize(km[sub[m] + m], res.pb, maj, min));
				tem.push_back(km[sub[m] + m]);
			}
			Ve_D pj;
			pj.resize(reg - wid + 1, 1.0 / ((reg - wid + 1) * mod + 1));
			for (int m = 0; m < mod; m++)
				res.pj.push_back(pj);
			res.entropy.resize(mod, 0);
			//res.last = test.loglike = 0;
			res = EM_tune(res, fa, mod, reg, eps);
			Ofs out(Out.data(), ios_base::app);
			out << "Initial K-mer:" << endl;
			for (int m = 0; m < mod; m++)
				out << "\t" << tem[m];
			out << endl;
			out << "Features of Translation initiation site (TIS) signal:" << endl
				<< "Final likelihood score of this model:\t" << res.loglike << endl
				<< "Probability of background nucleotide bias:" << endl
				<< "\tA\tC\tG\tT" << endl;
			for (int r = 0; r < Base.length(); r++)
				out << "\t" << res.pb[Base[r]];
			out << endl << "Probabilities of TIS start position:" << endl
				<< "No signal:" << endl << "\t" << res.p0 << endl;
			for (int r = 0; r < res.pj.size(); r++) {
				out << "Signal" << r + 1 << ":" << endl;
				for (int q = 0; q < res.pj[r].size(); q++)
					out << "\t" << res.pj[r][q];
				out << endl;
			}
			for (int r = 0; r < res.rwm.size(); r++) {
				Str Motif;
				out << "Position weight matrices (PWM) of Motif" << r + 1 << ":" << endl;
				//<< "\tA\tC\tG\tT" << endl;
				for (int q = 0; q < res.rwm[r].size(); q++) {
					//out << q + 1;
					double last = 0;
					int base = 0;
					for (int s = 0; s < Base.length(); s++) {
						if (log2(res.rwm[r][q][Base[s]]) * res.rwm[r][q][Base[s]] * res.pb[Base[s]] > last) {
							last = log2(res.rwm[r][q][Base[s]]) * res.rwm[r][q][Base[s]] * res.pb[Base[s]];
							base = s;
						}
					}
					Motif += Base[base];
					for (int s = 0; s < Base.length(); s++)
						out << "\t" << res.rwm[r][q][Base[s]] * res.pb[Base[s]];
					out << endl;
				}
				out << "Consensus sequenece of Motif" << r + 1 << ":\t" << Motif << endl;
				log << "\t" << Motif << "\t" << res.entropy[r];
				out << "Relative entropy of Motif" << r + 1 << ":\t" << res.entropy[r] << endl;
			}
			log << "\t" << res.trace[res.trace.size() - 1];
			for (int j = 0; j < res.trace.size(); j++)
				log << "\t" << res.trace[j];
			log << endl;
			log.close();
			out.close();
		}
	}
	return "TIS signals of " + Mod + " types have been modeled!";
}

Str tune_record(Str Out, Str Que, Ma_Str_Str fa, rec red, motif mot, Ma_Str_Ve_Str dis) {
	Out += red.Genome + ".tis.rec.dat";
	Ofs out(Out.data());
	int utr = 0;
	out << "Genomic Accession\tStrand\tLocus\tFeature\tTSS\tTSS Type";
	for (int m = 0; m < mot.rwm.size(); m++)
		out << "\tSignal" << m << " Type\tSignal" << m << " score\tSignal" << m << " Start\tSignal" << m << " Sequence";
	out << "\tHighest Probability" << endl;
	int fin = red.order.size() - 1;
	Str Last = red.order[fin];
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		if (!red.record[Synonym].tss.empty()) {
			out << red.record[Synonym].Replicon << "\t" << red.record[Synonym].Strand << "\t" << Synonym << "\t" << red.record[Synonym].Feature;
			bool LDL = false;
			bool LDD = false;
			bool INT = false;
            bool DIS = false;
            bool PRX = false;
            bool EXT = false;
			for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
				if (iter == red.record[Synonym].tss.begin())
					out << "\t" << iter->second.utr;
				else
					out << ";" << iter->second.utr;
				if (iter->second.utr < -5 && iter->second.utr >= -500) {
					if (iter->second.utr < -30)
						LDD = true;
					else if (iter->second.utr < -20)
						DIS = true;
					else if (iter->second.utr < -10)
						PRX = true;
					else
						EXT = true;
				}
				//if (abs(iter->second.utr) <= 7)
				if (iter->second.utr <= 0 && iter->second.utr >= -5)
					LDL = true;
                if (iter->second.utr > 0 && iter->second.utr <= 5)
                    INT = true;
			}
			if (LDL)
				out << "\tLeaderless";
            else if (EXT)
                out << "\tExtended Leaderless";
            else if(PRX)
                out << "\tProximal Scanning Region";
            else if(DIS)
                out << "\tDistal Scanning Region";
			else if (LDD)
				out << "\tLeadered";
            else if(INT)
                out << "\tInternal";
			else
				out << "\tOther";
			Ve_D score;
			for (int m = 0; m < mot.rwm.size(); m++) {
				if (!fa[Synonym].empty()) {
					double lsc = 0;
					int lst = 0;
					double pb = 0;
                    double pos = 0;
					for (int j = 0; j < mot.pj[m].size(); j++) {
						double sc = 1;
                        //pos += mot.pj[m][j];
						for (int k = 0; k < mot.rwm[m].size(); k++)
							sc *= mot.rwm[m][k][fa[Synonym][j + k]] / mot.pb[fa[Synonym][j + k]];
						if (sc > lsc) {
							lsc = sc;
							lst = j;
						}         
						sc *= mot.pj[m][j];
						pb += sc;
					}
                    //pb *= pos;
					int w = mot.rwm[m].size();
					int d = mot.pj[m].size();
					int s = lst - w - d + 1;
					score.push_back(pb);
					//score.push_back(lsc);
					out << "\t" << dis[Que][m] << "\t" << pb << "\t" << s << "\t" << fa[Synonym].substr(lst, w);
					//out << "\t" << dis[Que][m] << "\t" << lsc << "\t" << s << "\t" << fa[Synonym].substr(lst, w);
				}
				else
					out << "\tNA\tNA\tNA\tNA";
			}
			double best = 0;
			int which = 0;
			for (int m = 0; m < score.size(); m++) {
				if (score[m] > best) {
					which = m;
					best = score[m];
				}
			}
			out << "\t" << dis[Que][which] << endl;
		}
		Last = Synonym;
	}
	out.close();
	return "All information about " + red.Genome + " have been recorded!";
}

