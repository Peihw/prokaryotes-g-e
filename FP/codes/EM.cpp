#include "EM.h"
using namespace std;

bool check_g(Str In, Str Type) {
	bool RES = false;
	if (Type == "AG") {
		Str SD = "AAGGAGGTGA";
		for (int j = 0; j < In.length() - 4; j++) {
			Str Chk = In.substr(j, 5);
			int pos = SD.find(Chk);
			if (pos != Str::npos)
				RES = true;
		}
	}
	for (int j = 0; j < In.length() - 5; j++) {
		int a, c, g;
		a = c = g = 0;
		int ag, cg;
		ag = cg = 0;
		int bag, bcg;
		bag = bcg = 0;
		for (int k = 0; k < 6; k++) {
			if (In[j + k] == 'G') {
				g++;
				ag++;
				cg++;
			}
			else if (In[j + k] == 'A') {
				a++;
				ag++;
				if (cg > bcg)
					bcg = cg;
				cg = 0;
			}
			else if (In[j + k] == 'C') {
				c++;
				cg++;
				if (ag > bag)
					bag = ag;
				ag = 0;
			}
			else {
				if (ag > bag)
					bag = ag;
				if (cg > bcg)
					bcg = cg;
				ag = cg = 0;
			}
		}
		if (ag > bag)
			bag = ag;
		if (cg > bcg)
			bcg = cg;
		if ((g + a > 4 || bag > 3) && g > 0 && Type == "AG")
			RES = true;
		if ((g + c > 4 || bcg > 3) && g > 0 && Type == "CG")
			RES = true;
	}
	int ag, cg, g;
	ag = cg = g = 0;
	for (int i = 0; i < In.length(); i++) {
		if (In[i] == 'A')
			ag++;
		else if (In[i] == 'C')
			cg++;
		else if (In[i] == 'G') {
			ag++;
			cg++;
			g++;
		}
		if (ag > 5 && g > 0 && Type == "AG")
			RES = true;
		if (cg > 5 && g > 0 && Type == "CG")
			RES = true;
	}
	if (In[2] == 'T' && In[3] == 'A' && In[7] == 'T' && In[4] != 'A' && In[5] != 'T' && In[6] != 'T')
		RES = false;
	return RES;
}

bool check_promoter(Str In, Str Type) {
	bool RES = false;
	Str End = In.substr(4, 4);
	//TTTG signal
	if (In[0] == 'T' && In[1] == 'A' && In[2] != 'G' && In[3] != 'G' && In[4] == 'T' && In[5] == 'T' && In[6] == 'T' && In[7] != 'T' && Type == "TG")
		RES = true;
	//TTTG signal shift
	if (Type == "TG shift") {
		if (In[1] == 'T' && In[2] == 'A' && In[3] != 'G' && In[4] != 'G' && In[5] == 'T' && In[6] == 'T' && In[7] == 'T')
			RES = true;
		else if (In[2] == 'T' && In[3] == 'A' && In[4] != 'G' && In[5] != 'G' && In[6] == 'T' && In[7] == 'T')
			RES = true;
		int pos = In.find("TTTG");
		if (pos < 4 && pos != Str::npos)
			RES = true;
	}
	//Standard TA signal
	if (In[2] == 'T' && In[3] == 'A' && In[4] != 'A' && In[5] != 'T' && In[6] != 'T' && In[7] == 'T' && Type == "STA")
		RES = true;
	//TA signal
	if (Type == "TA") {
		for (int i = 0; i < In.length() - 5; i++) {
			if (In[i] == 'T' && In[i + 1] == 'A' && In[i + 2] != 'A' && In[i + 3] != 'T' && In[i + 4] != 'T' && In[i + 5] == 'T') {
				RES = true;
				break;
			}
		}
	}
	//TA signal shift
	if (Type == "TA shift") {
		if (In[0] == 'A' && In[1] != 'A' && In[2] != 'T' && In[3] != 'T' && In[4] == 'T')
			RES = true;
		else if (In[3] == 'T' && In[4] == 'A' && In[5] != 'A' && In[6] != 'T' && In[7] != 'T')
			RES = true;
		else if (In[4] == 'T' && In[5] == 'A' && In[6] != 'A' && In[7] != 'T')
			RES = true;
	}
	return RES;
}

align dynamic_program(Str Ref, Str Que) {
	align res;
	res.score = res.cont = 0;
	Ve_Ve_P mat;
	Pa_I_I fin;
	//cout << Ref << "\t" << Que << endl;
	for (int i = 0; i < Ref.length(); i++) {
		Ve_P col;
		for (int j = 0; j < Que.length(); j++) {
			dp tmp;
			if (!i || !j) {
				if (Ref[i] == Que[j])
					tmp.continuous = tmp.best = tmp.score = 1;
				else
					tmp.continuous = tmp.best = tmp.score = 0;
				tmp.trace = make_pair(-1, -1);
			}
			else {
				if (Ref[i] == Que[j]) {
					tmp.score = mat[i - 1][j - 1].score + 1;
					tmp.continuous = mat[i - 1][j - 1].continuous + 1;
					tmp.trace = make_pair(i - 1, j - 1);
					if (tmp.continuous > mat[i - 1][j - 1].best)
						tmp.best = tmp.continuous;
					else
						tmp.best = mat[i - 1][j - 1].best;
				}
				else {
					tmp.score = mat[i - 1][j - 1].score - 0.5;
					tmp.trace = make_pair(i - 1, j - 1);
					tmp.continuous = 0;
					tmp.best = mat[i - 1][j - 1].best;
				}
			}
			if (j == Que.length() - 1 || i == Ref.length() - 1) {
				if (tmp.score > res.score || (tmp.score == res.score && tmp.best > res.cont)) {
					res.score = tmp.score;
					res.cont = tmp.best;
					fin = make_pair(i, j);
				}
			}
			col.push_back(tmp);
		}
		mat.push_back(col);
	}
	int cont = 0;
	while (fin.first >= 0 && fin.second >= 0) {
		if (mat[fin.first][fin.second].continuous > cont)
			cont = mat[fin.first][fin.second].continuous;
		res.trace.push_back(fin);
		fin = mat[fin.first][fin.second].trace;
	}
	res.cont = cont;
	//res.ccorr = content(Ref, Que);
	return res;
}

double composition_correlation(Str Ref, Str Que) {
	double res = 0;
	Ma_C_D ref, que;
	for (int i = 0; i < Ref.length(); i++)
		ref[Ref[i]] ++;
	for (int i = 0; i < Que.length(); i++)
		que[Que[i]] ++;
	double p, r, q;
	p = r = q = 0;
	for (int i = 0; i < Base.length(); i++) {
		p += ref[Base[i]] * que[Base[i]];
		r += ref[Base[i]] * ref[Base[i]];
		q += que[Base[i]] * que[Base[i]];
	}
	res = p / sqrt(r) / sqrt(q);
	return res;
}

Ve_Str alignment(Ma_Str_Str fa, int wid, int lim, int mod, Ma_C_D pb) {
	Ve_Str res;
	Ma_Str_I count;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		for (int i = 0; i < iter->second.length() - wid + 1; i++) {
			Str Kmer = iter->second.substr(i, wid);
			count[Kmer] ++;
		}
	}
	Ve_Pa_Str_I rank(count.begin(), count.end());
	sort(rank.begin(), rank.end(), CmpByValue());
	int r = rank.size();
	r *= 0.05;
	int last = 0;
	int all = 0;
	bool AG, CG, TA, TG, ST, SG;
	AG = CG = TA = TG = ST = SG = false;
	Str FTA, FTG;
	int ta = 0;
	int tg = 0;
	for (int i = 0; i < rank.size(); i++) {
		if (i <= r) {
			all += rank[i].second;
			last = rank[i].second;
			if (check_promoter(rank[i].first, "TA")) {
				cout << rank[i].first << " TA found!" << endl;
				ta++;
				if (check_promoter(rank[i].first, "STA") && !ST) {
					FTA = rank[i].first;
					ST = true;
				}
				else if (!TA) {
					FTA = rank[i].first;
					TA = true;
				}
			}
			if (check_promoter(rank[i].first, "TG")) {
				cout << rank[i].first << " TTTG found!" << endl;
				tg++;
				if (check_promoter(rank[i].first, "STG") && !SG) {
					FTG = rank[i].first;
					SG = true;
				}
				else if (!TG) {
					FTG = rank[i].first;
					TG = true;
				}
			}
		}
		else if (i > r) {
			r = i;
			if (last != rank[i].second)
				break;
			else {
				all += rank[i].second;
				if (check_promoter(rank[i].first, "TA")) {
					cout << rank[i].first << " TA found!" << endl;
					ta++;
					if (check_promoter(rank[i].first, "STA") && !ST) {
						FTA = rank[i].first;
						ST = true;
					}
					else if (!TA) {
						FTA = rank[i].first;
						TA = true;
					}
				}
				if (check_promoter(rank[i].first, "TG")) {
					cout << rank[i].first << " TTTG found!" << endl;
					tg++;
					if (check_promoter(rank[i].first, "STG") && !SG) {
						FTG = rank[i].first;
						SG = true;
					}
					else if (!TG) {
						FTG = rank[i].first;
						TG = true;
					}
				}
			}
		}
	}
	//Distance Matrix
	Ma_Str_Ma_Str_A mat;
	for (int i = 0; i < r - 1; i++) {
		for (int j = i; j < r; j++) {
			align tp = dynamic_program(rank[i].first, rank[j].first);
			mat[rank[i].first][rank[j].first] = tp;
			mat[rank[j].first][rank[i].first] = tp;
		}
	}

	//Initialize candidates
	Ve_Str cand;
	double thr = 6.5;
	Pa_D_D bpd;
	bpd.first = bpd.second = 0;
	Ve_Pa_Str_I range;
	Ma_Str_Se_Str grp;
	for (int i = 0; i < r; i++)
		range.push_back(rank[i]);
	while (thr > 4) {
		int sum = 0;
		Ve_Str tmpc;
		Se_Str filter;
		difference best;
		Ma_Str_Se_Str tmpg;
		Ve_Pa_Str_I tmp = range;
		best.Cent = rank[0].first;
		AG = CG = false;
		while (tmpc.size() < mod && !best.Cent.empty()) {
			best.dis = best.size = best.cont = best.size = best.pos = 0;
			best.mem.clear();
			best.Cent.clear();
			//Candidate selection
			for (int i = 0; i < tmp.size(); i++) {
				difference diff;
				diff.Cent = tmp[i].first;
				diff.cont = diff.corr = diff.dis = 0;
				diff.size = tmp[i].second;
				diff.pos = i;
				bool A = check_g(tmp[i].first, "AG");
				bool C = check_g(tmp[i].first, "CG");
				if ((!AG || (AG && !A)) && (!CG || (CG && !C))) {
					for (int j = 0; j < tmp.size(); j++) {
						align tp = mat[tmp[i].first][tmp[j].first];
						diff.dis += tp.score * tmp[j].second;
						diff.cont += tp.cont * tmp[j].second;
					}
					if (diff.dis > best.dis)
						best = diff;
					else if (diff.dis == best.dis && diff.cont > best.cont)
						best = diff;
					else if (diff.dis == best.dis && diff.cont == best.cont && diff.size > best.size)
						best = diff;
				}
			}
			if (!best.Cent.empty()) {
				if (!AG)
					AG = check_g(best.Cent, "AG");
				if (!CG)
					CG = check_g(best.Cent, "CG");
				tmpc.push_back(best.Cent);
				Ve_Pa_Str_I tmpn;
				best.mem.insert(best.Cent);
				for (int i = 0; i < tmp.size(); i++) {
					align tp = mat[tmp[i].first][best.Cent];
					if (tp.score > thr) {
						best.mem.insert(tmp[i].first);
						filter.insert(tmp[i].first);
						sum += tmp[i].second;
					}
					else
						tmpn.push_back(tmp[i]);
				}
				tmp.clear();
				for (int i = 0; i < tmpn.size(); i++) {
					bool FIND = false;
					for (Se_Str::const_iterator iter = best.mem.begin(); iter != best.mem.end(); iter++) {
						align tp = mat[tmpn[i].first][*iter];
						if (tp.score > 6) {
							FIND = true;
							break;
						}
					}
					if (!FIND)
						tmp.push_back(tmpn[i]);
					else {
						best.mem.insert(tmpn[i].first);
						filter.insert(tmpn[i].first);
						sum += tmpn[i].second;
					}
				}
				tmpg[best.Cent] = best.mem;
			}
		}
		if (tmpc.size() * sum * 1.0 / all > bpd.second) {
			bpd.first = thr;
			bpd.second = tmpc.size() * sum * 1.0 / all;
			cand = tmpc;
			grp = tmpg;
		}
		thr -= 0.5;
	}
	bool TAF = false;
	bool TGF = false;
	for (int k = 0; k < cand.size(); k++) {
		cout << cand[k] << endl;
		Str Cand = cand[k];
		if (check_promoter(cand[k], "TA") && !check_promoter(cand[k], "STA")) {
			Pa_Str_I sta;
			sta.first = cand[k];
			sta.second = 0;
			bool STA = false;
			for (Se_Str::const_iterator itr = grp[cand[k]].begin(); itr != grp[cand[k]].end(); itr++) {
				if (check_promoter(*itr, "STA") && count[*itr] > sta.second) {
					sta.first = *itr;
					sta.second = count[*itr];
					STA = true;
				}
			}
			if (STA)
				Cand = sta.first;
			cout << "Standard TA found!\t" << Cand << endl;
		}
		bool TACR = false;
		bool TASFT = true;
		int tggt = cand[k].find("TATGGT");
		int tgct = cand[k].find("TATGCT");
		int tgcta = cand[k].find("TGCTA");
		if (((!check_promoter(cand[k], "TA") && check_promoter(cand[k], "TA shift")) || tggt != Str::npos || tgct != Str::npos || tgcta != Str::npos) && TA) {
			TASFT = false;
			Pa_Str_I sta;
			sta.first = cand[k];
			sta.second = 0;
			for (Se_Str::const_iterator itr = grp[cand[k]].begin(); itr != grp[cand[k]].end(); itr++) {
				Str Tmp = *itr;
				tggt = Tmp.find("TATGGT");
				tgct = Tmp.find("TATGCT");
				tgcta = Tmp.find("TGCTA");
				if (check_promoter(*itr, "TA") && tggt == Str::npos && tgct == Str::npos && tgcta == Str::npos) {
					if ((check_promoter(sta.first, "STA") == check_promoter(*itr, "STA") && count[*itr] > sta.second)
						|| (check_promoter(*itr, "STA") && !check_promoter(sta.first, "STA"))) {
						sta.first = *itr;
						sta.second = count[*itr];
					}
					TACR = true;
				}
			}
			if (TACR)
				Cand = sta.first;
			cout << "TA shift corrected!\t" << Cand << endl;
		}
		bool TGCR = false;
		bool TGSFT = true;
		if (check_promoter(cand[k], "TG shift") && !check_promoter(cand[k], "TG") && TG) {
			TGSFT = false;
			Pa_Str_I sta;
			sta.first = "";
			sta.second = 0;
			for (Se_Str::const_iterator itr = grp[cand[k]].begin(); itr != grp[cand[k]].end(); itr++) {
				if (check_promoter(*itr, "TG")) {
					if (count[*itr] > sta.second) {
						sta.first = *itr;
						sta.second = count[*itr];
					}
					TGCR = true;
				}
			}
			if (TGCR)
				Cand = sta.first;
		}
		if (cand[k] == "AGATTCCG") {
			Pa_Str_I sta;
			sta.first = "";
			sta.second = 0;
			bool FIND = false;
			for (Se_Str::const_iterator itr = grp[cand[k]].begin(); itr != grp[cand[k]].end(); itr++) {
				Str Tmp;
				if (Tmp[0] == 'T' && Tmp[1] == 'A' && Tmp[2] == 'G') {
					if (count[*itr] > sta.second) {
						sta.first = *itr;
						sta.second = count[*itr];
					}
					FIND = true;
				}
				if (FIND)
					Cand = sta.first;
			}
		}
		if (grp[cand[k]].size() > 1 && (TGSFT || TGCR) && (TASFT || TACR)) {
			res.push_back(Cand);
			if (check_promoter(Cand, "TA"))
				TAF = true;
			if (check_promoter(Cand, "TG"))
				TGF = true;
		}
	}
	if (!TAF && !FTA.empty() && ta > 1) 
		res.push_back(FTA);
	if (!TGF && !FTG.empty() && tg > 1) 
		res.push_back(FTG);
	if (res.size() < mod) {
		for (int i = 0; i < lim; i++)
			res.push_back(rank[i].first);
	}
	return res;
}

Ve_Ma_C_D pwm_assign(Str Seq, Ma_C_D pb, double maj, double min) {
	Ve_Ma_C_D res;
	for (int i = 0; i < Seq.size(); i++) {
		Ma_C_D temp;
		for (int j = 0; j < Base.size(); j++) {
			if (Seq[i] == Base[j])
				temp[Base[j]] = maj / pb[Base[j]];
			else
				temp[Base[j]] = min / pb[Base[j]];
		}
		res.push_back(temp);
	}
	return res;
}

Ma_C_D background(Ma_Str_Str fa) {
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
	//for (int j = 0; j < Base.length(); j++)
	//res[Base[j]] /= sum;
	return res;
}

Ma_Str_D background_score(Ma_Str_Str fa, Ma_C_D pb) {
	Ma_Str_D res;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		double b = 0;
		for (int i = 0; i < iter->second.size(); i++)
			b += log2(pb[iter->second[i]]);
		res[iter->first] = b;
	}
	return res;
}

motif EM(motif test, Ma_Str_Str fa, int mod, int reg, double eps, int step) {
	motif res = test;
	Ma_Str_Ma_C_D bct;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		for (int i = 0; i < iter->second.length(); i++)
			bct[iter->first][iter->second[i]] ++;
	}
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
		//Ma_Str_D bs = background_score(fa, res.pb);
		int k = 0;
		double p0 = 0;
		for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
			double sumk = res.p0;
			for (int j = 0; j < reg - res.rwm[0].size() + 1; j++) {
				for (int m = 0; m < mod; m++) {
					lam[m][k][j] = res.pj[m][j];
					for (int i = 0; i < res.rwm[m].size(); i++) {
						if (iter->second[i + j] == 'A' || iter->second[i + j] == 'C' || iter->second[i + j] == 'G' || iter->second[i + j] == 'T')
							lam[m][k][j] *= res.rwm[m][i][iter->second[i + j]];
					}
					sumk += lam[m][k][j];
				}
			}
			for (int v = 0; v < Base.length(); v++)
				res.loglike += bct[iter->first][Base[v]] * log2(res.pb[Base[v]]);
			res.loglike += log2(sumk - res.p0);
			//res.loglike += log2(sumk);
			//cout << n << "\t" << iter->first << "\t" << bs[iter->first] << "\t" << res.loglike << endl;
			p0 += res.p0 / sumk;
			//for (int v = 0; v < Base.length(); v++)
			//	B[Base[v]] += p0 * bct[iter->first][Base[v]];
			for (int j = 0; j < reg - res.rwm[0].size() + 1; j++) {
				for (int m = 0; m < mod; m++) {
					lam[m][k][j] /= sumk;
					C[m][j] += lam[m][k][j];
					for (int i = 0; i < res.rwm[m].size(); i++) {
						D[m][i][iter->second[i + j]] += lam[m][k][j];
						B[iter->second[i + j]] -= lam[m][k][j];
					}
					//cout << B['C'] << "\t" << lam[m][k][j] << endl;
				}
			}
			k++;
			//cout << "Background log likelihood:\t" << bs[iter->first] << "\tLog likelihood:\t" << res.loglike << endl;
		}
		if (abs(res.loglike - res.last) > eps * abs(res.loglike)) {
			//cout << res.loglike << endl;
			res.trace.push_back(res.loglike);
			res.last = res.loglike;
			res.loglike = 0;
			double sumb = 0;
			for (int v = 0; v < Base.length(); v++) {
				B[Base[v]] += res.pc[Base[v]];
				if (B[Base[v]] < 0) {
					cout << "Error in iterating background probability of " << Base[v] << "!" << endl;
					B[Base[v]] = 1.0e-5;
				}
				sumb += B[Base[v]];
			}
			//pb estimation
			for (int v = 0; v < Base.length(); v++) {
				res.pb[Base[v]] = B[Base[v]] / sumb;
				//if (res.pb[Base[v]] < 0)
				//	res.pb[Base[v]] = fabs(res.pb[Base[v]]);
				//cout << Base[v] << "\t" << res.pc[Base[v]] << "\t" << B[Base[v]] << "\t" << res.pb[Base[v]] << endl;
			}
			double sumc = 0;
			for (int m = 0; m < mod; m++) {
				//cout << endl << m << endl;
				for (int i = 0; i < res.rwm[m].size(); i++) {
					double sumd = 0;
					//cout << i;
					for (int v = 0; v < Base.length(); v++)
						sumd += D[m][i][Base[v]];
					//RWM estimation
					for (int v = 0; v < Base.length(); v++) {
						D[m][i][Base[v]] /= sumd;
						res.rwm[m][i][Base[v]] = D[m][i][Base[v]] / res.pb[Base[v]];
						//cout << "\t" << res.rwm[m][i][Base[v]];
					}
					//cout << endl;
				}
				for (int j = 0; j < C[m].size(); j++)
					sumc += C[m][j];
			}
			for (int m = 0; m < mod; m++) {
				//pj estimation
				//cout << m << ":";
				for (int j = 0; j < C[m].size(); j++) {
					res.pj[m][j] = C[m][j] / (sumc + p0);
					//cout << "\t" << res.pj[m][j];
				}
				//cout << endl;
			}
			//p0 estimation
			res.p0 = p0 / (sumc + p0);
			//cout << res.p0 << endl;
		}
		else {
			res.trace.push_back(res.loglike);
			res.last = res.loglike;
			res.entropy.clear();
			for (int m = 0; m < mod; m++) {
				double entropy = 0;
				for (int i = 0; i < res.rwm[m].size(); i++) {
					for (int v = 0; v < Base.length(); v++) {
						if (res.rwm[m][i][Base[v]] == 0)
							entropy += 0;
						else
							entropy += res.pb[Base[v]] * res.rwm[m][i][Base[v]] * log2(res.rwm[m][i][Base[v]]);
					}
				}
				res.entropy.push_back(entropy);
			}
			break;
		}
	}
	return res;
}

Str TIS_motif(Str Out, Str GenoID, Ma_Str_Str fa, double thr, int mod, int wid, double maj, double min, double eps, int trial, Ma_C_D pb) {
	Str Mod;
	convertFromNumber(Mod, mod);
	Out += GenoID + ".sig.dat";
	int reg = 0;
	Ma_Str_Str::const_iterator ir = fa.begin();
	reg = ir->second.length();
	//Ve_Str km = alignment(Out, fa, wid, thr, mod);
	//Ma_C_D npb = background_probability(fa);
	Ve_Str km = alignment(fa, wid, thr, mod, pb);
	Ofs out(Out.data());
	out << "Selected K-mers:" << endl;
	for (int l = 0; l < km.size(); l++) 
		out << km[l] << endl;
	cout << "Selected K-mer types: " << km.size() << endl;
	long long int max = 1;
	Ve_I limit;
	motif res;
	res.last = -1.0e16;
	res.p0 = 1;
	for (int m = 0; m < mod; m++) {
		limit.push_back(km.size() - mod + 1);
		max *= km.size() - mod + 1;
	}
	int ins = 0;
	Ve_M cand;
	Ve_Str fin;
	for (int k = 0; k < max; k++) {
		De_I sub = subscript(k, limit);
		bool CHECK = true;
		for (int m = 1; m < mod; m++) {
			//cout << sub[m - 1] + m - 1 << "\t" << sub[m] + m << endl;
			if (sub[m - 1] + m - 1 >= sub[m] + m)
				CHECK = false;
		}
		if (CHECK) {
			ins++;
			cout << "Initialization seed: " << ins << endl;
			motif temp;
			temp.loglike = 0;
			//temp.pb = pb;
			temp.pc = background(fa);
			double sum = 0;
			for (int i = 0; i < Base.size(); i++)
				sum += temp.pc[Base[i]];
			for (int i = 0; i < Base.size(); i++)
				temp.pb[Base[i]] = temp.pc[Base[i]] / sum;
			temp.p0 = 1.0 / ((reg - wid + 1) * mod + 1);
			cout << "Initial Sequences:";
			Ve_Str ini;
			for (int m = 0; m < mod; m++) {
				Ve_Ma_C_D mat = pwm_assign(km[sub[m] + m], temp.pb, maj, min);
				temp.rwm.push_back(mat);
				//temp.rwm.push_back(pwm_assign(km[sub[m] + m], temp.pb, maj, min));
				temp.cons.push_back(km[sub[m] + m]);
				cout << "\t" << km[sub[m] + m];
				ini.push_back(km[sub[m] + m]);
			}
			Ve_D pj;
			pj.resize(reg - wid + 1, 1.0 / ((reg - wid + 1) * mod + 1));
			for (int m = 0; m < mod; m++)
				temp.pj.push_back(pj);
			temp.entropy.resize(mod, 0);
			//temp = EM(temp, fa, mod, reg, eps, trial);
			temp = EM(temp, fa, mod, reg, eps);
			cout << "\t" << temp.trace.size() << "\t" << temp.p0 << "\t" << temp.last << endl;
			if (temp.last > res.last) {
				res = temp;
				fin = ini;
			}
		}
	}
	cout << "Finally chosen K-mers for initialization:" << endl;
	for (int k = 0; k < fin.size(); k++)
		cout << "\t" << fin[k];
	cout << endl;
	res = EM(res, fa, mod, reg, eps);
	out << "Initial k-mer:";
	for (int m = 0; m < mod; m++)
		out << "\t" << res.cons[m];
	out << endl;
	out << "Genome:" << endl
		<< GenoID << endl
		<< "Features of translation initiation site (TIS) signal:" << endl
		<< "Final likelihood score of this model:\t" << res.loglike << endl
		<< "Probability of background nucleotide bias:" << endl
		<< "\tA\tC\tG\tT" << endl;
	for (int r = 0; r < Base.length(); r++)
		out << "\t" << res.pb[Base[r]];
	out << endl << "Probabilities of TIS start position:" << endl
		<< "No signal:" << endl << "\t" << res.p0 << endl;
	for (int r = 0; r < res.pj.size(); r++) {
		out << "Motif" << r + 1 << ":" << endl;
		for (int q = 0; q < res.pj[r].size(); q++)
			out << "\t" << res.pj[r][q];
		out << endl;
	}
	bool EXC = false;
	for (int s = 0; s < Base.length(); s++) {
		if (res.pb[Base[s]] > 0.4) {
			cout << Base[s] << " as " << res.pb[Base[s]] << " exceeds 0.4!" << endl;
			EXC = true;
		}
	}
	for (int r = 0; r < res.rwm.size(); r++) {
		Str Motif;
		out << "Position weight matrix (PWM) of Motif" << r + 1 << ":" << endl;
		//<< "\tA\tC\tG\tT" << endl;
		for (int q = 0; q < res.rwm[r].size(); q++) {
			//out << q + 1;
			double last = 0;
			int base = 0;
			for (int s = 0; s < Base.length(); s++) {
				if (log2(res.rwm[r][q][Base[s]]) * res.pb[Base[s]] > last) {
					last = log2(res.rwm[r][q][Base[s]]) * res.pb[Base[s]];
					base = s;
				}
			}
			Motif += Base[base];
			for (int s = 0; s < Base.length(); s++)
				out << "\t" << res.rwm[r][q][Base[s]] * res.pb[Base[s]];
			out << endl;
		}
		out << "Consensus sequenece of Motif" << r + 1 << ":\t" << Motif << endl;
		out << "Relative entropy of Motif" << r + 1 << ":\t" << res.entropy[r] << endl;
	}
	out.close();
	return "TIS signals of " + Mod + " types have been modeled!";
}

