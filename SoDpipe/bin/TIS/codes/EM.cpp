#include "EM.h"
using namespace std;

bool check_cg(Str In) {
	bool RES = false;
	for (int j = 0; j < In.length() - 5; j++) {
		int g = 0;
		int c = 0;
		for (int k = 0; k < 6; k++) {
			if (In[j + k] == 'G')
				g++;
			else if (In[j + k] == 'C' || In[j + k] == 'S')
				c++;
		}
		if (c + g > 4 && g > 0 && c > 0)
			RES = true;
	}
	int g = 0;
	int c = 0;
	for (int i = 0; i < In.length(); i++) {
		if (In[i] == 'C')
			c++;
		if (In[i] == 'G')
			g++;
	}
	if (c + g > 5 && c > 0 && g > 0)
		RES = true;
	int ta = In.find("TAGGC");
	int at = In.find("AGGCT");
	if (ta != Str::npos || at != Str::npos)
		RES = false;
	return RES;
}

bool check_ag(Str In) {
	bool RES = false;
	Str SD = "AAGGAGGTGA";
	for (int j = 0; j < In.length() - 4; j++) {
		Str Chk = In.substr(j, 5);
		int pos = SD.find(Chk);
		if (pos != Str::npos)
			RES = true;
	}
	for (int j = 0; j < In.length() - 5; j++) {
		int a = 0;
		int g = 0;
		int cont = 0;
		int bestcont = 0;
		bool C = false;
		for (int k = 0; k < 6; k++) {
			if (In[j + k] == 'A' || In[j + k] == 'R') {
				a++;
				cont++;
			}
			else if (In[j + k] == 'G') {
				g++;
				cont++;
			}
			else {
				if (In[j + k] == 'C')
					C = true;
				if (cont > bestcont)
					bestcont = cont;
				cont = 0;
			}
		}
		if (cont > bestcont)
			bestcont = cont;
		if (((a + g > 4 && !C) || bestcont > 3) && a > 0 && g > 0)
			//if ((bestcont > 4 && a > 0 && g > 0) || pos != Str::npos)
			RES = true;
	}
	/*int g = 0;
	int a = 0;
	for (int i = 0; i < In.length(); i++) {
	if (In[i] == 'A')
	a++;
	if (In[i] == 'G')
	g++;
	}
	if (a + g > 5 && a > 0 && g > 0)
	RES = true;
	//int tac = In.find("TAGGC");
	//int agc = In.find("AGGCT");
	//int ata = In.find("ATAAT");*/
	int taa = In.find("TAGGAT");
	int tac = In.find("TAGGCT");
	int tag = In.find("TAGGGT");
	int aga = In.find("TAGAAT");
	if (taa != Str::npos || tag != Str::npos || tac != Str::npos || aga != Str::npos)
		RES = false;
	return RES;
}

align dynamic_program(Str Ref, Str Que) {
	align res;
	res.score = 1000;
	Ve_Ve_P mat;
	Pa_I_I fin;
	//cout << Ref << "\t" << Que << endl;
	for (int i = 0; i < Ref.length(); i++) {
		Ve_P col;
		for (int j = 0; j < Que.length(); j++) {
			dp tmp;
			//cout << Ref[i] << "\t" << Que[j] << "\t" << i << "\t" << j << endl;
			if (!i || !j) {
				if (Ref[i] == Que[j])
					tmp.score = 0;
				else
					tmp.score = 1;
				if (max(i, j) < 3)
					tmp.score += max(i, j) * 0.5;
				else
					tmp.score += max(i, j) - 1;
				tmp.trace = make_pair(-1, -1);
			}
			else {
				if (Ref[i] == Que[j]) {
					tmp.score = mat[i - 1][j - 1].score;
					tmp.trace = make_pair(i - 1, j - 1);
				}
				else {
					tmp.score = mat[i - 1][j - 1].score + 1;
					tmp.trace = make_pair(i - 1, j - 1);
				}
			}
			//cout << Ref[i] << "\t" << Que[j] << "\t" << i << "\t" << j
			//	<< "\t" << tmp.score << "\t" << tmp.trace.first << "\t" << tmp.trace.second << endl;
			if (tmp.score < res.score && (j == Que.length() - 1 || i == Ref.length() - 1)) {
				res.score = tmp.score;
				fin = make_pair(i, j);
			}
			col.push_back(tmp);
		}
		mat.push_back(col);
	}
	while (fin.first >= 0 && fin.second >= 0) {
		res.trace.push_back(fin);
		fin = mat[fin.first][fin.second].trace;
	}
	return res;
}

double comp_correlation(Str Ref, Str Que) {
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

align similar(Str Ref, Str Que) {
	align res;
	res.score = 0;
	Ve_Ve_P mat;
	Pa_I_I fin;
	//cout << Ref << "\t" << Que << endl;
	for (int i = 0; i < Ref.length(); i++) {
		Ve_P col;
		for (int j = 0; j < Que.length(); j++) {
			dp tmp;
			if (!i || !j) {
				if (Ref[i] == Que[j])
					tmp.continuous = tmp.score = 1;
				else
					tmp.continuous = tmp.score = 0;
				tmp.trace = make_pair(-1, -1);
			}
			else {
				if (Ref[i] == Que[j]) {
					tmp.score = mat[i - 1][j - 1].score + 1;
					tmp.continuous = mat[i - 1][j - 1].continuous + 1;
					tmp.trace = make_pair(i - 1, j - 1);
				}
				else {
					tmp.score = mat[i - 1][j - 1].score;
					tmp.trace = make_pair(i - 1, j - 1);
					tmp.continuous = 0;
				}
			}
			if (tmp.score > res.score && (j == Que.length() - 1 || i == Ref.length() - 1)) {
				res.score = tmp.score;
				fin = make_pair(i, j);
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
	return res;
}

asb pwm_transfer(Pa_Str_I ini) {
	asb res;
	for (int i = 0; i < ini.first.length(); i++) {
		Ma_C_D tmpc;
		tmpc[ini.first[i]] = ini.second;
		res.mat.push_back(tmpc);
	}
	res.stt = 0;
	return res;
}

asb assemble(asb last, align alg, Pa_Str_I curr) {
	asb res;
	int n = alg.trace.size() - 1;
	//Reference Sequence Protrudes
	//cout << alg.trace[n].first << "\t" << alg.trace[n].second << "\t" << last.stt << "\t" << curr.first << "\t" << curr.second << endl;
	if (alg.trace[n].first >= alg.trace[n].second) {
		res = last;
		for (int i = 0; i < curr.first.size(); i++) {
			//cout << last.stt + alg.trace[n].second + i << "\t" << last.mat.size() << endl;
			if (last.stt + alg.trace[n].first + i < last.mat.size()) {
				//if(curr.second > res.mat[last.stt + alg.trace[n].second + i][curr.first[i]])
				res.mat[last.stt + alg.trace[n].first + i][curr.first[i]] += curr.second;
			}
			else {
				Ma_C_D tmpc;
				tmpc[curr.first[i]] = curr.second;
				res.mat.push_back(tmpc);
			}
		}
	}
	//Query Sequence Protrudes
	else {
		if (last.stt >= alg.trace[n].second) {
			res = last;
			for (int i = 0; i < curr.first.size(); i++) {
				//if (curr.second > res.mat[last.stt - alg.trace[n].second + i][curr.first[i]])
				res.mat[last.stt - alg.trace[n].second + i][curr.first[i]] += curr.second;
			}
		}
		else {
			for (int i = 0; i < curr.first.size(); i++) {
				Ma_C_D tmpc;
				tmpc[curr.first[i]] = curr.second;
				res.mat.push_back(tmpc);
			}
			for (int i = 0; i < last.mat.size(); i++) {
				if (alg.trace[n].second - last.stt + i < res.mat.size()) {
					for (int j = 0; j < Base.size(); j++) {
						//if (last.mat[i][Base[j]] > res.mat[alg.trace[n].second - last.stt + i][Base[j]])
						res.mat[alg.trace[n].second - last.stt + i][Base[j]] += last.mat[i][Base[j]];
					}
				}
				else {
					Ma_C_D tmpc = last.mat[i];
					tmpc[curr.first[i]] = curr.second;
					res.mat.push_back(tmpc);
				}
			}
			res.stt = alg.trace[n].second;
		}
	}
	return res;
}

char nucleotide_code(Str In, int sub) {
	char res = 'N';
	if (!sub)
		res = In[sub];
	else if (sub == 1) {
		Str Sub = In.substr(0, sub + 1);
		if (Sub == "AC")
			res = 'M';
		else if (Sub == "AG")
			res = 'R';
		else if (Sub == "AT" || Sub == "AU")
			res = 'W';
		else if (Sub == "CG")
			res = 'S';
		else if (Sub == "CT" || Sub == "CU")
			res = 'Y';
		else if (Sub == "GT" || Sub == "GU")
			res = 'K';
	}
	else if (sub == 2) {
		Str Sub = In.substr(0, sub + 1);
		if (Sub == "ACG")
			res = 'V';
		else if (Sub == "ACT" || Sub == "ACU")
			res = 'H';
		else if (Sub == "AGT" || Sub == "AGU")
			res = 'D';
		else if (Sub == "CGT" || Sub == "CGU")
			res = 'B';
	}
	return res;
}

Ve_S alignment(Ma_Str_Str fa, int wid, int cand, int mod, Ma_C_D pb, double maj, double min) {
	Ve_S res;
	Ma_Str_I count;
	Ma_Str_Ve_D pos;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		for (int i = 0; i < iter->second.length() - wid + 1; i++) {
			Str Kmer = iter->second.substr(i, wid);
			count[Kmer] ++;
			if (pos[Kmer].empty())
				pos[Kmer].resize(iter->second.length() - wid + 1, 0);
			pos[Kmer][i] ++;
		}
	}
	Ve_Pa_Str_I rank(count.begin(), count.end());
	sort(rank.begin(), rank.end(), CmpByValue());
	int r = rank.size();
	r *= 0.05;
	int last = 0;
	int all = 0;
	bool TA = false;
	bool TG = false;
	for (int i = 0; i < rank.size(); i++) {
		//log << rank[i].first << "\t" << rank[i].second << endl;
		all += rank[i].second;
		if (!TA) {
			int posta = rank[i].first.find("TATAAT");
			int postagc = rank[i].first.find("TAGGCT");
			int postact = rank[i].first.find("TACACT");
			if (posta != Str::npos || postagc != Str::npos || postact != Str::npos)
				TA = true;
		}
		if (!TG) {
			int postg = rank[i].first.find("TTTG");
			if (postg == 4 && rank[i].first[0] == 'T')
				TG = true;
		}
		if (i <= r)
			last = rank[i].second;
		else if (i > r) {
			r = i;
			if (last != rank[i].second)
				break;
			else
				all += rank[i].second;
		}
	}
	Ma_Str_Ma_Str_A mat;
	for (int i = 0; i < r - 1; i++) {
		for (int j = i; j < r; j++) {
			mat[rank[i].first][rank[j].first] = similar(rank[i].first, rank[j].first);
			mat[rank[j].first][rank[i].first] = similar(rank[j].first, rank[i].first);
			mat[rank[j].first][rank[i].first].ccorr = mat[rank[i].first][rank[j].first].ccorr = comp_correlation(rank[i].first, rank[j].first);
			mat[rank[i].first][rank[j].first].dis = mat[rank[i].first][rank[j].first].score + mat[rank[i].first][rank[j].first].ccorr;
			mat[rank[j].first][rank[i].first].dis = mat[rank[j].first][rank[i].first].score + mat[rank[j].first][rank[i].first].ccorr;//  + mat[rank[i].first][rank[j].first].pcorr +mat[rank[i].first][rank[j].first].cont;
		}
	}
	bool AG = check_ag(rank[0].first);
	bool CG = check_cg(rank[0].first);
	//Filteration
	Ma_Str_I grz;
	Ma_Str_Se_Str grp;
	//grp[rank[0].first].insert(rank[0].first);
	double thr = wid / 2.0;
	Pa_D_D bpd;
	bpd.first = bpd.second = 0;
	while (thr <= wid) {
		Ve_Str tmpo;
		Ma_Str_I tmpz;
		Ma_Str_Se_Str tmpg;
		tmpg[rank[0].first].insert(rank[0].first);
		for (int i = 1; i < r; i++) {
			bool RED = false;
			if ((AG && check_ag(rank[i].first)) || (CG && check_cg(rank[i].first)))
				RED = true;
			for (Ma_Str_Se_Str::const_iterator iter = tmpg.begin(); iter != tmpg.end(); iter++) {
				if (mat[iter->first][rank[i].first].dis >= thr)
					RED = true;
			}
			if (!RED && tmpg.size() < cand) {
				tmpg[rank[i].first].insert(rank[i].first);
				tmpo.push_back(rank[i].first);
				if (!AG)
					AG = check_ag(rank[i].first);
				if (!CG)
					CG = check_cg(rank[i].first);
			}
			if (tmpg.size() == cand)
				break;
		}
		//Classification
		for (int i = 0; i < r; i++) {
			Pa_Str_D best;
			best.second = 0;
			for (int j = 0; j < tmpo.size(); j++) {
				AG = check_ag(tmpo[j]);
				CG = check_cg(tmpo[j]);
				if (mat[tmpo[j]][rank[i].first].dis > best.second && AG == check_ag(rank[i].first) && CG == check_cg(rank[i].first)) {
					best.first = tmpo[j];
					best.second = mat[tmpo[j]][rank[i].first].dis;
				}
			}
			if (best.second > thr)
				tmpg[best.first].insert(rank[i].first);
		}
		//Count group size
		double sum = 0;
		cout << "Group center:" << endl;
		for (int j = 0; j < tmpo.size(); j ++) {
			cout << "\t" << tmpo[j];
			tmpz[tmpo[j]] = 0;
			for (Se_Str::const_iterator ite = tmpg[tmpo[j]].begin(); ite != tmpg[tmpo[j]].end(); ite++) {
				tmpz[tmpo[j]] += count[*ite];
				sum += count[*ite];
			}
		}
		cout << endl;
		int m = tmpg.size();
		if (bpd.first < sum / all * m) {
			bpd.first = sum / all * m;
			bpd.second = thr;
			grp = tmpg;
			grz = tmpz;
		}
		cout << "Scoring:\t" << sum << "\t" << all << "\t" << sum / all << "\t" << m << "\t" << thr 
			<< "\t" << bpd.first << "\t" << bpd.second << "\t" << grp.size() << endl;
		thr += 0.5;
	}
	Ve_Pa_Str_I grank(grz.begin(), grz.end());
	sort(grank.begin(), grank.end(), CmpByValue());
	/*tot << "Center\t" << GenoID;
	for (int i = 0; i < grank.size(); i++) {
		if (i < mod)
			tot << "\t" << grank[i].first;
		log << ">" << grank[i].first << "\t" << grank[i].second << endl;
		for (Se_Str::const_iterator ite = grp[grank[i].first].begin(); ite != grp[grank[i].first].end(); ite++)
			log << "\t" << *ite;
		log << endl;
	}
	tot << endl;
	tot << "After Assembled\t" << GenoID;*/
	int app = 0;
	AG = CG = false;
	Ve_Str candidates;
	for (int i = 0; i < grank.size(); i++) {
		if (grp[grank[i].first].size() < 2)
			continue;
		Pa_Str_I tmp = make_pair(grank[i].first, count[grank[i].first]);
		asb pwm = pwm_transfer(tmp);
		pwm.stt = 0;
		//cout << ">" << grank[i].first << endl;
		for (Se_Str::const_iterator it = grp[grank[i].first].begin(); it != grp[grank[i].first].end(); it++) {
			//cout << "\t" << *it;
			tmp = make_pair(*it, count[*it]);
			if (*it != grank[i].first)
				pwm = assemble(pwm, mat[grank[i].first][*it], tmp);
		}
		//cout << endl;
		Str Csn(pwm.mat.size(), 'A');
		//log << endl << "\tA\tC\tG\tT\tA\tC\tG\tT" << endl;
		for (int j = 0; j < pwm.mat.size(); j++) {
			double sum = 0;
			for (int l = 0; l < Base.size(); l++) {
				//log << "\t" << pwm.mat[j][Base[l]];
				sum += pwm.mat[j][Base[l]];
			}
			double max = 0;
			Str Tmp = Base;
			int sub = 0;
			for (int l = 0; l < Base.size(); l++) {
				double perc = (pwm.mat[j][Base[l]] + 0.5) / (sum + 2);
				//log << "\t" << perc;
				if (perc > max) {
					max = perc;
					sub = 0;
					Tmp[sub] = Base[l];
				}
				else if (perc == max) {
					sub++;
					Tmp[sub] = Base[l];
				}
			}
			Csn[j] = nucleotide_code(Tmp, sub);
			//log << endl;
		}
		Str Fin = Csn.substr(pwm.stt, wid);
		//log << "Final Consensus:\t" << Fin << "\tStart Site:\t" << pwm.stt << endl;
		int start = pwm.stt;
		cout << Csn << "\t" << Fin << endl;
		if (!check_ag(Fin)) {
			if (TA) {
				int ta = Csn.find("TATAA");
				int tat = Csn.find("TATAAT");
				int act = Csn.find("TATACT");
				int tac = Csn.find("TACAAT");
				int cac = Csn.find("TACACT");
				int tag = Csn.find("TAGAAT");
				int tgc = Csn.find("TAGGCT");
				int gcc = Csn.find("TAGCCT");
				int gcg = Csn.find("TAGCGT");
				int cag = Csn.find("TACAGT");
				bool COR = false;
				if (tat != Str::npos) {
					start = tat - 2;
					COR = true;
				}
				else if (tgc != Str::npos) {
					start = tgc - 2;
					COR = true;
				}
				else if (act != Str::npos) {
					start = act - 2;
					COR = true;
				}
				else if (tac != Str::npos) {
					start = tac - 2;
					COR = true;
				}
				else if (cac != Str::npos) {
					start = cac - 2;
					COR = true;
				}
				else if (tag != Str::npos) {
					start = tag - 2;
					COR = true;
				}
				else if (gcc != Str::npos) {
					start = gcc - 2;
					COR = true;
				}
				else if (gcg != Str::npos) {
					start = gcg - 2;
					COR = true;
				}
				else if (cag != Str::npos) {
					start = cag - 2;
					COR = true;
				}
				else if (ta != Str::npos) {
					start = ta - 2;
					COR = true;
				}
				if (start < 0)
					start = 0;
				Fin = Csn.substr(start, wid);
				if (COR) {
					cout << "TATA Check!" << endl;
					if (Fin.length() != wid) {
						cout << "Correction failed!" << endl;
						continue;
					}
					else
						cout << "Corrected as " << Fin << endl;
				}
			}
			if (TG) {
				int ttg = Csn.find("TTTG");
				int tta = Csn.find("TATTTTTA");
				int tct = Csn.find("TATCTTT");
				bool COR = false;
				bool SUC = false;
				if (ttg != Str::npos)
					COR = true;
				while (ttg != Str::npos) {
					if (ttg - 4 >= 0) {
						if (Csn[ttg - 4] == 'T' && Csn[ttg - 3] == 'A') {
							SUC = true;
							start = ttg - 4;
							break;
						}
					}
					int next = ttg + 1;
					ttg = Csn.find("TTTG", next);
				}
				if (!SUC && (tta != Str::npos || tct != Str::npos)) {
					if (tta != Str::npos)
						start = tta;
					if (tct != Str::npos)
						start = tct;
					COR = true;
					SUC = true;
				}
				Fin = Csn.substr(start, wid);
				if (COR) {
					cout << "TTTG Check:" << endl;
					if (SUC) {
						if (Fin.length() < wid) {
							cout << "Correction failed!" << endl;
							continue;
						}
						else 
							cout << "Corrected as " << Fin << endl;
					}
					else {
						cout << "Correction failed!" << endl;
						continue;
					}
				}
			}
		}
		asb fin;
		for (int j = start; j < start + wid; j++) {
			int hit = 0;
			int bsc = 0;
			double sum = 0;
			for (int l = 0; l < Base.length(); l++) {
				//cout << "\t" << pwm.mat[j][Base[l]];
				if (pwm.mat[j][Base[l]] > bsc) {
					bsc = pwm.mat[j][Base[l]];
					hit = l;
				}
				if (pwm.mat[j][Base[l]] > 1.0e-5) {
					hit++;
					sum += pwm.mat[j][Base[l]];
				}
			}
			Ma_C_D mat;
			//cout << "\t" << sum << "\t" << hit << endl;
			/*for (int l = 0; l < Base.length(); l++) {
				if (hit < 4) {
					if (pwm.mat[j][Base[l]] == bsc)
						mat[Base[l]] = maj / pb[Base[l]];
					else
						mat[Base[l]] = min / pb[Base[l]];
				}
				else
					mat[Base[l]] = (pwm.mat[j][Base[l]] + 0.5) / (sum + 2) / pb[Base[l]];
				//cout << "\t" << mat[Base[l]] << "\t" << pb[Base[l]];
			}*/
			//cout << endl;
			//double maj = 0.52;
			//double min = 0.16;
			for (int l = 0; l < Base.length(); l++) {
				if (pwm.mat[j][Base[l]] > 1.0e-5) {
					if (hit == 1)
						mat[Base[l]] = maj / pb[Base[l]];
					else if (hit == 2)
						mat[Base[l]] = 0.8 * pwm.mat[j][Base[l]] / sum / pb[Base[l]];
					else
						mat[Base[l]] = (pwm.mat[j][Base[l]] + 0.5) / (sum + 2) / pb[Base[l]];
				}
				else {
					if (hit == 1)
						mat[Base[l]] = min / pb[Base[l]];
					else if (hit == 2)
						mat[Base[l]] = 0.1 / pb[Base[l]];
					else
						mat[Base[l]] = (pwm.mat[j][Base[l]] + 0.5) / (sum + 2) / pb[Base[l]];
				}
				//log << "\t" << mat[Base[l]];
			}
			//log << endl;
			fin.mat.push_back(mat);
		}
		bool FIND = false;
		for (int j = 0; j < candidates.size(); j++) {
			if (Fin == candidates[j])
				FIND = true;
		}
		if (!FIND) {
			fin.Consensus = Fin;
			if ((!AG || (AG && !check_ag(Fin))) && (!CG || (CG && !check_cg(Fin)))) {
				//log << "Candidate approved!" << endl;
				if (res.size() < cand)
					res.push_back(fin);
				else
					break;
					//tot << "\t" << Fin;
				if (!AG)
					AG = check_ag(Fin);
				if (!CG)
					CG = check_cg(Fin);
			}
			//else if (AG && check_ag(Fin))
				//log << "AG rich!" << endl;
			//else if (CG && check_cg(Fin))
				//log << "CG rich!" << endl;
		}
	}
	cout << res.size() << endl;
	if (res.size() < mod) {
		res.clear();
		for (int i = 0; i < cand; i++) {
			asb fin;
			fin.mat = pwm_assign(rank[i].first, pb, maj, min);
			fin.Consensus = rank[i].first;
			res.push_back(fin);
		}
	}
	return res;
}

Ve_Ma_C_D pwm_assign(Str Seq, Ma_C_D pb, double maj, double min) {
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
			//double sumk = res.p0;
			double sumk = 0;
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
			res.loglike += log2(sumk);
			//cout << n << "\t" << iter->first << "\t" << bs[iter->first] << "\t" << res.loglike << endl;
			//p0 += res.p0 / sumk;
			p0 += res.p0 / (sumk + res.p0);
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
	Ve_S km = alignment(fa, wid, thr, mod, pb, maj, min);
	Ofs out(Out.data());
	out << "Selected K-mers:" << endl;
	for (int l = 0; l < km.size(); l++) 
		out << km[l].Consensus << endl;
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
		bool Check = true;
		for (int m = 1; m < mod; m++) {
			//cout << sub[m - 1] + m - 1 << "\t" << sub[m] + m << endl;
			if (sub[m - 1] + m - 1 >= sub[m] + m)
				Check = false;
		}
		if (Check) {
			ins++;
			cout << "Initialization seed: " << ins << endl;
			motif temp;
			temp.loglike = 0;
			temp.pb = pb;
			temp.pc = background(fa);
			temp.p0 = 1.0 / ((reg - wid + 1) * mod + 1);
			cout << "Initial Sequences:";
			Ve_Str ini;
			for (int m = 0; m < mod; m++) {
				temp.rwm.push_back(km[sub[m] + m].mat);
				//temp.rwm.push_back(pwm_assign(km[sub[m] + m], temp.pb, maj, min));
				temp.cons.push_back(km[sub[m] + m].Consensus);
				cout << "\t" << km[sub[m] + m].Consensus;
				ini.push_back(km[sub[m] + m].Consensus);
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

