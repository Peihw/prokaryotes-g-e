#include "common.h"
using namespace std;

//Commonly Used Functions

Ve_Str string_parse(con_Str Line, con_Str Key) {
	Ve_Str result;
	int posB = 0;
	int posE = (int)Line.find(Key);
	while (posE != Str::npos) {
		Str Temp = Line.substr(posB, posE - posB);
		posB = posE + (int)Key.size();
		posE = (int)Line.find(Key, posB);
		result.push_back(Temp);
		//cout<<Temp<<"\t";
	}
	Str Temp = Line.substr(posB, posE - posB);
	//cout<<Temp<<endl;
	result.push_back(Temp);
	return result;
}

De_I subscript(long long int i, Ve_I limit) {
	De_I result;
	int final = limit.size() - 1;
	long long int prod = 1;
	for (int j = final; j >= 0; j--) {
		if (!j)
			result.push_front(i / prod);
		else {
			if (limit[j]) {
				result.push_front((i / prod) % limit[j]);
				prod *= limit[j];
			}
			else
				result.push_front(0);
			//cout<<i<<"\t"<<((i / prod) % limit[j])<<"\t"<<limit.size()<<endl;
		}
	}
	return result;
}

//Data Loading Function
ftt read_in_ftt(Str In) {
	ftt result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str temp = string_parse(Line, "\t");
				if (temp[0] != "gene") {
					feature cur;
					cur.Feature = temp[0];
					cur.Class = temp[1];
					result.AssemID = temp[2];
					cur.Unit = temp[3];
					cur.Seq_type = temp[4];
					int pos = temp[5].find(".");
					cur.Chromosome = temp[5].substr(0, pos);
					pos = temp[6].find(".");
					cur.GenoID = temp[6].substr(0, pos);
					result.Chr.insert(cur.GenoID);
					//cout << cur.GenoID << endl;
					convertFromString(cur.start, temp[7]);
					convertFromString(cur.end, temp[8]);
					cur.Strand = temp[9];
					cur.ProdID = temp[10];
					cur.Nr = temp[11];
					cur.RelID = temp[12];
					cur.Product = temp[13];
					cur.Symbol = temp[14];
					cur.GeneID = temp[15];
					cur.Synonym = temp[16];
					convertFromString(cur.interval, temp[17]);
					convertFromString(cur.product, temp[18]);
					cur.Attributes = temp[19];
					result.S2F[cur.Synonym] = cur;
					result.order.push_back(cur.Synonym);
				}
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Str read_in_fas(Str In) {
	Ma_Str_Str result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Label, Seq;
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				if (Line[0] == '>') {
					if (!Seq.empty() && Seq != "Error!" && !Label.empty())
						result[Label] = Seq;
					//cout << Label << "\t" << Seq << endl;
					Seq.clear();
					Ve_Str temp = string_parse(Line, "|");
					Label = temp[2];
				}
				else
					Seq += Line;
			}
		}
		if (!Seq.empty() && Seq != "Error!") 
			result[Label] = Seq;
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Str read_in_fa(Str In) {
	Ma_Str_Str result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Label, Seq;
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				if (Line[0] == '>') {
					int pos = Seq.find("N");
					if (!Seq.empty() && Seq != "Error!" && !Label.empty() && pos == Str::npos) 
						result[Label] = Seq;
					Seq.clear();
					Label = Line.substr(1, Line.length() - 1);
				}
				else
					Seq += Line;
			}
		}
		if (!Seq.empty() && Seq != "Error!")
			result[Label] = Seq;
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_M read_in_motif(Str In) {
	Ma_Str_M result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		bool TIS = false;
		bool PWM = false;
		Ve_Ma_C_D wm;
		motif curr;
		Str Type;
		while (!in.eof()) {
			Str Line;
			getline(in, Line);
			if (!Line.empty()) {
                int post = Line.find("Signal type");
				int posg = Line.find("Genome");
				int posb = Line.find("background");
				int posn = Line.find("No signal");
				int poss = Line.find("TIS start");
				int posw = Line.find("PWM");
				int posm = Line.find("Motif");
				if (post != Str::npos || posg != Str::npos) {
					if(!Type.empty())
						result[Type] = curr;
					getline(in, Type);
					curr.pb.clear();
					curr.pj.clear();
					curr.rwm.clear();
					curr.cons.clear();
					curr.entropy.clear();
				}
				else if (posb != Str::npos) {
					getline(in, Line);
					double a, c, g, t;
					a = c = g = t = 0;
					in >> a >> c >> g >> t;
					getline(in, Line);
					curr.pb['A'] = a;
					curr.pb['C'] = c;
					curr.pb['G'] = g;
					curr.pb['T'] = t;
				}
				else if (posn != Str::npos) {
					in >> curr.p0;
					getline(in, Line);
				}
				else if (poss != Str::npos) 
					TIS = true;
				else if (posw != Str::npos) {
					PWM = true;
					TIS = false;
					if (!wm.empty()) {
						curr.rwm.push_back(wm);
						wm.clear();
					}
				}
				else if (TIS && posm != Str::npos) {
					getline(in, Line);
					Ve_Str temp = string_parse(Line,"\t");
					Ve_D tem;
					for (int i = 1; i < temp.size(); i++) {
						double pj = 0;
						convertFromString(pj, temp[i]);
						tem.push_back(pj);
					}
					curr.pj.push_back(tem);
				}
				else if (PWM) {
					//getline(in, Line);
					int pose = Line.find("entropy");
					int posc = Line.find("Consensus");
					if (pose != Str::npos) {
						double entropy = 0;
						Ve_Str tmp = string_parse(Line, "\t");
						if (tmp[1] != "-nan(ind)")
							convertFromString(entropy, tmp[1]);
						else
							entropy = -100;
						curr.entropy.push_back(entropy);
                        PWM = false;
					}
					else if (posc != Str::npos) {
						curr.rwm.push_back(wm);
						wm.clear();
                        Ve_Str tmp = string_parse(Line, "\t");
						curr.cons.push_back(tmp[1]);
					}
					else {
						Ve_Str temp = string_parse(Line, "\t");
						Ma_C_D tem;
						for (int i = 1; i < temp.size(); i++) {
							double wi = 0;
							convertFromString(wi, temp[i]);
							tem[Base[i - 1]] = wi;
							//cout << "\t" << wi << "\t" << tem.size();
						}
						wm.push_back(tem);
						//cout << "\t" << wm.size() << endl;
					}
				}
			}
		}
        result[Type] = curr;
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Ve_Str read_in_distance(Str In) {
	Ma_Str_Ve_Str result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		while (!in.eof()) {
			Str Line;
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str tmp = string_parse(Line, "\t");
				result[tmp[0]].push_back(tmp[3]);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Ve_L read_in_blast(Str In) {
	Ma_Str_Ve_L result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		return result;
	}
	else {
		while (!in.eof()) {
			Str Line, Temp;
			in >> Temp;
			if (Temp == "#")
				getline(in, Line);
			else {
				if (!Temp.empty()) {
					blast temp;
					in >> temp.SID >> temp.id >> temp.len >> temp.mis >> temp.gap >> temp.qstart >> temp.qend
						>> temp.sstart >> temp.send >> temp.evalue >> temp.bit;
					temp.QID = Temp;
					result[temp.QID].push_back(temp);
				}
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Ma_Str_O read_in_ortholog(Str In) {
	Ma_Str_Ma_Str_O result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			ortholog temp;
			in >> temp.OID >> temp.QID >> temp.o2qi >> temp.q2oi >> temp.olen >> temp.qlen >> temp.o2qe >> temp.q2oe;
			getline(in, Line);
			if (!temp.OID.empty()) {
				Ve_Str tmo = string_parse(temp.OID, "|");
				Ve_Str tmq = string_parse(temp.QID, "|");
				result[tmq[1]][tmo[1]] = temp;
				//cout << tmq[1] << "\t" << tmo[1] << "\t" << result[tmq[1]][tmo[1]].OID << endl;
				//cout << temp.QID << "\t" << temp.OID << endl;
				//if (tmq[1] == "JK_RS00005")
				//	cout << tmq[1] << "|\t" << tmo[1] << "|\t" << result[tmq[1]][tmo[1]].OID << endl;
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

/*Ma_Str_R read_in_record_v2(Str In) {
	Ma_Str_R result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		while (!in.eof()) {
			Str Line;
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str temp = string_parse(Line, "\t");
				record rec;
				rec.Strand = temp[6];
				rec.Synonym = temp[1];
				convertFromString(rec.start, temp[4]);
				convertFromString(rec.end, temp[5]);
				if (temp[7] == "TA")
					rec.Type = "Leaderless";
				else if (temp[7] == "NA")
					rec.Type = "None";
				else if (temp[7] == "SD")
					rec.Type = "Cannonical";
				else if (temp[7] == "AT")
					rec.Type = "Atypical";
				convertFromString(rec.prob, temp[10]);
				convertFromString(rec.tis, temp[9]);
				rec.Seq = temp[8];
				rec.Produt = temp[11];
				result[temp[1]] = rec;
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}*/

Ma_C_D background_probability(Str In) {
	Ma_C_D res;
	Str Seq;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty() && Line[0] != '>') 
				Seq += Line;
		}
	}
	for (int i = 0; i < Seq.size(); i++) 
		res[Seq[i]] ++;
	res['C'] = res['G'] = (res['C'] + res['G']) / Seq.length() / 2.0;
	res['A'] = res['T'] = (res['A'] + res['T']) / Seq.length() / 2.0;
	//cout << "A:\t" << res['A'] << "\tC:\t" << res['C'] << "\tG:\t" << res['G'] << "\tT:\t" << res['T'] << endl;
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

rec read_in_record(Str In, Str Geno) {
	rec res;
	res.Genome = Geno;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str tmp = string_parse(Line, "\t");
				record red;
				red.opr = 0;
				red.Replicon = tmp[0];
				red.Strand = tmp[1];
				red.Name = tmp[3];
				red.Product = tmp[4];
				convertFromString(red.start, tmp[5]);
				convertFromString(red.end, tmp[6]);
				red.TISRef = tmp[7];
				if (tmp[8] != "NA")
					convertFromString(red.shift, tmp[8]);
				else
					red.shift = 0;
				red.Feature = tmp[9];
				if (tmp.size() > 10) {
					//TSS
					if (tmp[10] != "NA") {
						Ve_Str tp = string_parse(tmp[10], ";");
						for (int i = 0; i < tp.size(); i++) {
							Ve_Str tem = string_parse(tp[i], "|");
							int tss = 0;
							convertFromString(tss, tem[0]);
							convertFromString(red.tss[tss].utr, tem[1]);
							Ve_Str te = string_parse(tem[2], "&");
							for (int j = 0; j < te.size(); j++)
								red.tss[tss].phase.insert(te[j]);
							te = string_parse(tem[3], "&");
							for (int j = 0; j < te.size(); j++)
								red.tss[tss].ref.insert(te[j]);
						}
					}
					if (tmp[11] != "NA") {
						Ve_Str tp = string_parse(tmp[11], ";");
						for (int i = 0; i < tp.size(); i++) {
							Ve_Str tem = string_parse(tp[i], "|");
							int tts = 0;
							convertFromString(tts, tem[0]);
							Ve_Str te = string_parse(tem[1], "&");
							for (int j = 0; j < te.size(); j++)
								red.tts[tts].phase.insert(te[j]);
							te = string_parse(tem[2], "&");
							for (int j = 0; j < te.size(); j++)
								red.tts[tts].ref.insert(te[j]);
						}
					}
				}
				if(tmp.size() > 12){
					if (tmp[12] != "NA") {
						Ve_Str tp = string_parse(tmp[12], ";");
						for (int i = 0; i < tp.size(); i++) {
							Ve_Str te = string_parse(tp[i], "|");
							int pos, rel;
							pos = rel = 0;
							convertFromString(pos, te[0]);
							convertFromString(rel, te[1]);
							red.ldl[pos] = rel;
						}
					}
					if (tmp[13] != "NA")
						convertFromString(red.opr, tmp[13]);
					else
						red.opr = 0;
					red.TSS = tmp[14];
					red.TTS = tmp[15];
				}
				res.record[tmp[2]] = red;
				res.order.push_back(tmp[2]);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ma_Str_Str read_in_cog(Str In) {
	Ma_Str_Str res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Label, COG, Line;
		in >> Label >> Line >> COG;
		getline(in, Line);
		if (!Label.empty()) {
			Str Anch = "cds_";
			int posb = Label.find(Anch);
			if (posb != Str::npos) {
				int pose = Label.find("_", posb + Anch.length());
				Str Real = Label.substr(posb + Anch.length(), pose - posb - Anch.length());
				res[Real] = COG;
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}


Ma_Str_T read_in_exp_tss(Str In) {
	Ma_Str_T res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			//tss_exp tmp;
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str tmp = string_parse(Line, "\t");
				res[tmp[2]].Replicon = tmp[0];
				res[tmp[2]].Strand = tmp[1];
				res[tmp[2]].Synonym = tmp[2];
				int a = 0;
				int b = 0;
				convertFromString(a, tmp[3]);
				convertFromString(b, tmp[4]);
				res[tmp[2]].tss.push_back(make_pair(a, b));
				Ve_Str tp = string_parse(tmp[6], ";");
				for(int i = 0; i < tp.size(); i++)
					res[tmp[2]].tssref.insert(tp[i]);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ma_Str_Pa_Str_Str read_in_rps1(Str In) {
	Ma_Str_Pa_Str_Str res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str tmp = string_parse(Line, "\t");
				res[tmp[0]] = make_pair(tmp[1], tmp[3]);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ma_Str_X read_in_taxonomy(Str In) {
	Ma_Str_X res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			taxonomy tmp;
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str tp = string_parse(Line, "\t");
				tmp.AssemID = tp[0];
				tmp.Domain = tp[1];
				int posF = tp[2].find("Firmicutes");
				int posP = tp[2].find("roteobacteria");
				if (posF != Str::npos) {
					tmp.Phylum = "Firmicutes";
					if (tp[3] == "Streptococcus" || tp[3] == "Lactococcus")
						tmp.Sub = "Streptococcaceae";
					else if (tp[2] == "Firmicutes - Bacilli")
						tmp.Sub = "Bacilli-Others";
					else if (tp[2] == "Firmicutes - Clostridia")
						tmp.Sub = "Clostridia";
					else
						tmp.Sub = "Firmicutes-Others";
				}
				else if (posP != Str::npos) {
					tmp.Phylum = "Proteobacteria";
					if (tp[3] == "Anaplasma" || tp[3] == "Ehrlichia" || tp[3] == "Midichloria" || tp[3] == "Neorickettsia"
						|| tp[3] == "Orientia" || tp[3] == "Rickettsia" || tp[3] == "Wolbachia")
						tmp.Sub = "Rickettsiales";
					else if (tp[2] == "Alphaproteobacteria")
						tmp.Sub = "Alphaproteobacteria-Others";
					else if (tp[3] == "Beggiatoa" || tp[3] == "Cycloclasticus" || tp[3] == "Francisella"
						|| tp[3] == "Methylophaga" || tp[3] == "Piscirickettsia" || tp[3] == "Thiomicrospira"
						|| tp[3] == "Thioploca")
						tmp.Sub = "Thiotrichales";
					else if (tp[3] == "Buchnera")
						tmp.Sub = "Buchnera";
					else if (tp[2] == "Gammaproteobacteria - Others")
						tmp.Sub = "Gammaproteobacteria-Others";
					else if (tp[2] == "Gammaproteobacteria - Enterobacteria")
						tmp.Sub = "Enterobacteriales-Others";
					else if (tp[2] == "Unclassified Proteobacteria")
						tmp.Sub = "Proteobacteria-Unclassified";
					else
						tmp.Sub = tp[2];
				}
				else if (tp[2] == "Unclassified Terrabacteria group") {
					tmp.Phylum = "Terrabacteria-Unclassified";
					tmp.Sub = "NA";
				}
				else if (tp[2] == "Unclassified Bacteria") {
					tmp.Phylum = "Bacteria-Unclassified";
					tmp.Sub = "NA";
				}
				else if (tp[2] == "Unclassified Archaea") {
					tmp.Phylum = "Archaea-Unclassified";
					tmp.Sub = "NA";
				}
				else {
					tmp.Phylum = tp[2];
					tmp.Sub = "NA";
				}
				tmp.Genus = tp[3];
				tmp.Species = tp[4];
				res[tmp.AssemID] = tmp;
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}
