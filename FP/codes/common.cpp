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
				tmp.Gram = tp[2];
				tmp.Phylum = tp[3];
				tmp.Class = tp[4];
				tmp.Order = tp[5];
				tmp.Family = tp[6];
				tmp.Genus = tp[7];
				tmp.Species = tp[8];
				res[tmp.AssemID] = tmp;
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
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
					if (!Type.empty())
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
					Ve_Str temp = string_parse(Line, "\t");
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
				if (tmp.size() > 12) {
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
				if (tmp.size() > 16) {
					red.Signal = tmp[16];
					red.Type = tmp[17];
					if (tmp[18] != "NA")
						convertFromString(red.prob, tmp[18]);
					else
						red.prob = -1;
					if (tmp[19] != "NA")
						convertFromString(red.sig, tmp[19]);
					else
						red.sig = 1;
					red.Seq = tmp[20];
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

Ma_Str_Ve_Str read_in_distance(Str In) {
	Ma_Str_Ve_Str res;
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
				res[tmp[0]].push_back(tmp[12]);
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
		while (!in.eof()) {
			Str Synonym, COG;
			in >> Synonym >> COG;
			if (!Synonym.empty())
				res[Synonym] = COG;
		}
	}
	cout << In + " has been loaded!" << endl;
	return res;
}

