#include "common.h"
using namespace std;

//Commonly used functions

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

char complement(char S, bool NUC) {//NUC true means DNA, otherwise means RNA.
	if ((S == 'A' || S == 'a') && NUC)
		return 'T';
	else if ((S == 'A' || S == 'a') && !NUC)
		return 'U';
	else if (S == 'C' || S == 'c')
		return 'G';
	else if (S == 'G' || S == 'g')
		return 'C';
	else if (S == 'T' || S == 'U' || S == 't' || S == 'u')
		return 'A';
	else
		return 'N';
}

Str antisense(Str Forward, bool NUC) {
	Str Res = Forward;
	for (int i = 0; i < Forward.size(); i++) {
		int l = Forward.length() - 1;
		Res[i] = complement(Forward[l - i], NUC);
	}
	return Res;
}

Str transcribe(Str DNA) {
	Str Res = DNA;
	for (int i = 0; i < DNA.size(); i++) {
		if (Res[i] == 'T')
			Res[i] = 'U';
	}
	return Res;
}

Str cleave(Str Seq, Str Strand, int start, int end, bool NUC) {//By default two ends are included!
	Str Res;
	int len = Seq.length();
	//cout << Seq.length() + start - 1 << "\t" << end - start << "\t" << start << "\t" << end << "\t" << Seq.length() << "\t" << (start > len) << endl;
	if (start > end) {
		//cout << 1 << endl;
		Res = Seq.substr(start - 1, Seq.length() - start + 1);
		Res += Seq.substr(0, end);
	}
	else if (start > len) {
		//cout << 2 << endl;
		Res = Seq.substr(start - Seq.length() - 1, end - Seq.length());
	}
	else if (end <= 0) {
		//cout << 3 << endl;
		Res = Seq.substr(Seq.length() + start - 1, end - start + 1);
	}
	else if (start <= 0) {
		//cout << 4 << endl;
		Res = Seq.substr(Seq.length() + start - 1, 1 - start);
		Res += Seq.substr(0, end);
	}
	else if (end > len) {
		//cout << 5 << endl;
		Res = Seq.substr(start - 1, Seq.length() - start + 1);
		Res += Seq.substr(0, end - Seq.length());
	}
	else {
		//cout << 6 << endl;
		Res = Seq.substr(start - 1, end - start + 1);
	}
	//cout<< Res << endl;
	if (Strand == "+" || Strand == "forward" || Strand == "Forward") {
		if (!NUC)
			Res = transcribe(Res);
		return Res;
	}
	else if (Strand == "-" || Strand == "reverse" || Strand == "Reverse")
		return antisense(Res, NUC);
	else
		return "Error!";
}

Pa_I_I locate(Ve_Str order, int loc, Ma_Str_R red, Str Strand, int len, bool UTR) {
	Pa_I_I res;
	res.first = order.size() / 2;
	int upper = 0;
	int lower = order.size() - 1;
	Str Synonym = order[res.first];
	while (abs(upper - lower) > 1) {
		Synonym = order[res.first];
		if ((Strand == "+" && UTR) || (Strand == "-" && !UTR)) {
			//cout << loc << "\t" << red[Synonym].start << endl;
			if (loc - red[Synonym].start == 0)
				break;
			else {
				if (loc - red[order[res.first]].start > 0)
					upper = res.first;
				else
					lower = res.first;
				res.first = (upper + lower) / 2;
			}
		}
		else if ((Strand == "+" && !UTR) || (Strand == "-" && UTR)) {
			//cout << loc << "\t" << ano[Synonym].start << endl;
			if (loc - red[Synonym].end == 0)
				break;
			else {
				if (red[Synonym].end - loc > 0)
					lower = res.first;
				else
					upper = res.first;
				res.first = (upper + lower + 1) / 2;
			}
		}
	}
	Synonym = order[res.first];
	if ((Strand == "+" && UTR) || (Strand == "-" && !UTR)) {
		double cov = abs(loc - red[order[res.first]].end) * 1.0 / (red[order[res.first]].end - red[order[res.first]].start);
		if ((loc - red[order[res.first]].start > 0 && cov < 0.5 && abs(loc - red[order[res.first]].start) > abs(loc - red[order[lower]].start))
			|| loc - red[order[res.first]].end > 0)
			res.first = lower;
		res.second = loc - red[order[res.first]].start;
		if (res.first == order.size() - 1 && abs(loc - len - red[order[0]].start) < abs(res.second)) {
			res.first = 0;
			res.second = loc - len - red[order[0]].start;
			cout << loc << "\t" << len << "\t" << red[order[0]].start << "\t" << res.second << endl;
		}
	}
	if ((Strand == "+" && !UTR) || (Strand == "-" && UTR)) {
		int n = order.size() - 1;
		double cov = abs(loc - red[order[res.first]].start) * 1.0 / (red[order[res.first]].end - red[order[res.first]].start);		
		if ((red[order[res.first]].end - loc > 0 && cov < 0.5 && abs(red[order[res.first]].end - loc) > abs(red[order[upper]].end - loc))
			|| red[order[res.first]].start - loc > 0)
			res.first = upper;
		res.second = red[order[res.first]].end - loc;
		if (res.first == 0 && abs(loc + len - red[order[n]].end) < abs(res.second)) {
			res.first = n;
			res.second = red[order[n]].end - len - loc;
		}
	}
	//cout << Strand << "\t" << Synonym << "\t" << upper << "\t" << lower << "\t" << res.first << "\t" << loc << "\t" << ano[order[res]].end << endl;
	return res;
}

//Data loading functions

Ma_Str_Str read_in_codon(Str In) {
	Ma_Str_Str result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			Str Key, Cont;
			in >> Key >> Cont;
			getline(in, Line);
			Cont += Line;
			if (!Key.empty())
				result[Key] = Cont;
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Str read_in_fna(Str In) {
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
					if (!Label.empty())
						result[Label] = Seq;
					Seq.clear();
					Label.clear();
					//int posb = Line.find("|ref|");
					//int pos = Line.find("lasmid");
					//if (pos == Str::npos) {
					int pose = Line.find(".");
					Label = Line.substr(1, pose - 1);
					//cout << Label << endl;
					//}
				}
				else
					Seq += Line;
			}
		}
		if (!Label.empty())
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
					if (!Label.empty())
						result[Label] = Seq;
					Seq.clear();
					Label = Line.substr(1, Line.length() - 1);
				}
				else
					Seq += Line;
			}
		}
		if (!Label.empty())
			result[Label] = Seq;
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

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
				//for (int i = 0; i < temp.size(); i++)
					//cout << "\t" << temp[i];
				//cout << endl;
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
					/*if (cur.Strand == "+")
						result.forward[cur.GenoID].push_back(cur.Synonym);
					else if (cur.Strand == "-")
						result.reverse[cur.GenoID].push_back(cur.Synonym);*/
				}
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_K read_in_kegg(Str In) {
	Ma_Str_K result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Domain, Kingdom, Phylum, Genus;
		//bool UN = false;
		while (!in.eof()) {
			Str Probe, Line;
			in >> Probe;
			if (Probe == "A") {
				in >> Domain;
				getline(in, Line);
			}
			else if (Probe == "B") {
				in >> Kingdom;
				getline(in, Line);
			}
			else if (Probe == "C") {
				in >> Phylum;
				getline(in, Line);
				int pos = Line.find("(");
				Phylum += Line.substr(0, pos - 1);
			}
			else if (Probe == "D") {
				in >> Line;
				if (Genus != "unclassified" || Genus != "Unclassified")
					Genus = Line;
				result[Genus].Domain = Domain;
				result[Genus].Kingdom = Kingdom;
				result[Genus].Phylum = Phylum;
				result[Genus].Genus = Genus;
				getline(in, Line);
			}
			else
				getline(in, Line);
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ve_Str read_in_column(Str In) {
	Ve_Str res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty())
				res.push_back(Line);
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ma_I_M read_in_med(Str In) {
	Ma_I_M result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		return result;
	}
	else {
		Str Line;
		while (!in.eof()) {
			med tmp;
			double pts = 0;
			in >> tmp.p5 >> tmp.p3 >> tmp.Strand;
			getline(in, Line);
			Ve_Str temp = string_parse(Line, "       ");
			if (temp.size() > 2)
				convertFromString(pts, temp[2]);
			if (!tmp.Strand.empty() && pts >= 0.70) {
				if (tmp.Strand == "+")
					result[tmp.p3] = tmp;
				else {
					int a = tmp.p3;
					tmp.p3 = tmp.p5;
					tmp.p5 = a;
					result[tmp.p3] = tmp;
				}
			}
			//else
			//	cout << tmp.p5 << "\t" << tmp.p3 << "\t" << tmp.p5 - tmp.p3 << "\t" << tmp.Strand << "\t" << pts << endl;
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Str read_in_trans(Str In) {
	Ma_Str_Str res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			Str Fna, Ftt;
			in >> Fna >> Ftt;
			getline(in,Line);
			if (!Fna.empty()) {
				res[Fna] = Ftt;
				res[Ftt] = Fna;
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

Ma_Str_L read_in_blast(Str In) {
	Ma_Str_L res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Probe;
		while (!in.eof()) {
			in >> Probe;
			if (!Probe.empty()) {
				if (Probe[0] == '#')
					getline(in, Line);
				else {
					blast tmp;
					in >> tmp.Ref >> tmp.id >> tmp.len >> tmp.mis >> tmp.gap >> tmp.qst >> tmp.qed >> tmp.rst >> tmp.red >> tmp.eval >> tmp.bit;
					tmp.Que = Probe;
					Str CDS = "_cds_";
					int posb = tmp.Que.find(CDS); 
					if (posb != Str::npos) {
						int pose = tmp.Que.rfind("_");
						Str Label = tmp.Que.substr(posb + CDS.length(), pose - posb - CDS.length());
						//cout << Label << endl;
						if(res[Label].Que.empty() || res[Label].id < tmp.id)
							res[Label] = tmp;
					}
					else {
						Ve_Str tem = string_parse(tmp.Que, "|");
						if(res[tem[2]].Que.empty() || res[tem[2]].id < tmp.id)
							res[tem[2]] = tmp;
					}
				}
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ma_Str_Str read_in_cog_csv(Str In) {
	Ma_Str_Str res;
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
				Ve_Str tmp = string_parse(Line, ",");
				res[tmp[0]] = tmp[6];
			}
		}
	}
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
				if (red.Strand == "+")
					res.forward[red.Replicon].push_back(tmp[2]);
				else
					res.reverse[red.Replicon].push_back(tmp[2]);
				res.record[tmp[2]] = red;
				//cout << tmp[2] << "\t" << res.record[tmp[2]].Product << endl;
				res.order.push_back(tmp[2]);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ve_E read_in_exp(Str In) {
	Ve_E res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		Ve_Str label = string_parse(Line, "\t");
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				ex tp;
				Ve_Str tmp = string_parse(Line, "\t");
				tp.tts = -1;
				for (int i = 0; i < label.size(); i++) {
					if (label[i] == "Genomic Accession")
						tp.Replicon = tmp[i];
					else if (label[i] == "Strand")
						tp.Strand = tmp[i];
					else if (label[i] == "Phase") {
						Ve_Str tem = string_parse(tmp[i], ";");
						for (int j = 0; j < tem.size(); j++)
							tp.phase.insert(tem[j]);
					}
					else if (label[i] == "Reference"){
						Ve_Str tem = string_parse(tmp[i], ";");
						for (int j = 0; j < tem.size(); j++)
							tp.ref.insert(tem[j]);
					}
					else if (label[i] == "TSS") {
						int a = 0;
						Ve_Str tss = string_parse(tmp[i], ";");
						for (int i = 0; i < tss.size(); i++) {
							convertFromString(a, tss[i]);
							tp.tss.insert(a);
						}
					}
					else if (label[i] == "TTS" && !tmp[i].empty())
						convertFromString(tp.tts, tmp[i]);
				}
				res.push_back(tp);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ma_Str_T read_in_tu(Str In) {
	Ma_Str_T res;
	//First is operon ID, second is label for start and end.
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		int id = 0;
		Ma_Str_Se_Str opr;
		getline(in, Line);
		while (!in.eof()) {
			Str Strand, Genes;
			in >> Strand >> Genes;
			getline(in, Line);
			if (!Strand.empty()) {
				Str ID;
				Ve_Str genes = string_parse(Genes, ";");
				Se_Str tmp;
				for (int i = 0; i < genes.size(); i++) {
					if (!genes[i].empty()) {
						if (!res[genes[i]].ID.empty()) {
							for (Se_Str::const_iterator iter = opr[res[genes[i]].ID].begin(); iter != opr[res[genes[i]].ID].end(); iter++)
								tmp.insert(*iter);
						}
						else
							tmp.insert(genes[i]);
						if (i == 0 && Strand == "+" || i == genes.size() - 1 && Strand == "-" || i == 0 && Strand == "&") {
							if (res[genes[i]].TSS.empty())
								res[genes[i]].TSS = "Operonic TSS";
							else if (res[genes[i]].TSS != "Operonic TSS")
								res[genes[i]].TSS = "Sub-operonic TSS";
						}
						else {
							if (res[genes[i]].TSS.empty())
								res[genes[i]].TSS = "None";
							else if (res[genes[i]].TSS == "Operonic TSS")
								res[genes[i]].TSS = "Sub-operonic TSS";
						}
						if (i == 0 && Strand == "-" || i == genes.size() - 1 && Strand == "+" || i == genes.size() - 1 && Strand == "&") {
							if (res[genes[i]].TTS.empty())
								res[genes[i]].TTS = "Operonic TTS";
							else if (res[genes[i]].TTS != "Operonic TTS")
								res[genes[i]].TTS = "Sub-operonic TTS";
						}
						else {
							if (res[genes[i]].TTS.empty())
								res[genes[i]].TTS = "None";
							else if (res[genes[i]].TTS == "Operonic TTS")
								res[genes[i]].TTS = "Sub-operonic TTS";
						}
					}
				}
				id++;
				convertFromNumber(ID, id);
				for (Se_Str::const_iterator ite = tmp.begin(); ite != tmp.end(); ite++) {
					res[*ite].ID = ID;
					opr[ID] = tmp;
				}
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

/*Ma_Str_E read_in_trans_record(Str In) {
	Ma_Str_E res;
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
				int a, b, c;
				a = b = c = 0;
				convertFromString(a, tmp[3]);
				convertFromString(b, tmp[4]);
				convertFromString(c, tmp[5]);
				res[tmp[2]].tss.insert(a);
				res[tmp[2]].utr.insert(b);
				res[tmp[2]].tts.insert(c);
				Ve_Str tp = string_parse(tmp[7], ";");
				for (int i = 0; i < tp.size(); i++)
					res[tmp[2]].phase.insert(tp[i]);
				tp = string_parse(tmp[8], ";");
				for (int i = 0; i < tp.size(); i++)
					res[tmp[2]].ref.insert(tp[i]);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}*/

Ma_Str_Str read_in_exftt(Str In) {
	Ma_Str_Str res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Old1, Old2;
		getline(in, Line);
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				Ve_Str temp = string_parse(Line, "\t");
				if (temp[0] != "gene") {
					Str Strand = temp[9];
					Str Synonym = temp[16];
					if (!Old1.empty()) 
						res[Old1] = Synonym;
					if (!Old2.empty())
						res[Old2] = Synonym;
					res[Synonym] = Strand;
				}
				else {
					Old1 = temp[temp.size() - 1];
					Str Target = "old_locus_tag=";
					int pos = Old1.find(Target);
					if (pos != Str::npos) {
						Old1 = Old1.substr(pos + Target.length(), Old1.length() - pos - Target.length());
						Ve_Str tmp = string_parse(Old1, ",");
						Old1 = tmp[0];
						if(tmp.size() > 1)
							Old2 = tmp[1];
					}
				}
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ve_Str read_in_tu_genes(Str In) {
	Ve_Str res;
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
			if(!Line.empty())
				res.push_back(Line);
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ve_Pa_Str_I read_in_opr(Str In) {
	Ve_Pa_Str_I res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line;
		getline(in, Line);
		while (!in.eof()) {
			int op = 0;
			Str Gene;
			in >> op >> Gene;
			if (!Gene.empty())
				res.push_back(make_pair(Gene, op));
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}
