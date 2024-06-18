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
	}
	Str Temp = Line.substr(posB, posE - posB);
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
		}
	}
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

Str translate(Str Cds, Ma_Str_Str cod) {
	Str Res = "";
	int c = 3;
	for (int j = 0; j < Cds.length() / c; j++) {
		Str Codon = Cds.substr(j * c, c);
		if (!cod[Codon].empty()) {
			if (Codon == "UGA" && j != Cds.length() / c - 1)
				Res += "W";
			else if (j == Cds.length() / c - 1)
				continue;
			else
				Res += cod[Codon];
		}
		//else
		//	Res += "X";
	}
	Res[0] = 'M';
	return Res;
}

Str transverse(Str Forward) {
	Str Res = Forward;
	for (int i = 0; i < Forward.size(); i++) {
		int l = Forward.length() - 1;
		Res[i] = Forward[l - i];
	}
	return Res;
}

Str cleave(Str Seq, Str Strand, int start, int end, bool NUC) {//By default two ends are included!
	Str Res;
	int len = Seq.length();
	//cout << Seq.length() + start - 1 << "\t" << end - start << "\t" << start << "\t" << end << "\t" << Seq.length() << "\t" << (start > len) << endl;
	if (start <= 0)
		start += Seq.length();
	if (end <= 0)
		end += Seq.length();
	if (start > Seq.length())
		start -= Seq.length();
	if (end > Seq.length())
		end -= Seq.length();
	if (start > end) {
		//cout << 1 << endl;
		Res = Seq.substr(start - 1, Seq.length() - start + 1);
		Res += Seq.substr(0, end);
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

//File Loading Functions

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
					//int posb = Line.find("|ref|");
					int pose = Line.find(".");
					Label = Line.substr(1, pose - 1);
				}
				else
					Seq += Line;
			}
		}
		result[Label] = Seq;
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Str read_in_fa(Str In, int num) {
	Ma_Str_Str res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Label;
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				if (Line[0] == '>') {
					Label = Line.substr(1, Line.length() - 1);
					//Ve_Str temp = string_parse(Label, "|");
					//Label = temp[num - 1];
				}
				else
					res[Label] = Line;
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
			//cout << Line << endl;
			if (!Line.empty()) {
				record red;
				Ve_Str tmp = string_parse(Line, "\t");
				red.Replicon = tmp[0];
				red.Strand = tmp[1];
				red.Name = tmp[3];
				red.Product = tmp[4];
				convertFromString(red.start, tmp[5]);
				convertFromString(red.end, tmp[6]);
				red.TISRef = tmp[7];
				if (red.Strand == "+")
					red.tis = red.start;
				else
					red.tis = red.end;
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
				if(tmp.size() > 16){
					red.Signal = tmp[16];
					red.Type = tmp[17];
					if (tmp[18] != "NA")
						convertFromString(red.prob, tmp[18]);
					else
						red.prob = 0;
					if (tmp[19] != "NA")
						convertFromString(red.sig, tmp[19]);
					else
						red.sig = 0;
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
				//cout << Line << "\t" << temp.size() << endl;
				if (temp[0] != "gene" && !temp[0].empty()) {
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

Ma_Str_Ve_Str read_in_keywords(Str In) {
	Ma_Str_Ve_Str result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Label;
		while (!in.eof()) {
			getline(in, Line);
			if (!Line.empty()) {
				if (Line[0] == '>')
					Label = Line.substr(1, Line.length() - 1);
				else
					result[Label].push_back(Line);
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

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

Ma_Str_Str read_in_label(Str In) {
	Ma_Str_Str result;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Label;
		getline(in, Line);
		while (!in.eof()) {
			in >> Label >> Line;
			//cout<<Label<<"\t"<<Line<<endl;
			if (!Label.empty()) {
				result[Label] = Line;
				result[Line] = Label;
			}
			getline(in, Line);
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return result;
}

Ma_Str_Str read_in_rfa(Str In, Str Geno) {
	Ma_Str_Str res;
	Ifs in(In.data());
	if (!in) {
		cout << "Can't open " + In + "!" << endl;
		exit(0);
	}
	else {
		Str Line, Label;
		while (!in.eof()) {
			getline(in, Line);
			if (Line[0] == '>') {
				int pos = Line.find(Geno);
				if(pos != Str::npos)
					Label = Line.substr(1, Line.length() - 1);
				else
					Label.clear();
			}
			else {
				if(!Label.empty())
					res[Line] = Label;
			}
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Se_Str read_in_list(Str In) {
	Se_Str res;
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
				res.insert(Line);
		}
	}
	in.close();
	cout << In + " has been loaded!" << endl;
	return res;
}

Ma_Str_T read_in_taxonomy(Str In) {
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
						tmp.Sub = "Bacilli-Streptococcaceae";
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
						tmp.Sub = "Alphaproteobacteria-Rickettsiales";
					else if (tp[2] == "Alphaproteobacteria")
						tmp.Sub = "Alphaproteobacteria-Others";
					else if (tp[3] == "Beggiatoa" || tp[3] == "Cycloclasticus" || tp[3] == "Francisella"
						|| tp[3] == "Methylophaga" || tp[3] == "Piscirickettsia" || tp[3] == "Thiomicrospira"
						|| tp[3] == "Thioploca")
						tmp.Sub = "Gammaproteobacteria-Thiotrichales";
					//else if (tp[3] == "Buchnera")
					//	tmp.Sub = "Buchnera";
					else if (tp[2] == "Gammaproteobacteria - Others")
						tmp.Sub = "Gammaproteobacteria-Others";
					else if (tp[2] == "Gammaproteobacteria - Enterobacteria")
						tmp.Sub = "Gammaproteobacteria-Enterobacteriales";
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

double gc_content(Ma_Str_Str fna) {
	double res = 0;
	Str All;
	for (Ma_Str_Str::const_iterator iter = fna.begin(); iter != fna.end(); iter++)
		All += iter->second;
	for (int i = 0; i < All.length(); i++) {
		if (All[i] == 'C' || All[i] == 'G')
			res += 1;
	}
	res /= All.length();
	return res;
}