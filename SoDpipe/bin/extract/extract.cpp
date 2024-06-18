#include "extract.h"
using namespace std;

typedef std::map<Str, std::set<Str> > Ma_Str_Se_Str;

Str gene(Str Out, Str Geno, Ma_Str_Ve_Str keyword, ftt feat, Ma_Str_Str fa, Ma_Str_Str cod, Ma_Str_Str label, int reg) {
	Out += Geno + "." + keyword["Range"][0] + ".fa";
	Ofs out(Out.data());
	Str Log = Out + ".log";
	Ofs log(Log.data(), ios_base::app);
	//cout << gene.size() << gene["Type"].size() << endl;
	bool EXT = false;
	for (int i = 0; i < feat.order.size(); i++) {
		bool FIND = true;
		bool FILT = false;
		int pos = feat.S2F[feat.order[i]].Feature.find("rotein_coding");
		if (feat.S2F[feat.order[i]].Feature == keyword["Type"][0] || ((keyword["Type"][0] == "CDS" || keyword["Type"][0] == "AA")
			&& (pos != Str::npos || feat.S2F[feat.order[i]].Feature == "CDS"))) {
			int count = 0;
			if (!keyword["Positive"].empty() || !keyword["Negative"].empty()) {
				for (int k = 0; k < keyword["Positive"].size(); k++) {
					if (feat.S2F[feat.order[i]].Product.find(keyword["Positive"][k]) == Str::npos)
						count++;
				}
				if (count * 1.0 / keyword["Positive"].size() > 0.5)
					FIND = false;
				for (int k = 0; k < keyword["Negative"].size(); k++) {
					if (feat.S2F[feat.order[i]].Product.find(keyword["Negative"][k]) != Str::npos)
						FILT = true;
				}
			}
		}
		else
			FIND = false;
		if (FIND && !FILT) {
			EXT = true;
			log << Geno << "\t" << feat.S2F[feat.order[i]].GenoID << "\t" << feat.order[i] << "\t" << feat.S2F[feat.order[i]].Feature
				<< "\t" << feat.S2F[feat.order[i]].Strand << "\t" << feat.S2F[feat.order[i]].start << "\t" << feat.S2F[feat.order[i]].end
				<< "\t" << feat.S2F[feat.order[i]].Product << endl;
			//cout << feat.S2F[feat.order[i]].GenoID << "\t" << fa[feat.S2F[feat.order[i]].GenoID] << endl;
			Str DNA = fa[feat.S2F[feat.order[i]].GenoID];
			if (DNA.empty())
				DNA = fa[label[feat.S2F[feat.order[i]].GenoID]];
			//cout << label.size() << endl;
			//cout << feat.S2F[feat.order[i]].GenoID << "\t" << label[feat.S2F[feat.order[i]].GenoID] << "\t" << DNA.length() << endl;
			if (!DNA.empty() && !feat.S2F[feat.order[i]].Strand.empty()) {
				out << ">" << Geno << "|" << feat.S2F[feat.order[i]].GenoID << "|" << feat.order[i] << "|" << keyword["Range"][0] << endl;
				int start, end;
				start = end = 0;
				if (keyword["Range"][0] == "Content") {
					Str Seq = cleave(DNA, feat.S2F[feat.order[i]].Strand, feat.S2F[feat.order[i]].start, feat.S2F[feat.order[i]].end, true);
					if (keyword["Type"][0] == "AA")
						out << translate(Seq, cod) << endl;
					else
						out << Seq << endl;
				}
				else {
					if (keyword["Range"][0] == "5'UTR") {
						if (feat.S2F[feat.order[i]].Strand == "+") {
							start = feat.S2F[feat.order[i]].start - reg;
							end = feat.S2F[feat.order[i]].start - 1;
						}
						else {
							start = feat.S2F[feat.order[i]].end + 1;
							end = feat.S2F[feat.order[i]].end + reg;
						}
					}
					else {
						if (feat.S2F[feat.order[i]].Strand == "+") {
							start = feat.S2F[feat.order[i]].end + 1;
							end = feat.S2F[feat.order[i]].end + reg;
						}
						else {
							start = feat.S2F[feat.order[i]].start - reg;
							end = feat.S2F[feat.order[i]].start - 1;

						}
					}
					//cout << start << "\t" << end << endl;
					Str Seq = cleave(DNA, feat.S2F[feat.order[i]].Strand, start, end, true);
					out << Seq << endl;
				}
			}
		}
	}
	if (!EXT)
		log << feat.AssemID << "\tNo Possible Result!" << endl;
	log.close();
	out.close();
	return "Selected sequences have been all extracted!";
}

Str utr(Str Out, Ma_Str_Ve_Str gene, ftt feat, Ma_Str_Str fna, int len, bool PRIMER, Ma_Str_Str label) {//true for 5' end and false for 3' end 
	Out += ".utr.fa";
	Ofs out(Out.data());
	for (int i = 0; i < feat.order.size(); i++) {
		bool FIND = true;
		bool FILT = false;
		for (Ma_Str_Ve_Str::const_iterator iter = gene.begin(); iter != gene.end(); iter++) {
			if (iter->first == "Positive") {
				for (int j = 0; j < iter->second.size(); j++) {
					if (feat.S2F[feat.order[i]].Product.find(iter->second[j]) == Str::npos)
						FIND = false;
				}
			}
			else {
				for (int j = 0; j < iter->second.size(); j++) {
					if (feat.S2F[feat.order[i]].Product.find(iter->second[j]) != Str::npos)
						FILT = true;
				}
			}
		}
		if (FIND && !FILT) {
			out << ">" << feat.AssemID << "|" << feat.S2F[feat.order[i]].GenoID << "|" << feat.order[i] << "|" << feat.S2F[feat.order[i]].Strand
				<< "|" << feat.S2F[feat.order[i]].start << ".." << feat.S2F[feat.order[i]].end << endl;
			int start, end;
			if (PRIMER) {
				start = feat.S2F[feat.order[i]].start - len;
				end = feat.S2F[feat.order[i]].start - 1;
			}
			else {
				start = feat.S2F[feat.order[i]].end + 1;
				end = feat.S2F[feat.order[i]].end + len;
			}
			Str Chr = fna[feat.S2F[feat.order[i]].GenoID];
			if (Chr.empty())
				Chr = fna[label[feat.S2F[feat.order[i]].GenoID]];
			Str Seq = cleave(fna[feat.S2F[feat.order[i]].GenoID], feat.S2F[feat.order[i]].Strand, start, end, true);
			out << Seq << endl;
		}
	}
	out.close();
	return "Untranslated region have been all extracted!";
}

Str regions(Str Out, Str Opt, Str Geno, rec red, Ma_Str_Str fna, Ma_Str_Str cod, int up, int down) {
	Str Output = Out + Geno + "." + Opt + ".fa";
	Str Exp = Out + Geno + "." + Opt + ".exp.fa";
	Str Ldl = Out + Geno + "." + Opt + ".ldl.fa";
	Ofs out(Output.data());
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		int sta = 0;
		int end = 0;
		Ve_Pa_I_I pos;
		if (red.record[Synonym].Strand == "+") {
			//cout << Synonym << "\t" << red.record[Synonym].tss.size() << "\t" << red.record[Synonym].Type << endl;
			if (Opt == "tis") {
				sta = red.record[Synonym].start - up;
				end = red.record[Synonym].start - 1;
				Pa_I_I tmp = make_pair(sta, end);
				pos.push_back(tmp);
			}
			else if (Opt == "stc") {
				sta = red.record[Synonym].start;
				end = red.record[Synonym].start + 2;
				Pa_I_I tmp = make_pair(sta, end);
				pos.push_back(tmp);
			}
			else {
				if (Opt == "mrna")
					end = red.record[Synonym].start + down - 1;
				else if (Opt == "urbs")
					end = red.record[Synonym].start + red.record[Synonym].sig - 1;
				if (red.record[Synonym].tss.empty()) {
					if (red.record[Synonym].Type != "TA-like" && red.record[Synonym].Type != "TANNTTTG-like")
						sta = red.record[Synonym].start - up;
					else
						sta = red.record[Synonym].start;
					Pa_I_I tmp = make_pair(sta, end);
					pos.push_back(tmp);
				}
				else {
					int min = -100;
					if (red.record[Synonym].Type != "TA-like" && red.record[Synonym].Type != "TANNTTTG-like") {
						for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
							if (iter->second.utr > min && iter->second.utr < red.record[Synonym].sig) {
								sta = iter->first;
								Pa_I_I tmp = make_pair(sta, end);
								pos.push_back(tmp);
							}
						}
						if (pos.empty()) {
							sta = red.record[Synonym].start - up;
							Pa_I_I tmp = make_pair(sta, end);
							pos.push_back(tmp);
						}
					}
					else {
						min = -5;
						for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
							if (iter->second.utr >= min && iter->second.utr <= 0) {
								sta = iter->first;
								Pa_I_I tmp = make_pair(sta, end);
								pos.push_back(tmp);
							}
						}
						if (pos.empty()) {
							sta = red.record[Synonym].start;
							Pa_I_I tmp = make_pair(sta, end);
							pos.push_back(tmp);
						}
					}
				}
			}
		}
		else {
			if (Opt == "tis") {
				sta = red.record[Synonym].end + 1;
				end = red.record[Synonym].end + up;
				//cout << Synonym << "\t" << sta << "\t" << end << "\t" << fna[red.record[Synonym].Replicon].length() << endl;
				Pa_I_I tmp = make_pair(sta, end);
				pos.push_back(tmp);
			}
			else if (Opt == "stc") {
				sta = red.record[Synonym].end - 2;
				end = red.record[Synonym].end;
				Pa_I_I tmp = make_pair(sta, end);
				pos.push_back(tmp);
			}
			else {
				if (Opt == "mrna")
					sta = red.record[Synonym].end - down + 1;
				else if (Opt == "urbs")
					sta = red.record[Synonym].end - red.record[Synonym].sig + 1;
				if (red.record[Synonym].tss.empty()) {
					if (red.record[Synonym].Type != "TA-like" && red.record[Synonym].Type != "TANNTTTG-like")
						end = red.record[Synonym].end + up;
					else
						end = red.record[Synonym].end;
					Pa_I_I tmp = make_pair(sta, end);
					pos.push_back(tmp);
				}
				else {
					int min = -100;
					if (red.record[Synonym].Type != "TA-like" && red.record[Synonym].Type != "TANNTTTG-like") {
						for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
							if (iter->second.utr > min && iter->second.utr < red.record[Synonym].sig) {
								end = iter->first;
								Pa_I_I tmp = make_pair(sta, end);
								pos.push_back(tmp);
							}
						}
						if (pos.empty()) {
							end = red.record[Synonym].end + up;
							Pa_I_I tmp = make_pair(sta, end);
							pos.push_back(tmp);
						}
					}
					else {
						min = -5;
						for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
							if (iter->second.utr > min && iter->second.utr < 0) {
								end = iter->first;
								Pa_I_I tmp = make_pair(sta, end);
								pos.push_back(tmp);
							}
						}
						if (pos.empty()) {
							end = red.record[Synonym].end;
							Pa_I_I tmp = make_pair(sta, end);
							pos.push_back(tmp);
						}
					}
				}
			}
		}
		if (Opt == "cds" || Opt == "gene") {
			sta = red.record[Synonym].start;
			end = red.record[Synonym].end;
			Pa_I_I tmp = make_pair(sta, end);
			pos.push_back(tmp);
		}
		Se_Str::iterator iter;
		Str Seq;
		for (int j = 0; j < pos.size(); j++) {
			//cout << "Position:\t" << red.record[Synonym].Strand << "\t" << Synonym << "\t" << pos[j].first << "\t" << pos[j].second;
			if (Opt == "tis" || Opt == "cds" || Opt == "gene" || Opt == "stc")
				Seq = cleave(fna[red.record[Synonym].Replicon], red.record[Synonym].Strand, pos[j].first, pos[j].second, true);
			else if ((Opt == "mrna" || (Opt == "urbs" && red.record[Synonym].Type != "TA-like" && red.record[Synonym].Type != "TANNTTTG-like")) 
				&& red.record[Synonym].Type != "NA" && red.record[Synonym].Type != "Atypical")
				Seq = cleave(fna[red.record[Synonym].Replicon], red.record[Synonym].Strand, pos[j].first, pos[j].second, false);
			if (Opt == "cds")
				Seq = translate(Seq, cod);
			if (Seq != "Error!" && !Seq.empty() && (Opt == "gene" || ((Opt == "mrna" || Opt == "urbs" || Opt == "tis" || Opt == "cds" || Opt == "stc") && red.record[Synonym].Feature == "CDS"))) {
				if (Opt == "tis") {
					if (!red.record[Synonym].tss.empty()) {
						bool LDL = false;
						bool TRS = false;
						for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
							//cout << Synonym << "\t" << iter->first << "\t" << iter->second.utr << endl;
							if (abs(iter->second.utr) <= 5)
								LDL = true;
							if (iter->second.utr <= 5 && iter->second.utr >= -500)
								TRS = true;
						}
						if (LDL) {
							Ofs ldl(Ldl.data(), ios_base::app);
							ldl << ">" << red.Genome << "|" << red.record[Synonym].Replicon << "|" << Synonym;
							if (red.record[Synonym].Strand == "+")
								ldl << "|" << pos[j].first << "|" << Opt << endl << Seq << endl;
							else
								ldl << "|" << pos[j].second << "|" << Opt << endl << Seq << endl;
							ldl.close();
						}
						if (TRS) {
							Ofs exp(Exp.data(), ios_base::app);
							exp << ">" << red.Genome << "|" << red.record[Synonym].Replicon << "|" << Synonym;
							if (red.record[Synonym].Strand == "+")
								exp << "|" << pos[j].first << "|" << Opt << endl << Seq << endl;
							else
								exp << "|" << pos[j].second << "|" << Opt << endl << Seq << endl;
							exp.close();
						}
					}
					else if (!red.record[Synonym].TSS.empty() && red.record[Synonym].TSS != "NA") {
						Ofs exp(Exp.data(), ios_base::app);
						exp << ">" << red.Genome << "|" << red.record[Synonym].Replicon << "|" << Synonym;
						if (red.record[Synonym].Strand == "+")
							exp << "|" << pos[j].first << "|" << Opt << endl << Seq << endl;
						else
							exp << "|" << pos[j].second << "|" << Opt << endl << Seq << endl;
						exp.close();
					}
				}
				out << ">" << red.Genome << "|" << red.record[Synonym].Replicon << "|" << Synonym;
				if (red.record[Synonym].Strand == "+")
					out << "|" << pos[j].first << "|" << Opt << endl << Seq << endl;
				else
					out << "|" << pos[j].second << "|" << Opt << endl << Seq << endl;
			}
		}
	}
	out.close();
	return "All " + Opt +  " sequences of " + Geno + " have been extracted!";
}


Str ips(Str Out, Str Geno, Ma_Str_Str fa, Ma_Str_Str fna) {
	Out += Geno + ".ips";
	double gc = gc_content(fna);
	Ofs out(Out.data());
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		if (iter == fa.begin()) 
			out << iter->second.length() << "\t" << fa.size() << "\t" << gc << endl;
		if(iter->second != "Error!")
			out << iter->second << endl;
	}
	out.close();
	return "Ips file for " + Geno + " have been prepared!";
}

Str nonredundant(Str Out, Ma_Str_Str rfa) {
	Ofs out(Out.data(), ios_base::app);
	for (Ma_Str_Str::const_iterator iter = rfa.begin(); iter != rfa.end(); iter++) {
		out << ">" << iter->second << endl << iter->first << endl;
	}
	out.close();
	return "Redundent sequences have been removed!";
}

Str rRNA(Str Out, Ma_Str_Str fa, Ma_Str_T tax, int len) {
	Ma_Str_Se_Str rrna;
	for (Ma_Str_Str::const_iterator iter = fa.begin(); iter != fa.end(); iter++) {
		Ve_Str tmp = string_parse(iter->first, "|");
		Str Rev = iter->second.substr(iter->second.length() - len, len);
		Str For = antisense(Rev,true);
		rrna[tmp[0]].insert(For);

		/*Str Output;
		if (tax[tmp[0]].Sub != "NA")
			Output = Out + tax[tmp[0]].Sub + ".fa";
		else
			Output = Out + tax[tmp[0]].Phylum + ".fa";
		Ofs out(Output.data(), ios_base::app);*/
		//out << ">" << iter->first << endl << iter->second.substr(iter->second.length() - len, len) << endl;
		//out.close();
	}
	for (Ma_Str_Se_Str::const_iterator ite = rrna.begin(); ite != rrna.end(); ite++) {
		Str Output;
		if (tax[ite->first].Sub != "NA")
			Output = Out + tax[ite->first].Sub + ".fa";
		else if (tax[ite->first].Phylum == "Proteobacteria-Unclassified")
			Output = Out + "Gammaproteobacteria-Others.fa";
		else if (ite->first == "GCA_001729625.1_ASM172962v1")
			Output = Out + "Gammaproteobacteria-Enterobacteriales.fa";
		else if (tax[ite->first].Phylum == "Chloroflexi" || tax[ite->first].Phylum == "Peregrinibacteria" || tax[ite->first].Phylum == "Saccharibacteria" ||
			tax[ite->first].Phylum == "Terra-Unclassified")
			Output = Out + "Chloroflexi.fa";
		else if (ite->first == "GCA_001007975.1_ASM100797v1" || ite->first == "GCA_001029715.1_ASM102971v1" || ite->first == "GCA_001029635.1_ASM102963v1" ||
			ite->first == "GCA_001029675.1_ASM102967v1" || ite->first == "GCA_001029735.1_ASM102973v1" || ite->first == "GCA_001029755.1_ASM102975v1" ||
			ite->first == "GCA_001029775.1_ASM102977v1" || ite->first == "GCA_001029795.1_ASM102979v1" || ite->first == "GCF_000503875.1_ASM50387v1")
			Output = Out + "Chloroflexi.fa";
		else if (tax[ite->first].Phylum == "Aquificae" || tax[ite->first].Phylum == "Caldiserica" || tax[ite->first].Phylum == "Dictyoglomi" ||
			tax[ite->first].Phylum == "Thermotogae" || tax[ite->first].Phylum == "Thermodesulfobacteria" || tax[ite->first].Phylum == "Synergistetes")
			Output = Out + "Thermophiles.fa";
		else if (tax[ite->first].Phylum == "Bacteroidetes" || tax[ite->first].Phylum == "Calditrichaetes" || tax[ite->first].Phylum == "Chlorobi" ||
			tax[ite->first].Phylum == "Deinococcus-Thermus" || tax[ite->first].Phylum == "Nitrospirae")
			Output = Out + "CB_Group.fa";
		else if (ite->first == "GCF_000091165.1_ASM9116v1" || ite->first == "GCF_000146065.2_ASM14606v1")
			Output = Out + "CB_Group.fa";
		else if (tax[ite->first].Phylum == "Chlamydiae" || tax[ite->first].Phylum == "Fibrobacteres" || tax[ite->first].Phylum == "Gemmatimonadetes" ||
			tax[ite->first].Phylum == "Kiritimatiellaeota" || tax[ite->first].Phylum == "Planctomycetes" || tax[ite->first].Phylum == "Verrucomicrobia")
			Output = Out + "PVC_Group.fa";
		else if (ite->first == "GCF_000503835.1_ASM50383v1" || ite->first == "GCF_001443605.1_ASM144360v1")
			Output = Out + "PVC_Group.fa";
		else if (tax[ite->first].Phylum == "Acidobacteria" || tax[ite->first].Phylum == "Armatimonadetes" || tax[ite->first].Phylum == "Chrysiogenetes" ||
			tax[ite->first].Phylum == "Deferribacteres" || tax[ite->first].Phylum == "Elusimicrobia" || tax[ite->first].Phylum == "Fusobacteria" ||
			tax[ite->first].Phylum == "Spirochaetes")
			Output = Out + "Minorities.fa";
		else if (tax[ite->first].Domain == "Archaea")
			Output = Out + tax[ite->first].Domain + ".fa";
		else
			Output = Out + tax[ite->first].Phylum + ".fa";
		Ofs out(Output.data(), ios_base::app);
		for (Se_Str::const_iterator it = ite->second.begin(); it != ite->second.end(); it++)
			out << ">" << ite->first << endl << *it << endl;
		out.close();
	}
	return "16S rRNA tails have been extracted!";
}
