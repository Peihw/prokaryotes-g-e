#include "prokaryotes.h"
using namespace std;

Ma_Str_Str identify(Ve_Str id, Str Key) {
	Ma_Str_Str res;
	for (int i = 0; i < id.size(); i++) {
		int posb = id[i].find("_");
		int pose = id[i].find(Key);
		Str ID = id[i].substr(posb + 1, pose - posb - 1);
		res[ID] = id[i];
	}
	return res;
}

Str correlate(Str Out, Ve_Str fna, Ve_Str ftt) {
	Ma_Str_Str fn = identify(fna, "_genomic.fna");
	Ma_Str_Str ft = identify(ftt, "_feature_table");
	Ofs out(Out.data());
	out << "Assembly ID\tFna file\tFtt file" << endl;
	for (Ma_Str_Str::const_iterator iter = fn.begin(); iter != fn.end(); iter++) {
		if (ft[iter->first].empty())
			out << iter->first << "\t" << iter->second << "\tNA" << endl;
		else
			out << iter->first << "\t" << iter->second << "\t" << ft[iter->first] << endl;
	}
	out.close();
	return "FNA ID and FTT ID are correlated!";
}

Str prokaryotes(Str Fna, Str Ftt, Str Out, Ma_Str_Str fna, Ma_Str_K keg, ftt ano) {
	/*Ve_Str genoid;
	for (Se_Str::const_iterator ite = ano.Chr.begin(); ite != ano.Chr.end(); ite++) {
		genoid.push_back(*ite);
		cout << *ite << endl;
	}*/
	bool INI = false;
	Ifs in(Out.data());
	if (!in)
		INI = true;
	in.close();
	Ofs out(Out.data(), ios_base::app);
	//int num = 0;
	if (INI)
		out << "Fna assembly ID\tFtt assembly ID\tFna chromosome ID\tFtt chromosome ID\tType\tDomain\tPhylum\tGenus\tGenome full name" << endl;
	for (Ma_Str_Str::const_iterator iter = fna.begin(); iter != fna.end(); iter++) {
		Str NC, Type, Genome, Genus;
		bool GENO = false;
		bool GENUS = false;
		Genome = "";
		Type = "Chromosome";
		int posp = iter->first.find("lasmid");
		int posh = iter->first.find("phage");
		int posc = iter->first.find("contig");
		int posC = iter->first.find("Contig");
		int posd = iter->first.find("draft");
		int poss = iter->first.find("scaffold");
		int posS = iter->first.find("Scaffold");
		int posD = iter->first.find("DRAFT");
		if (posp != Str::npos)
			Type = "Plasmid";
		else if (posh != Str::npos)
			Type = "Phage";
		else if (posc != Str::npos || posd != Str::npos || poss != Str::npos || posD != Str::npos || posC != Str::npos || posS != Str::npos)
			Type = "Draft";
		Ve_Str temp = string_parse(iter->first, " ");
		for (int i = 0; i < temp.size(); i++) {
			if (!i) {
				int pos = temp[i].find(".");
				NC = iter->first.substr(0, pos);
				GENUS = true;
			}
			else if (temp[i][temp[i].length() - 1] == ',') {
				Genome += " " + temp[i].substr(0, temp[i].length() - 1);
				GENO = true;
			}
			else if (GENUS) {
				if (temp[i] != "Candidatus" && temp[i] != "UNVERIFIED:" && temp[i] != "Chlamydophila")
					Genus = temp[i];
				else if (temp[i] == "Chlamydophila")
					Genus = "Chlamydia";
				else
					Genus = temp[i + 1];
				if (Genus[0] == '[' || Genus[0] == '\'')
					Genus = Genus.substr(1, Genus.length() - 1);
				if (Genus[Genus.length() - 1] == ']')
					Genus = Genus.substr(0, Genus.length() - 1);
				Genome = temp[i];
				GENUS = false;
			}
			else if (!GENO)
				Genome += " " + temp[i];
		}
		//if (num < genoid.size())
		if(ano.Chr.find(NC) != ano.Chr.end())
			out << Fna << "\t" << Ftt << "\t" << NC << "\t" << NC << "\t" << Type;
		else {
			out << Fna << "\t" << Ftt << "\t" << NC << "\t";
			for (Se_Str::const_iterator iter = ano.Chr.begin(); iter != ano.Chr.end(); iter++)
				out << *iter << ";";
			out << "\t" << Type;
		}
		if (!keg[Genus].Domain.empty())
			out << "\t" << keg[Genus].Kingdom << "\t" << keg[Genus].Phylum << "\t" << Genus;
		else
			out << "\tNA\tNA\t" << Genus;
		out << "\t" << Genome << endl;
		//num++;
	}
	out.close();
	return "All prokaryotes have been listed!";
}

Str nc_fasta(Str Out, Ma_Str_Str fna) {
	for (Ma_Str_Str::const_iterator iter = fna.begin(); iter != fna.end(); iter++) {
		int pos = iter->first.find(".");
		Str NC = iter->first.substr(0, pos);
		Str Output = Out + NC + ".fna";
		Ofs out(Output.data());
		out << ">" << iter->first << endl << iter->second << endl;
		out.close();
	}
	return "Sequences have been outputed by chromosome!";
}

Str ptt_output(Str Out, Str NC, ftt ano) {
	Out += NC + ".med";
	Ofs out(Out.data());
	out << "Location\tStrand\tGene" << endl;
	for (int i = 0; i < ano.order.size(); i++) {
		Str Gene = ano.order[i];
		if (ano.S2F[Gene].GenoID == NC && ano.S2F[Gene].Feature == "CDS")
			out << ano.S2F[Gene].start << ".." << ano.S2F[Gene].end << "\t" << ano.S2F[Gene].Strand << "\t" << Gene << endl;
	}
	out.close();
	return "Annotations have been reformed in med format!";
}

Str tritisa_record(Str Out, Str Chr, Ma_I_M tisa, ftt ano, Ma_Str_Str trans) {
	Out += ".tritisa.rec.dat";
	bool INI = false;
	Ifs in(Out.data());
	if (!in)
		INI = true;
	in.close();
	Ofs out(Out.data(),ios_base::app);
	if (INI)
		out << "Genomic Accession\tStrand\tLocus\tGene Name\tProduct\tGene Start\tGene End\tAnnotation Reference\tStart Codon Shift(5\'-,3\'+)\tFeature" << endl;
	for (int i = 0; i < ano.order.size(); i++) {
		Str Synonym = ano.order[i];
		if (ano.S2F[Synonym].GenoID == Chr){ //&& ano.S2F[Synonym].Feature == "CDS") {
			if (trans[ano.S2F[Synonym].GenoID].empty())
				out << ano.S2F[Synonym].GenoID;
			else
				out << trans[ano.S2F[Synonym].GenoID];
			out << "\t" << ano.S2F[Synonym].Strand << "\t" << Synonym << "\t" << ano.S2F[Synonym].Symbol << "\t" << ano.S2F[Synonym].Product;
			int p3, p5;
			p3 = p5 = 0;
			if (ano.S2F[Synonym].Strand == "+") {
				p3 = ano.S2F[Synonym].end;
				p5 = ano.S2F[Synonym].start;
			}
			else {
				p3 = ano.S2F[Synonym].start;
				p5 = ano.S2F[Synonym].end;
			}
			if (tisa[p3].Strand.empty()) {
				int pos = ano.AssemID.find("GCA");
				if (pos != Str::npos)
					out << "\t" << ano.S2F[Synonym].start << "\t" << ano.S2F[Synonym].end << "\tGenbank\tNA";
				else
					out << "\t" << ano.S2F[Synonym].start << "\t" << ano.S2F[Synonym].end << "\tRefSeq\tNA";
			}
			else {
				if (p5 - tisa[p3].p5 == 0 || abs(p5 - tisa[p3].p5) > 36) {
					int pos = ano.AssemID.find("GCA");
					if (pos != Str::npos)
						out << "\t" << ano.S2F[Synonym].start << "\t" << ano.S2F[Synonym].end << "\tGenbank\tNA";
					else
						out << "\t" << ano.S2F[Synonym].start << "\t" << ano.S2F[Synonym].end << "\tRefSeq\tNA";
				}
				else {
					if (tisa[p3].Strand == "+")
						out << "\t" << tisa[p3].p5 << "\t" << tisa[p3].p3 << "\tTriTISA\t" << tisa[p3].p5 - p5;
					else
						out << "\t" << tisa[p3].p3 << "\t" << tisa[p3].p5 << "\tTriTISA\t" << p5 - tisa[p3].p5;
				}
			}
			out << "\t" << ano.S2F[Synonym].Feature << endl;
		}
	}
	out.close();
	return "All Tritisa data have been recorded!";
}

Str blast2cog(Str Out, Str Geno, ftt ano, Ma_Str_L bla, Ma_Str_Str csv) {
	Out += Geno + ".cog";
	Ofs out(Out.data());
	for (int i = 0; i < ano.order.size(); i ++) {
		Str Gene = ano.order[i];
		if (!bla[Gene].Que.empty()) {
			Ve_Str tmp = string_parse(bla[Gene].Ref, "|");
			out << Gene << "\t" << csv[tmp[1]] << endl;
		}
		else if (!bla[ano.S2F[Gene].ProdID].Que.empty()) {
			Ve_Str tmp = string_parse(bla[ano.S2F[Gene].ProdID].Ref, "|");
			out << Gene << "\t" << csv[tmp[1]] << endl;
		}
	}
	out.close();
	return "COG annotation of " + Geno + " have been done!";
}

Str reannotate(Str Out, Str Geno, Ma_Str_Str fna, ftt ano) {
	Out += Geno + "_feature_table.txt";
	Ofs out(Out.data());
	out << "# feature\tclass\tassembly\tassembly_unit\tseq_type\tchromosome\tgenomic_accession\tstart\tend\tstrand\tproduct_accession\tnon-redundant_refseq\trelated_accession\tname\tsymbol\tGeneID\tlocus_tag\tfeature_interval_length\tproduct_length\tattributes" << endl;
	for (int i = 0; i < ano.order.size(); i++) {
		Str Synonym = ano.order[i];
		out << ano.S2F[Synonym].Feature << "\t" << ano.S2F[Synonym].Class << "\t" << ano.AssemID << "\t" << ano.S2F[Synonym].Unit
			<< "\t" << ano.S2F[Synonym].Seq_type << "\t" << ano.S2F[Synonym].Chromosome << "\t" << ano.S2F[Synonym].GenoID
			<< "\t" << ano.S2F[Synonym].start << "\t" << ano.S2F[Synonym].end;
		if (ano.S2F[Synonym].Strand != "+" && ano.S2F[Synonym].Strand != "-") {
			Str Seq = cleave(fna[ano.S2F[Synonym].GenoID], "+", ano.S2F[Synonym].start, ano.S2F[Synonym].end, true);
			Str Stop = Seq.substr(Seq.length() - 3, 3);
			cout << "+:\t" << Seq.substr(0, 3) << "\t" << Seq.substr(Seq.length() - 3, 3) << endl;
			if (Stop == "TAA" || Stop == "TGA" || Stop == "TAG")
				ano.S2F[Synonym].Strand = "+";
			Seq = cleave(fna[ano.S2F[Synonym].GenoID], "-", ano.S2F[Synonym].start, ano.S2F[Synonym].end, true);
			Stop = Seq.substr(Seq.length() - 3, 3);
			cout << "-:\t" << Seq.substr(0, 3) << "\t" << Seq.substr(Seq.length() - 3, 3) << endl;
			if (Stop == "TAA" || Stop == "TGA" || Stop == "TAG")
				ano.S2F[Synonym].Strand = "-";
			out << "\t" << ano.S2F[Synonym].Strand;
		}
		else
			out << "\t" << ano.S2F[Synonym].Strand;
		out << "\t" << ano.S2F[Synonym].ProdID << "\t" << ano.S2F[Synonym].Nr << "\t" << ano.S2F[Synonym].RelID
			<< "\t" << ano.S2F[Synonym].Product << "\t" << ano.S2F[Synonym].Symbol << "\t" << ano.S2F[Synonym].GeneID << "\t" << Synonym 
			<< "\t" << ano.S2F[Synonym].interval << "\t" << ano.S2F[Synonym].product << "\t" << ano.S2F[Synonym].Attributes << endl;
	}
	out.close();
	return "Reannotaion of " + Geno + " have been done!";
}

Str exp_map(Str Out, Str Geno, rec red, Ma_Str_Str fna, Ve_E ep) {
	Str Output = Out + Geno + ".exp.rec.dat";
	Str Ope = Out + Geno + ".ope.dat";
	Str Ph = Out + Geno + ".phase.dat";
	Ma_Str_E exp;
	Ma_Str_Str op;
	Ma_Str_Se_Str ph;
	for (int i = 0; i < ep.size(); i++) {
		Pa_I_I tss, tts;
		Str Stt, End;
		Se_Str temp;
		if (ep[i].Strand == "+") {
			if (ep[i].tts > 0) {
				tts = locate(red.forward[ep[i].Replicon], ep[i].tts, red.record, "+", fna[ep[i].Replicon].length(), false);
				End = red.forward[ep[i].Replicon][tts.first];
			}
			for (Se_I::const_iterator ite = ep[i].tss.begin(); ite != ep[i].tss.end(); ite++) {
				tss = locate(red.forward[ep[i].Replicon], *ite, red.record, "+", fna[ep[i].Replicon].length(), true);
				Stt = red.forward[ep[i].Replicon][tss.first];
				//cout << i << "\t" << ep[i].Strand << "\t" << *ite << "\t" << Stt << endl;
				exp[Stt].ts[*ite].utr = tss.second;
				for (Se_Str::const_iterator iter = ep[i].phase.begin(); iter != ep[i].phase.end(); iter++) {
					exp[Stt].ts[*ite].phase.insert(*iter);
					Str Phase = *iter;
					if (!Phase.empty())
						temp.insert(Phase);
				}
				for (Se_Str::const_iterator iter = ep[i].ref.begin(); iter != ep[i].ref.end(); iter++)
					exp[Stt].ts[*ite].ref.insert(*iter);
			}
		}
		else {
			if (ep[i].tts > 0) {
				tts = locate(red.reverse[ep[i].Replicon], ep[i].tts, red.record, "-", fna[ep[i].Replicon].length(), false);
				End = red.reverse[ep[i].Replicon][tts.first];
			}
			for (Se_I::const_iterator ite = ep[i].tss.begin(); ite != ep[i].tss.end(); ite++) {
				tss = locate(red.reverse[ep[i].Replicon], *ite, red.record, "-", fna[ep[i].Replicon].length(), true);
				Stt = red.reverse[ep[i].Replicon][tss.first];
				exp[Stt].ts[*ite].utr = tss.second;
				for (Se_Str::const_iterator iter = ep[i].phase.begin(); iter != ep[i].phase.end(); iter++) {
					exp[Stt].ts[*ite].phase.insert(*iter);
					Str Phase = *iter;
					if(!Phase.empty())
						temp.insert(Phase);
				}
				for(Se_Str::const_iterator iter = ep[i].ref.begin(); iter != ep[i].ref.end(); iter++)
					exp[Stt].ts[*ite].ref.insert(*iter);
			}
		}
		if (ep[i].tts > 0) {
			exp[End].tt[ep[i].tts].utr = tts.second;
			for (Se_Str::const_iterator iter = ep[i].phase.begin(); iter != ep[i].phase.end(); iter++) {
				exp[End].tt[ep[i].tts].phase.insert(*iter);
				Str Phase = *iter;
				if (!Phase.empty())
					temp.insert(Phase);
			}
			for (Se_Str::const_iterator iter = ep[i].ref.begin(); iter != ep[i].ref.end(); iter++)
				exp[End].tt[ep[i].tts].ref.insert(*iter);
			Str Ope;
			if (ep[i].Strand == "+"){
				double ratio = (tss.first - tts.first) * 1.0 / red.forward[ep[i].Replicon].size();
				if (tss.first <= tts.first) {
					if (tts.first - tss.first < 100) {
						for (int j = tss.first; j <= tts.first; j++) {
							if (j == tss.first)
								Ope = red.forward[ep[i].Replicon][j];
							else
								Ope += ";" + red.forward[ep[i].Replicon][j];
							for(Se_Str::const_iterator ite = temp.begin(); ite != temp.end(); ite ++)
								ph[red.forward[ep[i].Replicon][j]].insert(*ite);
						}
					}
				}
				else if (ratio > 0.5) {
					for (int j = tss.first; j < red.forward.size(); j++) {
						if (j == tss.first)
							Ope = red.forward[ep[i].Replicon][j];
						else
							Ope += ";" + red.forward[ep[i].Replicon][j];
						for (Se_Str::const_iterator ite = temp.begin(); ite != temp.end(); ite++)
							ph[red.forward[ep[i].Replicon][j]].insert(*ite);
					}
					for (int j = 0; j <= tts.first; j++) {
						Ope += ";" + red.forward[ep[i].Replicon][j];
						for (Se_Str::const_iterator ite = temp.begin(); ite != temp.end(); ite++)
							ph[red.forward[ep[i].Replicon][j]].insert(*ite);
					}
				}
				//cout << Ope << endl;
			}
			else {
				double ratio = (tts.first - tss.first) * 1.0 / red.reverse[ep[i].Replicon].size();
				if (tss.first >= tts.first) {
					if (tss.first - tts.first < 100) {
						for (int j = tts.first; j <= tss.first; j++) {
							if (j == tts.first)
								Ope = red.reverse[ep[i].Replicon][j];
							else
								Ope += ";" + red.reverse[ep[i].Replicon][j];
							for (Se_Str::const_iterator ite = temp.begin(); ite != temp.end(); ite++)
								ph[red.reverse[ep[i].Replicon][j]].insert(*ite);
						}
					}
				}
				else if (ratio > 0.5) {
					for (int j = tts.first; j < red.reverse.size(); j++) {
						if (j == tts.first)
							Ope = red.reverse[ep[i].Replicon][j];
						else
							Ope += ";" + red.reverse[ep[i].Replicon][j];
						for (Se_Str::const_iterator ite = temp.begin(); ite != temp.end(); ite++)
							ph[red.reverse[ep[i].Replicon][j]].insert(*ite);
					}
					for (int j = 0; j <= tss.first; j++) {
						Ope += ";" + red.reverse[ep[i].Replicon][j];
						for (Se_Str::const_iterator ite = temp.begin(); ite != temp.end(); ite++)
							ph[red.reverse[ep[i].Replicon][j]].insert(*ite);
					}
				}
			}
			op[Ope] = ep[i].Strand;
		}
	}
	if (!op.empty()) {
		Ofs ope(Ope.data());
		ope << "Strand\tGenes" << endl;
		for (Ma_Str_Str::const_iterator iter = op.begin(); iter != op.end(); iter++) {
			if(!iter->first.empty())
				ope << iter->second << "\t" << iter->first << endl;
		}
		ope.close();
	}
	if (!ph.empty()) {
		Ofs pha(Ph.data());
		pha << "Locus\tPhase" << endl;
		for (Ma_Str_Se_Str::const_iterator iter = ph.begin(); iter != ph.end(); iter++) {
			pha << iter->first;
			for (Se_Str::const_iterator ite = iter->second.begin(); ite != iter->second.end(); ite++) {
				if (ite == iter->second.begin())
					pha << "\t" << *ite;
				else
					pha << ";" << *ite;
			}
			pha << endl;
		}
		pha.close();
	}
	Ofs out(Output.data());
	out << "Genomic Accession\tStrand\tLocus\tGene Name\tProduct\tGene Start\tGene End\tTIS Annotation Reference\tTIS shift (5\'-,3\'+)\tFeature"
		<< "\tExperimentally Verified TSS (Position|5\'UTR Length|Phase|Reference)\tExperimentally Verified TTS (Position|Phase|Reference)" << endl;
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		out << red.record[Synonym].Replicon << "\t" << red.record[Synonym].Strand << "\t" << Synonym << "\t" << red.record[Synonym].Name
			<< "\t" << red.record[Synonym].Product << "\t" << red.record[Synonym].start << "\t" << red.record[Synonym].end
			<< "\t" << red.record[Synonym].TISRef;
		if (!red.record[Synonym].shift)
			out << "\tNA";
		else
			out << "\t" << red.record[Synonym].shift;
		out << "\t" << red.record[Synonym].Feature;
		if (!exp[Synonym].ts.empty()) {
			for (Ma_I_S::const_iterator iter = exp[Synonym].ts.begin(); iter != exp[Synonym].ts.end(); iter++) {
				if (iter == exp[Synonym].ts.begin())
					out << "\t" << iter->first << "|" << iter->second.utr;
				else
					out << ";" << iter->first << "|" << iter->second.utr;
				if (!iter->second.phase.empty()) {
					for (Se_Str::const_iterator ite = iter->second.phase.begin(); ite != iter->second.phase.end(); ite++) {
						if (ite == iter->second.phase.begin())
							out << "|" << *ite;
						else
							out << "&" << *ite;
					}
				}
				else
					out << "|NA";
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
		if(!exp[Synonym].tt.empty()){
			for (Ma_I_S::const_iterator iter = exp[Synonym].tt.begin(); iter != exp[Synonym].tt.end(); iter++) {
				if (iter == exp[Synonym].tt.begin())
					out << "\t" << iter->first;
				else
					out << ";" << iter->first;
				if (!iter->second.phase.empty()) {
					for (Se_Str::const_iterator ite = iter->second.phase.begin(); ite != iter->second.phase.end(); ite++) {
						if (ite == iter->second.phase.begin())
							out << "|" << *ite;
						else
							out << "&" << *ite;
					}
				}
				else
					out << "|NA";
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
		out << endl;
	}
	out.close();
	return "Experimentally verified TUs or TSSs of " + Geno + " have been mapped!";
}

Str transcript_map(Str Out, Str Geno, Ma_Str_T tu, rec red) {
	Str Output = Out + Geno + ".trans.rec.dat";
	Ofs out(Output.data());
	out<< "Genomic Accession\tStrand\tLocus\tGene Name\tProduct\tGene Start\tGene End\tTIS Annotation Reference\tTIS shift (5\'-,3\'+)\tFeature"
		<< "\tExperimentally Verified TSS (Position|5\'UTR Length|Phase|Reference)\tExperimentally Verified TTS (Position|Phase|Reference)"
		<< "\tExperimentally Verified Leaderless TIS (Position, 5\'End)\tOperon ID\tRole of TSSs\tRole of TTSs"<< endl;
	int id = 0;
	Str Last;
	Ma_Str_I corr;
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		out << red.record[Synonym].Replicon << "\t" << red.record[Synonym].Strand << "\t" << Synonym << "\t" << red.record[Synonym].Name
			<< "\t" << red.record[Synonym].Product << "\t" << red.record[Synonym].start << "\t" << red.record[Synonym].end
			<< "\t" << red.record[Synonym].TISRef;
		if (!red.record[Synonym].shift)
			out << "\tNA";
		else
			out << "\t" << red.record[Synonym].shift;
		out << "\t" << red.record[Synonym].Feature;
		Ma_I_I ldl;
		if (!red.record[Synonym].tss.empty()) {
			for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
				if(abs(iter->second.utr) <= 5)
					ldl[iter->first] = iter->second.utr;
				if (iter == red.record[Synonym].tss.begin())
					out << "\t" << iter->first << "|" << iter->second.utr;
				else
					out << ";" << iter->first << "|" << iter->second.utr;
				if (!iter->second.phase.empty()) {
					for (Se_Str::const_iterator ite = iter->second.phase.begin(); ite != iter->second.phase.end(); ite++) {
						if (ite == iter->second.phase.begin())
							out << "|" << *ite;
						else
							out << "&" << *ite;
					}
				}
				else
					out << "|NA";
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
		if (!ldl.empty()) {
			for (Ma_I_I::const_iterator ite = ldl.begin(); ite != ldl.end(); ite++) {
				if (ite == ldl.begin())
					out << "\t" << ite->first << "|" << ite->second;
				else
					out << ";" << ite->first << "|" << ite->second;
			}
		}
		else
			out << "\tNA";
		if (!tu[Synonym].ID.empty()) {
			if (Last != tu[Synonym].ID && corr[tu[Synonym].ID] <= 0) {
				id++;
				out << "\t" << id;// << "\t" << tu[Synonym].ID;
				corr[tu[Synonym].ID] = id;
			}
			else
				out << "\t" << corr[tu[Synonym].ID];// << "\t" << tu[Synonym].ID;
			Last = tu[Synonym].ID;
		}
		else
			out << "\tNA";
		if (!tu[Synonym].TSS.empty()) {
			if (tu[Synonym].TSS == "None" && !red.record[Synonym].tss.empty())
				out << "\t" << "Sub-operonic TSS";
			else
				out << "\t" << tu[Synonym].TSS;
		}
		else if (!red.record[Synonym].tss.empty()) 
			out << "\tVerified TSS";
		else
			out << "\tNA";
		if (!tu[Synonym].TTS.empty()) {
			if (tu[Synonym].TTS == "NA" && !red.record[Synonym].tts.empty())
				out << "\t" << "Sub-operonic TTS";
			else
				out << "\t" << tu[Synonym].TTS;
		}
		else
			out << "\tNA";
		out << endl;
	}
	out.close();
	return "All transcripts of " + Geno + " have been mapped!";
}


Str operon(Str Out, Str Geno, Ve_Str genes, Ve_Pa_Str_I opr, Ma_Str_Str exp) {
	Str Output = Out + Geno + ".ope.dat";
	Ofs out(Output.data());
	out << "Strand\tGenes" << endl;
	if (!genes.empty()) {
		for (int i = 0; i < genes.size(); i++) {
			Ma_Str_I count;
			Ve_Str tmp = string_parse(genes[i], ";");
			Str Genes = "";
			for (int j = 0; j < tmp.size(); j++) {
				Str Synonym = exp[tmp[j]];
				if (!Synonym.empty()) {
					if (Genes.empty())
						Genes = Synonym;
					else
						Genes += ";" + Synonym;
					count[exp[Synonym]] ++;
				}
			}
			if (count["+"] > count["-"])
				out << "+\t" << Genes << endl;
			else if (count["+"] < count["-"])
				out << "-\t" << Genes << endl;
			else
				cout << "Error at " << genes[i] << endl;
		}
	}
	if (!opr.empty()) {
		int last = -1;
		Str Genes, Strand;
		Ma_Str_I count;
		for (int i = 0; i < opr.size(); i++) {
			Str Synonym = exp[opr[i].first];
			if (Synonym == "+" || Synonym == "-")
				Synonym = opr[i].first;
			if (!Synonym.empty()) {
				if (last == opr[i].second) {
					Genes += ";" + Synonym;
					count[exp[Synonym]]++;
				}
				else {
					if (!Genes.empty()) {
						if (count["+"] > count["-"])
							out << "+\t" << Genes << endl;
						else if (count["+"] < count["-"])
							out << "-\t" << Genes << endl;
						else
							cout << "Error at " << Genes << endl;
					}
					Genes = Synonym;
					last = opr[i].second;
					count.clear();
					count[exp[Synonym]]++;
				}
			}
		}
	}
	out.close();
	return "All operons have been remapped!";
}

