#include "scoring.h"
using namespace std;

Str RPS1(Str Out, Str GenoID, rec red, Ma_Str_X tax) {
	Ofs out(Out.data(), ios_base::app);
	bool FOUND = false;
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		Ve_Str tmp = string_parse(red.record[Synonym].Product, " ");
		bool R, P, S;
		R = P = S = false;
		for (int j = 0; j < tmp.size(); j++) {
			if (tmp[j] == "ribosomal" || tmp[j] == "Ribosomal" || tmp[j] == "Ribosome" || tmp[j] == "ribosome")
				R = true;
			if (tmp[j] == "protein" || tmp[j] == "Protein")
				P = true;
			if (tmp[j] == "S1")
				S = true;
		}
		if (R & P & S) {
			out << GenoID << "\t" << Synonym << "\t" << red.record[Synonym].Product << "\t" << tax[GenoID].Phylum << endl;
			FOUND = true;
			break;
		}
	}
	if (!FOUND)
		out << GenoID << "\tNA\tNA\t" << tax[GenoID].Phylum << endl;
	out.close();
	return "Detection of RPS1 have been done for " + GenoID + "!";
}

Str dist(Str Out, Ma_Str_M ref, Ma_Str_M que, Str GenoID, int len, double tthr, double sthr, double pthr, Ma_Str_Pa_Str_Str check) {
	bool IN = true;
	Ifs in(Out.data());
	if (!in)
		IN = false;
	in.close();
	Ofs out(Out.data(), ios_base::app);
	bool RBS = false;
	if (!IN) {
		out << "Genome Accession\tMotif\tConsensus\tEntropy\tType";
		for (Ma_Str_M::const_iterator it = ref.begin(); it != ref.end(); it++) {
			out << "\t Distance to " << it->first;
			if (it->first == "Canonical Leaderless")
				out << "\tDistance to Canonical RBS";
		}
		out << endl;
	}
	for (int i = 0; i < que[GenoID].rwm.size(); i++) {
		Str SD = "AAGGAGGTGA";
        bool EXC = false;
        for(int l = 0; l < Base.size(); l++){
            if(que[GenoID].pb[Base[l]] > 0.4){
                EXC = true;
                cout << GenoID << " with " << Base[l] << " as " << que[GenoID].pb[Base[l]] << " exceeds 0.4!" << endl;
            }
        }
		double tmp = SD_distance(SD, que[GenoID].rwm[i], que[GenoID].pb, len, EXC);
		out << GenoID << "\tSignal" << i + 1 << "\t" << que[GenoID].cons[i] << "\t" << que[GenoID].entropy[i];
		Ma_Str_D score;
		score["Canonical RBS"] = tmp;
		double dis = 1000;
		Str Type;
		for (Ma_Str_M::const_iterator iter = ref.begin(); iter != ref.end(); iter++) {
			for (int j = 0; j < iter->second.rwm.size(); j++) {
				tmp = TA_distance(iter->second.rwm[j],iter->second.pb, que[GenoID].rwm[i], que[GenoID].pb, len + 1, EXC);
				score[iter->first] = tmp;
				if (tmp < dis) {
					dis = tmp;
					Type = iter->first;
				}
			}
		}
		if (score["Canonical RBS"] < sthr)
			out << "\tCanonical RBS";
		else if (score["Canonical Leaderless"] < tthr)
			out << "\tCanonical Leaderless";
		else
			out << "\t" << Type;
		for (Ma_Str_D::const_iterator ite = score.begin(); ite != score.end(); ite++)
			out << "\t" << ite->second;
		out << endl;
		/*double cldl = 0;
		if (ref.rwm.size() > 1) {
			Ma_C_D bak;
			if (ref.rwm[0].size() > 7) {
				bak['A'] = 0.522086;
				bak['C'] = 0.139728;
				bak['G'] = 0.120138;
				bak['T'] = 0.218047;
			}
			else {
				bak['A'] = 0.27987;
				bak['C'] = 0.210288;
				bak['G'] = 0.210519;
				bak['T'] = 0.299323;
			}
			cldl = TA_distance(ref.rwm[0], bak, que[0].rwm[i], que[0].pb, len + 1, EXC);
		}
		else
			cldl = TA_distance(ref.rwm[0], ref.pb, que[0].rwm[i], que[0].pb, len + 1, EXC);
		double nrbs = 0;
		double nldl = 0;
		double atg = 0;
		if (ref.rwm.size() > 1) 
			nrbs = TA_distance(ref.rwm[1], ref.pb, que[0].rwm[i], que[0].pb, len + 1, EXC);
		if(ref.rwm.size() > 2)
			nldl = TA_distance(ref.rwm[2], ref.pb, que[0].rwm[i], que[0].pb, len + 1, EXC);
		//Str ATG = "ATG";
		//atg = SD_distance(ATG, que.rwm[i], que.pb, ATG.length(), EXC);
		out << GenoID << "\tMotif" << i + 1 << "\t" << que.cons[i];
		if (ref.rwm.size() <= 2) {
			if (crbs < sthr && cldl < tthr)
				out << "\tAmbiguous";
			else if (crbs < sthr && cldl >= tthr) 
				out << "\tCanonical RBS";
			else if (crbs >= sthr && cldl < tthr)
				out << "\tCanonical Leaderless";
			else {
				if (nrbs < pthr && check[GenoID].first != "NA")
					out << "\tNon-canonical RBS";
				//else if (atg < 4)
					//out << "\tPutative uAUG";
				else
					out << "\tAtypical";
			}
		}
		else {
			if (crbs < pthr)
				out << "\tCanonical RBS";
			else if (nrbs < sthr && nldl < tthr) {
				if (nrbs > nldl)
					out << "\tNon-canonical Leaderless";
				else
					out << "\tNon-canonical RBS";
			}
			else if (nrbs < sthr && nldl >= tthr)
				out << "\tNon-canonical RBS";
			else if (nrbs >= sthr && nldl < tthr)
				out << "\tNon-canonical Leaderless";
			else {
				if (cldl < 5)
					out << "\tCanonical Leaderless";
				//else if (atg < 4)
					//out << "\tPutative uAUG";
				else
					out << "\tAtypical";
			}
		}
		out << "\t" << que.entropy[i] << "\t" << crbs << "\t" << cldl << "\t" << nrbs << "\t" << nldl << endl;// << "\t" << atg << endl;*/
	}
	out.close();
	return "The distances to SD and TATA standard signal have been calculated!";
}

double SD_distance(Str SD, Ve_Ma_C_D que, Ma_C_D pb, int len, bool EXC) {
	double result = 1000;
	for (int i = 0; i < que.size() - len + 1; i++) {
		for (int j = 0; j < SD.length() - len + 1; j++) {
			double score = 0;
			for (int k = 0; k < len; k++) {
                /*double sum = 0;
				if (!EXC) {
					for (int l = 0; l < Base.length(); l++)
						sum += que[i + k][Base[l]] / pb[Base[l]];
					score -= log2(que[i + k][SD[j + k]] / pb[SD[j + k]] / sum);
				}
				else
					score -= log2(que[i + k][SD[j + k]]);*/
			    score -= log2(que[i + k][SD[j+k]]);
				//cout << "\t" << que[i + k][SD[j + k]];
			}
			if (score < result)
				result = score;
			//cout << i << "\t" << SD.substr(j, len) << "\t" << score << "\t" << result << endl;
		}
	}
	return result;
}

double TA_distance(Ve_Ma_C_D ref, Ma_C_D rpb, Ve_Ma_C_D que, Ma_C_D qpb, int len, bool EXC) {
	double result = 1000;
	for (int i = 0; i < ref.size() - len + 1; i++) {
		for (int j = 0; j < que.size() - len + 1; j++) {
			double score = 0;
			for (int k = 0; k < len; k++) {
				//Normalize
				double dis = 0;
				/*double sumq, sumr, dot;
				sumq = sumr = dot = 0;
				for (int l = 0; l < Base.length(); l++) {
					sumq += que[j + k][Base[l]] / qpb[Base[l]];
					sumr += ref[i + k][Base[l]] / rpb[Base[l]];
				}*/
				for (int l = 0; l < Base.length(); l++) {
                    double a, b;
                    //if(!EXC){
					a = que[j + k][Base[l]] / qpb[Base[l]];// / sumq;
					b = ref[i + k][Base[l]] / rpb[Base[l]];// / sumr;
                    /*}
                    else{
                        a = que[j + k][Base[l]];
                        b = ref[i + k][Base[l]];
                    }*/
					//dot += a * b;
                    //score += b * log2(b / a);
					dis += (a - b)*(a - b);
				}
				dis = sqrt(dis);
				score += dis;
			}
			if (score < result) {
				result = score;
				//cout << "Ref: " << i << "\tQue: " << j << "\t" << result << endl;
			}
            //result = sqrt(result);
		}
	}
	return result;
}

Str tis_record(Str Out, Str Que, Ma_Str_Str fa, rec red, motif mot, Ma_Str_Ve_Str dis, int len) {
	Out += red.Genome + ".tis.rec.dat";
	Ofs out(Out.data());
	int utr = 0;
	Ma_Str_S sinfo = prob(fa, mot, dis[Que]);
	out << "Genomic Accession\tStrand\tLocus\tGene Name\tProduct\tGene Start\tGene End\tTIS Annotation Reference\tTIS shift (5\'-,3\'+)\tFeature"
		<< "\tExperimentally Verified TSS (Position|5\'UTR Length|Phase|Reference)\tExperimentally Verified TTS (Position|Phase|Reference)"
		<< "\tExperimentally Verified Leaderless TIS (Position|5\'End Length)\tOperon ID\tRole of TSSs\tRole of TTSs"
		<< "\tTranscription Initiation Signal\tTranscription Signal Type\tTranscription Signal Probability"
		<< "\tTranscription Signal Start\tTranscription Signal Sequence" << endl;
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		//cout << Synonym << "\t";
		//if (!fa[Synonym].empty() && fa[Synonym] != "Error!")
		out << red.record[Synonym].Replicon << "\t" << red.record[Synonym].Strand << "\t" << Synonym
			<< "\t" << red.record[Synonym].Name << "\t" << red.record[Synonym].Product
			<< "\t" << red.record[Synonym].start << "\t" << red.record[Synonym].end
			<< "\t" << red.record[Synonym].TISRef;
		if (red.record[Synonym].shift)
			out << "\t" << red.record[Synonym].shift;
		else
			out << "\tNA";
		out << "\t" << red.record[Synonym].Feature;
		if (!red.record[Synonym].tss.empty()) {
			for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
				if (iter == red.record[Synonym].tss.begin())
					out << "\t" << iter->first << "|" << iter->second.utr;
				else
					out << ";" << iter->first << "|" << iter->second.utr;
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
		if (!red.record[Synonym].ldl.empty()) {
			for (Ma_I_I::const_iterator ite = red.record[Synonym].ldl.begin(); ite != red.record[Synonym].ldl.end(); ite++) {
				if (ite == red.record[Synonym].ldl.begin())
					out << "\t" << ite->first << "|" << ite->second;
				else
					out << ";" << ite->first << "|" << ite->second;
			}
		}
		else
			out << "\tNA";
		if (red.record[Synonym].opr)
			out << "\t" << red.record[Synonym].opr;
		else
			out << "\tNA";
		if (!red.record[Synonym].TSS.empty())
			out << "\t" << red.record[Synonym].TSS;
		else if (!red.record[Synonym].tss.empty())
			out << "\tVerified TSS";
		else
			out << "\tNA";
		if (!red.record[Synonym].TTS.empty())
			out << "\t" << red.record[Synonym].TTS;
		else
			out << "\tNA";
		if (!sinfo[Synonym].Name.empty()) {
			out << "\t" << sinfo[Synonym].Name << "\t" << sinfo[Synonym].Type << "\t" << sinfo[Synonym].pb
				<< "\t" << sinfo[Synonym].start << "\t" << sinfo[Synonym].Seq << endl;
		}
		else
			out << "\tNA\tNA\tNA\tNA\tNA" << endl;
	}
	out.close();
	return "All information about " + red.Genome + " have been recorded!";
}

Ma_Str_S prob(Ma_Str_Str fa, motif mot, Ve_Str dis) {
	Ma_Str_S result;
	for (Ma_Str_Str::iterator iter = fa.begin(); iter != fa.end(); iter++) {
		scan temp;
		temp.pb = temp.sc = 0;
		double sum, lsc;
		int lst = 0;
		sum = lsc = 0;
		for (int i = 0; i < mot.rwm.size(); i++) {
			double pb = 0;
			for (int j = 0; j < mot.pj[i].size(); j++) {
				double sc = 1;
				sc *= mot.pj[i][j];
				//cout << sc << "\t" << i << "\t" << j << endl;
				for (int k = 0; k < mot.rwm[i].size(); k++)
					sc *= mot.rwm[i][k][iter->second[j + k]] / mot.pb[iter->second[j + k]];
				//cout << sc << "\t" << i << "\t" << j << endl;
				if (sc > lsc) {
					lsc = sc;
					lst = j;
				}
				pb += sc;
			}
			sum += pb;
			if (pb > temp.pb) {
				temp.pb = pb;
				temp.sc = lsc;
				int w = mot.rwm[i].size();
				int d = mot.pj[i].size();
				temp.start = lst - w - d + 1;
				temp.Seq = iter->second.substr(lst, w);
				temp.Type = dis[i];
				int a = i + 1;
				Str N;
				convertFromNumber(N, a);
				temp.Name = "Signal" + N;
			}
		}
		//cout << temp.pb << "\t" << temp.sc << endl;
		sum += mot.p0;
		if (temp.pb > mot.p0)
			temp.pb /= sum;
		else {
			temp.pb = mot.p0 / sum;
			temp.sc = 0;
			temp.Seq = "NA";
			temp.start = 0;
			temp.Name = "None";
			temp.Type = "NA";
		}
		result[iter->first] = temp;
	}
	return result;
}

Str ortholog_statistics(Str Out, Str Ori, Str Que, Ma_Str_R ori, Ma_Str_R que, Ma_Str_Ma_Str_O orth) {
	//cout << ori.size() << "\t" << que.size() << "\t" << orth.size() << endl;
	Ofs out(Out.data(), ios_base::app);//TATA, SDTA, SDSD, TASD;
	Ifs in(Out.data());
	if (!in)
		out << "Origin Genome\tQuery Genome\tCannonical to Cannonical\tCannonical to Leaderless\tLeaderless to Cannonical\tLeaderless to Leaderless\tDistance (log2)\tDistance (ln)\tDistance (lg)" << endl;
	Ma_Str_I counts;
	for (Ma_Str_R::const_iterator iter = que.begin(); iter != que.end(); iter++) {
		for (Ma_Str_R::const_iterator ite = ori.begin(); ite != ori.end(); ite++) {
			//cout << iter->first << "\t" << ite->first << endl;
			if (!orth[iter->first][ite->first].OID.empty()) {
				//cout << iter->first << "\t" << ite->first << "\t" << iter->second.Type << "\t" << ite->second.Type << endl;
				counts[iter->second.Type + ite->second.Type] ++;
			}
		}
	}
	double deno = counts["CannonicalCannonical"] * counts["LeaderlessLeaderless"] - counts["CannonicalLeaderless"] * counts["LeaderlessCannonical"];
	double bsd = counts["CannonicalCannonical"] + counts["CannonicalLeaderless"];
	double esd = counts["CannonicalCannonical"] + counts["LeaderlessCannonical"];
	double eta = counts["CannonicalLeaderless"] + counts["LeaderlessLeaderless"];
	double bta = counts["LeaderlessCannonical"] + counts["LeaderlessLeaderless"];
	double dist2, diste, dist10;
	dist2 = diste = dist10 = 0;
	//cout << deno << "\t" << bsd << "\t" << esd << "\t" << eta << "\t" << bta << endl;
	if (fabs(deno * bsd * esd * eta * bta) > 1.0e-3) {
		dist2 = -0.25 * (2 * log2(abs(deno)) - log2(bsd) - log2(esd) - log2(eta) - log2(bta));
		diste = -0.25 * (2 * log2(abs(deno)) - log2(bsd) - log2(esd) - log2(eta) - log(bta));
		dist10 = -0.25 * (2 * log10(abs(deno)) - log10(bsd) - log10(esd) - log10(eta) - log10(bta));
	}
	out << Ori << "\t" << Que << "\t" << counts["CannonicalCannonical"] << "\t" << counts["CannonicalLeaderless"] << "\t" << counts["LeaderlessCannonical"]
		<< "\t" << counts["LeaderlessLeaderless"] << "\t" << dist2 << "\t" << diste << "\t" << dist10 << endl;
	out.close();
	return "Statistics of orthologous genes have been finished!";
}

Str signal_statistics(Str Out, Str Ori, Ma_Str_R ori) {
	Ofs out(Out.data());
	out << "Genome ID\tCannonical TIS Genes\tLeaderless TIS Genes\tNon-TIS Genes\tAtypical TIS Genes\tTotal Genes" << endl << Ori;
	Ma_Str_I counts;
	for (Ma_Str_R::const_iterator iter = ori.begin(); iter != ori.end(); iter++)
		counts[iter->second.Type] ++;
	for (Ma_Str_I::const_iterator ite = counts.begin(); ite != counts.end(); ite++) {
		if (counts["Cannonical"] > 0)
			out << counts["Cannonical"];
		else
			out << 0;
		if (counts["Leaderless"] > 0)
			out << "\t" << counts["Leaderless"];
		else
			out << "\t" << 0;
		if (counts["None"] > 0)
			out << "\t" << counts["None"];
		else
			out << "\t" << 0;
		if (counts["Atypical"] > 0)
			out << "\t" << counts["Atypical"];
		else
			out << "\t" << 0;
		out << "\t" << ori.size() << endl;
	}
	out.close();
	return "Statistics of cannonical and leaderless signals have been finished!";
}

Str motif2fasta(Str Out, Str GenoID, motif mot, int depth) {
	Str Log = Out + "WebLOGO.sh";
	Ofs log(Log.data(), ios_base::app);
	for (int i = 1; i <= mot.rwm.size(); i++) {
		Str Num;
		convertFromNumber(Num, i);
		Str Output = Out + GenoID + ".Signal" + Num + ".fa";
		log << "./weblogo\t-F png\t--fontsize 10\t-c classic\t--number-interval 1\t--resolution 1200\t--composition \"{";
		for (int l = 0; l < Base.length(); l++) {
			int pb = mot.pb[Base[l]] * 10000 + 0.5;
			if (!l)
				log << "\'" << Base[l] << "\':" << pb / 100.0;
			else
				log << ", \'" << Base[l] << "\':" << pb / 100.0;
		}
		log << "}\"\t< " << Output << " >\t" << Output << ".png" << endl;
		Ofs out(Output.data());
		for (int j = 0; j < depth; j++) {
			out << ">" << j + 1 << endl;
			for (int k = 0; k < mot.rwm[i - 1].size(); k++) {
				double thr = 0;
				for(int l = 0; l < Base.length(); l ++){
					thr += mot.rwm[i - 1][k][Base[l]];
					if (j < depth * thr || l == Base.length() - 1) {
						//cout << j << "\t" << k << "\t" << mot.rwm[i - 1][k][Base[l]] << "\t" << thr << "\t" << thr * depth << endl;
						out << Base[l];
						break;
					}
				}
			}
			out << endl;
		}
		out.close();
	}
	log.close();
	return "Motifs of " + GenoID + " have been transfered to fasta!";
}

Str max_entropy(Str Out, motif ref, Str GenoID) {
	Ofs out(Out.data(), ios_base::app);
	out << GenoID;
	double sum = 0;
	for (int i = 0; i < ref.entropy.size(); i++) {
		out << "\t" << ref.entropy[i];
		sum += ref.entropy[i];
	}
	out << "\t" << sum << "\t" << sum / ref.entropy.size() << endl;
	out.close();
	return "Motifs of " + GenoID + " have been listed!";
}
