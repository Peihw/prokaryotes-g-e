#include "estimate.h"
using namespace std;

Str estimate(Str Out, Str GenoID, rec red) {
	bool INI = false;
	Ifs in(Out.data());
	if (!in)
		INI = true;
	in.close();
	Ofs out(Out.data(), ios_base::app);
	if (INI)
		out << "Assembly ID\tVerified Leadered gene\tVerified Leaderless gene\tVerified Internal TSS\tVerified Remote TSS"
		<< "\tPredicted Leadered gene\tLeadered TP\tLeadered FP\tLeadered TN\tLeadered FN\tLeadered Recall\tLeadered Precision"
		<< "\tPredicted Leaderless gene\tLeaderless TP\tLeaderless FP\tLeaderless TN\tLeaderless FN\tLeaderless Recall\tLeaderless Precision" << endl;
	Ma_Str_I counts;
	for (int i = 0; i < red.order.size(); i++) {
		Str Synonym = red.order[i];
		if (!red.record[Synonym].tss.empty() && red.record[Synonym].Feature == "CDS") {
			bool LDD, LDL, INT, RMT;
			LDD = LDL = INT = RMT = false;
			Str Exp = "";
			for (Ma_I_S::const_iterator iter = red.record[Synonym].tss.begin(); iter != red.record[Synonym].tss.end(); iter++) {
				//cout << Synonym << "\t" << iter->first << "\t" << iter->second.utr << endl;
				if (iter->second.utr < -500)
					RMT = true;
				else if (iter->second.utr < -5 && iter->second.utr >= -500)
					LDD = true;
				else if (iter->second.utr < 5 && iter->second.utr >= -5)
					LDL = true;
				else
					INT = true;
			}
			if (LDL)
				Exp = "Leaderless";
			else if (LDD)
				Exp = "Leadered";
			else if (RMT)
				Exp = "Remote";
			else
				Exp = "Internal";
			counts[Exp] ++;
			if (Exp == "Leadered" || Exp == "Leaderless") {
				Str Pred = "";
				int ldl = red.record[Synonym].Type.find("TA");
				if (ldl != Str::npos)
					Pred = "Leaderless";
				//int ldd = red.record[Synonym].Type.find("RBS");
				//if (ldd != Str::npos)
				else
					Pred = "Leadered";
				if (red.record[Synonym].Type == "Atypical")
					Pred = "Atypical";
				if (red.record[Synonym].Type == "NA")
					Pred = "NA";
				counts["Pred" + Pred] ++;
				if (Exp == "Leadered") {
					if (Pred == "Leadered") {
						counts["LDDTP"] ++;
						counts["LDLTN"] ++;
					}
					else {
						counts["LDDFN"] ++;
						if (Pred == "Leaderless")
							counts["LDLFP"] ++;
						else
							counts["LDLTN"] ++;
					}
				}
				if (Exp == "Leaderless") {
					if (Pred == "Leaderless") {
						counts["LDLTP"] ++;
						counts["LDDTN"] ++;
					}
					else {
						counts["LDLFN"] ++;
						if (Pred == "Leadered")
							counts["LDDFP"] ++;
						else
							counts["LDDTN"] ++;
					}
				}
			}
		}
	}
	double lddrc, lddps, ldlps, ldlrc;
	lddps = lddrc = ldlps = ldlrc = -1;
	if (counts["LDDTP"] + counts["LDDFP"] > 0)
		lddps = counts["LDDTP"] * 1.0 / (counts["LDDTP"] + counts["LDDFP"]);
	if (counts["LDDTP"] + counts["LDDFN"] > 0)
		lddrc = counts["LDDTP"] * 1.0 / (counts["LDDTP"] + counts["LDDFN"]);
	if (counts["LDLTP"] + counts["LDLFN"] > 0)
		ldlrc = counts["LDLTP"] * 1.0 / (counts["LDLTP"] + counts["LDLFN"]);
	if (counts["LDLTP"] + counts["LDLFP"] > 0)
		ldlps = counts["LDLTP"] * 1.0 / (counts["LDLTP"] + counts["LDLFP"]);
	out << GenoID << "\t" << counts["Leadered"] << "\t" << counts["Leaderless"] << "\t" << counts["Internal"] << "\t" << counts["Remote"] << "\t"
		<< counts["PredLeadered"] << "\t" << counts["LDDTP"] << "\t" << counts["LDDFP"] << "\t" << counts["LDDTN"] << "\t" << counts["LDDFN"] << "\t" << lddrc << "\t" << lddps << "\t"
		<< counts["PredLeaderless"] << "\t" << counts["LDLTP"] << "\t" << counts["LDLFP"] << "\t" << counts["LDLTN"] << "\t" << counts["LDLFN"] << "\t" << ldlrc << "\t" << ldlps << "\t" << endl;
	out.close();
	return "Estimation of " + GenoID + " have been conducted!";
}