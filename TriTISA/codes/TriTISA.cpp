#include"GenomeInfo.h"
#include"OftenUsedOperatLib.h"
#include"SequenceTransform.h"   
#include"time.h"

//Paramters to be calculated
M_D trueTIS, upFalse, downFalse;                              
int maxCandidateN = 0;
double pFU, pFD, pT;      
Ve_D qi;                                                                                                                                                                                                           
//Settings
int UPBP, DOWNBP, MAXORDER, ORDER;
Str INPUTFORMAT;
Set_I STARTS, STOPS;                                                                                     
//Initiate parameters
void initiateSettings();
void initiatePWMs( const GenomeInfo& genome, const Ve_Location allORFs, 
				  M_D& trueTIS, M_D& upFalse, 
				  M_D& downFalse);//
//Make prediction
void revise( GenomeInfo& genome, Ve_Location& test);
Ve_Location scoreAllCandidate( GenomeInfo& genome, Ve_Location locs );
//load and save results                                                                                            
void resultToFile( std::vector<Location_T>& result, Str resultFilename );                                                 
Ve_Location getRvsLocation( Str fileName );                                                                               
void paramtersOut();

int seqLen = 0;  //size of genome sequence 
int maximalInterationNum = 20;//The iteration will stop if more than 
                              //99.9% of the TISs remain unchanged.
                              //In rare cases, it will not converge 
                              //in such strong condition and the maximal
							  //iteration number is set as 20.
Str predictions;// = "tis_predictions.txt";
Str allCandidates;// = "all_candidates.txt";
Str paramters;// = "paramters.txt";
                                                                   
int main( int argc, char* args[] )                                                                                                       
{
	cout<<"----------------------------------------------------------------\n";
	cout<<"-------------------- WELCOME TO TriTISA ------------------------\n";
	cout<<"----------------------------------------------------------------\n";
	if( argc != 4)                                                                                                    
	{                                                                                                                 
		cout<<"usage: TriTISA.exe FASTA_genome_sequence initial_annotations "<<endl;                                                                         
		return -1;                                                                                                   
	}                                                                                                           
	initiateSettings();
	//Read genome sequence                                                                      
	Str seq, fna(args[1]);
	predictions = allCandidates = paramters = args[3];
	predictions += ".tritisa.dat";
	allCandidates += ".cand.dat";
	paramters += ".para.dat";
	SequenceTransform_T::char2FileDigitalSeq( fna, seq );                                                    
	seqLen = seq.size();                                                                                              
	GenomeSeq geneSeq( seq );  
	GenomeInfo genome(geneSeq);
	cout<<"----------------------------------------------------------------\n";
	cout<<"INPUT INFORMATION...\n";
	cout<<"  GENOME SEQUENCE: "<<args[1]<<endl;
	cout<<"    GENOME SIZE (BP): "<<seqLen<<endl;
	cout<<"    GENOMIC GC CONTENT (%): "<<std::setprecision(3)<<(100*2*genome.background[1])<<endl;
	fileFormatConversion(Str(args[2]), predictions, INPUTFORMAT );
	Ve_Location test = 	getRvsLocation(predictions.data()); 
	Ve_Location::iterator it = test.begin();
	for( ; it != test.end(); ++it ){
		it->fileForm2Location( seqLen );
	}
	cout<<"  NUMBER OF GENE: "<<test.size()<<endl;
	cout<<"----------------------------------------------------------------\n";
	cout<<"CASCADE COMBINATION OF DIFFERENT MARKOV MODELS...\n";
	for( ORDER = 0; ORDER < MAXORDER+0.1; ++ORDER ){
		if( ORDER < 1 ){
			cout<<"  POST-PROCESS ORIGNAL INPUT TISs BY A 0-TH ORDER MARKOV MODEL:\n";
		}else{
			cout<<"  POST-PROCESS OUTPUT TISs FROM A "<<(ORDER-1)<<"-TH ORDER MARKOV MODEL BY \n"
				  "                                A "<<ORDER<<"-TH ORDER MARKOV MODEL:\n";
		}
		revise( genome, test );
	}
	--ORDER;
	cout<<"----------------------------------------------------------------\n";
	cout<<"OUTPUTS...\n";
	cout<<"  BAYESIAN SCORES FOR ALL CANDIDATES: all_candidates.txt\n";
	Ve_Location all = scoreAllCandidate( genome, test );	
	resultToFile( all, allCandidates);
	it = test.begin();
	for( ; it != test.end(); ++it ){
		it->location2FileForm( seqLen );
	}
	cout<<"  TIS PREDICTION FOR EACH GENE: tis_predictions.txt\n";
	resultToFile( test, predictions );		
	cout<<"  PARAMTERS LERANED: paramters.txt\n";	
	paramtersOut();	
	cout<<"----------------------------------------------------------------\n";
	cout<<"--------------- BY hugangqing@ctb.pku.edu.cn -------------------\n";
	cout<<"---------------- OPEN UNDER GNU GPL LICENCE --------------------\n";	
	cout<<"----------------------------------------------------------------\n";	
	return 1;
}

void revise( GenomeInfo& genome, Ve_Location& test){
	cout<<"    TRAIN PARAMTERS VIA ITERATION...\n";
	cout<<std::setw(15)<<"    % TIS UNCHANGED:\n";
	for( int p = 0; p < maximalInterationNum; ++p ){
		//Initiate PWMs
		int sz = pow(4,ORDER+1);
		upFalse = M_D(sz, UPBP+3+DOWNBP,1);
		downFalse = M_D(sz, UPBP+3+DOWNBP,1);
		trueTIS = M_D( sz, UPBP+3+DOWNBP,1);
		initiatePWMs( genome, test, trueTIS, upFalse, downFalse );
		//Relocate TIS
		Ve_Location::iterator it = test.begin();
		for( ; it != test.end(); ++it ){
			const char* seq = it->isPositive ? genome.positiveSeq.data() : genome.negtiveSeq.data();
			it->origLocation = it->location;	
			x2LongestORF(seq,it->location);
			if( it->location.first < UPBP || it->location.second + DOWNBP + 1 > seqLen ){
				continue;
			}
			double maxScore = -1;
			int hint = it->location.first ;
			int indexc = 0;
			for( ; hint < it->location.second -30 && indexc < maxCandidateN; hint += 3 ){
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if( STARTS.find(subStr) != STARTS.end() ){
					++indexc;
					int beg = hint - UPBP;
					double t = 0, d = 0, u = 0;
					int index = 0;
					int f = 0;
					for( ; f < ORDER ; ++f ){
						index += pow(4,ORDER-f)*(seq[beg+f]-'0');
					}
					for( f = 0; f < 4; ++f ){
						t += exp(trueTIS(index+f,0)); 
						d += exp(downFalse(index+f,0));
						u += exp(upFalse(index+f,0));
					}
					t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
					int i = 0;
					for( ; i < trueTIS.getColum(); ++i ){
						int index = 0;
						int f = 0;
						for( ; f < ORDER ; ++f ){
							index += pow(4,ORDER-f)*(seq[beg+i+f]-'0');
						}
						double pt = 0, pt2 = 0, pd = 0, pu = 0;
						for( f = 0; f < 4; ++f ){
							pt += exp(trueTIS(index+f,i)); 
							pd += exp(downFalse(index+f,i));
							pu += exp(upFalse(index+f,i));
						}
						t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
						d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
						u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
					}
					t = exp(t), d = exp(d); u = exp(u);
					double score = (t)/(t+u+d);
					if( score > maxScore ){
						maxScore = it->trueScore = score;
						it->upFalseScore = (u)/(t+u+d);
						it->downFalseScore = (d)/(t+u+d);
						it->location.first = hint;
					}
				}
			}
		}                                            
		double q = 0, n = 0;
		for( int j = 0; j < test.size(); ++j ){
			{
				++n;
				if( test[j].origLocation == test[j].location )
					++q;
			}
		}
		cout<<"    "<<std::setw(15)<<std::setprecision(3)<<(q/n)*100<<"\n";
		if( q/n > 0.999 ){
			break;
		}                                                  
  	} 
}                                                                                                                         

Ve_Location scoreAllCandidate( GenomeInfo& genome, Ve_Location locs ){
	Ve_Location rst;
	Ve_Location::iterator it = locs.begin();
	for( it = locs.begin(); it != locs.end(); ++it ){
		const char* seq = it->isPositive ? genome.positiveSeq.data() : genome.negtiveSeq.data();
		x2LongestORF(seq,it->location);
		if( it->location.first < UPBP || it->location.second + DOWNBP + MAXORDER + 1 > seqLen ){
			continue;
		}
		int hint = it->location.first ;
		for( ; hint < it->location.second -30 ; hint += 3 ){
			int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
			if( STARTS.find(subStr) != STARTS.end() ){
				int beg = hint - UPBP;
				double t = 0, d = 0, u = 0;
				int index = 0;
				int f = 0;
				for( ; f < ORDER ; ++f ){
					index += pow(4,ORDER-f)*(seq[beg+f]-'0');
				}
				for( f = 0; f < 4; ++f ){
					t += exp(trueTIS(index+f,0)); 
					d += exp(downFalse(index+f,0));
					u += exp(upFalse(index+f,0));
				}
				t = log(t) + log(pT),  d = log(d) + log(pFD), u = log(u) + log(pFU);
				int i = 0;
				for( ; i < trueTIS.getColum(); ++i ){
					int index = 0;
					int f = 0;
					for( ; f < ORDER ; ++f ){
						index += pow(4,ORDER-f)*(seq[beg+i+f]-'0');
					}
					double pt = 0, pt2 = 0, pd = 0, pu = 0;
					for( f = 0; f < 4; ++f ){
						pt += exp(trueTIS(index+f,i)); 
						pd += exp(downFalse(index+f,i));
						pu += exp(upFalse(index+f,i));
					}
					t += trueTIS(index+(seq[beg+i+ORDER]-'0'),i) - log(pt);
					d += downFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pd);
					u += upFalse(index+(seq[beg+i+ORDER]-'0'),i) - log(pu);
				}
				t = exp(t), d = exp(d); u = exp(u);
				double score = (t)/(t+u+d);
				Location_T loc = *it;
				loc.location.first = hint;
				loc.trueScore = (t)/(t+u+d);
				loc.downFalseScore = (d)/(t+u+d);
				loc.upFalseScore = (u)/(t+u+d);
				loc.location2FileForm( seqLen );
				rst.push_back(loc);
			}
		}  
	}
	return rst;
}//*/

void initiatePWMs( const GenomeInfo& genome, const Ve_Location allORFs, 
				  M_D& trueTIS, M_D& upFalse, 
				  M_D& downFalse){
	Ve_D bkg(4,0);
	int C = 0;
	double sdN = 0, nsdN = 0;
	Ve_Location::const_iterator it = allORFs.begin();
	for( ; it != allORFs.end(); ++it ){
		if( it->location.second - it->location.first > 300 ){//To avoid the impact from false positive gene
			const Str& seq = it->isPositive ? genome.positiveSeq : genome.negtiveSeq;
			//up false TIS, regardless of frame
			if( it->location.first + DOWNBP + MAXORDER + 1< seqLen ){
				int hint = MAX(0,it->location.first - 200);
				for( ; hint < it->location.first - 3; ++hint ){
					int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
					if( subStr == 14 || subStr == 46 || subStr == 62 || subStr == 30){//)
						int beg = hint - UPBP;
						if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1 < seqLen  ){
							int i = 0;
							for( ; i < upFalse.getColum(); ++i ){
								int index = 0;
								for( int f = 0; f < ORDER + 1; ++f ){
									index += pow(4,ORDER-f)*(seq[beg+i+f]-'0');
								}
								++upFalse(index,i);							
								++bkg[seq[beg+i]-'0'];
								++C;
							}
						}
					}
				}
			}//*/
			//down false TIS
			int hint = it->location.first + 3;
			for( ; hint < it->location.second; hint += 3 ){
				int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
				if(STARTS.find(subStr) != STARTS.end()){//)
					int beg = hint - UPBP;
					if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1 < seqLen  ){
						int i = 0;
						for( ; i < upFalse.getColum(); ++i ){
							int index = 0;
							for( int f = 0; f < ORDER + 1; ++f ){
								index += pow(4,ORDER-f)*(seq[beg+i+f]-'0');
							}
							++downFalse(index,i);							
						}
					}
				}
			}
			//true TIS
			int beg = it->location.first - UPBP;
			if( beg > 0 && beg + UPBP + 3 + DOWNBP + 1< seqLen  ){
				int i = 0;
				for( ; i < trueTIS.getColum(); ++i ){
					int index = 0;
					for( int f = 0; f < ORDER + 1; ++f ){
						index += pow(4,ORDER-f)*(seq[beg+i+f]-'0');
					}
					++trueTIS(index,i);
				}
			}//*/
		}
	}
	//Normalize PWM
	trueTIS.toBeAveraged();
	upFalse.toBeAveraged();
	downFalse.toBeAveraged();
	trueTIS.toBeLoged();
	upFalse.toBeLoged();
	downFalse.toBeLoged();

	int i = 0;
	for( ; i < bkg.size(); ++i ){
		bkg[i] /= C;
	}
	qi = Ve_D();
	double P = 0;
	double ufN = 0;
	double a = bkg[0], c = bkg[1], g = bkg[2], t = bkg[3];
	for( i = 0; ; ++i ){
		double p = (a*t*g+g*t*g+t*t*g)/(a*t*g+g*t*g+t*t*g+t*g*a+t*a*a+t*a*g);
		double q = (t*g*a+t*a*a+t*a*g)/(a*t*g+g*t*g+t*t*g+t*g*a+t*a*a+t*a*g);
		P += q*pow(p,i);
		qi.push_back(q*pow(p,i));
		if( P > 0.9999 ){
			maxCandidateN = i+1;
			break;
		}else{
			ufN += i * q*pow(p,i);
		}
	}

	pT = 1;
	pFU = ufN;
	pFD = maxCandidateN - 1 - pFU;
}

void initiateSettings(){
	std::ifstream in("settings");
	if( !in.good() ){
		cout<<"Settings.txt is missing in this folder."<<endl;
		exit(0);
	}
	cout<<"SETTINGS ..."<<endl;
	while( !in.eof() ){
		Str line;
		std::getline(in,line);
		if( !line.empty() && (line.find("#") == std::string::npos) ){
			int pos = line.find("=");
			Str key = line.substr(0,pos), value = line.substr(pos+1);
			if( key == "MAXORDER" ){
				MAXORDER = atoi( value.data() );
				cout<<"  MAX ORDER OF MARKOV MODEL: "<<MAXORDER<<"\n";
			}//*/
			if( key == "UPSTREAM" ){
				UPBP = atoi( value.data() );
				cout<<"  UPSTREAM WINDOW LENGTH: "<<UPBP<<endl;
			}
			if( key == "DOWNSTREAM" ){
				DOWNBP = atoi( value.data() );
				cout<<"  DOWNSTREAM WINDOW LENGTH: "<<DOWNBP<<endl;
			}
			if( key == "START" ){
				cout<<"  START CODON: "<<value<<endl;
				value = SequenceTransform_T::char2DigitalSeq(value);
				int c = ((value[0] - 48)<<4) + ((value[1] - 48)<<2) + value[2] - 48;
				STARTS.insert( c );
			}
			if( key == "STOP" ){
				cout<<"  STOP CODON: "<<value<<endl;
				value = SequenceTransform_T::char2DigitalSeq(value);
				int c = ((value[0] - 48)<<4) + ((value[1] - 48)<<2) + value[2] - 48;
				STOPS.insert( c );
			}
			if( key == "INPUT" ){
				cout<<"  INPUT FORMAT: "<<value<<endl;
				INPUTFORMAT = value;
			}			
		}
	}
}

Ve_Location getRvsLocation( Str fileName ){                                                                               
	Ve_Location results;                                                                                              
	std::ifstream in( fileName.data() );   
	if(!in.good() ){
		cout<<"Open "<<fileName<<" errors"<<endl;
		exit(1);
	}                                                                                          
	while( !in.eof() )                                                                                                
	{                                                                                                                 
		Location_T tmp;                                                                                           
		Str c;                                                                                                    
		in>>tmp.location.first>>tmp.location.second>>c;                                                     
		if( c.empty() )                                                                                           
			break;                                                                                            
		if( c == "+" )                                                                                            
			tmp.isPositive = true;                                                                            
		else                                                                                                      
			tmp.isPositive = false;                                                                                                                                                                           
		std::getline( in, c );                                                                                  
 
		//genome is restricted to linear chrimosome 
		int len = abs(tmp.location.first-tmp.location.second);
		if( tmp.location.first > 0 && tmp.location.second > 0                                                     
			&& tmp.location.first < seqLen && tmp.location.second < seqLen 
			&& len < seqLen/2.
			&& (len % 3 == 2)
			)                                  
			results.push_back( tmp );                                                                         
	}            
	in.close();
	return results;                                                                                                   
}   
                                                                                                                       
void resultToFile( std::vector<Location_T>& result, Str resultFilename )                                                  
{                                                                                                                         
	ResultLSortRule_T sortRule;                                                                                       
	std::sort( result.begin(), result.end(), sortRule );                                                              
	std::ofstream outResult( (resultFilename).data() );                                                                                                                                                                                                      
	int i = 0;
	for(; i < result.size(); ++i )                                                                               
	{  
		outResult<<std::setw(10)<<(result[i].location.first)                                     
				<<std::setw(10)<<result[i].location.second                                        
				<<" "<<(result[i].isPositive ? '+' : '-' ) 
				<<std::setw(15)<<(result[i].location.second-result[i].location.first-2)/3
				<<std::setw(15)<<result[i].upFalseScore
				<<std::setw(15)<<result[i].trueScore
				<<std::setw(15)<<result[i].downFalseScore;
		outResult<<endl;  
	}                                                                                                                 
	outResult.close();                                                                                               
}                                                                                                                         
     
void paramtersOut(){
	ofstream out(paramters.data());
	
	out<<MAXORDER<<"-TH MARKOV MODEL FOR"<<endl;
	out<<"UPSTREAM FALSE TIS"<<endl; 
	upFalse.toBeExp();
	matrixOut(out,upFalse);
	out<<"TRUE TIS"<<endl; 
	trueTIS.toBeExp();
	matrixOut(out,trueTIS);
	out<<"DOWNSTREAM FALSE TIS"<<endl; 
	downFalse.toBeExp();
	matrixOut(out,downFalse);
	
	out<<"PROBABILITIES FOR A TOTAL OF I CANDIDATES UPSTREAM TO TRUE:"<<endl;
	int i = 0;
	for( ; i < qi.size(); ++i ){
		out<<i<<" "<<qi[i]<<endl;
	}
	out<<"THE FOLLOWING PARAMTERS ARE CALCUALTED BASED ON THE ABOVE FIGURES.\n";
	out<<"MAXIMAL NUMBER OF CANDIDATES CONSIDERED: "<<endl;
	out<<maxCandidateN<<endl;	
	out<<"PRIOR PROBABILITIES FOR UPSTREAM FALSE, TRUE, AND DOWNSTREAM FALSE TISS:"<<endl;
	double sum = pFU + pFD + 1;
	out<<pFU/sum<<" "<<1/sum<<" "<<pFD/sum<<endl;
	out.close();
}
                                                                                                                                                                                                                                                                                                                                          
