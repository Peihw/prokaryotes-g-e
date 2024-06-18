#include"OftenUsedOperatLib.h"
//math function
double GetAver( Ve_D& v )
{
	double aver = 0;
	int i = 0;
	for( ; i < v.size(); ++i )
		aver += v[i];
	return (aver / v.size());
}

double GetDerivate( Ve_D& tmp, bool perio  )
{
	Ve_D v = tmp;
	double aver = GetAver( v );
	double der = 0;
	int i = 0;
	for( ; i < v.size(); ++i )
		der += ( v[i] - aver ) * ( v[i] - aver );
	return sqrt( der / v.size() );
}
//triplet XTG or ORF
void x2LongestORF( const char* seq, Pa_I_I& location )
{
	int STPPos = location.second;
	int ATGPos = -1;
	int tmp = ((seq[STPPos] - 48)<<4) + ((seq[STPPos+1] - 48)<<2) + seq[STPPos + 2] - 48;

	for( STPPos -= 3; ; STPPos -= 3 )
	{
		if( STPPos < 0 )
		{
			location.first = ATGPos;
			return;
		}
		int subStr = ((seq[STPPos] - 48)<<4) + ((seq[STPPos+1] - 48)<<2) + seq[STPPos + 2] - 48;
		//			TAA-300			TGA-320			TAG-302
		if( STOPS.find(subStr) != STOPS.end() ){
			location.first = ATGPos;
			return;
		}
		//			ATG-032			GTG-232          TTG-332        CTG-132
		if( STARTS.find(subStr) != STARTS.end() )
			ATGPos = STPPos;
	}
	location.first = ATGPos;
}

//read and output matrix
std::ostream& matrixOut(std::ostream & out,const M_D& matrix){
	out<<std::setw(15)<<matrix.getLine()<<std::setw(15)<<matrix.getColum()<<std::endl;
	int col=0;
	for(;col<matrix.getColum();++col)
	{
		int lin=0;
		for(;lin<matrix.getLine();++lin)
			out<<matrix(lin,col)<<' ';
		out<<std::endl;
	}
	return out;
}

std::istream& matrixIn(std::istream & in,M_D& matrix)
{
	in>>matrix.line>>matrix.colum;
	matrix.data.resize(matrix.line*matrix.colum);
	int col=0;
	for(;col<matrix.colum;++col){
		int lin=0;
		for(;lin<matrix.line;++lin)
			in>>matrix(lin,col);
	}
	return in;
}


void fileFormatConversion(Str inf, Str outf, Str tag ){
	std::ifstream in( inf.data() );
	if( !in.good() ){
		cout<<"File "<<inf<<" not found"<<endl;
		exit(0);
	}
	std::ofstream out( outf.data() );
	if( tag == "REFSEQ" || tag == "GENBANK"){
		Str line;
		while( line.find("Location") == std::string::npos 
			&& !in.eof() ){
			std::getline( in, line );
		}
		while( !in.eof() ){
			Str a_b;
			Str c;
			in>>a_b>>c;
			if( !a_b.empty() ){
				int pos = a_b.find("..");
				int a = atoi(a_b.substr(0,pos).data());
				int b = atoi(a_b.substr(pos+2).data());
				out<<a<<" "<<b<<" "<<c<<endl;
			}
			std::getline(in,c);
		}
		in.close();
	}
	if( tag == "MED" ){
		while( !in.eof() ){
			Str line;
			std::getline( in,line );
			out<<line<<endl;
		}
	}
	if( tag == "GLIMMER3" ){
		Str line;
		std::getline( in, line );
		while( !in.eof() ){
			Str a, b, s;
			in>>a>>a>>b>>s;
			if( !a.empty() ){
				if( s.find("+") == std::string::npos ){
					out<<atoi(b.data())<<" "<<atoi(a.data())<<"  -"<<endl;
				}else{
					out<<atoi(a.data())<<" "<<atoi(b.data())<<"  +"<<endl;
				}
			}
			std::getline( in,s);
		}
	}
	if( tag == "GENEMARKS" ){
		Str line;
		while( line.find("#") == std::string::npos && !in.eof() ){
			std::getline( in, line );
		}
		while( !in.eof() ){
			Str a, b, s;
			in>>s>>s>>a>>b;
			if( !a.empty() ){
				out<<atoi(a.data())<<" "<<atoi(b.data())<<" "<<s<<endl;
			}
			std::getline( in,s);
		}
	}
	out.close();
}
