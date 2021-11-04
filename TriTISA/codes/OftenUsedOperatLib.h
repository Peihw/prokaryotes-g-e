//A set of functions commonly used
#ifndef OFTENUSEDOPERATLIB_H
#define OFTENUSEDOPERATLIB_H

#include"TypeDef.h"
//math function
inline int MIN( int a, int b ){
	return a > b ? b:a;
}
inline int MAX( int a, int b ){
	return a < b ? b:a;
}
inline double MAX( double a, double b ){
	return a < b ? b:a;
}
double GetDerivate( Ve_D& v, bool l = false );
double GetAver( Ve_D& v );
//triplet XTG or ORF
void x2LongestORF( const char* seq,  Pa_I_I& location );
//read and output matrix
std::istream& matrixIn(std::istream & in,M_D& matrix);
std::ostream& matrixOut(std::ostream & out,const M_D& matrix);
//File format conversion
void fileFormatConversion(Str in, Str out, Str tag);

#endif