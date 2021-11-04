/*
Defining the abbreviation form of STL data structures.
Module of Defining Data Structures 
Aurthor: Longshu Yang(D.E.Rommel)
Latest update: 2016/02/05
*/

#ifndef TYPEDEFBASE_H
#define TYPEDEFBASE_H

//IO streams
#include<fstream>
#include<iostream>
#include<sstream>

//STL class
#include<string>
#include<vector>
#include<utility>//pair class
#include<deque>
#include<map>
#include<set>

//STL function
#include<iomanip>
#include"assert.h"
#include<algorithm>
#include<cstdio>
#include<cmath>
#include<limits>

//IO streams
typedef std::ifstream Ifs;
typedef std::ofstream Ofs;
typedef std::ostringstream Os;
typedef std::stringstream Ss;

//String
typedef std::string Str;
typedef const std::string con_Str;

//Vector definition
typedef std::vector<int> Ve_I;
typedef std::vector<bool> Ve_B;
typedef std::vector<double> Ve_D;
typedef std::vector<std::string> Ve_Str;
typedef std::vector<long long int> Ve_LLI;
typedef std::vector<std::vector<int> > Ve_Ve_I;
typedef std::vector<std::vector<double> > Ve_Ve_D;
typedef std::vector<std::vector<std::string> > Ve_Ve_Str;
typedef std::vector<std::vector<std::vector<double> > > Ve_Ve_Ve_D;

//Pair definition
typedef std::pair<int, int> Pa_I_I;
typedef std::pair<int, bool> Pa_I_B;
typedef std::pair<char, char> Pa_C_C;
typedef std::pair<int, double> Pa_I_D;
typedef std::pair<bool, double> Pa_B_D;
typedef std::pair<double, double> Pa_D_D;
typedef std::pair<std::string, int> Pa_Str_I;
typedef std::pair<std::string, bool> Pa_Str_B;
typedef std::pair<std::string, double> Pa_Str_D;
typedef std::pair<std::string, std::string> Pa_Str_Str;
typedef std::pair<std::pair<double, double>, std::pair<double, double> > Pa_Pa_D_D;

//Vector mixed with pair
typedef std::vector<std::pair<int,int> > Ve_Pa_I_I;
typedef std::vector<std::pair<char, int> > Ve_Pa_C_I;
typedef std::vector<std::pair<bool, bool> > Ve_Pa_B_B;
typedef std::vector<std::pair<double, int> > Ve_Pa_D_I;
typedef std::vector<std::pair<double, double> > Ve_Pa_D_D;
typedef std::vector<std::pair<std::string, int> > Ve_Pa_Str_I;
typedef std::vector<std::pair<std::string, bool> > Ve_Pa_Str_B;
typedef std::vector<std::pair<std::string, double > > Ve_Pa_Str_D;
typedef std::vector<std::pair<std::string, std::string> > Ve_Pa_Str_Str;
typedef std::pair<std::vector<double>, std::vector<double> > Pa_Ve_D;
typedef std::pair<std::vector<double>, std::vector<int> > Pa_Ve_D_Ve_I;

//Deque definition
typedef std::deque<int> De_I; 
typedef std::deque<std::string> De_Str;
typedef std::deque<std::pair<std::string, std::string> > De_Pa_Str_Str;

//Map definition
typedef std::map<int, int> Ma_I_I;
typedef std::map<char, double> Ma_C_D;
typedef std::map<int, std::string> Ma_I_Str;
typedef std::map<std::string, int> Ma_Str_I;
typedef std::map<std::string, double> Ma_Str_D;
typedef std::map<std::string, std::string> Ma_Str_Str;
typedef std::vector<std::map<char, double> > Ve_Ma_C_D;
typedef std::map<int, std::map<int, std::string> > Ma_I_Ma_I_Str;
typedef std::map<std::string, std::pair<int, int> > Ma_Str_Pa_I_I;

//Set definition
typedef std::set<int> Se_I;
typedef std::set<double> Se_D;
typedef std::set<std::string> Se_Str;



typedef std::vector<std::vector< std::vector<int> > > Ve_Ve_Ve_I;


//typedef std::map<int, std::map<int, Str> > Map_I_Map_I_Str;

#endif //TYPEDEFBASE_H