#ifndef TYPEDEFBASE_H
#define TYPEDEFBASE_H

#include<vector>
#include<utility>
#include<fstream>
#include<iostream>
#include<string>
#include<set>
#include<map>
#include<iomanip>
#include<algorithm>
#include<list>
#include<map>
#include"assert.h"
#include"math.h"

typedef std::vector<int> Ve_I;
typedef std::vector<std::string> Ve_Str;
typedef std::vector<double> Ve_D;
typedef std::vector< std::vector<double> > Ve_Ve_D;
typedef std::pair<int, int> Pa_I_I;
typedef std::pair<double, double> Pa_D_D;
using std::make_pair;

typedef std::ifstream ifstream;
typedef std::ofstream ofstream;
using std::cout;

typedef std::string Str;
typedef const std::string con_Str;

typedef std::set<int> Set_I;
typedef std::set<Str> Set_Str;
typedef std::set<double> Set_D;


typedef std::map<int,int> Map_I_I;

using std::endl;
using std::setw;


//Settings
extern int UPBP, DOWNBP, MAXORDER;
extern Str INPUTFORMAT;
extern Set_I STARTS, STOPS;
#endif //TYPEDEFBASE_H