#ifndef MAXTRIX_H
#define MAXTRIX_H

#include"TypeDefBase.h"
				
class M_D 
{
public:
	M_D(unsigned rows, unsigned cols,double init=0.):
	  line(rows),colum(cols),data(rows*cols, init ){};
	M_D(){}
	inline double& operator()(unsigned i,unsigned j)							//operator () 
	{
//		assert(i<line&&j<colum);
		return data[i*colum+j];
	}
	inline double operator()(unsigned i,unsigned j)const					//operator() 
	{			
//		assert(i<line&&j<colum);
		return data[i*colum+j];
	}
	M_D& toBeAveraged();
	M_D& toBeLoged();
	M_D& toBeExp();

	inline int getLine()const {return line;}	
	inline int getColum()const {return colum;}
	M_D(std::vector<double>& d,unsigned rows, unsigned cols)
		:line(rows),colum(cols),data(d){};
	std::vector<double> data;
    unsigned line, colum;
};

#endif//difine MAXTRIX_H.



