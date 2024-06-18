#include "Matrix.h"

M_D& M_D::toBeLoged()
{
	int col=0;
	for(;col<colum;col++)
	{
		int lin=0;
		for(;lin<line;lin++)
			if(data[ lin * colum + col ] != 0 )
				data[ lin * colum + col ] = log(data[ lin * colum + col ]);
			else
				data[ lin * colum + col ] = log(0.00000001);
	}
	return *this;
}


M_D& M_D::toBeExp()
{
	int col=0;
	for(;col<colum;col++)
	{
		int lin=0;
		for(;lin<line;lin++)
			if(data[ lin * colum + col ] != 0 )
				data[ lin * colum + col ] = exp(data[ lin * colum + col ]);
	}
	return *this;
}


M_D& M_D::toBeAveraged()									
{
	int col=0;
	for(;col<colum;++col)
	{
		double sum=0;
		int lin=0;
		for(;lin<line;lin++)
			sum+=data[ lin * colum + col ];
		for(lin=0;lin<line;++lin)
			data[ lin * colum + col ]=data[ lin * colum + col ]/sum;
	}
	return *this;
}