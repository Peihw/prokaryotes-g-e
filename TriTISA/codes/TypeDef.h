#ifndef TYPEDEF_H
#define TYPEDEF_H

#include"TypeDefBase.h"

#include"Matrix.h"
typedef std::vector<M_D > Ve_M_D;
typedef std::pair< Ve_D, M_D > Pa_Ve_D_M_D;

class Location_T
{
public:
	Location_T() : isPositive( false ),dis2PreSTP(0)
	{
	}

	bool operator <( const Location_T& rhs )const
	{
		return location < rhs.location;
	}

	bool operator ==( const Location_T& rhs )const
	{
		if( isPositive != rhs.isPositive )
			return false;
		return isPositive ? location.second == rhs.location.second
			: location.first == rhs.location.first;
	}

	void location2FileForm( int seqLen )
	{
		if( isPositive )
		{
			location.first += 1;
			location.second += 3;
		}
		else
		{
			int tmp = location.first;
			location.first = abs( seqLen - location.second ) - 2;
			location.second = abs( seqLen - tmp );
		}
	}
	void fileForm2Location( int seqLen )
	{
		if( isPositive )
		{
			location.first -= 1;
			location.second -= 3;
		}
		else
		{
			int tmp = location.first;
			location.first = abs( seqLen - location.second );
			location.second = abs( seqLen - tmp ) - 2;
		}
	}
	void setOrigLocation(Pa_I_I location){
		origLocation = location;
	}
	//ORF features
	Pa_I_I location;//|A|TG........|T|TA. Start and stop codon, start may be revised later.
	Pa_I_I origLocation;//a temporal status for "location".
	bool isPositive;//is in positive strand
	Str seq;
	//start codon features
	int dis2PreSTP;
	//scores
	double trueScore, upFalseScore, downFalseScore;
};

class ResultLSortRule_T
{
public:
	bool operator () ( Location_T lhs, Location_T rhs )
	{
		if( lhs.location.second == rhs.location.second)
			return lhs.location.first < rhs.location.first;
		else
			return lhs.location.second > rhs.location.second;
	}
};

typedef std::vector< Location_T > Ve_Location;

#endif //TYPEDEF_H
