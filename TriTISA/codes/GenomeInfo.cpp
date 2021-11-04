#include "GenomeInfo.h"

GenomeInfo::GenomeInfo( const GenomeSeq& genomeSeq )
:background(4,0)
{
	positiveSeq = genomeSeq.positiveSeq;
	negtiveSeq = genomeSeq.negtiveSeq;
	seqLen = genomeSeq.seqLen;
	double s = 0;
	int i = 0;
	for( ; i < seqLen; ++i ){
		++background[positiveSeq[i]-'0'];
		++s;
	}
	for( i = 0; i < background.size(); ++i )
		background[i] /= s;
}

int GenomeInfo::getNextTISPosition( const char* seq, int hint )
{
	for( ; ; hint += 1 )
	{
		if( hint > seqLen - 3 )
			return - 1;
		int subStr = ((seq[hint] - 48)<<4) + ((seq[hint+1] - 48)<<2) + seq[hint + 2] - 48;
		//			ATG-032			GTG-232          TTG-332        CTG-132
		if( subStr == 14 || subStr == 46 || subStr == 62 )//|| subStr == 30)
			return hint;
	}
}

int GenomeInfo::getNextPhaseTIS( const char* seq, int hint )
{
	int TIS = hint;
	for( ; ; )
	{
		TIS = getNextTISPosition( seq, TIS ) + 1;
		if( TIS == 0 )
			return -1;
		if( ( TIS - hint - 1 ) % 3 == 0 )
			return TIS - 1;
	}
}

