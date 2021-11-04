#ifndef GenomeInfo_H
#define GenomeInfo_H
#include"TypeDef.h"
#include"GenomeSeq.h"
//genome sequence.
class GenomeInfo
{
public:
	GenomeInfo( const GenomeSeq& GenomeSeq );
	int getNextPhaseTIS( const char* seq, int hint = 0 );
	Str negtiveSeq;	
	Str positiveSeq;
	int seqLen;
	Ve_D background;//Genomic GC content
private:
	int getNextTISPosition( const char* seq, int hint = 0 );
	enum { boundOfORF = 90 };
};

#endif// GenomeInfo_H