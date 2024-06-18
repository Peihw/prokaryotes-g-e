#include "GenomeSeq.h"
#include"SequenceTransform.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

GenomeSeq::GenomeSeq( const con_Str& seq1 )
:positiveSeq( seq1 ), negtiveSeq( seq1 ), seqLen( seq1.size() )
{
	std::reverse( negtiveSeq.begin(), negtiveSeq.end() );	
	std::for_each( negtiveSeq.begin(), negtiveSeq.end(),	
		SequenceTransform_T::ToOppRule_T() );
}