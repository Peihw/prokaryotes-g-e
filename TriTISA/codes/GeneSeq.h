
#if !defined(AFX_GENESEQ_H__75C88F65_DF9E_4F91_873F_EF99EB36E7FA__INCLUDED_)
#define AFX_GENESEQ_H__75C88F65_DF9E_4F91_873F_EF99EB36E7FA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include"TypeDef.h"
//contain genome sequence, no operation function
class GeneSeq  
{
public:
	GeneSeq(const con_Str& seq);
	const char* seq;
	Str negtiveSeq;								
	Str positiveSeq;						
	int seqLen;								
};

#endif // !defined(AFX_GENESEQ_H__75C88F65_DF9E_4F91_873F_EF99EB36E7FA__INCLUDED_)
