#pragma once
#include<string>
#include<map>
#include<vector>
#include"CBasis.h"
using namespace std;
class CParam
{
public:
	void Init(string ParamFilePath);

	string m_strInputType;
	FastaType m_fastaType;
	string m_strIdentResultDirectoryPath;
	string m_strProteinSequenceFilePath;
	string m_strPIPLFAQResultPath;

	string m_strExprimentDesignPath;
	bool m_bIfExistContamProtein;
	string m_strContaminantProteinPrefix; //The prefix of Contaminant protein ID;

	int m_iAllowMinPeptideLength;
	int m_iAllowMaxPeptideLength;

	int mf_LoadParams(string ParamFilePath);
	void mf_SetDataDependParams(const vector<CPeptide>&Peptides);
	void mf_GetBasicParmas(CParam param);
	bool mf_GetResultPath(string ParamFilePath);
};

