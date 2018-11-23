
#include"Parameters.h"
#include<map>

#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif


class CMaxQuantIO
{
private:
	map<string, string> m_mapExperimentNameAndPeptideIntensityName;
	int m_iNumberOfExperiments;

public:
	void mf_GetAttributesName(string ExperimentDesignPath);
	// 读入肽段信息，如果bGroupUnique为真，只读入group unique 肽段；如果bGroupUnique为假，只读入group shared 肽段。
	bool mf_LoadPeptides(const CParam &param, string PeptideFilePath, vector<CPeptide>&cpeptides, bool bGroupUnique = true);
	//读入peptides.txt 中所有非污染和非反库肽段
	void mf_LoadAllPeptides(const CParam &param, string PeptideFilePath, vector<CPeptide>&cpeptides);
};

class CProteinIO
{
public:
	// 根据鉴定蛋白的ID提取蛋白的序列；
	void mf_LoadIdentifiedProteins(const CParam &param, const vector<CPeptide>&peptides, vector<CProtein> & S, map<string, vector<int>> &mapIdentProteinsAndPeptideIndex);
	bool mf_SaveProteinSequences(string Path, vector<CProtein> cproteins);
	bool mf_SaveProteinQuantInfo(string Path, vector<CProtein> cproteins);
	CProtein m_cprotein;

};

class CDataIO
{
public:
	void mf_LoadPeptidesFromMaxQuant(const CParam &param, vector<CPeptide>& Peptides);
	void mf_LoadProteins(const CParam &param, vector<CProtein> &cproteins, \
		const vector<CPeptide>& cpeptides);


	void mf_SaveIdentifiedPeptides(string strPath, const vector<CPeptide>&cPeptides);
	void mf_SaveIdentProteinGroups(string strPath, const vector<CProtein>&cProteins);

	bool mf_SaveProteinQuantInfo(string Path, vector<CProtein> cproteins);

};

class CCleavSiteIO
{
public:
	CCleavSiteIO(void);
	~CCleavSiteIO();
	bool SaveNmcSite(const string StrFileName, vector<CCleavSite> cUnmcSites);
	bool SaveMcSite(const string StrFileName, vector<CCleavSite> cMcSites);
	bool SaveRMcSite(const string StrFileName, vector<CCleavSite> cMcSites);
	bool SaveKMcSite(const string StrFileName, vector<CCleavSite> cMcSites);
};