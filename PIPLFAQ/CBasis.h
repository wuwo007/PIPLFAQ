#pragma once
#include <map>
#include <string>
#include <vector>
#include<regex>

using namespace std;
#define BUFFERLENGTH 200000 
const int coniFlankingRegionWidth = 15;
enum DataType{ MaxQuantTpye = 0, pfindType = 1, PANDAType = 2, mzQuantMLType = 3 };
typedef std::regex FastaType;
#define MAX_PEPTIDE_LENGTH 128
#define NEAR_NUM 4
const int Obs = 4;
const double Threhold = 0.4;
class CExperimentInfo
{
public:
	vector<string> vecstrExperimentNames;
	int iExperimentNumbers;
	void GetExperimentInfo();
};
enum CleaveSiteMode
{
	CleaveSiteR = 0, CleaveSiteK = 1, CleavageSiteNon = 2,
};
class CNineMer      //9连子类
{
public:
	void AssignNineMer(string s);  //用长度为9的字符串对新构建的9连子类进行赋值
	string nine;                //9连子的内容
	void Show(void);            //输出9连子，以检验其正确性；
	void Clear(void);          //9连子清空
	void Copy(CNineMer ninemer);  //用原9连子对新类进行赋值
};


class CProtein        //蛋白质类
{
public:
	string m_strFastaHeader;   //记录来自于fasta文件的蛋白头
	string m_strProteinName;   //蛋白质名称
	string m_strSequence;      //蛋白质序列

	vector<string> m_vUniquePeptideSequneces;
	vector<string> m_vSharedPeptideSequneces;
	vector<double> m_vUniquePeptideIntensities;
	vector<double> m_vSharedPeptideIntensities;
	vector<int> m_vUniquePeptideLeftLocations;
	vector<int> m_vUniquePeptideRightLocations;
	vector<int> m_vSharedPeptideLeftLocations;
	vector<int> m_vSharedPeptideRightLocations;

	bool m_bIfContainUniquePep;

	string m_strAgentUniquePepSequence;
	double m_dAgentUniquePepIntensity;
	int m_iPepNumInMaxSet;

	//vector<string> m_vPeptidesSequences;  //蛋白对应的鉴定肽段
	//vector<double> m_vIdentPeptideIntensities; //蛋白对应的鉴定肽段的Intensities；
	//vector<bool> m_vbIfPeptidesShared; // 蛋白对应的鉴定肽段是否是共享肽；若一个肽段对应超过一个蛋白，在将该肽段定义为共享肽；
	//vector<int> m_vPeptideLeftLocations; // 鉴定肽段左端在蛋白序列中的位置，从0开始；
	//vector<int> m_vPeptideRightLocations; // 鉴定肽段右端在蛋白序列中的位置
	
	map<string, int> m_mapIndentPeptidesSC;       //蛋白对应的鉴定肽段的SC



	int m_iGroupUniquePeptidesNum;
	int m_iGroupSharedPeptidesNum;
	double m_dSequenceCoverage;

	void GetNines(int i, CNineMer& nine);  //根据位点的位置找到相应9连子；
	void Clear(void);  // 清空蛋白对应的属性信息，并调用CProteinEnzymeInfo的clear函数
	// 获得该肽段在该蛋白中所在位置左右coniFlankingRegionWidth个氨基酸的序列信息
	string mf_GetPeptidesAdjacentSequence(string PeptideSequence)const;
	void Show(void);
};
class ProteinInfer
{
public:
	//以蛋白为单位整理鉴定结果，确定每个肽段在蛋白序列中出现的位置
	//（统计下在同一个蛋白序列出现多次的肽段的数目，目前对这种肽段，暂取第一次匹配的位置）,
	// 计算唯一肽的intensity CV。
	void m_fGetPepLocationInProteins(vector<CProtein>& cproteins);


	// 对每一个蛋白，将该蛋白关联的所有唯一肽段归并为一个肽段，
	//肽段强度取蛋白质序列上某一个位置上出现的肽段的信号强度最大值，
	//肽段序列取该位置上的肽段序列的并集。
	void m_fMergeOverlapUniquePeptides(vector<CProtein>& cproteins);

};
class CPeptideEnzyme
{
public:
	//考虑一个肽段在蛋白中多次出现的情况；
	size_t m_iLocNum;             //肽段在相应蛋白质中出现的次数；
	vector<int> m_vecLocations;    //肽段在蛋白质序列中的首位氨基酸（酶切位点）在相应蛋白质中的位置；

	vector<CNineMer> m_vecRNineMer;      //肽段在蛋白质序列中的末位氨基酸（酶切位点）对应的9连子
	vector<CNineMer> m_vecLNineMer;      //肽段在蛋白质序列中的前一位氨基酸（酶切位点）对应的9连子
	vector<double> m_vecLPros;        //肽段在蛋白质序列中末位氨基酸对应的酶切概率；
	vector<double> m_vecRPros;        //肽段在蛋白质序列中的前一位氨基酸（酶切位点）对应的酶切概率；
	vector<bool> m_vecbBeforeIfKRs;        //肽段在蛋白质序列中的前一位氨基酸是否为K或者R,正常为True
	vector<bool> m_vecbBeforceIfExist; //肽段在蛋白质序列中的前一位氨基酸是否存在或者前一个氨基酸是蛋白序列的第一个M的情况，正常为true
	bool m_bIfEndKR;         //肽段在蛋白质序列中的末位氨基酸是否为K或者R，正常为true

	int m_iMCAANum;    //肽段序列中的漏切位点个数
	vector<int> m_vecMCLocations;   //肽段在蛋白质序列中的漏切的氨基酸位点在相应蛋白质中的位置；
	map<int, vector<double>> m_mapMCSitePros;   //肽段在蛋白质序列中的漏切的氨基酸位点对应的9连子对应的酶切概率
	map<int, vector<CNineMer>> m_mapMCNineMers; //肽段在蛋白质序列中的漏切的氨基酸位点对应的9连子
	void Init(void);   // 初始化
	void Clear(void); //清空肽段所携带的酶切信息


};
class CPeptide           //肽段类
{
public:
	CPeptide();

	string m_strPeptideID; //  肽段的ID
	string m_strPeptideSeq; //肽段序列
	string m_strPep2ProteinName;  //肽段对应的蛋白质名称
	vector<string> m_vstrProteinsOfPeptide; //肽段对应的所有鉴定蛋白名
	vector<string> m_vstrAgentProteinsOfSharedPeptide; //shared肽段对应的所有鉴定蛋白group的代表蛋白；
	double m_dPeptideIntesity;
	bool m_bProteinsShared;

	double m_dDigestionPro;    // 肽段酶切效率
	double m_dDectectionPro;  // 肽段可检测性；
	CPeptideEnzyme m_cPeptideEnzyme;  // 肽段对应的酶切信息；
	bool m_IfFindProtein;      //标注是否找到了相应的蛋白质，对于找不到相应蛋白质的肽段在后面便不再考虑；
	int m_iSCNumber;   // 肽段对应的Spectrum Count
	bool m_bIdentified; //是否被鉴定到了
	bool m_bIfInTrainSet;
	char m_cPrecedAA;
	char m_cPostAA;

	void AssignMCNum(void);       //对已知序列的肽段的漏切位点个数进行赋值
	bool Location(CProtein &P); //确定肽段末位氨基酸（酶切位点）对应的9连子和在相应蛋白质中的位置；
	void CalculPro(); // 由位点的酶切概率计算肽段的酶切概率

	bool IfmcSiteEqual(CPeptide cpeptide);
	void Show(void);
	bool Save(string Path);
	void Clear(void);  //清空肽段的信息并调用CPeptideEnzyme的Clear函数；

};
enum CleaveSiteType{ Missed, Cut, None };
class CCleavSite         //酶切位点类
{
public:
	CCleavSite()
		:m_strProteinName(""),
		m_iLocation(0),
		m_iLPeptideNum(0),
		m_iRPeptideNum(0),
		m_cleavesitemod(CleavageSiteNon),
		m_CleaveType(None)
	{}
	CleaveSiteType m_CleaveType; // 标识位点归为了酶切位点还是漏切位点
	string m_strProteinName;  //酶切位点所在的蛋白质名称
	size_t m_iLocation;       //酶切位点在相应蛋白质中的位置；
	CleaveSiteMode m_cleavesitemod; //酶切位点的类型K/R
	CNineMer m_cNineMer;       //酶切位点对应的9连子
	size_t m_iLPeptideNum;     //酶切位点左边肽段出现的次数
	size_t m_iRPeptideNum;         //酶切位点右边肽段出现的次数

	int m_iMCSiteNum;        //酶切位点如果发生漏切，发生的漏切次数；
	int m_iMcSiteNumByHMDPeptides; // 由高可检测性但是未被检测到的肽段引起的漏切次数；

	// 以下四个数据结构仅对漏切位点有效
	vector<string> m_vsMcLeftPeptides; // 漏切位点鉴定肽段对应的左半部分肽段
	vector<string> m_vsMcRightPeptides; //漏切位点鉴定肽段对应的右半部分肽段
	vector<double> m_vdMcLeftPepDetect; //漏切位点鉴定肽段对应的左半部分肽段的肽段可检测性
	vector<double> m_vdMcRightPepDetect;//漏切位点鉴定肽段对应的右半部分肽段的肽段可检测性

	// 以下两个数据结构仅对酶切位点有效
	vector<string> m_vsNoMCLeftPeptides;     //酶切位点左端出现的鉴定肽段
	vector<string> m_vsNoMCRightPeptides;    //酶切位点右端出现的鉴定肽段
	vector<string> m_vsNoMCConcatedPeptides; //酶切位点的左右鉴定肽段的串联肽段；
	vector<double> m_vdNoMCContaedPepDetec;  //酶切位点的左右鉴定肽段的串联肽段的肽段可检测性；

	bool IfmcSiteEqual(CCleavSite candidate);  //判断两个漏切酶切位点类是否相同（冗余）
	void Clear(void);
};//酶切位点类

double CalculateCV(const vector<double>&s);

// 从大到小排序
void DescendSortAndGetIndex(vector<double> &v, vector<int>&indexTemp, int left, int right);
void DescendSortAndGetIndex(double* v, vector<int>&indexTemp, int left, int right);

void AscendSortAndGetIndex(double* v, vector<int>&indexTemp, int left, int right);
void AscendSortAndGetIndex(vector<double>&v, vector<int>&indexTemp, int left, int right);

double AUC(double *test_targets, double* output, int length);
string Unicode2Multibyte(_TCHAR *chars);
bool is_dir(const char* path);
bool is_file(const char* FILENAME);
bool fStringToBool(string str, bool &bl);
string fInt2String(int i);
string fDouble2String(double d);

//自动去掉行尾的换行符和sep符号。
vector<string> Split(string str, char sep);

void GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns);
void GetAttributesFromFirstRow(string str, map<string, int>& mapAttributesAndColumns);

map<string, string> GetParametersFromFile(const string &ParamFilePath);
void SaveVector(const string &strPath, const vector<string>& vVec);
void SavePeptidesEnzymeInfo(string strPath, const vector<CPeptide> cpeptides);
