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
class CNineMer      //9������
{
public:
	void AssignNineMer(string s);  //�ó���Ϊ9���ַ������¹�����9��������и�ֵ
	string nine;                //9���ӵ�����
	void Show(void);            //���9���ӣ��Լ�������ȷ�ԣ�
	void Clear(void);          //9�������
	void Copy(CNineMer ninemer);  //��ԭ9���Ӷ�������и�ֵ
};


class CProtein        //��������
{
public:
	string m_strFastaHeader;   //��¼������fasta�ļ��ĵ���ͷ
	string m_strProteinName;   //����������
	string m_strSequence;      //����������

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

	//vector<string> m_vPeptidesSequences;  //���׶�Ӧ�ļ����Ķ�
	//vector<double> m_vIdentPeptideIntensities; //���׶�Ӧ�ļ����Ķε�Intensities��
	//vector<bool> m_vbIfPeptidesShared; // ���׶�Ӧ�ļ����Ķ��Ƿ��ǹ����ģ���һ���Ķζ�Ӧ����һ�����ף��ڽ����Ķζ���Ϊ�����ģ�
	//vector<int> m_vPeptideLeftLocations; // �����Ķ�����ڵ��������е�λ�ã���0��ʼ��
	//vector<int> m_vPeptideRightLocations; // �����Ķ��Ҷ��ڵ��������е�λ��
	
	map<string, int> m_mapIndentPeptidesSC;       //���׶�Ӧ�ļ����Ķε�SC



	int m_iGroupUniquePeptidesNum;
	int m_iGroupSharedPeptidesNum;
	double m_dSequenceCoverage;

	void GetNines(int i, CNineMer& nine);  //����λ���λ���ҵ���Ӧ9���ӣ�
	void Clear(void);  // ��յ��׶�Ӧ��������Ϣ��������CProteinEnzymeInfo��clear����
	// ��ø��Ķ��ڸõ���������λ������coniFlankingRegionWidth���������������Ϣ
	string mf_GetPeptidesAdjacentSequence(string PeptideSequence)const;
	void Show(void);
};
class ProteinInfer
{
public:
	//�Ե���Ϊ��λ������������ȷ��ÿ���Ķ��ڵ��������г��ֵ�λ��
	//��ͳ������ͬһ���������г��ֶ�ε��Ķε���Ŀ��Ŀǰ�������ĶΣ���ȡ��һ��ƥ���λ�ã�,
	// ����Ψһ�ĵ�intensity CV��
	void m_fGetPepLocationInProteins(vector<CProtein>& cproteins);


	// ��ÿһ�����ף����õ��׹���������Ψһ�Ķι鲢Ϊһ���ĶΣ�
	//�Ķ�ǿ��ȡ������������ĳһ��λ���ϳ��ֵ��Ķε��ź�ǿ�����ֵ��
	//�Ķ�����ȡ��λ���ϵ��Ķ����еĲ�����
	void m_fMergeOverlapUniquePeptides(vector<CProtein>& cproteins);

};
class CPeptideEnzyme
{
public:
	//����һ���Ķ��ڵ����ж�γ��ֵ������
	size_t m_iLocNum;             //�Ķ�����Ӧ�������г��ֵĴ�����
	vector<int> m_vecLocations;    //�Ķ��ڵ����������е���λ�����ᣨø��λ�㣩����Ӧ�������е�λ�ã�

	vector<CNineMer> m_vecRNineMer;      //�Ķ��ڵ����������е�ĩλ�����ᣨø��λ�㣩��Ӧ��9����
	vector<CNineMer> m_vecLNineMer;      //�Ķ��ڵ����������е�ǰһλ�����ᣨø��λ�㣩��Ӧ��9����
	vector<double> m_vecLPros;        //�Ķ��ڵ�����������ĩλ�������Ӧ��ø�и��ʣ�
	vector<double> m_vecRPros;        //�Ķ��ڵ����������е�ǰһλ�����ᣨø��λ�㣩��Ӧ��ø�и��ʣ�
	vector<bool> m_vecbBeforeIfKRs;        //�Ķ��ڵ����������е�ǰһλ�������Ƿ�ΪK����R,����ΪTrue
	vector<bool> m_vecbBeforceIfExist; //�Ķ��ڵ����������е�ǰһλ�������Ƿ���ڻ���ǰһ���������ǵ������еĵ�һ��M�����������Ϊtrue
	bool m_bIfEndKR;         //�Ķ��ڵ����������е�ĩλ�������Ƿ�ΪK����R������Ϊtrue

	int m_iMCAANum;    //�Ķ������е�©��λ�����
	vector<int> m_vecMCLocations;   //�Ķ��ڵ����������е�©�еİ�����λ������Ӧ�������е�λ�ã�
	map<int, vector<double>> m_mapMCSitePros;   //�Ķ��ڵ����������е�©�еİ�����λ���Ӧ��9���Ӷ�Ӧ��ø�и���
	map<int, vector<CNineMer>> m_mapMCNineMers; //�Ķ��ڵ����������е�©�еİ�����λ���Ӧ��9����
	void Init(void);   // ��ʼ��
	void Clear(void); //����Ķ���Я����ø����Ϣ


};
class CPeptide           //�Ķ���
{
public:
	CPeptide();

	string m_strPeptideID; //  �Ķε�ID
	string m_strPeptideSeq; //�Ķ�����
	string m_strPep2ProteinName;  //�Ķζ�Ӧ�ĵ���������
	vector<string> m_vstrProteinsOfPeptide; //�Ķζ�Ӧ�����м���������
	vector<string> m_vstrAgentProteinsOfSharedPeptide; //shared�Ķζ�Ӧ�����м�������group�Ĵ����ף�
	double m_dPeptideIntesity;
	bool m_bProteinsShared;

	double m_dDigestionPro;    // �Ķ�ø��Ч��
	double m_dDectectionPro;  // �Ķοɼ���ԣ�
	CPeptideEnzyme m_cPeptideEnzyme;  // �Ķζ�Ӧ��ø����Ϣ��
	bool m_IfFindProtein;      //��ע�Ƿ��ҵ�����Ӧ�ĵ����ʣ������Ҳ�����Ӧ�����ʵ��Ķ��ں���㲻�ٿ��ǣ�
	int m_iSCNumber;   // �Ķζ�Ӧ��Spectrum Count
	bool m_bIdentified; //�Ƿ񱻼�������
	bool m_bIfInTrainSet;
	char m_cPrecedAA;
	char m_cPostAA;

	void AssignMCNum(void);       //����֪���е��Ķε�©��λ��������и�ֵ
	bool Location(CProtein &P); //ȷ���Ķ�ĩλ�����ᣨø��λ�㣩��Ӧ��9���Ӻ�����Ӧ�������е�λ�ã�
	void CalculPro(); // ��λ���ø�и��ʼ����Ķε�ø�и���

	bool IfmcSiteEqual(CPeptide cpeptide);
	void Show(void);
	bool Save(string Path);
	void Clear(void);  //����Ķε���Ϣ������CPeptideEnzyme��Clear������

};
enum CleaveSiteType{ Missed, Cut, None };
class CCleavSite         //ø��λ����
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
	CleaveSiteType m_CleaveType; // ��ʶλ���Ϊ��ø��λ�㻹��©��λ��
	string m_strProteinName;  //ø��λ�����ڵĵ���������
	size_t m_iLocation;       //ø��λ������Ӧ�������е�λ�ã�
	CleaveSiteMode m_cleavesitemod; //ø��λ�������K/R
	CNineMer m_cNineMer;       //ø��λ���Ӧ��9����
	size_t m_iLPeptideNum;     //ø��λ������Ķγ��ֵĴ���
	size_t m_iRPeptideNum;         //ø��λ���ұ��Ķγ��ֵĴ���

	int m_iMCSiteNum;        //ø��λ���������©�У�������©�д�����
	int m_iMcSiteNumByHMDPeptides; // �ɸ߿ɼ���Ե���δ����⵽���Ķ������©�д�����

	// �����ĸ����ݽṹ����©��λ����Ч
	vector<string> m_vsMcLeftPeptides; // ©��λ������Ķζ�Ӧ����벿���Ķ�
	vector<string> m_vsMcRightPeptides; //©��λ������Ķζ�Ӧ���Ұ벿���Ķ�
	vector<double> m_vdMcLeftPepDetect; //©��λ������Ķζ�Ӧ����벿���Ķε��Ķοɼ����
	vector<double> m_vdMcRightPepDetect;//©��λ������Ķζ�Ӧ���Ұ벿���Ķε��Ķοɼ����

	// �����������ݽṹ����ø��λ����Ч
	vector<string> m_vsNoMCLeftPeptides;     //ø��λ����˳��ֵļ����Ķ�
	vector<string> m_vsNoMCRightPeptides;    //ø��λ���Ҷ˳��ֵļ����Ķ�
	vector<string> m_vsNoMCConcatedPeptides; //ø��λ������Ҽ����ĶεĴ����ĶΣ�
	vector<double> m_vdNoMCContaedPepDetec;  //ø��λ������Ҽ����ĶεĴ����Ķε��Ķοɼ���ԣ�

	bool IfmcSiteEqual(CCleavSite candidate);  //�ж�����©��ø��λ�����Ƿ���ͬ�����ࣩ
	void Clear(void);
};//ø��λ����

double CalculateCV(const vector<double>&s);

// �Ӵ�С����
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

//�Զ�ȥ����β�Ļ��з���sep���š�
vector<string> Split(string str, char sep);

void GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns);
void GetAttributesFromFirstRow(string str, map<string, int>& mapAttributesAndColumns);

map<string, string> GetParametersFromFile(const string &ParamFilePath);
void SaveVector(const string &strPath, const vector<string>& vVec);
void SavePeptidesEnzymeInfo(string strPath, const vector<CPeptide> cpeptides);
