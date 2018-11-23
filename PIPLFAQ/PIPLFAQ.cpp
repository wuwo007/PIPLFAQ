// PIPLFAQ.cpp : �������̨Ӧ�ó������ڵ㡣
/*

*/

#include "stdafx.h"
#include"CDataIO.h"
#include<time.h>
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	string ParamFilePath = Unicode2Multibyte(argv[1]);
	CParam params;
	params.Init(ParamFilePath);

	clock_t t_cleavageModelStart = clock();
	CDataIO dataIO;

	vector<CProtein> cproteins;
	vector<CPeptide> cpeptides;

	//1����MaxQuant�Ľ���ļ� peptides.txt �ж�ȡ���м��������ĶΣ��Լ����п��ܵ��ף�
	dataIO.mf_LoadPeptidesFromMaxQuant(params, cpeptides);
	//2����fasta�ļ��ж�ȡ��Щ���ܵ��׵ĵ������У�
	dataIO.mf_LoadProteins(params, cproteins, cpeptides);
	// 3���Ե���Ϊ��λ������������ȷ��ÿ���Ķ��ڵ��������г��ֵ�λ�ã�ͳ������ͬһ���������г��ֶ�ε��Ķε���Ŀ��Ŀǰ�������ĶΣ���ȡ��һ��ƥ���λ�ã���
	ProteinInfer proteininfer;
	proteininfer.m_fGetPepLocationInProteins(cproteins);

	//4����ÿһ�����ף����õ��׹���������Ψһ�Ķι鲢Ϊһ���ĶΣ�
	//�Ķ�ǿ��ȡ������������ĳһ��λ���ϳ��ֵ��Ķε��ź�ǿ�����ֵ���Ķ�����ȡ��λ���ϵ��Ķ����еĲ�����
	proteininfer.m_fMergeOverlapUniquePeptides(cproteins);
	
	string strProteinQuantInfoPath = params.m_strPIPLFAQResultPath + "//ProteinQuantInfo.txt";
	dataIO.mf_SaveProteinQuantInfo(strProteinQuantInfoPath, cproteins);

	return 0;
}

