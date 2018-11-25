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
	dataIO.LoadPeptidesFromMaxQuant(params, cpeptides);
	//2����fasta�ļ��ж�ȡ��Щ���ܵ��׵ĵ������У�
	dataIO.LoadProteins(params, cproteins, cpeptides);
	// 3���Ե���Ϊ��λ������������ȷ��ÿ���Ķ��ڵ��������г��ֵ�λ�ã�ͳ������ͬһ���������г��ֶ�ε��Ķε���Ŀ��Ŀǰ�������ĶΣ���ȡ��һ��ƥ���λ�ã���
	ProteinInfer proteininfer;
	proteininfer.GetPepLocationInProteins(cproteins);

	//��ÿһ�����ף���Ψһ�Ķ���Ϊ�ڵ㣬����Ķ�֮�����ص����������Ķνڵ�֮���бߡ�
	//�ϲ���ͼ�����е���ͨ��ͼ����ÿһ����ͨ��ͼ�鲢Ϊһ���ĶΣ��Ķ�ǿ��ȡ����ͨ��ͼ�����Ķ�Ȩ�ؼ�������Ķνڵ�Ȩ��֮�͵����ֵ��
	//������ͨ��ͼ��Ӧ�ĵ�����������ĳһ��λ���ϳ��ֵ��Ķε��ź�ǿ�����ֵ���Ķ�����ȡ��λ���ϵ��Ķ����еĲ�����
	proteininfer.MergeOverlapUniquePeptides(cproteins);
	
	proteininfer.CalculateUniquePepIntensitiesCV(cproteins, false);
	proteininfer.CalculateUniquePepIntensitiesCV(cproteins, true);

	string strProteinQuantInfoPath = params.m_strPIPLFAQResultPath + "//ProteinQuantInfo.txt";
	dataIO.mf_SaveProteinQuantInfo(strProteinQuantInfoPath, cproteins);

	return 0;
}

