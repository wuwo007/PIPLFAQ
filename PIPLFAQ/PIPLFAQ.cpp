// PIPLFAQ.cpp : 定义控制台应用程序的入口点。
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

	//1）从MaxQuant的结果文件 peptides.txt 中读取所有鉴定到的肽段，以及所有可能蛋白；
	dataIO.LoadPeptidesFromMaxQuant(params, cpeptides);
	//2）从fasta文件中读取这些可能蛋白的蛋白序列；
	dataIO.LoadProteins(params, cproteins, cpeptides);
	// 3）以蛋白为单位整理鉴定结果，确定每个肽段在蛋白序列中出现的位置（统计下在同一个蛋白序列出现多次的肽段的数目，目前对这种肽段，暂取第一次匹配的位置）。
	ProteinInfer proteininfer;
	proteininfer.GetPepLocationInProteins(cproteins);

	//对每一个蛋白，以唯一肽段作为节点，如果肽段之间有重叠，则两个肽段节点之间有边。
	//合并该图中所有的连通子图，将每一个连通子图归并为一个肽段，肽段强度取该连通子图中中肽段权重及其关联肽段节点权重之和的最大值，
	//即该连通子图对应的蛋白质序列上某一个位置上出现的肽段的信号强度最大值，肽段序列取该位置上的肽段序列的并集；
	proteininfer.MergeOverlapUniquePeptides(cproteins);
	
	proteininfer.CalculateUniquePepIntensitiesCV(cproteins, false);
	proteininfer.CalculateUniquePepIntensitiesCV(cproteins, true);

	string strProteinQuantInfoPath = params.m_strPIPLFAQResultPath + "//ProteinQuantInfo.txt";
	dataIO.mf_SaveProteinQuantInfo(strProteinQuantInfoPath, cproteins);

	return 0;
}

