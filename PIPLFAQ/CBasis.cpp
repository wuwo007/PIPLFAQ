
#include "stdafx.h"
#include"CBasis.h"
#include <sys/stat.h> 
using namespace std;

void CNineMer::AssignNineMer(string s)
{
	nine = s;
}
void CNineMer::Copy(CNineMer ninemer)
{
	nine = ninemer.nine;
}
void CNineMer::Show(void)
{
	cout << nine;
}
void CNineMer::Clear(void)
{
	nine.clear();
}


CPeptide::CPeptide()
:m_dDectectionPro(0.0),
m_dDigestionPro(0.0),
m_iSCNumber(0),
m_bIfInTrainSet(false),
m_cPostAA(' '),
m_cPrecedAA(' '),
m_dPeptideIntesity(0.0)
{
	m_cPeptideEnzyme.Init();
}
void CPeptideEnzyme::Init(void)
{
	m_iLocNum = 0;
	m_vecLocations.clear();
	m_vecRNineMer.clear();
	m_vecLNineMer.clear();
	m_vecLPros.clear();
	m_vecRPros.clear();
	m_vecbBeforceIfExist.clear();
	m_vecbBeforeIfKRs.clear();
	m_bIfEndKR = true;
	m_iMCAANum = 0;
	m_vecMCLocations.clear();
	m_mapMCSitePros.clear();
}
bool CPeptide::IfmcSiteEqual(CPeptide cpeptide)
{
	return  (m_strPeptideSeq == cpeptide.m_strPeptideSeq) && (m_strPep2ProteinName == cpeptide.m_strPep2ProteinName);
}

void CPeptide::AssignMCNum(void)
{
	m_cPeptideEnzyme.m_iMCAANum = 0;
	for (size_t i = 0; i<m_strPeptideSeq.size() - 1; i++)
	{
		if (m_strPeptideSeq[i] == 'K' || m_strPeptideSeq[i] == 'R')
		{
			m_cPeptideEnzyme.m_iMCAANum++;
		}
	}
	if (m_strPeptideSeq[m_strPeptideSeq.size() - 1] != 'K'&&m_strPeptideSeq[m_strPeptideSeq.size() - 1] != 'R')
		m_cPeptideEnzyme.m_bIfEndKR = false;          //肽段末位氨基酸不是K或者R
}


void CPeptide::CalculPro()
{

	m_dDigestionPro = m_cPeptideEnzyme.m_vecLPros[0] * m_cPeptideEnzyme.m_vecRPros[0];//  m_LeftCleavagePro*m_RightCleavagePro;

	for (size_t j = 0; j<m_cPeptideEnzyme.m_iMCAANum; j++)
	{
		m_dDigestionPro = m_dDigestionPro*(1 - m_cPeptideEnzyme.m_mapMCSitePros[0][j]);
	}

}
bool CPeptide::Save(string path)
{
	ofstream ofFile;
	ofFile.open(path.c_str(), ios::app);
	if (!ofFile) {
		cout << "Error when writing \"" << "\"." << endl;
		return false;
	}

	ofFile << "m_strPeptideSeq\t" << m_strPeptideSeq << "m_strPep2ProteinName\t" << m_strPep2ProteinName << endl;
	ofFile.close();
	return true;

}

void CPeptide::Show(void)
{
	cout << "m_strPeptideSeq: " << m_strPeptideSeq << " m_strPep2protein: " << m_strPep2ProteinName << endl;
	cout << "m_NineMer: ";
	cout << " m_iMCAANum: " << m_cPeptideEnzyme.m_iMCAANum << endl;
	cout << " m_iLocNum: " << m_cPeptideEnzyme.m_iLocNum << endl;
}
//清空肽段所携带的酶切信息
void CPeptideEnzyme::Clear(void)
{
	m_iLocNum = 0;
	m_vecLocations.clear();
	m_vecRNineMer.clear();
	m_vecLNineMer.clear();
	m_vecLPros.clear();
	m_vecRPros.clear();
	m_vecbBeforeIfKRs.clear();
	m_vecbBeforceIfExist.clear();
	m_bIfEndKR = true;
	m_iMCAANum = 0;
	m_vecMCLocations.clear();
	m_mapMCSitePros.clear();
	m_mapMCNineMers.clear();
}

void CPeptide::Clear(void)
{
	m_bProteinsShared = false;
	m_strPeptideSeq.clear();
	m_strPep2ProteinName.clear();
	m_dDigestionPro = 0.0;
	m_dDectectionPro = 0.0;
	m_IfFindProtein = false;
	m_iSCNumber = 0;
	//m_dPeptideIntensity = 0.0;
	m_cPeptideEnzyme.Clear();
	m_vstrProteinsOfPeptide.clear();
	m_vstrAgentProteinsOfSharedPeptide.clear();
	m_bIdentified = false;
	m_bIfInTrainSet = false;
	m_cPrecedAA = ' ';
	m_cPostAA = ' ';
}

void CProtein::GetNines(int i, CNineMer& nine)
{
	string Tempstr;
	nine.Clear();
	if ((i >= NEAR_NUM) && (i<m_strSequence.size() - NEAR_NUM))            //这个地方一个假设：所有的蛋白长度至少为4；  注意&&号与||号！！！！！！
	{
		Tempstr = Tempstr + m_strSequence.substr(i - NEAR_NUM, 2 * NEAR_NUM + 1);
		nine.AssignNineMer(Tempstr);

	}
	else if (i<NEAR_NUM) //已经考虑了i=-1的情况；
	{		//一个左端酶切位点;
		for (int j = 0; j<NEAR_NUM - i; j++)
			Tempstr = Tempstr + 'Z';                           //对两端不足的，以Z补齐；
		Tempstr = Tempstr + m_strSequence.substr(0, i + NEAR_NUM + 1);

		nine.AssignNineMer(Tempstr);
	}
	else if (i >= m_strSequence.size() - NEAR_NUM)
	{		//一个右端酶切位点; 
		Tempstr = Tempstr + m_strSequence.substr(i - NEAR_NUM, m_strSequence.size() - i + NEAR_NUM);
		for (size_t j = 0; j<NEAR_NUM - (m_strSequence.size() - i - 1); j++)
			Tempstr = Tempstr + 'Z';
		//cout<<Tempstr<<endl;                    /////////////////////
		nine.AssignNineMer(Tempstr);
		//m_csite.m_ninemers.push_back(TempNinemer);
	}
	//} //end if
}

void CProtein::Show(void)
{
	cout << "m_strProteinName: " << m_strProteinName << endl;
	cout << "m_strSequence: " << m_strSequence << endl;

}

void CProtein::Clear(void)
{
	m_strProteinName.clear();
	m_strSequence.clear();
	m_iGroupUniquePeptidesNum = 0;
	m_iGroupSharedPeptidesNum = 0;

	m_vUniquePeptideSequneces.clear();
	m_vSharedPeptideSequneces.clear();
	m_vUniquePeptideIntensities.clear();
	m_vSharedPeptideIntensities.clear();
	m_vUniquePeptideLeftLocations.clear();
	m_vUniquePeptideRightLocations.clear();

	m_bIfContainUniquePep = false;

	m_vMergedUniquePeptideSequneces.clear();
	m_vMergedUniquePeptideIntensities.clear();
	m_vMergedUniquePeptideLeftLocations.clear();
	m_vMergedUniquePeptideRightLocations.clear();

}

void ProteinInfer::GetPepLocationInProteins(vector<CProtein>& cproteins)
{
	int iLoc;
	int iLeftLoc = 0, iRightLoc = 0;
	int iMultiLocPeptideNum = 0;
	for (int i = 0; i < cproteins.size(); i++)
	{
		// unique peptides
		for (int j = 0; j < cproteins[i].m_vUniquePeptideSequneces.size(); j++)
		{
			iLoc = cproteins[i].m_strSequence.find(cproteins[i].m_vUniquePeptideSequneces[j]);
			if (iLoc == cproteins[i].m_strSequence.npos)
			{
				cout << "Cannot find peptide " << cproteins[i].m_vUniquePeptideSequneces[j] << " in protein " << cproteins[i].m_strSequence << endl;
				flog.mf_Input("Cannot find peptide " + cproteins[i].m_vUniquePeptideSequneces[j] + " in protein " + cproteins[i].m_strSequence + "\n");
				flog.mf_Destroy();
				exit(-1);
			}
			cproteins[i].m_vUniquePeptideLeftLocations.push_back(iLoc);
			cproteins[i].m_vUniquePeptideRightLocations.push_back(iLoc + cproteins[i].m_vUniquePeptideSequneces[j].size() - 1);

			iLoc = cproteins[i].m_strSequence.find(cproteins[i].m_vUniquePeptideSequneces[j], iLoc + 1);
			if (iLoc != cproteins[i].m_strSequence.npos)
			{
				iMultiLocPeptideNum++;
			}
		}

		// shared peptides
		for (int j = 0; j < cproteins[i].m_vSharedPeptideSequneces.size(); j++)
		{
			iLoc = cproteins[i].m_strSequence.find(cproteins[i].m_vSharedPeptideSequneces[j]);
			if (iLoc == cproteins[i].m_strSequence.npos)
			{
				cout << "Cannot find peptide " << cproteins[i].m_vSharedPeptideSequneces[j] << " in protein " << cproteins[i].m_strSequence << endl;
				flog.mf_Input("Cannot find peptide " + cproteins[i].m_vSharedPeptideSequneces[j] + " in protein " + cproteins[i].m_strSequence + "\n");
				flog.mf_Destroy();
				exit(-1);
			}
			cproteins[i].m_vSharedPeptideLeftLocations.push_back(iLoc);
			cproteins[i].m_vSharedPeptideRightLocations.push_back(iLoc + cproteins[i].m_vSharedPeptideSequneces[j].size() - 1);

			iLoc = cproteins[i].m_strSequence.find(cproteins[i].m_vSharedPeptideSequneces[j], iLoc + 1);
			if (iLoc != cproteins[i].m_strSequence.npos)
			{
				iMultiLocPeptideNum++;
			}
		}
	}
	cout << "The number of peptides which appear multiple times in one protein sequence is  " << iMultiLocPeptideNum << endl;
	flog.mf_Input("The number of peptides which appear multiple times in one protein sequence is  " + fInt2String(iMultiLocPeptideNum) + "\n");

}

bool ProteinInfer::m_fIfPeptideOverlap(int b1, int e1, int b2, int e2)
{
	if (b1 >= b2&&b1 <= e2)
		return true;
	else if (e1 >= b2&&e1 <= e2)
		return true;
	else if (b2 >= b1&&b2 <= e1)
		return true;
	else if (e2 >= b1&&e2 <= e1)
		return true;
	else
		return  false;
}
void ProteinInfer::CreatAdjTableGraph(const CProtein& cprotein, AdjTableGraph & adjtable)
{
	if (cprotein.m_vUniquePeptideSequneces.size() != cprotein.m_vUniquePeptideLeftLocations.size() ||
		cprotein.m_vUniquePeptideSequneces.size() != cprotein.m_vUniquePeptideRightLocations.size())
	{
		cout << "The dimensions of m_vUniquePeptideSequneces,m_vUniquePeptideLeftLocations and m_vUniquePeptideRightLocations of protein " 
			<< cprotein.m_strProteinName << " are not same.\n"; 
		flog.mf_Input("The dimensions of m_vUniquePeptideSequneces,m_vUniquePeptideLeftLocations and m_vUniquePeptideRightLocations of protein "
			+ cprotein.m_strProteinName + " are not same.\n");
		flog.mf_Destroy();
		exit(-1);
	}
	adjtable.adjTable.clear();
	adjtable.adjTableWeight.clear();
	for (int i = 0; i < cprotein.m_vUniquePeptideSequneces.size(); i++)
	{
		adjtable.adjTable.push_back(vector<int>{});
		adjtable.adjTableWeight.push_back(vector<double>{});
	}
	for (int i = 0; i < cprotein.m_vUniquePeptideSequneces.size(); i++)
	{
		for (int j = 0; j < cprotein.m_vUniquePeptideSequneces.size(); j++)
		{
			if (j!=i&&m_fIfPeptideOverlap(cprotein.m_vUniquePeptideLeftLocations[i], cprotein.m_vUniquePeptideRightLocations[i],
				cprotein.m_vUniquePeptideLeftLocations[j],cprotein.m_vUniquePeptideRightLocations[j]))
			{
				adjtable.adjTable[i].push_back(j);
				adjtable.adjTableWeight[i].push_back(cprotein.m_vUniquePeptideIntensities[j]);
			}
		}
	}

}

// 邻接表深度优先搜索算法迭代
vector<int> ProteinInfer::m_fAdjTableDFS(const AdjTableGraph& graph, int startNode)
{
	int vertexNum = graph.adjTable.size();
	vector<int> visited(vertexNum, 0);
	vector<int> visitOrder;
	stack<int> trace;
	trace.push(startNode);
	visited[startNode] = 1;

	while (!trace.empty())
	{
		int currentNode = trace.top();
		trace.pop();
		visitOrder.push_back(currentNode);

		if (graph.adjTable[currentNode].size() > 0)
		{
			for (size_t i = 0; i < graph.adjTable[currentNode].size(); ++i)
			{
				if (visited[graph.adjTable[currentNode][i]] == 0)
				{
					trace.push(graph.adjTable[currentNode][i]);
					visited[graph.adjTable[currentNode][i]] = 1;
				}
			}
		}
	}

	return visitOrder;
}

void ProteinInfer::GetSubGraphsByDFS(const AdjTableGraph& adjtable, vector<vector<int>>& vSubSets)
{
	vSubSets.clear();
	int iNodeNum = adjtable.adjTable.size();
	vector<bool> vbVisited(iNodeNum, false);
	for (int i = 0; i < iNodeNum; i++)
	{
		if (!vbVisited[i])
		{
			vector<int> visitOrder = m_fAdjTableDFS(adjtable, i);
			vSubSets.push_back(visitOrder);
			for (int j = 0; j < visitOrder.size(); j++)
			{
				vbVisited[visitOrder[j]] = true;
				//cout << visitOrder[j] << " ";
			}
			//cout << endl;
		}
	}
}

void ProteinInfer::CalculateNodeIntensity(const AdjTableGraph& adjtable, const vector<double>& vSelfIntensitiesOfNodes, vector<double>& vNodeIntensities)
{
	int iNodeNum = adjtable.adjTableWeight.size();
	vNodeIntensities.clear();
	vNodeIntensities.resize(iNodeNum, 0.0);

	if (vSelfIntensitiesOfNodes.size() != iNodeNum)
	{
		cout << "Error:\t the number peptide intensities is wrong.\n";
		flog.mf_Input("Error:\t the number peptide intensities is wrong.\n");
		flog.mf_Destroy();
		exit(-1);
	}

	for (int i = 0; i < iNodeNum; i++)
	{
		vNodeIntensities[i] += vSelfIntensitiesOfNodes[i];
		for (int j = 0; j < adjtable.adjTableWeight[i].size(); j++)
		{
			vNodeIntensities[i] += adjtable.adjTableWeight[i][j];
		}
	}
}
double ProteinInfer::MergeSubgraphIntensity(const vector<int>& vSubSet, const AdjTableGraph& adjtable, const vector<double>& vNodeIntensities)
{
	if (vSubSet.size() == 1)
		return vNodeIntensities[vSubSet[0]];

	double dMaxIntensity=0.0;
	for (int i = 0; i < vSubSet.size(); i++)
	{
		if (dMaxIntensity < vNodeIntensities[vSubSet[i]])
			dMaxIntensity = vNodeIntensities[vSubSet[i]];
	}
	return dMaxIntensity;
}
string ProteinInfer::MergeSubgraphsequence(const vector<int >& vSubSet, const string& proteinsequence, const vector<int>& LeftlLoc,
	const vector<int>& RightLoc, int& iLeft, int& iRight)
{
	if (vSubSet.size() == 1)
	{
		iLeft = LeftlLoc[vSubSet[0]];
		iRight = RightLoc[vSubSet[0]];
		return proteinsequence.substr(iLeft, iRight - iLeft + 1);
	}

	iLeft = proteinsequence.size(), iRight = 0;
	if (LeftlLoc.size() != RightLoc.size())
	{
		cout << "The numbe of left locations and right locations are not same.\n";
		flog.mf_Input("The numbe of left locations and right locations are not same.\n");
		flog.mf_Destroy();
		exit(-1);
	}
	for (int i = 0; i < vSubSet.size(); i++)
	{
		if (LeftlLoc[vSubSet[i]] < iLeft)
			iLeft = LeftlLoc[vSubSet[i]];
		if (RightLoc[vSubSet[i]]>iRight)
			iRight = RightLoc[vSubSet[i]];
	}

	if (iLeft >= proteinsequence.size() || iRight >= proteinsequence.size())
	{
		cout << "The location index is larger than the protein length.\n";
		flog.mf_Input("The location index is larger than the protein length.\n");
		flog.mf_Destroy();
		exit(-1);
	}
	return proteinsequence.substr(iLeft, iRight - iLeft + 1);
}


void ProteinInfer::MergeOverlapUniquePeptides(vector<CProtein>& cproteins)
{
	AdjTableGraph adjtable;
	vector<vector<int>> vSubSets;
	double dPepIntensityTemp;
	string strPepSequenceTemp;
	int iLeftLoc;
	int iRightLoc;
	for (int i = 0; i < cproteins.size(); i++)
	{
		if (cproteins[i].m_bIfContainUniquePep)
		{
			// 根据蛋白的unique肽段构建邻接表
			CreatAdjTableGraph(cproteins[i], adjtable);
			//cout << cproteins[i].m_strProteinName << endl;
			// 根据图的邻接表，求其所有连通子图
			GetSubGraphsByDFS(adjtable, vSubSets);
			vector<double> vNodeIntensities;
			CalculateNodeIntensity(adjtable, cproteins[i].m_vUniquePeptideIntensities, vNodeIntensities);
			cproteins[i].m_vMergedUniquePeptideSequneces.clear();
			cproteins[i].m_vMergedUniquePeptideLeftLocations.clear();
			cproteins[i].m_vMergedUniquePeptideRightLocations.clear();
			cproteins[i].m_vMergedUniquePeptideIntensities.clear();

			for (int j = 0; j < vSubSets.size(); j++)
			{
				dPepIntensityTemp = MergeSubgraphIntensity(vSubSets[j], adjtable, vNodeIntensities);
				strPepSequenceTemp = MergeSubgraphsequence(vSubSets[j], cproteins[i].m_strSequence,
					cproteins[i].m_vUniquePeptideLeftLocations, cproteins[i].m_vUniquePeptideRightLocations, iLeftLoc, iRightLoc);
				cproteins[i].m_vMergedUniquePeptideSequneces.push_back(strPepSequenceTemp);
				cproteins[i].m_vMergedUniquePeptideLeftLocations.push_back(iLeftLoc);
				cproteins[i].m_vMergedUniquePeptideRightLocations.push_back(iRightLoc);
				cproteins[i].m_vMergedUniquePeptideIntensities.push_back(dPepIntensityTemp);
			}
		}		
	}
}

void ProteinInfer::CalculateUniquePepIntensitiesCV(vector<CProtein>& cproteins, bool bIfMerged)
{
	if (bIfMerged)
	{
		for (int i = 0; i < cproteins.size(); i++)
		{
			cproteins[i].m_dCVOfMergedUniquePeptideIntensities = CalculateCV(cproteins[i].m_vMergedUniquePeptideIntensities);
		}
	}
	else
	{
		for (int i = 0; i < cproteins.size(); i++)
		{
			cproteins[i].m_dCVOfUniquePeptideIntensities = CalculateCV(cproteins[i].m_vUniquePeptideIntensities);

		}
	}
}

string CProtein::mf_GetPeptidesAdjacentSequence(string PeptideSequence)const
{
	string::size_type peptideLocation = 0, begin = 0, Getwidth = 0, end;
	end = PeptideSequence.find("{");
	PeptideSequence = PeptideSequence.substr(0, end);
	peptideLocation = m_strSequence.find(PeptideSequence);
	if (peptideLocation != m_strSequence.npos)
	{

		if ((peptideLocation >= coniFlankingRegionWidth) && (peptideLocation <= m_strSequence.size() - PeptideSequence.size() - coniFlankingRegionWidth))
		{
			begin = peptideLocation - coniFlankingRegionWidth;
			Getwidth = PeptideSequence.size() + 2 * coniFlankingRegionWidth;
		}
		else if (peptideLocation < coniFlankingRegionWidth)
		{
			begin = 0;
			Getwidth = peptideLocation + PeptideSequence.size() + coniFlankingRegionWidth;
		}
		else if (peptideLocation > m_strSequence.size() - PeptideSequence.size() - coniFlankingRegionWidth)
		{
			begin = peptideLocation - coniFlankingRegionWidth;
			Getwidth = m_strSequence.size() - peptideLocation - 1 + PeptideSequence.size() + coniFlankingRegionWidth;
		}
		return m_strSequence.substr(begin, Getwidth);
	}
	else
	{
		cout << "Cannot find peptide " << PeptideSequence << " in protein " << m_strSequence << endl;
		cout << "The protein sequence is " << m_strSequence << endl;
		return "NULL";
	}

}


void CCleavSite::Clear(void)
{
	m_strProteinName.clear();
	m_iLocation = 0;
	m_cNineMer.Clear();
	m_iLPeptideNum = 0;
	m_iMCSiteNum = 0;
	//m_iProteinIndex=0;
	m_iRPeptideNum = 0;
	m_CleaveType = None;
	m_vsMcLeftPeptides.clear();
	m_vsMcRightPeptides.clear();
	m_vdMcLeftPepDetect.clear();
	m_vdMcRightPepDetect.clear();

	m_vsNoMCLeftPeptides.clear();
	m_vsNoMCRightPeptides.clear();
	m_vsNoMCConcatedPeptides.clear();
	m_vdNoMCContaedPepDetec.clear();


}
bool CCleavSite::IfmcSiteEqual(CCleavSite candidate)
{
	bool Bequal = true;
	return(candidate.m_strProteinName == m_strProteinName&&candidate.m_iLocation == m_iLocation);
	//for(int i=0;i<2*NEAR_NUM+1;i++)
	//	if(candidate.m_cNineMer.nine[i]!=m_cNineMer.nine[i])
	//           Bequal=false;
	//return(Bequal);
}

// 从大到小排序
void DescendSortAndGetIndex(vector<double> &v, vector<int>&indexTemp, int left, int right)
{
	//cout << right - left << endl;
	if (left < right)
	{
		double key = v[left];
		int indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] <= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] > key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		DescendSortAndGetIndex(v, indexTemp, left, low - 1);
		DescendSortAndGetIndex(v, indexTemp, low + 1, right);
	}
}

void DescendSortAndGetIndex(double* v, vector<int>&indexTemp, int left, int right)
{
	//cout << right - left << endl;
	if (left < right)
	{
		double key = v[left];
		int indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] <= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] > key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		DescendSortAndGetIndex(v, indexTemp, left, low - 1);
		DescendSortAndGetIndex(v, indexTemp, low + 1, right);
	}
}
void AscendSortAndGetIndex(double* v, vector<int>&indexTemp, int left, int right)
{
	if (left < right)
	{
		double key = v[left];
		int indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] >= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] < key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		AscendSortAndGetIndex(v, indexTemp, left, low - 1);
		AscendSortAndGetIndex(v, indexTemp, low + 1, right);
	}
}
void AscendSortAndGetIndex(vector<double>&v, vector<int>&indexTemp, int left, int right)
{
	if (left < right)
	{
		double key = v[left];
		int indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] >= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] < key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		AscendSortAndGetIndex(v, indexTemp, left, low - 1);
		AscendSortAndGetIndex(v, indexTemp, low + 1, right);
	}

}

double AUC(double *test_targets, double* output, int length)
{
	ofstream ofile("TestAUC.txt");
	for (int i = 0; i < length; i++)
	{
		ofile << test_targets[i] << "\t" << output[i] << endl;
	}
	ofile.close();
	vector<int> vecindexTemp;
	for (int i = 0; i < length; i++)
	{
		vecindexTemp.push_back(i);
	}
	AscendSortAndGetIndex(output, vecindexTemp, 0, length - 1);

	int M = 0, N = 0;
	double sigma = 0.0;
	for (int i = 0; i < length; i++)
	{
		if (test_targets[i] == 1.0)
		{
			M++;
		}
		else
		{
			N++;
		}
	}
	for (int j = M + N; j >= 1; j--)
	{
		if (test_targets[vecindexTemp.at(j - 1)] == 1.0)
			sigma = sigma + j;
	}
	return ((sigma - (M + 1)*M / 2) / (M*N));
}

string Unicode2Multibyte(_TCHAR * chars)
{
	DWORD dwNum = WideCharToMultiByte(CP_OEMCP, NULL, chars, -1, NULL, 0, NULL, FALSE);
	char *psText;
	psText = new char[dwNum];
	if (!psText)
	{
		delete[]psText;
	}
	WideCharToMultiByte(CP_OEMCP, NULL, chars, -1, psText, dwNum, NULL, FALSE);
	string strTemp = psText;
	delete[]psText;
	psText = NULL;
	return strTemp;
}

bool is_dir(const char* path) {
	struct _stat buf = { 0 };
	_stat(path, &buf);
	return buf.st_mode & _S_IFDIR;
}

bool is_file(const char* FILENAME)
{
	fstream _file;
	_file.open(FILENAME, ios::in);
	if (!_file)
	{
		return false;
	}
	else
	{
		return true;
	}

}

bool fStringToBool(string str, bool &bl)
{
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	if (str == "true")
		bl = true;
	else if (str == "false")
		bl = false;
	else
	{
		return false;
	}
	return true;
}
string fInt2String(int i)
{
	char buf[10];
	sprintf(buf, "%d", i);
	string str = buf;
	return str;
}
string fDouble2String(double d)
{
	char buf[10];
	sprintf(buf, "%f", d);
	string str = buf;
	return str;
}

vector<string> Split(string str, char sep)
{
	vector<string> vecStrTemp;
	vecStrTemp.clear();
	int iBegin = 0, iEnd = 0;
	iEnd = str.find(sep, iBegin);
	while (iEnd != str.npos)
	{
		vecStrTemp.push_back(str.substr(iBegin, iEnd - iBegin));
		iBegin = iEnd + 1;
		iEnd = str.find(sep, iBegin);
	}
	if (str != "")
	{
		if (str[str.size() - 1] == sep || str[str.size() - 1] == '\n')
		{
			vecStrTemp.push_back(str.substr(iBegin, str.size() - iBegin - 1));
		}
		else
		{
			vecStrTemp.push_back(str.substr(iBegin, str.size() - iBegin));
		}
	}
	return vecStrTemp;
}

map<string, string> GetParametersFromFile(const string &ParamFilePath)
{
	map<string, string> mapAttributeNameAndValues;
	ifstream fin(ParamFilePath.c_str());
	if (!fin)
	{
		cout << "Error:\tCannot open parameter file: " << ParamFilePath << endl;
	}
	string strAttributeName, strAtributeValue;
	string strTemp = "";
	int iTemp1 = 0, iTemp2 = 0;
	while (getline(fin, strTemp))
	{
		if (strTemp == "")
		{
			continue;
		}
		iTemp2 = strTemp.find("=", 0);
		strAttributeName = strTemp.substr(0, iTemp2);
		iTemp1 = strTemp.find("\"");
		if (iTemp1 == strTemp.npos)
		{
			cout << "Error:\tCannot read parameter from " << strTemp << endl;
			return mapAttributeNameAndValues;
		}
		iTemp2 = strTemp.find("\"", iTemp1 + 1);
		if (iTemp2 == strTemp.npos)
		{
			cout << "Error:\tCannot read parameter from " << strTemp << endl;
			return mapAttributeNameAndValues;
		}
		strAtributeValue = strTemp.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		mapAttributeNameAndValues.insert(pair<string, string>(strAttributeName, strAtributeValue));
	}
	fin.close();
	return mapAttributeNameAndValues;
}

void SaveVector(const string &strPath, const vector<string>& vVec)
{
	ofstream oFile(strPath);
	if (!oFile)
	{
		cout << "Cannot open file " << strPath << endl;
		flog.mf_Input("Cannot open file " + strPath + "\n");
		flog.mf_Destroy();
		exit(-1);
	}
	for (int i = 0; i < vVec.size(); i++)
	{
		oFile << vVec[i] << endl;
	}
	oFile.close();
}
void SavePeptidesEnzymeInfo(string strPath, const vector<CPeptide> cpeptides)
{
	ofstream oFile(strPath);
	if (!oFile){
		cout << "Cannot open file " << strPath << endl;
		flog.mf_Input("Cannot open file " + strPath + "\n");
		flog.mf_Destroy();
		exit(-1);
	}
	else{
		cout << "Save digested peptides into " << strPath << endl;
		flog.mf_Input("Save digested peptides into " + strPath + "\n");
	}
	oFile << "PeptideSequence\tProteinName\tIfIdentified\tRNineMer\tLNineMer\tMCNum\tMCNineMers\tPepLength\n";
	map<int, vector<CNineMer>>::const_iterator mapMCNineMersIter;
	for (int i = 0; i < cpeptides.size(); i++)
	{
		oFile << cpeptides[i].m_strPeptideSeq << "\t";
		oFile << cpeptides[i].m_strPep2ProteinName << "\t";
		oFile << cpeptides[i].m_bIdentified << "\t";
		int j = 0;
		for (j = 0; j < cpeptides[i].m_cPeptideEnzyme.m_vecRNineMer.size(); j++)
		{
			oFile << cpeptides[i].m_cPeptideEnzyme.m_vecRNineMer[j].nine << ";";
		}
		oFile << "\t";
		for (j = 0; j < cpeptides[i].m_cPeptideEnzyme.m_vecLNineMer.size(); j++)
		{
			oFile << cpeptides[i].m_cPeptideEnzyme.m_vecLNineMer[j].nine << ";";
		}
		oFile << "\t";
		if (cpeptides[i].m_cPeptideEnzyme.m_mapMCNineMers.size() == 0){
			oFile << 0 << "\t";
			oFile << "\t";
		}
		else {
			mapMCNineMersIter = cpeptides[i].m_cPeptideEnzyme.m_mapMCNineMers.begin();
			oFile << mapMCNineMersIter->second.size() << "\t";
			for (int k = 0; k < mapMCNineMersIter->second.size(); k++)
			{
				oFile << mapMCNineMersIter->second[k].nine << ";";
			}
		}
		oFile << cpeptides[i].m_strPeptideSeq.size() << "\n";
	}
	oFile.close();
}

void GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns)
{
	char *pstr1;
	int icolumns = 0;

	pstr1 = strstr(pstr, "\t");
	while (pstr1 != NULL)
	{
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			cout << "Error:\tThe format of proteinGroups.txt  is wrong.\n";
			flog.mf_Input("Error:\tThe format of proteinGroups.txt  is wrong.\n");
			flog.mf_Destroy();
			exit(-11);
		}
		mapAttributesAndColumns.insert(pair<string, int>(pstr, icolumns));
		icolumns++;
		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
	}
	pstr1 = strstr(pstr, "\n");
	if (pstr1 != NULL)
	{
		*pstr1 = '\0';
	}
	else
	{
		cout << "Error:\tThe format of proteinGroups.txt  is wrong.\n";
		flog.mf_Input("Error:\tThe format of proteinGroups.txt  is wrong.\n");
		flog.mf_Destroy();
		exit(1);
	}
	mapAttributesAndColumns.insert(pair<string, int>(pstr, icolumns));
}

void GetAttributesFromFirstRow(string str, map<string, int>& mapAttributesAndColumns)
{
	int iBegin = 0, iEnd = 0;
	int icolumns = 0;
	string strSub;
	mapAttributesAndColumns.clear();
	iEnd = str.find("\t", iBegin);
	while (iEnd != str.npos)
	{
		strSub = str.substr(iBegin, iEnd - iBegin);
		mapAttributesAndColumns.insert(pair<string, int>(strSub, icolumns));
		icolumns++;
		if (iEnd == str.size() - 1)
		{
			break;
		}
		iBegin = iEnd + 1;
		iEnd = str.find("\t", iBegin);
	}
	strSub = str.substr(iBegin, str.size() - iBegin - 1);
	mapAttributesAndColumns.insert(pair<string, int>(strSub, icolumns));

}

double CalculateCV(const vector<double>&s)
{  
	vector<double> nonZeroS;
	nonZeroS.clear();
	for (int i = 0; i < s.size(); i++)
	{
		if (s[i] != 0.0)
		{
			nonZeroS.push_back(s[i]);
		}
	}

	double mean = 0.0;
	double sd = 0.0, sd2 = 0.0;
	int len = nonZeroS.size();


	for (int i = 0; i < len; i++)
	{
		mean += (nonZeroS[i]);
	}

	if (len <= 1)
		return 0.0;
	else
	{
		mean = mean / len;
		for (int i = 0; i < len; i++)
		{
			sd += ((nonZeroS[i]) - mean)*((nonZeroS[i]) - mean);
		}
	}

	sd = sqrt(sd / len);

	if (mean <= 0)
		return -1;
	else	return sd / mean;
}
