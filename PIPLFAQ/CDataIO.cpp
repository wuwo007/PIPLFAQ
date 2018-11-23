#include"stdafx.h"
#include"CDataIO.h"
#include<sys/stat.h>
#include<regex>
using namespace std;
using std::fclose;

void CDataIO::mf_SaveIdentifiedPeptides(string strPath, const vector<CPeptide>&cPeptides)
{
	ofstream oFile(strPath.c_str());
	if (!oFile)
	{
		cout << "Cannot open " << strPath << endl;
		flog.mf_Input("Cannot open " + strPath + "\n");
		exit(-1);
	}
	else
	{
		cout << "Save identified peptides into " << strPath << endl;
		flog.mf_Input("Save identified peptides into " + strPath + "\n");
	}
	oFile << "Peptide Sequence\tProtein Id\tSC\n";
	int iPepSize = cPeptides.size();
	for (int i = 0; i <iPepSize; i++)
	{
		oFile << cPeptides.at(i).m_strPeptideSeq << "\t";
		oFile << cPeptides[i].m_strPep2ProteinName << "\t";
		oFile << cPeptides[i].m_iSCNumber << endl;
	}

	oFile.close();
}
void CDataIO::mf_SaveIdentProteinGroups(string strPath, const vector<CProtein>&cProteins)
{
	ofstream oFile(strPath.c_str());
	if (!oFile){
		cout << "Cannot open " << strPath << endl;
		flog.mf_Input("Cannot open " + strPath + "\n");
		exit(-1);
	}
	else{
		cout << "Save identified protein groups into " << strPath << endl;
		flog.mf_Input("Save identified protein groups into " + strPath + "\n");
	}

	oFile << "Protein ID\tNumber of Group unique peptides\tNumber of Group shared peptides\tSequence Coverage\n";
	for (size_t i = 0; i < cProteins.size(); i++){
		oFile << cProteins[i].m_strProteinName << "\t";
		oFile << cProteins[i].m_iGroupUniquePeptidesNum << "\t";
		oFile << cProteins[i].m_iGroupSharedPeptidesNum << "\t";
		oFile << cProteins[i].m_dSequenceCoverage << "\n";
	}

	oFile.close();
}
/*mf_GetAttributesName
Because some attribute name change with the experiment design,
so we need determine the attribute name according to the experiment design file.
*/
void CMaxQuantIO::mf_GetAttributesName(string ExperimentDesignPath)
{
	map<string, int> mapExperiments;
	map<string, int>::iterator mapExperimentsIter;
	ifstream fin(ExperimentDesignPath.c_str());
	if (!fin)
	{
		cout << "Cannot open file " << ExperimentDesignPath << endl;
		flog.mf_Input("Error:\tCannot open file " + ExperimentDesignPath + "\n");
		flog.mf_Destroy();
		exit(2);
	}
	string strTemp1;
	m_iNumberOfExperiments = 0;
	getline(fin, strTemp1);//jump the first row.
	vector<string> vStrTemps;
	vStrTemps = Split(strTemp1, '\t');
	int iExperimentColumn = 0;
	int i = 0;
	int iExpNum = vStrTemps.size();
	for (i = 0; i < iExpNum; i++)
	{
		if (vStrTemps[i] == "Experiment")
		{
			iExperimentColumn = i;
			break;
		}
	}
	if (i == vStrTemps.size())
	{
		cout << "Cannot find Experiment column in the experimentDesignTempate.txt\n";
		flog.mf_Input("Cannot find Experiment column in the experimentDesignTempate.txt\n");
		flog.mf_Destroy();
		exit(-1);
	}

	string strTemp2;
	while (getline(fin, strTemp1))
	{
		vStrTemps = Split(strTemp1, '\t');
		strTemp2 = vStrTemps[iExperimentColumn];
		mapExperiments.insert(pair<string, int>(strTemp2, 0));
	}
	m_iNumberOfExperiments = mapExperiments.size();

	string strPeptideIntensityName;
	string striBAQ_IntensityName;
	string strLFQ_IntensityName;


	mapExperimentsIter = mapExperiments.begin();
	for (; mapExperimentsIter != mapExperiments.end(); mapExperimentsIter++)
	{
		strPeptideIntensityName = "Intensity " + mapExperimentsIter->first;
		m_mapExperimentNameAndPeptideIntensityName.insert(pair<string, string>(mapExperimentsIter->first, strPeptideIntensityName));
		//striBAQ_IntensityName = "iBAQ " + mapExperimentsIter->first;
		//m_mapExperimentNameAndiBAQIntensityName.insert(pair<string, string>(mapExperimentsIter->first, striBAQ_IntensityName));
		//if (mapExperimentsIter->first != "")
		//{
		//	strLFQ_IntensityName = "LFQ intensity " + mapExperimentsIter->first;
		//	m_mapExperimentNameAndLFQIntensityName.insert(pair<string, string>(mapExperimentsIter->first, strLFQ_IntensityName));
		//}
	}

	fin.close();
}

//Load peptide sequence and its corresponding protein from peptides.txt
bool CMaxQuantIO::mf_LoadPeptides(const CParam &param, string strPeptideFilePath, vector<CPeptide>&cpeptides, bool bGroupUnique)
{
	if (bGroupUnique)
	{
		cout << "Load group unique peptides from " << strPeptideFilePath << "\n";
		flog.mf_Input("Load group unique peptides from " + strPeptideFilePath + "\n");
	}
	else
	{
		cout << "Load group shared peptides from " << strPeptideFilePath << "\n";
		flog.mf_Input("Load group shared peptides from " + strPeptideFilePath + "\n");
	}

	// open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, strPeptideFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Cannot open " << strPeptideFilePath << endl;
		flog.mf_Input("Error:\tCannot open " + strPeptideFilePath + ". The path of Peptides.txt is needed here.\n");
		flog.mf_Destroy();
		exit(1);
	}

	char Buffer[BUFFERLENGTH];
	char *pstr;
	int iSequencecolumnNum = 0, iIDcolumnNum = 0, iProteinsNamecolumnNum = 0;
	int iGroupUniqueColumnNum = 0;
	int iReverseColumnNum = 0;
	int iContaminantColumnNum = 0;
	int iSCcolumnNum = 0;
	vector<int> veciIntensitycolumnNum;
	int icolumns = 0, count = 0;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;

	//get the columns of sequence°¢ID°¢Proteins°¢Intensity  by analysising the first row;
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns);

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Sequence");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iSequencecolumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Cannot find the column \"Sequence\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Sequence\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("id");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iIDcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Cannot find the column \"id\" in the peptides.txt." << endl;
		flog.mf_Input("Error:\tCannot find the column \"id\" in the peptides.txt.\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Proteins");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinsNamecolumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Cannot find the column \"Proteins\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Proteins\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Unique (Groups)");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iGroupUniqueColumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Cannot find the column \"Unique (Groups)\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Unique (Groups)\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Reverse");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iReverseColumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Cannot find the column \"Reverse\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Reverse\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Potential contaminant");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iContaminantColumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Cannot find the column \"Potential contaminant\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Potential contaminant\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	map<string, string>::iterator mapExperimentAndPeptideIntensityIter;
	mapExperimentAndPeptideIntensityIter = m_mapExperimentNameAndPeptideIntensityName.begin();
	for (; mapExperimentAndPeptideIntensityIter != m_mapExperimentNameAndPeptideIntensityName.end(); \
		mapExperimentAndPeptideIntensityIter++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find(\
			mapExperimentAndPeptideIntensityIter->second);
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			veciIntensitycolumnNum.push_back(mapAttrtibuteAndcolumnsIter->second);
			count++;
		}
		else
		{
			cout << "Error:\tCannot find the column \"" << \
				mapExperimentAndPeptideIntensityIter->second << "\" in the peptides.txt.";
			flog.mf_Input("Error:\tCannot find the column \"" + \
				mapExperimentAndPeptideIntensityIter->second + "\" in the peptides.txt.");
			flog.mf_Destroy();
			exit(1);
		}
	}

	/*mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("MS/MS Count");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
	iSCcolumnNum = mapAttrtibuteAndcolumnsIter->second;
	count++;
	}
	else
	{
	cout << "Cannot find the column \"MS/MS Count\" in the peptides.txt" << endl;
	flog.mf_Input("Error:\tCannot find the column \"MS/MS Count\" in the peptides.txt\n");
	flog.mf_Destroy();
	exit(1);
	}*/
	// extract spectra number from MS/MS IDs 
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("MS/MS IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iSCcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		count++;
	}
	else
	{
		cout << "Cannot find the column \"MS/MS IDs\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"MS/MS IDs\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	CPeptide peptideTemp;
	vector<string> vecStrTemp;
	vector<string> vecProteinNameTemps;
	vector<string>::iterator vecProteinNameIter;
	//vector<double> vecdIntensityTemp;
	string strLine;
	string strSequenceTemp;
	string strProteinNameTemp;
	string strProteinNameSuffix;
	string strReverseTemp;
	string strContaminantTemp;
	string strGroupUniqueTemp;
	vector<string> vecMSMStemp;

	bool bIfDelete;
	string strIDTemp;
	fgets(Buffer, BUFFERLENGTH, pFile);
	strLine = Buffer;
	vecStrTemp = Split(strLine, '\t');
	while (!feof(pFile))
	{
		bIfDelete = false;
		strSequenceTemp = vecStrTemp.at(iSequencecolumnNum);
		strIDTemp = vecStrTemp.at(iIDcolumnNum);
		strProteinNameTemp = vecStrTemp.at(iProteinsNamecolumnNum);
		strReverseTemp = vecStrTemp.at(iReverseColumnNum);
		strContaminantTemp = vecStrTemp.at(iContaminantColumnNum);
		strGroupUniqueTemp = vecStrTemp.at(iGroupUniqueColumnNum);

		// delete the peptide, if it is a decoy peptide
		if (strReverseTemp != "" || strContaminantTemp != "")
		{
			bIfDelete = true;
			fgets(Buffer, BUFFERLENGTH, pFile);
			pstr = Buffer;
			strLine = Buffer;
			vecStrTemp = Split(strLine, '\t');
			continue;
		}
		vecProteinNameTemps = Split(strProteinNameTemp, ';');
		//delete the contaminitant proteins names
		if (param.m_bIfExistContamProtein)
		{
			vecProteinNameIter = vecProteinNameTemps.begin();
			for (; vecProteinNameIter != vecProteinNameTemps.end();)
			{
				strProteinNameSuffix = (*vecProteinNameIter).substr(0, param.m_strContaminantProteinPrefix.length());

				if (strProteinNameSuffix == param.m_strContaminantProteinPrefix)
				{
					vecProteinNameIter = vecProteinNameTemps.erase(vecProteinNameIter);
				}
				else
				{
					vecProteinNameIter++;
				}
			}
			if (vecProteinNameTemps.size() == 0)
			{
				bIfDelete = true;
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				strLine = Buffer;
				vecStrTemp = Split(strLine, '\t');
				continue;
			}
		}

		if (bGroupUnique&&strGroupUniqueTemp == "no")
		{
			bIfDelete = true;
			fgets(Buffer, BUFFERLENGTH, pFile);
			pstr = Buffer;
			strLine = Buffer;
			vecStrTemp = Split(strLine, '\t');
			continue;
		}

		if (!bGroupUnique&&strGroupUniqueTemp == "yes")
		{
			bIfDelete = true;
			fgets(Buffer, BUFFERLENGTH, pFile);
			pstr = Buffer;
			strLine = Buffer;
			vecStrTemp = Split(strLine, '\t');
			continue;
		}

		//vecdIntensityTemp.clear();
		//int iIntSize = veciIntensitycolumnNum.size();
		//for (int i = 0; i < iIntSize; i++)
		//{
		//	vecdIntensityTemp.push_back(atof((vecStrTemp.at(veciIntensitycolumnNum.at(i))).c_str()));
		//}

		if (!bIfDelete)
		{
			peptideTemp.Clear();
			peptideTemp.m_strPeptideSeq = strSequenceTemp;
			peptideTemp.m_strPeptideID = strIDTemp;
			vecProteinNameIter = vecProteinNameTemps.begin();
			if (vecProteinNameIter == vecProteinNameTemps.end())
			{
				cout << "Warning:\tPeptide " << strSequenceTemp << \
					" do not have corresponding proteins.\n";
				flog.mf_Input("Warning:\tPeptide " + strSequenceTemp + \
					" do not have corresponding proteins.\n");
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				strLine = Buffer;
				vecStrTemp = Split(strLine, '\t');
				continue;
			}

			peptideTemp.m_strPep2ProteinName = *vecProteinNameIter;
			for (; vecProteinNameIter != vecProteinNameTemps.end(); vecProteinNameIter++)
			{
				peptideTemp.m_vstrProteinsOfPeptide.push_back(*vecProteinNameIter);
			}
			//peptideTemp.m_dPeptideIntensity = vecdIntensityTemp[0];
			vecMSMStemp.clear();
			vecMSMStemp = Split(vecStrTemp[iSCcolumnNum], ';');
			peptideTemp.m_iSCNumber = vecMSMStemp.size();
			/*peptideTemp.m_iSCNumber = atoi(vecStrTemp[iSCcolumnNum].c_str());*/

			cpeptides.push_back(peptideTemp);

			//if (peptideTemp.m_dPeptideIntensity != 0.0)
			//{
			//	cpeptides.push_back(peptideTemp);
			//}
		}
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
		strLine = Buffer;
		vecStrTemp = Split(strLine, '\t');
	}
	fclose(pFile);
	if (bGroupUnique)
	{
		cout << "\tLoaded " << cpeptides.size() << " group unique peptides.\n";
		flog.mf_Input("\tLoaded " + fInt2String(cpeptides.size()) + " group unique peptides.\n");
	}
	else
	{
		cout << "\tLoaded " << cpeptides.size() << " group shared peptides.\n";
		flog.mf_Input("\tLoaded " + fInt2String(cpeptides.size()) + " group shared peptides.\n");
	}

	return 1;

}

void CMaxQuantIO::mf_LoadAllPeptides(const CParam &param, string strPeptideFilePath, vector<CPeptide>&cpeptides)
{
	// open the file
	FILE * pFile;
	errno_t err;
	err = fopen_s(&pFile, strPeptideFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Cannot open " << strPeptideFilePath << endl;
		flog.mf_Input("Error:\tCannot open " + strPeptideFilePath + ". The path of Peptides.txt is needed here.\n");
		flog.mf_Destroy();
		exit(1);
	}

	char Buffer[BUFFERLENGTH];
	char *pstr;
	int iSequencecolumnNum = 0, iProteinsNamecolumnNum = 0;
	int iReverseColumnNum = 0;
	int iContaminantColumnNum = 0;
	int iIntensityColumnNum = 0;
	int icolumns = 0;

	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;

	//get the columns of sequence°¢ID°¢Proteins°¢Intensity  by analysising the first row;
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns);

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Sequence");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iSequencecolumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Cannot find the column \"Sequence\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Sequence\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Proteins");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinsNamecolumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Cannot find the column \"Proteins\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Proteins\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Reverse");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iReverseColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Cannot find the column \"Reverse\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Reverse\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Potential contaminant");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iContaminantColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Cannot find the column \"Potential contaminant\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Potential contaminant\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Intensity");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iIntensityColumnNum = mapAttrtibuteAndcolumnsIter->second;
	}
	else
	{
		cout << "Cannot find the column \"Intensity\" in the peptides.txt" << endl;
		flog.mf_Input("Error:\tCannot find the column \"Intensity\" in the peptides.txt\n");
		flog.mf_Destroy();
		exit(1);
	}

	CPeptide peptideTemp;
	vector<string> vecStrTemp;
	vector<string> vecProteinNameTemps;
	vector<string>::iterator vecProteinNameIter;
	string strLine;
	string strSequenceTemp;
	string strProteinNameTemp;
	string strReverseTemp;
	string strContaminantTemp;
	string strProteinNameSuffix;
	double dPeptideIntensity;

	bool bIfDelete;
	fgets(Buffer, BUFFERLENGTH, pFile);
	strLine = Buffer;
	vecStrTemp = Split(strLine, '\t');
	while (!feof(pFile))
	{
		bIfDelete = false;
		strSequenceTemp = vecStrTemp.at(iSequencecolumnNum);
		strProteinNameTemp = vecStrTemp.at(iProteinsNamecolumnNum);
		strReverseTemp = vecStrTemp.at(iReverseColumnNum);
		strContaminantTemp = vecStrTemp.at(iContaminantColumnNum);
		dPeptideIntensity = atof(vecStrTemp[iIntensityColumnNum].c_str());

		if (dPeptideIntensity == 0.0)
			bIfDelete = true;
		if (strReverseTemp != "" || strContaminantTemp != "")
		{
			bIfDelete = true;
			fgets(Buffer, BUFFERLENGTH, pFile);
			pstr = Buffer;
			strLine = Buffer;
			vecStrTemp = Split(strLine, '\t');
			continue;
		}

		//delete the contaminitant proteins names
		vecProteinNameTemps = Split(strProteinNameTemp, ';');
		if (param.m_bIfExistContamProtein)
		{
			vecProteinNameIter = vecProteinNameTemps.begin();
			for (; vecProteinNameIter != vecProteinNameTemps.end();)
			{
				strProteinNameSuffix = (*vecProteinNameIter).substr(0, param.m_strContaminantProteinPrefix.length());

				if (strProteinNameSuffix == param.m_strContaminantProteinPrefix||(*vecProteinNameIter).find("UPS")==(*vecProteinNameIter).npos)
				{
					vecProteinNameIter = vecProteinNameTemps.erase(vecProteinNameIter);
				}
				else
				{
					vecProteinNameIter++;
				}
			}
			if (vecProteinNameTemps.size() == 0)
			{
				bIfDelete = true;
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				strLine = Buffer;
				vecStrTemp = Split(strLine, '\t');
				continue;
			}
		}

		if (!bIfDelete)
		{
			peptideTemp.Clear();
			peptideTemp.m_strPeptideSeq = strSequenceTemp;
			vecProteinNameIter = vecProteinNameTemps.begin();
			if (vecProteinNameIter == vecProteinNameTemps.end())
			{
				cout << "Warning:\tPeptide " << strSequenceTemp << \
					" do not have corresponding proteins.\n";
				flog.mf_Input("Warning:\tPeptide " + strSequenceTemp + \
					" do not have corresponding proteins.\n");
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				strLine = Buffer;
				vecStrTemp = Split(strLine, '\t');
				continue;
			}

			peptideTemp.m_strPep2ProteinName = *vecProteinNameIter;
			for (; vecProteinNameIter != vecProteinNameTemps.end(); vecProteinNameIter++)
			{
				peptideTemp.m_vstrProteinsOfPeptide.push_back(*vecProteinNameIter);
			}
			if (vecProteinNameTemps.size() > 1)
			{
				peptideTemp.m_bProteinsShared = true;
			}
			else
			{
				peptideTemp.m_bProteinsShared = false;
			}
			peptideTemp.m_dPeptideIntesity = dPeptideIntensity;
			cpeptides.push_back(peptideTemp);

		}
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
		strLine = Buffer;
		vecStrTemp = Split(strLine, '\t');

	}

	fclose(pFile);
	cout << "\tLoaded " << cpeptides.size() << " peptides.\n";
	flog.mf_Input("\tLoaded " + fInt2String(cpeptides.size()) + " peptides.\n");

}
void CProteinIO::mf_LoadIdentifiedProteins(const CParam &param, const vector<CPeptide>&peptides, vector<CProtein> & S, \
	map<string, vector<int>> &mapIdentProteinsAndPeptideIndex)
{

	char Buffer[BUFFERLENGTH];
	FILE * pFile;
	errno_t err;
	// open file
	err = fopen_s(&pFile, param.m_strProteinSequenceFilePath.c_str(), "r");
	if (err != 0)
	{
		cout << "Error:\tCannot open " << param.m_strProteinSequenceFilePath << endl;
		flog.mf_Input("Error:\tCannot open " + param.m_strProteinSequenceFilePath + "\n");
		flog.mf_Destroy();
		exit(1);
	}

	ofstream ofileOfMissingProteins(param.m_strPIPLFAQResultPath + "\\MissingProteins.txt");
	char *pstr;
	string strIDTemp;
	string strSequenceTemp;
	string strFastaHeaderTemp;
	string strTemp;

	// first line
	fgets(Buffer, BUFFERLENGTH, pFile);
	while ((Buffer[0] == '\n') && (!feof(pFile)))  //allowing empty lines
	{
		fgets(Buffer, BUFFERLENGTH, pFile);
	}
	map<string, vector<int>>::iterator mapIdentProteinsAndPepIndexIter;
	std::match_results<std::string::const_iterator> Matchresult;
	bool valid;

	S.clear();
	pstr = Buffer;
	if (Buffer[0] == '>')  //get the protein id according to the regular expression fastatype;
	{
		strFastaHeaderTemp = Buffer;
		valid = std::regex_search(strFastaHeaderTemp, Matchresult, param.m_fastaType);
		if (valid == true)
		{
			strIDTemp = Matchresult[1];
		}
		else
		{
			cout << "Error:\tCannot parse the fasta file by the regular expression.\n";
			flog.mf_Input("Error:\tCannot parse the fasta file by the regular expression.\n");
			flog.mf_Destroy();
			exit(1);
		}
	}
	// read line by line
	int iPepIndex;
	while (fgets(Buffer, BUFFERLENGTH, pFile))
	{
		if (Buffer[0] == '\0')   //allow the blank lines
			continue;
		pstr = Buffer;

		if (Buffer[0] == '>')  //get the protein id according to the regular expression fastatype;
		{
			mapIdentProteinsAndPepIndexIter = mapIdentProteinsAndPeptideIndex.find(strIDTemp);
			if (mapIdentProteinsAndPepIndexIter != mapIdentProteinsAndPeptideIndex.end())
			{
				m_cprotein.Clear();
				m_cprotein.m_strFastaHeader = strFastaHeaderTemp;
				m_cprotein.m_strProteinName = strIDTemp;
				m_cprotein.m_strSequence = strSequenceTemp;
				for (int i = 0; i < mapIdentProteinsAndPepIndexIter->second.size(); i++)
				{
					iPepIndex = mapIdentProteinsAndPepIndexIter->second[i];
					if (!peptides[iPepIndex].m_bProteinsShared)
					{ // unique peptide
						m_cprotein.m_bIfContainUniquePep = true;
						m_cprotein.m_vUniquePeptideSequneces.push_back(peptides[iPepIndex].m_strPeptideSeq);
						m_cprotein.m_vUniquePeptideIntensities.push_back(peptides[iPepIndex].m_dPeptideIntesity);
					}
					else
					{
						m_cprotein.m_vSharedPeptideSequneces.push_back(peptides[iPepIndex].m_strPeptideSeq);
						m_cprotein.m_vSharedPeptideIntensities.push_back(peptides[iPepIndex].m_dPeptideIntesity);
					}
				}
				S.push_back(m_cprotein);
			}
			else
			{
				ofileOfMissingProteins << strIDTemp << endl;

			}

			strSequenceTemp.clear();
			strFastaHeaderTemp = Buffer;
			valid = std::regex_search(strFastaHeaderTemp, Matchresult, param.m_fastaType);
			if (valid == true)
			{
				strIDTemp = Matchresult[1];
			}
			else
			{
				cout << "Error:\tCannot parse the fasta file by the regular expression.\n";
				flog.mf_Input("Error:\tCannot parse the fasta file by the regular expression.\n");
				flog.mf_Destroy();
				exit(1);
			}
		}
		else // get the sequence of the protein id
		{
			strTemp.clear();
			strTemp = pstr;
			strSequenceTemp = strSequenceTemp + strTemp.substr(0, strTemp.size() - 1);
		}
	}
	mapIdentProteinsAndPepIndexIter = mapIdentProteinsAndPeptideIndex.find(strIDTemp);
	if (mapIdentProteinsAndPepIndexIter != mapIdentProteinsAndPeptideIndex.end())
	{
		m_cprotein.Clear();
		m_cprotein.m_strFastaHeader = strFastaHeaderTemp;
		m_cprotein.m_strProteinName = strIDTemp;
		m_cprotein.m_strSequence = strSequenceTemp;
		for (int i = 0; i < mapIdentProteinsAndPepIndexIter->second.size(); i++)
		{
			iPepIndex = mapIdentProteinsAndPepIndexIter->second[i];
			if (!peptides[iPepIndex].m_bProteinsShared)
			{ // unique peptide
				m_cprotein.m_bIfContainUniquePep = true;
				m_cprotein.m_vUniquePeptideSequneces.push_back(peptides[iPepIndex].m_strPeptideSeq);
				m_cprotein.m_vUniquePeptideIntensities.push_back(peptides[iPepIndex].m_dPeptideIntesity);
			}
			else
			{
				m_cprotein.m_vSharedPeptideSequneces.push_back(peptides[iPepIndex].m_strPeptideSeq);
				m_cprotein.m_vSharedPeptideIntensities.push_back(peptides[iPepIndex].m_dPeptideIntesity);
			}
		}
		S.push_back(m_cprotein);
	}
	else
	{
		ofileOfMissingProteins << strIDTemp << endl;
	}

	fclose(pFile);
	ofileOfMissingProteins.close();
}
bool CProteinIO::mf_SaveProteinSequences(string Path, vector<CProtein> cproteins)
{

	ofstream ofFile;
	ofFile.open(Path.c_str(), ios::out);
	if (!ofFile)
	{
		cout << "Cannot open \"" << Path << "\"." << endl;
		flog.mf_Input("Cannot open " + Path + "\n");
		return false;
	}
	int iProSize = cproteins.size();
	int iProLen = 0;
	for (int i = 0; i<iProSize; i++)
	{
		ofFile << cproteins[i].m_strFastaHeader;
		iProLen = cproteins[i].m_strSequence.size();
		for (int j = 0; j <iProLen; j = j + 60)
		{
			if (j + 60 < iProLen)
			{ // sequence is shorter than 60
				ofFile << cproteins[i].m_strSequence.substr(j, 60) << endl;
			}
			else
			{
				ofFile << cproteins[i].m_strSequence.substr(j, cproteins[i].m_strSequence.size() - j) << endl;
			}
		}

	}

	cout << "Saved " << cproteins.size() << " protein groups into " << Path << ".\n";
	flog.mf_Input("Load " + fInt2String(cproteins.size()) + " protein groups into " + Path + ".\n");

	ofFile.close();
	return true;

}
bool CProteinIO::mf_SaveProteinQuantInfo(string Path, vector<CProtein> cproteins)
{
	ofstream ofFile;
	ofFile.open(Path.c_str(), ios::out);
	if (!ofFile)
	{
		cout << "Error when writing " <<Path << endl;
		return false;
	}
	else
	{
		cout << "Save quantitative information of proteins into " << Path << "\n";
		flog.mf_Input("Save quantitative information of proteins into " + Path + "\n");
	}
	ofFile << "ProteinName\tShared peptides number\tUnique peptides number\tShared peptides sequences\tUnique peptides sequenes\
\tShared peptides intensities\tUnique peptides intensities\tPepNumInMaxSet\tAgent peptide sequence\tIntensity of agent peptide\n";
	map<string, int>::iterator mapPeptideSCIter;
	map<string, bool>::iterator mapbIfPeptidesSharedIter;
	for (size_t i = 0; i<cproteins.size(); i++)
	{
		ofFile << cproteins[i].m_strProteinName << "\t";
		ofFile << cproteins[i].m_vSharedPeptideSequneces.size() << "\t";
		ofFile << cproteins[i].m_vUniquePeptideSequneces.size()<<"\t";
		for (int j = 0; j < cproteins[i].m_vSharedPeptideSequneces.size(); j++)
		{
			ofFile << cproteins[i].m_vSharedPeptideSequneces[j] << ";";
		}
		ofFile << "\t";
		for (int j = 0; j < cproteins[i].m_vUniquePeptideSequneces.size(); j++)
		{
			ofFile << cproteins[i].m_vUniquePeptideSequneces[j] << ";";
		}
		ofFile << "\t";
		for (int j = 0; j < cproteins[i].m_vSharedPeptideIntensities.size(); j++)
		{
			ofFile << cproteins[i].m_vSharedPeptideIntensities[j] << ";";
		}
		ofFile << "\t";
		for (int j = 0; j < cproteins[i].m_vUniquePeptideIntensities.size(); j++)
		{
			ofFile << cproteins[i].m_vUniquePeptideIntensities[j] << ";";
		}
		ofFile << "\t";
		ofFile << cproteins[i].m_iPepNumInMaxSet << "\t";
		ofFile << cproteins[i].m_strAgentUniquePepSequence << "\t";
		ofFile << cproteins[i].m_dAgentUniquePepIntensity << "\t";
		ofFile << endl;
	}

	ofFile.close();
	return true;

}

void CDataIO::mf_LoadPeptidesFromMaxQuant(const CParam &param, vector<CPeptide>& cPeptides)
{
	CMaxQuantIO maxquantio;
	maxquantio.mf_GetAttributesName(param.m_strExprimentDesignPath);
	string strPeptideFilePath = param.m_strIdentResultDirectoryPath + "\\Peptides.txt";
	//maxquantio.mf_LoadPeptides(param, strPeptideFilePath, cPeptides, true);
	maxquantio.mf_LoadAllPeptides(param, strPeptideFilePath, cPeptides);


}
//∂¡»Îµ∞∞◊÷ –Ú¡–
void CDataIO::mf_LoadProteins(const CParam &param, vector<CProtein> &cproteins, \
	const vector<CPeptide>& cpeptides)
{
	CProteinIO cproteinio;
	map<string, vector<int>>mapIdentifiedProteinsAndPeptideIndex;
	map<string, vector<int>>::iterator mapIdentifiedProteinsAndPeptideIndexIter;
	map<string, vector<string>>::const_iterator mapProteinGroupsIter;

	//Get protein names which may exist iin the sample
	mapIdentifiedProteinsAndPeptideIndex.clear();
	int iGUSize = cpeptides.size();
	vector<int> vPeptideIndexTemp;
	for (int i = 0; i <iGUSize; i++)
	{
		for (int j = 0; j < cpeptides[i].m_vstrProteinsOfPeptide.size(); j++)
		{
			mapIdentifiedProteinsAndPeptideIndexIter = mapIdentifiedProteinsAndPeptideIndex.find(cpeptides[i].m_vstrProteinsOfPeptide[j]);
			if (mapIdentifiedProteinsAndPeptideIndexIter == mapIdentifiedProteinsAndPeptideIndex.end())
			{
				vPeptideIndexTemp.clear();
				vPeptideIndexTemp.push_back(i);
				mapIdentifiedProteinsAndPeptideIndex.insert(pair<string, vector<int>>(\
					cpeptides[i].m_vstrProteinsOfPeptide[j], vPeptideIndexTemp));
			}
			else
			{
				mapIdentifiedProteinsAndPeptideIndexIter->second.push_back(i);
			}

		}

	}

	cout << "Load " << mapIdentifiedProteinsAndPeptideIndex.size() << " proteins.\n";
	flog.mf_Input("Load " + fInt2String(mapIdentifiedProteinsAndPeptideIndex.size()) + " proteins.\n");

	cproteinio.mf_LoadIdentifiedProteins(param, cpeptides, cproteins, mapIdentifiedProteinsAndPeptideIndex);

#ifdef TEST
	string	strIdentifiedProteinsPath = param.m_strPIPLFAQResultPath + "\\IdentifiedProteins.fasta";
	cproteinio.mf_SaveProteinSequences(strIdentifiedProteinsPath, cproteins);
#endif
	cout << "There are " << cproteins.size() << " proteins which have sequence in the fasta file." << endl;
	flog.mf_Input("There are " + fInt2String(cproteins.size()) + " proteins which have sequence in the fasta file.\n");

}


bool CDataIO::mf_SaveProteinQuantInfo(string Path, vector<CProtein> cproteins)
{
	CProteinIO proteinio;
	return proteinio.mf_SaveProteinQuantInfo(Path, cproteins);

}

bool CCleavSiteIO::SaveNmcSite(const string StrFileName, vector<CCleavSite> cUnmcSites)
{
	ofstream ofFile;
	ofFile.open(StrFileName.c_str(), ios::out);

	if (!ofFile) {
		cout << "Error when writing \"" << StrFileName << "\"." << endl;
		return false;
	}
	ofFile << "m_cNineMer\tm_strProteinName\tm_iLocation\tMissedCutNum" << endl;
	for (size_t i = 0; i < cUnmcSites.size(); ++i)
	{
		for (int t = 0; t<2 * NEAR_NUM + 1; t++)
			ofFile << cUnmcSites[i].m_cNineMer.nine[t];
		ofFile << "\t" << cUnmcSites[i].m_strProteinName << "\t" << cUnmcSites[i].m_iLocation << "\t" << cUnmcSites[i].m_iMCSiteNum << endl;

	}
	ofFile.close();
	return true;
}

bool CCleavSiteIO::SaveMcSite(const string StrFileName, vector<CCleavSite> cMcSites)
{
	ofstream ofFile;
	ofFile.open(StrFileName.c_str(), ios::out);

	if (!ofFile)
	{
		cout << "Error when writing \"" << StrFileName << "\"." << endl;
		flog.mf_Input("Error when writing \"" + StrFileName + "\".\n");
		return false;
	}
	else
	{
		cout << "Save cleavage site into " << StrFileName << endl;
		flog.mf_Input("Save cleavage site into " + StrFileName + "\n");
	}
	//ofFile<<"m_cNineMer\tm_strProteinName\tm_iLocation\tMissedCutNum"<<endl;
	ofFile << "m_cNineMer\tm_strProteinName\tm_iMCSiteNum\tm_iLPeptideNum\tm_iRPeptideNum\tIfCleavageSite" << endl;
	//ofFile<<"m_cNineMer\tm_strProteinName\tm_iLPeptideNum\tm_iRPeptideNum\tMissedCutNum"<<endl;

	for (size_t i = 0; i < cMcSites.size(); ++i)
	{
		/*for(int t=0;t<2*NEAR_NUM+1;t++)
		ofFile<<cMcSites[i].m_cNineMer.nine[t];*/
		ofFile << cMcSites[i].m_cNineMer.nine;
		ofFile << "\t" << cMcSites[i].m_strProteinName;
		ofFile << "\t" << cMcSites[i].m_iMCSiteNum;
		ofFile << "\t" << cMcSites[i].m_iLPeptideNum;
		ofFile << "\t" << cMcSites[i].m_iRPeptideNum;


		//ofFile<<"\t"<<cMcSites[i].m_strProteinName<<"\t"<<cMcSites[i].m_iLocation<<"\t"<<cMcSites[i].m_bIfMCut<<endl;
		//ofFile<<"\t"<<cMcSites[i].m_bIfMCut<<endl;
		if (cMcSites[i].m_CleaveType == Cut)
		{
			ofFile << "\t" << 1 << endl;
		}
		else if (cMcSites[i].m_CleaveType == Missed)
		{
			ofFile << "\t" << 0 << endl;
		}


	}
	ofFile.close();
	return true;
}
bool CCleavSiteIO::SaveRMcSite(const string StrFileName, vector<CCleavSite> cMcSites)
{
	ofstream ofFile;
	ofFile.open(StrFileName.c_str(), ios::out);

	if (!ofFile)
	{
		cout << "Error when writing \"" << StrFileName << "\"." << endl;
		flog.mf_Input("Error when writing \"" + StrFileName + "\".\n");
		return false;
	}
	else
	{
		cout << "Save cleavage site into " << StrFileName << endl;
		flog.mf_Input("Save cleavage site into " + StrFileName + "\n");
	}
	//ofFile<<"m_cNineMer\tm_strProteinName\tm_iLocation\tMissedCutNum"<<endl;
	ofFile << "m_cNineMer\tm_strProteinName\tm_iMCSiteNum\tm_iLPeptideNum\tm_iRPeptideNum\tIfCleavageSite" << endl;
	//ofFile<<"m_cNineMer\tm_strProteinName\tm_iLPeptideNum\tm_iRPeptideNum\tMissedCutNum"<<endl;

	for (size_t i = 0; i < cMcSites.size(); ++i)
	{
		/*for(int t=0;t<2*NEAR_NUM+1;t++)
		ofFile<<cMcSites[i].m_cNineMer.nine[t];*/
		if (cMcSites[i].m_cNineMer.nine[4] == 'R')
		{
			ofFile << cMcSites[i].m_cNineMer.nine;
			ofFile << "\t" << cMcSites[i].m_strProteinName;
			ofFile << "\t" << cMcSites[i].m_iMCSiteNum;
			ofFile << "\t" << cMcSites[i].m_iLPeptideNum;
			ofFile << "\t" << cMcSites[i].m_iRPeptideNum;


			//ofFile<<"\t"<<cMcSites[i].m_strProteinName<<"\t"<<cMcSites[i].m_iLocation<<"\t"<<cMcSites[i].m_bIfMCut<<endl;
			//ofFile<<"\t"<<cMcSites[i].m_bIfMCut<<endl;
			if (cMcSites[i].m_CleaveType == Cut)
			{
				ofFile << "\t" << 1 << endl;
			}
			else if (cMcSites[i].m_CleaveType == Missed)
			{
				ofFile << "\t" << 0 << endl;
			}
		}



	}
	ofFile.close();
	return true;
}
bool CCleavSiteIO::SaveKMcSite(const string StrFileName, vector<CCleavSite> cMcSites)
{
	ofstream ofFile;
	ofFile.open(StrFileName.c_str(), ios::out);

	if (!ofFile)
	{
		cout << "Error when writing \"" << StrFileName << "\"." << endl;
		flog.mf_Input("Error when writing \"" + StrFileName + "\".\n");
		return false;
	}
	else
	{
		cout << "Save cleavage site into " << StrFileName << endl;
		flog.mf_Input("Save cleavage site into " + StrFileName + "\n");
	}
	//ofFile<<"m_cNineMer\tm_strProteinName\tm_iLocation\tMissedCutNum"<<endl;
	ofFile << "m_cNineMer\tm_strProteinName\tm_iMCSiteNum\tm_iLPeptideNum\tm_iRPeptideNum\tIfCleavageSite" << endl;
	//ofFile<<"m_cNineMer\tm_strProteinName\tm_iLPeptideNum\tm_iRPeptideNum\tMissedCutNum"<<endl;

	for (size_t i = 0; i < cMcSites.size(); ++i)
	{
		/*for(int t=0;t<2*NEAR_NUM+1;t++)
		ofFile<<cMcSites[i].m_cNineMer.nine[t];*/
		if (cMcSites[i].m_cNineMer.nine[4] == 'K')
		{
			ofFile << cMcSites[i].m_cNineMer.nine;
			ofFile << "\t" << cMcSites[i].m_strProteinName;
			ofFile << "\t" << cMcSites[i].m_iMCSiteNum;
			ofFile << "\t" << cMcSites[i].m_iLPeptideNum;
			ofFile << "\t" << cMcSites[i].m_iRPeptideNum;


			//ofFile<<"\t"<<cMcSites[i].m_strProteinName<<"\t"<<cMcSites[i].m_iLocation<<"\t"<<cMcSites[i].m_bIfMCut<<endl;
			//ofFile<<"\t"<<cMcSites[i].m_bIfMCut<<endl;
			if (cMcSites[i].m_CleaveType == Cut)
			{
				ofFile << "\t" << 1 << endl;
			}
			else if (cMcSites[i].m_CleaveType == Missed)
			{
				ofFile << "\t" << 0 << endl;
			}
		}



	}
	ofFile.close();
	return true;
}

CCleavSiteIO::CCleavSiteIO()
{}
CCleavSiteIO::~CCleavSiteIO()
{}