#include"stdafx.h"
#include"Parameters.h"

int CParam::mf_LoadParams(string paramfilepath)
{
	cout << "\tSet basic parameters\n";
	flog.mf_Input("Set basic parameters\n");

	map <string, string>  mapArttributesInFile;
	map<string, string>::iterator mapAttributesInfileIter;
	mapArttributesInFile = GetParametersFromFile(paramfilepath);

	mapAttributesInfileIter = mapArttributesInFile.find("InputType");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strInputType = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter InputType is not correct" << endl;
		flog.mf_Input("Parameter InputType is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}



	mapAttributesInfileIter = mapArttributesInFile.find("IdentResultDirectoryPath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strIdentResultDirectoryPath = mapAttributesInfileIter->second;

	}
	else
	{
		cout << "Parameter IdentResultDirectoryPath is not correct" << endl;
		flog.mf_Input("Parameter IdentResultDirectoryPath is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("ProteinSequenceFilePath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strProteinSequenceFilePath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter ProteinSequenceFilePath is not correct" << endl;
		flog.mf_Input("Parameter ProteinSequenceFilePath is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}
	mapAttributesInfileIter = mapArttributesInFile.find("IdentifierParsingRule");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_fastaType = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter IdentifierParsingRule is not correct" << endl;
		flog.mf_Input("Parameter IdentifierParsingRule is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("PIPLFAQResultPath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strPIPLFAQResultPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter PIPLFAQResultPath is not correct" << endl;
		flog.mf_Input("Parameter PIPLFAQResultPath is not correct\n");
		flog.mf_Destroy();
		exit(1);
	}

	if (m_strInputType == "MaxQuant")
	{
		m_strExprimentDesignPath = m_strIdentResultDirectoryPath + "\\experimentalDesignTemplate.txt";
	}

	
	mapAttributesInfileIter = mapArttributesInFile.find("IfExistContaminantProtein");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		if (!fStringToBool(mapAttributesInfileIter->second, m_bIfExistContamProtein))
		{
			cout << "Cannot convert " << mapAttributesInfileIter->second << " to bool!" << endl;
			flog.mf_Input("Error:\tCannot convert " + mapAttributesInfileIter->second + " to bool!\n");
			flog.mf_Destroy();
		}
	}
	else
	{
		cout << "Error:\tThe format of parameter IfExistContaminantProtein is not correct." << endl;
		flog.mf_Input("Error:\tThe format of parameter IfExistContaminantProtein is not correct.\n");
		flog.mf_Destroy();
		exit(1);
	}
	if (m_bIfExistContamProtein == true)
	{
		mapAttributesInfileIter = mapArttributesInFile.find("PrefixOfContaminantProtein");
		if (mapAttributesInfileIter != mapArttributesInFile.end())
		{
			m_strContaminantProteinPrefix = mapAttributesInfileIter->second;
		}
		else
		{
			cout << "Error:\tThe format of parameter PrefixOfContaminantProtein is not correct." << endl;
			flog.mf_Input("Error:\tThe format of parameter PrefixOfContaminantProtein is not correct.\n");
			flog.mf_Destroy();
			exit(1);
		}
	}

	return 0;

}

void CParam::Init(string ParamFilePath)
{
	if (!mf_GetResultPath(ParamFilePath))
	{
		cout << "Cannot read the result directory path from the parameter file.\n";
		exit(2);
	}
	flog.mf_Init(m_strPIPLFAQResultPath + "\\Log.txt", ios_base::app);

	cout << "Start of TestCHelpD:" << endl;
	flog.mf_Input("Start of TestCHelpD:\n");
	//S1 读入基本参数
	mf_LoadParams(ParamFilePath);
}

void CParam::mf_SetDataDependParams(const vector<CPeptide>&Peptides)
{
	m_iAllowMaxPeptideLength = 0;
	m_iAllowMinPeptideLength = 1000;
	for (int i = 0; i < Peptides.size(); i++)
	{
		if (m_iAllowMinPeptideLength > Peptides[i].m_strPeptideSeq.size())
		{
			m_iAllowMinPeptideLength = Peptides[i].m_strPeptideSeq.size();
		}
		if (m_iAllowMaxPeptideLength < Peptides[i].m_strPeptideSeq.size())
		{
			m_iAllowMaxPeptideLength = Peptides[i].m_strPeptideSeq.size();
		}
	}

	cout << "The minimum length of identified peptides is " << m_iAllowMinPeptideLength << "\n";
	flog.mf_Input("The minimum length of identified peptides is " + fInt2String(m_iAllowMinPeptideLength) + "\n");
	cout << "The maximum length of identified peptides is " << m_iAllowMaxPeptideLength << "\n";
	flog.mf_Input("The maximum length of identified peptides is " + fInt2String(m_iAllowMaxPeptideLength) + "\n");

}
void CParam::mf_GetBasicParmas(CParam param)
{
	this->m_strInputType = param.m_strInputType;
	this->m_strIdentResultDirectoryPath = param.m_strIdentResultDirectoryPath;
	this->m_strProteinSequenceFilePath = param.m_strProteinSequenceFilePath;
	this->m_strPIPLFAQResultPath = param.m_strPIPLFAQResultPath;
	this->m_iAllowMinPeptideLength = param.m_iAllowMinPeptideLength;
	this->m_iAllowMaxPeptideLength = param.m_iAllowMaxPeptideLength;
}

bool CParam::mf_GetResultPath(string ParamFilePath)
{
	map<string, string> mapAttributeNameAndValues;
	map<string, string>::iterator mapAttributeNameAndValuesIter;
	mapAttributeNameAndValues = GetParametersFromFile(ParamFilePath);
	mapAttributeNameAndValuesIter = mapAttributeNameAndValues.find("PIPLFAQResultPath");
	string strAtributeValue;
	if (mapAttributeNameAndValuesIter != mapAttributeNameAndValues.end())
	{
		m_strPIPLFAQResultPath = mapAttributeNameAndValuesIter->second;
		return true;
	}
	else
	{
		return false;
	}
}
