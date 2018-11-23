#include"stdafx.h"

CLog flog;
bool CLog::mf_Init(std::string strLogPath, std::ios_base::openmode mode)
{
	//string strLogPath = strResultPath + "\\log.txt";
	fstr.open(strLogPath.c_str(), mode);
	if (!fstr)
		return false;
	else
		return true;
}
void CLog::mf_Input(std::string str)
{
	if (!fstr.is_open())
	{
		std::cout << "Please check the log file is opening." << std::endl;
		exit(-4);
	}
	else
	{
		fstr << str;
	}

}
void CLog::mf_Destroy()
{
	fstr.close();
}