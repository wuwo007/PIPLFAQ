#include<fstream>
#include<string>

class CLog
{
public:

	bool mf_Init(std::string strResultPath, std::ios_base::openmode = std::ios::out);
	void mf_Destroy();
	void mf_Input(std::string str);
private:
	std::ofstream fstr;
};

extern CLog flog;
