#ifndef UT_LOG_H
#define UT_LOG_H

#include <iostream>
#include <fstream>

#define LOG(msg) Logger::getInstance().file << msg << "\n"; Logger::getInstance().file.flush(); \
	std::cout << msg << "\n"	

//#define LOG(msg) 

using namespace std;
class Logger
{
public:	
	static Logger & 
	getInstance()
	{
		if(_instance == NULL)
		{
			_instance = new Logger();

		}
		return *_instance;
	}

	fstream file;
private:
	Logger()
	{
		file.open ("log.txt",  fstream::out);		
	}

	~Logger(){} //Not implemented
	Logger(Logger const&){}; //Not implemented
	Logger& operator=(Logger const&){};  //Not implemented
	static Logger * _instance;
};


Logger* Logger::_instance = NULL;




#endif
