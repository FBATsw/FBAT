// myerror.h
#pragma once
//#include <cstring>
#include <string.h>

class ERROR {
	char msg[256];
	
	public:
	ERROR(char *s) { strncpy(msg, s, 255); }
	void append(char *s) { strncat(msg, s, 255-strlen(msg)); }
	void rethrow() { throw ERROR(msg); }
	char* what() { return msg; }
};

class Data_Error : public ERROR { public: Data_Error(char *s) : ERROR(s) { }; };

class File_Error : public ERROR { public: File_Error(char *s) : ERROR(s) { }; };

class Mem_Error : public ERROR { public: Mem_Error(char *s) : ERROR(s) { }; };

class Param_Error:  public ERROR { public: Param_Error(char *s) : ERROR(s) { }; };