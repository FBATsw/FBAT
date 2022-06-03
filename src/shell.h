#pragma once
#include <iostream>
#include <cstring>
#include "mystream.h"

const int kMaxCMDLen = 64;
const int kMaxLineLen = 102400;
const int kMaxPromptLen = 8;
const int kMaxCharFileName = 64;

const char escape_char = '#';

class CMD_SHELL;
typedef void (CMD_SHELL::*funcp)(char*);

struct COMMAND {
	char	name[kMaxCMDLen];
	char	*help;
	funcp dofunc;
};
	
class COMMANDCHAIN : public COMMAND {
	public:
	COMMANDCHAIN	*next;
	COMMANDCHAIN(char *s, funcp call, char *usage=NULL) 
		{ strncpy(name, s, kMaxCMDLen-1); dofunc=call; help=usage; next=NULL; }
	void usage(mystream &os) 
		{ if (help!=NULL) MYOUTPUT(os, "  " << help << endl) 
			else MYOUTPUT(os, "help for " << name << " is unavailable" << endl) }
};


class CMD_SHELL {	
	
	char line[kMaxLineLen];
	char prompt[kMaxPromptLen];
	char log_file_name[kMaxCharFileName];	
	bool log_on;
	bool shell_done;
	ofstream *f_log;
	COMMANDCHAIN *cmd;
	
	public:
	mystream os;
	CMD_SHELL() { cmd=NULL; strcpy(prompt, ">> "); f_log=NULL; log_on=false; shell_done=false; }
	COMMANDCHAIN* parsecmd(char *s, bool do_command=true);
	void addcmd(COMMAND cmds[], int n);
	void addcmd(char *s, funcp dofunc, char *usage=NULL);
	void init();
	void loop();
	void runcmd (char *c);
	void runscript(char *s);
	void bye(char *s);
	void setprompt(char *s) { strncpy(prompt, s, kMaxPromptLen-1); }
	void help(char *s);
	void usage(char *s);
	void log(char *s);
	virtual void call(COMMANDCHAIN *cmd, char *param);
	char* ask(char *question);
};

