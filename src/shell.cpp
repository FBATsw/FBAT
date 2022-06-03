#include <string.h>
#include "myutil.h"
#include "shell.h"
#include "mystream.h"
#include "myerror.h"

void CMD_SHELL::call(COMMANDCHAIN *cmd, char *param) {
	
	try {
		(this->*(cmd->dofunc))(param);
	}
	catch (Param_Error &err) {
		MYOUTPUT(os, err.what() << endl)
		cmd->usage(os);
	}
	catch (ERROR &err) {
		MYOUTPUT(os, err.what() << endl)
	}
	catch (...) {
		MYOUTPUT(os, "sorry, an unspecified error occured\n")
	}
}

// insert commandchain and sort by command name
void CMD_SHELL::addcmd(char *s, funcp fp, char *usage) {
	
	if (cmd==NULL) {
		cmd = new COMMANDCHAIN(s, fp, usage);
		return;
	}
	
	COMMANDCHAIN *prev, *next;
	
	prev = NULL;
	for (next=cmd; next!=NULL && strcmp(s,next->name)>0; next=next->next)
		prev = next;
		
	if (prev) {
		prev->next = new COMMANDCHAIN(s, fp, usage);
		prev->next->next = next;
	} else {
		cmd = new COMMANDCHAIN(s, fp, usage);
		cmd->next = next;
	}
	
	
}

void CMD_SHELL::addcmd(COMMAND cmds[], int n) {

	try {
		for (int i=0; i<n; i++)
			addcmd(cmds[i].name, cmds[i].dofunc, cmds[i].help);
	}
	catch (...) {
		MYOUTPUT(os, "program initialization failed. program will exit." << endl)
		bye("");
	}

}
	
COMMANDCHAIN* CMD_SHELL::parsecmd(char *s, bool do_command) {
	
	COMMANDCHAIN *command, *target;
	int	i, maxmatch;
	bool	ambig;
	char	*c1, *c2;
	
	target = NULL;
	maxmatch = 0;
	ambig = true;
	while (*s==' ') s++;
	// 2/18/2004 add escape character for notes
	if (*s == escape_char)
		return NULL;
		
	for (command=cmd; command!=NULL; command=command->next) {
		for (i=0,c1=command->name,c2=s; *c1 && *c2 && *c1==*c2; c1++,c2++)
			i++;
		if (i && *c1==0 && (*c2==' ' || *c2==0 || *c2=='\t')) {	// exact match
			maxmatch = i;
			ambig = false;
			target = command;
			break;
		}
		if (*c2 && *c2!=' ' && *c2!='\t' && *c2!='\n' && *c2!='\r')	// not match
			continue;
		if (i>maxmatch) {
			maxmatch = i;
			target = command;
			ambig = false;
		} 
		else if (i==maxmatch)
			ambig = true;
	}
	
	if (maxmatch>0) {
		if (!ambig) {
			if (target!=NULL && target->dofunc!=NULL) {
				if (do_command)
					call(target, s);
				return target;
			}
		} else
			MYOUTPUT(os, "command too ambigous" << endl)
	} 
	else if (WordCount(s)>0)
		MYOUTPUT(os, "command unknown" << endl)
	
	return NULL;
}

void CMD_SHELL::init() {
	
	static COMMAND cmds[] = {
		{	"quit",
			"quit -- quit program",
			static_cast<funcp> (&CMD_SHELL::bye)
		},
		{	"?",
			"? [command]  -- print out usage for command",
			static_cast<funcp> (&CMD_SHELL::help)
		},
		{	"log",
			"log [log_file,on,off] -- logging to log_file, toggle logging",
			static_cast<funcp> (&CMD_SHELL::log)
		},
		{	"run",
			"run script_file -- run batch commands from script_file",
			static_cast<funcp> (&CMD_SHELL::runscript)
		}
	};
	
	addcmd(cmds, 4);

	os.setf(ios::fixed, ios::floatfield);
	os.setf(ios::left, ios::adjustfield);

}
	
void CMD_SHELL::loop() {
	
	while (!shell_done) {
		MYOUTPUT(os, prompt)
		*line = 0;
		cin.getline(line,kMaxLineLen);
		if (*line) {
			if (log_on && f_log)
				*f_log << line << endl;
			if (strlen(line) >= kMaxLineLen-1)
				MYOUTPUT(os, "command line too long\n")
			else
				parsecmd(line);
		}
	}
}

void CMD_SHELL::help(char *s) {
	
	if (WordCount(s)==1) {
		MYOUTPUT(os, "commands available:\n")
		for (COMMANDCHAIN *command=cmd; command!=NULL; command=command->next)
			if (command->help!=NULL)
				command->usage(os);
	}
	else {
		COMMANDCHAIN *command=parsecmd(GetWord(s,2,NULL), false);
		if (command!=NULL)
			command->usage(os);
	}

}

void CMD_SHELL::usage(char *s) {
	
	COMMANDCHAIN *command = parsecmd(s, false);
	if (command!=NULL)
		command->usage(os);

}

void CMD_SHELL::bye(char *s) {
	
	MYOUTPUT(os,  "\n\n                GOOD BYE......\n\n\n" << endl)
	
	shell_done = true;
	
}

void CMD_SHELL::log(char *s) {
	
	int wc = WordCount(s);
	
	if (wc>2) {
		usage(s);
		return;
	}
	if (wc==2) {
		char word[64];
		GetWord(s, 2, word);
		if (strcmp(word,"on")==0)
			log_on = true;
		else if (strcmp(word,"off")==0)
			log_on = false;
		else {	// specified a log file
			if (f_log!=NULL) {
				f_log->close();
				delete f_log;
				f_log = NULL;
			}
			strncpy(log_file_name, word, kMaxCharFileName);
			if ((f_log=new ofstream(log_file_name, ios::out|ios::trunc)) == NULL) {
				MYOUTPUT(os, "failed to create log file \"" << log_file_name << "\"\n")
				log_on = false;
			} else {
				log_on = true;
				os.setfs(f_log);
				os.setf(ios::fixed, ios::floatfield);
				os.setf(ios::left, ios::adjustfield);
			}
		}
	}

	if (log_on) {
		MYOUTPUT(os, "logging to file \"" << log_file_name << "\" is on\n");
		os.setfs(f_log);
	} else {
		os.setfs(NULL);
		if (f_log!=NULL)
			MYOUTPUT(os, "logging to file \"" << log_file_name << "\" is off\n")
		else
			MYOUTPUT(os, "logging if off, no log file has been set\n")
	}

}

char* CMD_SHELL::ask(char *question) {
		
		MYOUTPUT(os, question << "  ")
		*line = 0;
		cin.getline(line,kMaxLineLen);
		
		return line;

}

// run script
void CMD_SHELL::runscript(char *s) {
	
	char *fname;
	if (WordCount(s)!=2)
		throw Param_Error("wrong number of arguments");
	fname = GetWord(s, 2, NULL);
	fstream fs(fname, ios::in);
	if (fs.fail()) {
		MYOUTPUT(os, "opening file \"" << fname << "\" failed\n")
		return;
	}
	while (fs.getline(line,kMaxLineLen) && !fs.eof()) 
		runcmd(line);
}

void CMD_SHELL::runcmd(char *c) {
	
	if (c==NULL || *c==0)
		return;
	strncpy(line, c, kMaxLineLen);
	line[kMaxLineLen-1] = 0;	
	MYOUTPUT(os, prompt)
	if (*line) {
		if (log_on && f_log)
			*f_log << line << endl;
		parsecmd(line);
	}
}

