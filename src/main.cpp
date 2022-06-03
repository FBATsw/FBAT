#include "myutil.h"
#include "myerror.h"
#include "mystream.h"
#include "fbat.h"

//#include <profiler.h>

int main() {
	
//	int i = ProfilerInit(collectSummary, PPCTimeBase, 100, 10);
	
//	ProfilerSetStatus(0);
	
	FBAT xdt;
	
	try {
		xdt.CMD_SHELL::init();
		xdt.DATA::init();

		xdt.init();
	}
	catch (...) {
		cout << "insuficient memory, initialization failed\n";
		return -1;
	}
	
	xdt.loop();
	
//	i = ProfilerDump("\pfbat.prof");
	
//	ProfilerTerm();
		
	return 0;
	
}