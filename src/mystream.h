#pragma once

#define USE_OUTPUT_MACRO

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#ifdef USE_OUTPUT_MACRO
	#define MYOUTPUT(os,a) { if (os.ocs!=NULL) *(os.ocs) << a; if (os.ofs!=NULL) *(os.ofs) << a; }
	class mystream {
		public:
		ostream	*ocs;
		ofstream *ofs;
		mystream(ostream *s1=&cout, ofstream *s2=NULL) {	ocs=s1; ofs=s2; }
		long setf(ios_base::fmtflags bits, ios_base::fmtflags mask) { if (ofs!=NULL) ofs->setf(bits, mask); return ocs->setf(bits, mask); }
		void unsetf(ios_base::fmtflags bits) { ocs->unsetf(bits); if (ofs!=NULL) ofs->unsetf(bits); }
		void setfs(ofstream *s) { ofs=s; }
		void setcs(ostream *s) { ocs=s; }
	};
#else
	#define mystream ostream
	#define MYOUTPUT(os,s) { os << s; }
#endif

