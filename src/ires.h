// ires.h
//   wai-ki yip (3/4/11)
//   intermediate fbat statistics 
//     contribution from each family from each pedigree
//
#pragma once
#include <string>
#include <vector>
#include "fbat.h"

/* option to weigh the contribution from each SNP */
#define $UNW	0
#define $WT	1

struct IRESITEM{
	int sel;
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	double s;
	double es;
	double v;
};

struct WTITEM{
	int sel;
	double wt;
};

//class IRES : public DATA{ 
class IRES{
private:
	int count;	// total number of items added
	int wc;		// total number of weight items added
	vector <struct IRESITEM> ir;
	vector <struct WTITEM> weight;
	int fam_total;
	int ped_count;
	int old_min_size;
	class FBAT *fb;
	bool empirical;	// use empirical variance -e option
	PEDIGREELIST *peds;	// global pedigree list
public:
	IRES(FBAT *fbat, bool empirical, PEDIGREELIST *peds, int pedcount);
	~IRES();

	int getsize() {return (int)ir.size();};
	void additem(int sel, PEDIGREELIST *ped, NUCFAMILYLIST *flist, double s, double es, double v);
	void addwt(int sel, double wt);
	void print();
	void setfam_total(int tf) { this->fam_total = tf;};
	int  getfam_total() { return this->fam_total;};
	void save_minsize(int ms) { this->old_min_size = ms;};
	int  restore_minsize() { return this->old_min_size;};
	void rvtest(int option);
};


