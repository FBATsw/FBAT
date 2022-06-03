#pragma once
#include <string.h>
#include "mysimplepointerlist.h"
#include "mystream.h"
//#include "sufficient_stat.h"
#include "xfunc.h"

class Sufficient_Stat;

class PERSON;
class NUCFAMILY;

const int kMaxLocusNameLen = 255;
const int kMaxAlleles = 80;
const int kSexLinked = 1;

const int kMaxIDLen = 64;
const int kMaxPedNameLen = 64;
const int kMaxMarriages = 4;

const int flag_uninformative = 1<<15;
const int flag_processed = 1<<14;
const int flag_selected = 1<<10;
const int flag_phased = 1<<9;

const double kVeryLargeNumber = 1.0e38;
const double kEvenLargerNumber = 2.0e38;

enum PHASE { ph_pat=1, ph_mat=2, ph_unknown=4, ph_derived=8 };

enum Trait_type { tr_outcome=1, tr_covariate=2};

class Genotype {
	public:
	char a[2];
	
	Genotype() { a[0] = 0; a[1] = 0; }

	void set(char a1, char a2) { a[0] = a1; a[1] = a2;}
	void set(Genotype &gt) { a[0] = gt.a[0]; a[1] = gt.a[1]; }

	bool defined(bool xlink_male=false) { return a[0]>0 && (xlink_male?(a[1]==0 || a[1]==a[0]):a[1]>0); }	// partial (one allele) genotype is not allowed
	bool equal(char a1, char a2) { return (a[0]==a1 && a[1]==a2) || (a[0]==a2 && a[1]==a1); }
	bool equal(Genotype &gt) { return equal(gt.a[0], gt.a[1]); }
	bool homozygote() { return a[0]>0 && a[0]==a[1]; }

	void swap() { char c = a[0]; a[0] = a[1], a[1] = c; }

	int index() { return gtindex(a[0],a[1]); }
	int allele_cnt(char a1, bool xlink_male=false) { return (a[0]==a1?1:0)+(a[1]==a1 && !xlink_male?1:0); }
};

struct Locus_Info {
	char name[kMaxLocusNameLen];
	int nAlleles;
	int chrNum;
	int pos;
	int flag;
	double cm;
	double *frq;
	double *frq_gt;
	int alleleName[kMaxAlleles];
	
	Locus_Info() { chrNum=0; pos=0, cm=0, nAlleles=0; alleleName[0]=0; frq=NULL; frq_gt=NULL; flag=0;}
	char allele_id(int aName);
	bool sexlinked() { return flag&kSexLinked; }
	void sexlinked(bool link) { if (link!=sexlinked()) flag ^= kSexLinked; }
	void setflag(int f) { flag |= f; }
	void resetflag (int f) { flag &= ~f; }
	int testflag(int f) { return flag & f; }
	~Locus_Info() { if (frq!=NULL) delete[] frq; if (frq_gt!=NULL) delete[] frq_gt; frq=NULL; frq_gt=NULL; }
};

struct Trait_Info {
	char name[kMaxLocusNameLen];
	int nlevels;	// 0: quantitative traits
	int *lvalues;	// values of level
	int flag;
	double mean[4];	// trait mean of genotype[0..2], mean[3]=sample mean
	double sd;	// trait standard deviation (presuming equal for genotype[0..2] and equal sample sd (h^2 neglegible)
	double penetrance[3];	// for nlevels=2
	
	Trait_Info() { name[0]=0;  flag=0; lvalues=NULL; }
	void set(char *s, int nl=0) { if (s) strncpy(name, s, kMaxLocusNameLen-2); nlevels=nl; }
	void setflag(int f) { flag |= f; }
	void resetflag (int f) { flag &= ~f; }
	int testflag(int f) { return flag & f; }
	bool match_name(char *s) { return (s && strcmp(s,name)==0); }
	double prob_y_g(double y, int g); // P(Y=y|G=g)
	~Trait_Info() { if (lvalues!=NULL) delete[] lvalues; };
};

class PERSON {
	private:
	int	flag;	
	
	public:
	char id[kMaxIDLen];	
	char sex;
	char affection;
	Genotype *gt;
	double	*pt;
	
	NUCFAMILY *ascend;
	NUCFAMILY *descend[kMaxMarriages];
	
	PERSON(char *s="");
	//PERSON(int nMarkers, int nTraits);
	void sex_status(char status);
	void affection_status(char status);
	void ascends(NUCFAMILY *fam);
	void descends(NUCFAMILY *fam);
	void notrait(int which) { if (pt!=NULL) pt[which] = kEvenLargerNumber; }	// set trait unknown
	bool hastrait(int which) { return which<0?(affection>0):(pt?(pt[which]<kVeryLargeNumber):false); }
	bool match_id(char *s) { return (s!=NULL && strcmp(id,s)==0); }
	bool good_gt(int loc_index, bool xlink=false) { return gt!=NULL && gt[loc_index].defined(xlink && sex==1); }
	double trait(int which) { return (pt!=NULL && which>=0? pt[which] : kVeryLargeNumber); }

	void setflag(int f) { flag |= f; }
	void resetflag (int f) { flag &= ~f; }
	int testflag(int f) { return flag & f; }
	
	bool testroot() { return ascend==NULL; }

	~PERSON() { if (pt!=NULL) delete[] pt; if (gt!=NULL) delete[] gt; }
};

typedef MySimplePointerList<class PERSON> PERSONLIST;

class NUCFAMILY {
	private:
	PERSONLIST *members;
	int flags;

	public:
	Sufficient_Stat *stat;
	
	NUCFAMILY(PERSON *f=NULL, PERSON *m=NULL);
	~NUCFAMILY();
	PERSON* fa() { return members->node; }
	PERSON* mo() { return members->next->node; }
	PERSONLIST* sibs() { return  members->next->next; }
	PERSONLIST* memb() { return members; }
	void fa(PERSON *p) { members->node = p; if (p) p->descends(this); }
	void mo(PERSON *p) { members->next->node = p; if (p) p->descends(this); }
	void sib(PERSON *p) { if (p) { new PERSONLIST(p,members); p->ascends(this); } }
	PERSON* id_lookup(char *id);

	int parentsGenotyped(int loc_index, bool xlink=false) { 
		return (fa() && fa()->gt && fa()->gt[loc_index].defined(xlink)?1:0)
				+ (mo() && mo()->gt && mo()->gt[loc_index].defined()?1:0); }
	int sibsGenotyped(int loc_index, bool xlink=false);	// return number of sibs genotypes at loc_index
	int sib_count(); // return number of sibs
	bool mendelerr(int loc_index, bool xlink=false);
	bool compatible_mating(Genotype &fg, Genotype &mg, int locidx, bool xlink);

	void reset_genotype(int loc_index);	// zero the genotype in all members

	void setflag(int f) { flags |= f; }
	void resetflag (int f) { flags &= ~f; }
	int testflag(int f) { return flags & f; }

};


typedef MySimplePointerList<class NUCFAMILY> NUCFAMILYLIST;

class PEDIGREE {
	public:
	char name[kMaxPedNameLen];
	NUCFAMILYLIST *nucfams;
	PERSONLIST *members;
	
	PEDIGREE(char *id="") { strncpy(name, id, kMaxPedNameLen-1); nucfams = NULL; members=NULL; }
	~PEDIGREE() { if (nucfams) delete nucfams; if (members) delete members; nucfams=NULL, members=NULL; }
	NUCFAMILY* fam_lookup(char *fid, char *mid);
	PERSON* id_lookup(char *id);
	bool match_name(char *s) { return (s && strcmp(name,s)==0); }
	void addperson(PERSON *per) { if (per) { PERSONLIST *li = new PERSONLIST(per, members); if (!members) members = li;} }
	void addfamily(NUCFAMILY *fam) { if (fam) { NUCFAMILYLIST *li = new NUCFAMILYLIST(fam, nucfams); if (!nucfams) nucfams = li;} }
	void removeTraits();
	void print_marker(int idx, mystream &os, Locus_Info *ti, bool xlink);
};

typedef MySimplePointerList<class PEDIGREE> PEDIGREELIST;



