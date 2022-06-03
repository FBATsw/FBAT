/*
 *  sufficient_stat.h
 *  fbat
 *
 *  Created by Xin Xu on 12/14/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once

#include <iostream>

#include "mysimplepointerlist.h"
#include "mystream.h"
#include "pedigree.h"
#include "newmat.h"

const int kotherallele = 99;

struct Compatible_parent_genotype {
	int ngt;
	Genotype gt[8][2];	// <=4 possible mating genotypes; 0=fa, 1=mo
	
	Compatible_parent_genotype() { ngt = 0; };
};

class Haplotype_List {
	public:
	char *h[2];
	int nmrk;
	double p;	// em freq
	double q;	// temorary holder for p in em calculation
	Haplotype_List *next;
	
	Haplotype_List(int n);
	Haplotype_List(int n, char *h1, char *h2, Haplotype_List *nxt=NULL);
	bool equal(char *h1, char *h2);	// test if haplotype=h1/h2
	bool identical(char *h1, char *h2);	// test if h1==h2
	bool compatible(char *h1, char *h2); // test if haplotype=h1/h2, but allow kotherallele
	Haplotype_List* find(char *h1, char *h2);
	~Haplotype_List();
};

class Haploid {
	public:
	int nmrk;
	char *h;
	double p;
	double q; // D'
	
	Haploid(int n, char *ht, double f=0.0);
	bool equal(char *ht);
	~Haploid() { if (h!=NULL) delete[] h; h = NULL; }
};

typedef MySimplePointerList<class Haploid> HaploidList;

class HaploidTable {
	public:
	int size;
	int nmrk;
	Locus_Info **loc_ptr;
	HaploidList *hpl;
	
	HaploidTable() { size=0; hpl=NULL; loc_ptr=NULL;}
	HaploidTable(Locus_Info *loc, int n, int idx[], Haplotype_List *hl);
	void maketable(Haplotype_List *hl);
	void dprime();
	int allele(char *h);	// return allele index starting at 1 (sequential in *hpl). 0=unknown/missing
	int gid(Haplotype_List *hl);	// return genotype index
	double hapfreq(int a);
	Haploid* find(char *h);
	void print(mystream &os);
	~HaploidTable() { if (hpl!=NULL) delete hpl; hpl = NULL; if (loc_ptr!=NULL) delete[] loc_ptr; }
};

class Mating_Haplotype {
	public:
	int nmrks;
	int hash_val;
	char *h[4];
	double pg[4];	// P[comm_gt|mating_haplotype]
	double pset;	// P[observed_offspring_gt|mating_haplotype];
	double frq;		// mating type frequency given observed_offspring_gt, added for weighted hbat calculation 10/25/03
	double q;		// mating type frequency given observed_offspring_gt under H1, used for power calculation
	double es;		// E(S|M,Y) under H1, used for power calculation
	Haplotype_List *hl[2];	// added for haplotype frequency estimation using EM, only pointer to the emhl list
	Mating_Haplotype *next;
	
	Mating_Haplotype(int n, Haplotype_List *h1=NULL, Haplotype_List *h2=NULL, Mating_Haplotype *nxt=NULL);
	int hash();
	int len();		// number of mating types
	int n_heterozygotes(); // return # of heteerozygote parents
	bool equal(Mating_Haplotype *mh);	// test if two mating types are equal
	bool gt_compatible(Genotype *g);
	bool gt_compatible(Genotype *g, int sex, bool xlink);
	bool ht_compatible(Haplotype_List *hl) { return hl!=NULL && (hl->equal(h[0],h[2]) || hl->equal(h[0],h[3]) || hl->equal(h[1],h[2]) || hl->equal(h[1],h[3])); }
	bool unique(int idx); // return false if mo or fa haplotype contains kotherallele
	//void update_pg(Genotype *gt[], int ng);
	void update_pg(Genotype *gt[], int sex[], int ng, bool sexlink);
	double get_p_set(int cnt[4], bool is_observed=false);
	double offspring_ht_prob(Haplotype_List *hl) { return hl==NULL?0:(hl->equal(h[0],h[2])?0.25:0)+(hl->equal(h[0],h[3])?0.25:0)+(hl->equal(h[1],h[2])?0.25:0)+(hl->equal(h[1],h[3])?0.25:0); }
	double offspring_ht_prob(Haplotype_List *hl, int sex, bool sexlink);
	~Mating_Haplotype();
};

class Mating_mm_genotype {	// multimarger genotype
	public:
	bool xlink;
	int nMarkers;
	int *seg[2];
	Genotype *gt[2];
	Mating_mm_genotype *next;
	
	Mating_mm_genotype(int nmrk, Genotype *mg[2], bool sexlink=false);
	bool og_by_origin_order(Genotype &pg1, Genotype &pg2, Genotype &og);	// if yes, og.a[0], a[1] are from pg1, pg2, respectively
	bool resolve(int noffs, Genotype *ogt[]);
	bool resolve_x(int noffs, Genotype *ogt[], int *sex);
	bool compatible_mating_haplotype(Haplotype_List *h1, Haplotype_List *h2, Genotype *ogt, int sex); // if yes, h1 and h2 is compatible with ogt
	Haplotype_List* enumberate_parent_haplotype(bool father);
	Mating_Haplotype* enumerate_mating_haplo(int noffs, Genotype *ogt[], int *sex, int maxmh);	
	
	~Mating_mm_genotype();
};

class Off_Genotype_Pattern {
	public:
	int cnt[4];	// cnt[i] = count of offspring with genotpye=common_ogt[i]
	double	pg;
	Off_Genotype_Pattern *next;
	
	Off_Genotype_Pattern(int count[4], double p, Off_Genotype_Pattern *nxt);
	~Off_Genotype_Pattern() { if (next) delete next; }
	void print() { cout << cnt[0] << "\t" << cnt[1] << "\t" << cnt[2] << "\t" << cnt[3] << "\t" << pg << "\t" << permut_n() << "\n"; }
	double sum_log_tao(int n);
	double permut_n();
	double permut_n_x(int sex[]);	// for xlinked
};

class Sufficient_Stat {
	static const int kMaxMembers = 20;
	public:
	bool xlink;
	int nMarkers;
	int n_offspring;
	int sex[kMaxMembers];
	double *otr[kMaxMembers];
	Genotype *pgt[2];
	Genotype *ogt[kMaxMembers];
	int fa_bool; //!/
	int mo_bool; //!/
	int maxcmh;
	Mating_Haplotype *mhap;
	
	int n_common_ogt;
	int common_ogt_sex[4];
	Genotype *common_ogt[4];		// offspring genotypes shared by all compatible mating type
	Haplotype_List *common_oht[4];	// lists of offspring haplotypes compatible for common_ogt[] and mhap

	Off_Genotype_Pattern *ogpat; // list of offspring patterns with the same cpf as observed and constant frequency ratio over observed in mhap
	
	double pg[4];
	double pgg[4][4];
	
	double pg_x[2][4];	//pg[sex][]
	double pgg_x[3][4][4];	// pgg[sex_pair][][]; 0 = male-male, 1 = female-female; 2 = male-female
	
	// used for permute
	ColumnVector *EGX;
	Matrix *GX;
	
	Sufficient_Stat(NUCFAMILY *fam, int nmrk, int loc_index[], bool sexlink, int maxcmht=1000);
	Sufficient_Stat(Genotype *fa_gt, Genotype *mo_gt,  Genotype *off_gt[], int sx[], int nmrk, int n_off, bool sexlink, int maxcmht);
	
	int enumerate_genotype_from_allele(int idx, Genotype *ppgt);	// enumerate genotype from observed alleles
	bool is_compatible_mating(Genotype &fg, Genotype &mg, int idx);	// check if fg/mg compatible with ogt
	Compatible_parent_genotype* enumerate_compatible_mating_genotype(int idx);
	
	Mating_Haplotype* compatible_mating_haplotypes();
	int mhap_hash(int &hash_code);	// return length of mhap list
	
	int get_common_ogt();
	int which_comm_ogt(Genotype *g, int sx);
	void compatible_offspring_gt_set();	// Lstar*
	void permute_off_set_distrib(int gid[], int ng);	// permute the number of offspring in each genotype group given the gt set gid[]
	void check_and_add_offpat(int cnt[]); // test if cnt[] (# of common_ogt[]) has fixed p ratio (Lstar**), add to ogpat if true
	
	void get_common_oht();		// derive common_oht[n_common_ogt]
	void update_comm_oht_freq(Haplotype_List *emlist);
	
	void pgpgg2();	// calculate pg[], pgg[][]
	void pgpgg_x();	// calculate pg_x[], pgg_x[][]
	void analyze(Haplotype_List *emlist=NULL); // calculate until pgppg2
	
	void output_info(mystream &os, HaploidTable *htab, char *fam_id, bool informative_only=false); // informative_only: only print if this families's Var(S)>0
	int fbat_stat_observed(int trt_id, double offset, double &tsum, double &tssq, ColumnVector &S, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab);
	int fbat_stat_observed_x(int trt_id, double offset, double tsum[], double tssq[], int noff[], ColumnVector &S, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab);
	bool fbat_stat(ColumnVector &EX, Matrix &VX, Matrix &CX, double gx[][3], HaploidTable *htab=NULL);	// VX=Var(X_i), CX=COV(X_i,X_j), gx[3]=gt coding
	bool fbat_stat_x(ColumnVector &EXm, ColumnVector &EXf, Matrix &VXm, Matrix &VXf, Matrix &CXmm, Matrix &CXff, Matrix &CXmf, double gx[][3], HaploidTable *htab=NULL);
	bool fbat_stat_x(int a_idx, double &exm, double &exf, double &vxm, double &vxf, double &cxmm, double &cxff, double &cxmf, double gx[3], HaploidTable *htab=NULL);
	//bool fbat_stat2(int trt_id, double offset, ColumnVector &ES, ColumnVector &S, Matrix &V, double gx[], int flag, HaploidTable *htab=NULL);
	bool fbat_stat3(int trt_id, double offset, ColumnVector &ES, ColumnVector &S, Matrix &V, Matrix &A, Matrix &B, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab=NULL);
	bool fbat_stat3_x(int trt_id, double offset, ColumnVector &ES, ColumnVector &S, Matrix &V, Matrix &A, Matrix &B, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab);
	//bool fbat_stat4(int trt_id, Matrix &VX, Matrix &CX, double &tsum, double &tssq, int &noff, double gx[], int flag, HaploidTable *htab=NULL);

	void comm_ogt_X(ColumnVector &X, int which, double gx[][3], HaploidTable *htab=NULL); // get a vector observed xscore X for which common_ogt
	double comm_ogt_X(int a_idx, int which, double gx[3], HaploidTable *htab=NULL);
	bool fbat_ex(ColumnVector &EX, double gx[][3], HaploidTable *htab=NULL);
	
	bool informative() { return (mhap!=NULL && n_common_ogt>1 && ogpat!=NULL); }
	bool marker_prob(int a_idx, double mrk_prob[]);
	
	void init_permute_X(int na, double gscore[][3], HaploidTable *htab=NULL);	// initialize GX, rearrange ogt, and EGX
	void permute(ColumnVector &S, int trt_id, double offset, int flag); // permute offspring genotype, return S of the permuted data
	
	double ex(int a_idx, double gx[][3], HaploidTable *htab=NULL);
	void ex(double mx[2], int a_idx, double gx[][3], HaploidTable *htab=NULL);
	bool xdist(int a_idx, double &ex, double &vx, double &cx, double gx[][3], HaploidTable *htab=NULL);
	double fbat_stat1(int trt_id, int na, int a_idx, double offset, double gx[][3], int flag, double &s, HaploidTable *htab=NULL); // return var(s) for allele a_idx, s=S-ES for allele idx
	
	bool fbat_stat_mp(int ntr, int trt_id[], int a_idx, ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab=NULL); // S & V for allele a_idx and multiple phenotypes 
	bool fbat_stat_x_mp(int ntr, int trt_id[], int a_idx, ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab=NULL);
	bool fbat_gee_mp(int na, int ntr, int trt_id[], ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab=NULL);
	bool fbat_gee_x_mp(int na, int ntr, int trt_id[], ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab);
	double fbat_es_h1(Trait_Info *tr_info, int trt_id, double offset, int a_idx, double gx[], double pa); // return E(S|H1)-E(S|H0) for allele a_idx where H1 is defined in trt[trt_id]. Only for bi-allelic markers
	
	~Sufficient_Stat();

};

