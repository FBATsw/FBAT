/*
 *  breakpoint.h
 *  fbat
 *
 *  Created by Xin Xu on 1/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "pedigree.h"

const int kMaxBreakPoints = 24;
const int kMinMarkers_for_brk = 50;
const int kMaxSibs = 20;

struct ordered_index {
	int idx;
	int pos;
};

class breakpoint {
	Locus_Info *loc;
	int nloc;
	int n_ord;
	int *ord;	// pos-ordered index of loc with pos<=0 removed
	int nbrks[2];
	int bk_pos[2][kMaxBreakPoints][3];	// [fa|mp][...][left bound|obs|right bound of the crossover];
	
	int compare_sib_phase(int *p1, int *p2, int nmrk, int *brk_pos, int min_informative_flank_mrks=5);
	int crossover(int p1, int p2, int digits);	// return crossover# between pattern p1 & p2
	bool og_origin(Genotype &pg1, Genotype &pg2, Genotype &og);
	bool trio_phase(PERSON *per, int min_informative_flank_mrks);
	float geno_drop_rate(Genotype *g);
	
public:
	int n_crossover_sib[2][kMaxSibs];	// [fa|mo][]
	int n_crossover_sib_sim[2][kMaxSibs];
	float avg_crossover[2];

	breakpoint(Locus_Info *l, int nl) { loc=l; nloc=nl; n_ord=0; ord=NULL; nbrks[0]=0; nbrks[1]=0; avg_crossover[0]=-1; avg_crossover[1]=-1;}
	void order_marker(int chr);
	bool get(NUCFAMILY *fam, int min_informative_flank_mrks=5, float max_drop_rate=0.1, bool from_both_parents=true);
	bool get_phased(NUCFAMILY *fam, int min_informative_flank_mrks=5, float max_drop_rate=0.1);	// 1+ parents are phased
	void simulate(NUCFAMILY *fam, int nc_min, int nc_max); // min and max # of cross
	void get_simulate(mystream &os, NUCFAMILY *fam, int nc_min, int nc_max, int ncycles=100, int min_informative_flank_mrks=5, float max_drop_rate=0.1);
	~breakpoint() { delete[] ord; }
	
};


