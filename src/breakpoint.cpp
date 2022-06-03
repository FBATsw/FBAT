/*
 *  breakpoint.cpp
 *  fbat
 *
 *  Created by Xin Xu on 1/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <cstdlib>
#include "breakpoint.h"
#include "rand.h"

int compare_ord_idx (const void * a, const void * b)
{
	return ((ordered_index*)a)->pos - ((ordered_index*)b)->pos;
}

// set og.a[0] from pg1, og.a[1] from pg2, return true if phasing is unambiguous
// if clear_gt && phase is ambiguous, reset og to 0
bool breakpoint::og_origin(Genotype &pg1, Genotype &pg2, Genotype &og) {
	
	if (og.a[0]==0) return false;
	
	if (og.homozygote()) return true;
	
	int i, j;
	
	j = -1;	// og.a[j] from pg1
	for (i=0; i<2; i++) {
		if (pg1.a[0]>0 && pg1.a[1]>0 && pg1.a[0]!=og.a[i] && pg1.a[1]!=og.a[i]) {
			j = 1-i;
			break;
		}
		else if (pg2.a[0]>0 && pg2.a[1]>0 && pg2.a[0]!=og.a[i] && pg2.a[1]!=og.a[i]) {
			j = i;
			break;
		}
	}
	
	if (j==1)
		og.swap();
	
	return j>=0;
	
}

void breakpoint::order_marker(int chr) {
	
	ordered_index *od = new ordered_index[nloc];
	
	int i;
	for (n_ord=0,i=0; i<nloc; i++)
		if (loc[i].chrNum==chr && loc[i].pos>0) {
			od[n_ord].idx = i;
			od[n_ord].pos = loc[i].pos;
			n_ord++;
		}
	qsort(od, n_ord, sizeof(ordered_index), compare_ord_idx);
	
	if (ord!=NULL)
		delete[] ord;
	ord = new int[n_ord];
	for (i=0; i<n_ord; i++)
		ord[i] = od[i].idx;
	
	delete[] od;

}

float breakpoint::geno_drop_rate(Genotype *g) {
	
	if (g==NULL)
		return 1.0;
	
	int i, k;
	for (k=0,i=0; i<n_ord; i++)
		if (g[ord[i]].defined())
			k++;
	
	return 1.0 - (float)k/(float)n_ord;

}

// changed to do one parent too
bool breakpoint::get(NUCFAMILY *fam, int min_informative_flank_mrks, float max_drop_rate, bool from_both_parents) {
	
	nbrks[0]=0; 
	nbrks[1]=0;
	avg_crossover[0] = -1;
	avg_crossover[1] = -1;
	
	if (n_ord < kMinMarkers_for_brk || fam==NULL)
		return false;
	
	bool p_missing[2], p_fail[2];
	
	p_missing[0] = (fam->fa()==NULL || fam->fa()->gt==NULL);
	p_missing[1] = (fam->mo()==NULL || fam->mo()->gt==NULL);
	p_fail[0] = (p_missing[0] || geno_drop_rate(fam->fa()->gt) >= max_drop_rate);
	p_fail[1] = (p_missing[1] || geno_drop_rate(fam->mo()->gt) >= max_drop_rate);

	if (p_fail[0] && p_fail[1])
		return false;
	
	if (from_both_parents && (p_fail[0] || p_fail[1]))
		return false;
	
	int nsibs = fam->sib_count();
	if (nsibs<2) 
		return false;
	
	PERSONLIST *sl;
	Genotype *pg1, *pg2, *og;
	
	int i, j, k, l, m;
	float g_drop_rate[kMaxSibs];
	for (j=0,k=0,sl=fam->sibs(); sl!=NULL; sl=sl->next,k++) {
		g_drop_rate[k] = (sl->node!=NULL && (og=sl->node->gt)!=NULL? geno_drop_rate(og) : 1.0);
		if (g_drop_rate[k]<max_drop_rate)
			j++;
	}
	if (j<2)
		return false;
	
	int **phase[2];	//[fa|mo][nsibs][nloc], =i:parental origin = pg.a[i], =-1 if pg is homozygotes or uninformative
	for (i=0; i<2; i++) {
		phase[i] = new int*[nsibs];
		for (j=0; j<nsibs; j++) {
			phase[i][j] = new int[n_ord];
			for (k=0; k<n_ord; k++)
				phase[i][j][k] = -1;
		}
	}
	
	// set phase
	pg1 = (p_missing[0]?new Genotype[nloc] : fam->fa()->gt);
	pg2 = (p_missing[1]?new Genotype[nloc] : fam->mo()->gt);
	int info_cnt[2], max_info_cnt[2], sib_idx_max__info_cnt[2];	// the sib with maximum infomative phasing points
	for (i=0; i<2; i++) {
		max_info_cnt[i] = -1;
		sib_idx_max__info_cnt[i] = -1;
	}
	for (k=0,sl=fam->sibs(); sl!=NULL; sl=sl->next,k++) {
		if (sl->node!=NULL && (og=sl->node->gt)!=NULL) {
			info_cnt[0] = 0;
			info_cnt[1] = 0;
			for (i=0; i<n_ord; i++) {
				j = ord[i];
				if (og_origin(pg1[j], pg2[j], og[j])) {
					if (pg1[j].a[0]>0 && pg1[j].a[0]!=pg1[j].a[1]) {
						phase[0][k][i] = (pg1[j].a[0]==og[j].a[0]? 0 : 1);
						info_cnt[0]++;
					}
					if (pg2[j].a[0]>0 && pg2[j].a[0]!=pg2[j].a[1]) {
						phase[1][k][i] = (pg2[j].a[0]==og[j].a[1]? 0 : 1);
						info_cnt[1]++;
					}
				}
			}
			for (i=0; i<2; i++)
				if (info_cnt[i]>max_info_cnt[i])
					sib_idx_max__info_cnt[i] = k;
			if (g_drop_rate[k]<max_drop_rate)
				sl->node->setflag(flag_phased);
		}
	}
	
	int *nb = new int[nsibs];
	int **bpos = new int*[nsibs];
	for (i=0; i<nsibs; i++)
		bpos[i] = new int[n_ord];
	int nb_merge;
	int *bpos_merge = new int[n_ord];
	int *idx = new int[n_ord];
	int *pat = new int[nsibs];
	int min_cr_cnt, min_cr_pat;
	bool redaundancy;
	
	for (i=0; i<2; i++) {
		if (p_fail[i])
			continue;
		nb_merge = 0;
		for (k=0; k<nsibs; k++)	// use sib_idx_max__info_cnt as reference
			if (k==sib_idx_max__info_cnt[i])
				nb[k] = 0;
			else {
				nb[k] = compare_sib_phase(phase[i][sib_idx_max__info_cnt[i]], phase[i][k], n_ord, bpos[k], min_informative_flank_mrks);
				// merge breakpoint by location
				for (l=0,j=0; j<nb[k]; j++) {
					while (l<nb_merge && bpos[k][j]>bpos_merge[l])
						l++;
					if (l==nb_merge || bpos[k][j]<bpos_merge[l]) {
						for (m=nb_merge; m>l; m--)
							bpos_merge[m] = bpos_merge[m-1];
						bpos_merge[l] = bpos[k][j];
						l++;
						nb_merge++;
					}
				}
			}
		// check redaundancy with bpos_merge, if bpos_merge[l] & bpos_merge[l+1] is close, and no informative markers lie between in sibs with cross_over at bpos_merge[l+1]
		for (l=1; l<nb_merge; l++) {
			redaundancy = true;
			for (k=0; k<nsibs; k++) {
				for (j=0; j<nb[k] && bpos[k][j]!=bpos_merge[l]; j++)
					;
				if (j<nb[k]) {	// match, check informativeness of markers between =bpos_merge[l] & =bpos_merge[l-1]
					for (m=bpos_merge[l]-1; m>=bpos_merge[l-1] && (phase[i][k][m]<0 || phase[i][sib_idx_max__info_cnt[i]][m]<0); m--)
						;
					if (m>=bpos_merge[l-1])
						redaundancy = false;
					else
						bpos[k][j] = bpos_merge[l-1];
				}
			}
			if (redaundancy) {
				for (j=l; j<nb_merge-1; j++)
					bpos_merge[j] = bpos_merge[j+1];
				nb_merge--;
			}
		}
		
		if (nb_merge>kMaxBreakPoints)
			cout << "WARNING: too many breakpoints observed. nb_merge=" << nb_merge << "\n";
		// assign sibs' crossover pattern according to reference
		for (k=0; k<nsibs; k++) {
			// map bpos[k] into index for bpos_merge[]
			pat[k] = 0;
			for (l=0,j=0; j<nb[k]; j++) {
				while (l<nb_merge && bpos[k][j]!=bpos_merge[l])
					l++;
				if (l>=nb_merge)
					cout << "SORRY, INTERNAL ERROR OCCURED\n";
				idx[j] = l;
			}
			idx[nb[k]] = nb_merge;
			for (j=0; j<nb[k]; j+=2)
				for (l=idx[j]; l<idx[j+1]; l++)
					pat[k] |= (2<<l);
		}
		
		// permute ref pat, search for a ref pat with minimal # of crossover
		min_cr_cnt = -1;
		min_cr_pat = 0;
		for (l=0; l<(2<<(nb_merge)); l++) {
			for (m=0,j=0; j<nsibs; j++)
				m += crossover(l, pat[j], nb_merge);
			if (min_cr_cnt==-1 || m<min_cr_cnt) {
				min_cr_cnt = m;
				min_cr_pat = l;
			}
		}
		
		// update n_crossover_sib
		for (l=0,j=0,k=0; k<nsibs; k++) {
			n_crossover_sib[i][k] = crossover(min_cr_pat, pat[k], nb_merge);
			if (g_drop_rate[k]<max_drop_rate) {
				l += n_crossover_sib[i][k];
				j++;
			}
		}
		avg_crossover[i] = (float)l/(float)j;
		
		// debug
		/*
		cout << "\nget parent " << i <<  "\nref sib =" << sib_idx_max__info_cnt[i] << "\n";
		for (k=0; k<nsibs; k++) {
			cout << "P" << i << "S" << k << ": ";
			for (j=0; j<nb[k]; j++)
				cout << bpos[k][j] << " ";
			cout << "\n";
		}
		cout << "nb_merge=" << nb_merge << " : ";
		for (k=0; k<nb_merge; k++)
			cout << bpos_merge[k] << " ";
		cout << "\n";
		
		cout << "ref sib =" << sib_idx_max__info_cnt[i] << ", pat=\n";
		for (k=0; k<nsibs; k++) {
			cout << "P" << i << "S" << k << ": ";
			for (l=pat[k]; l>0; l=(l>>1))
				cout << (l&1) << " ";
			cout << "\n";
		}
		cout << "final parental pat = ";
		for (l=min_cr_pat; l>0; l=(l>>1))
			cout << (l&1) << " ";
		cout << "\n";
		for (k=0; k<nsibs; k++) {
			cout << "P" << i << "S" << k << "=" << n_crossover_sib[i][k] << ": ";
			for (l=pat[k]^min_cr_pat; l>0; l=(l>>1))
				cout << (l&1) << " ";
			cout << "\n";
		}
		*/
		// end debug
	}
	
	if (p_missing[0])
		delete[] pg1;
	if (p_missing[1])
		delete[] pg2;
	
	for (i=0; i<2; i++) {
		for (j=0; j<nsibs; j++)
			delete [] phase[i][j];
		delete[] phase[i];
	}
	
	delete[] nb;
	delete[] idx;
	delete[] pat;
	delete[] bpos_merge;
	for (i=0; i<nsibs; i++)
		delete[] bpos[i];
	delete[] bpos;
	
	return true;
}

// return crossover# between pattern p1 & p2
int breakpoint::crossover(int p1, int p2, int digits) {
	
	int p = (p1^p2);
	int cnt = (p>0 && (p&(1<<(digits)))==0?1:0);
	int last = (p&1);
	for (p=(p>>1); p>0; p=(p>>1)) {
		if ((p&1) != last) {
			cnt++;
			last = (p&1);
		}
	}
	
	return cnt;
}

// return # of breakpoint between 2 sibs' phases
int breakpoint::compare_sib_phase(int *p1, int *p2, int nmrk, int *brk_pos, int min_informative_flank_mrks) {
	
	int i, j, newp, oldp = -1;
	
	int *informative_cnt = new int[nmrk];
	for (i=0; i<nmrk; i++)
		informative_cnt[i] = 0;

	int nbrks = 0;
	for (i=0; i<nmrk; i++)
		if (p1[i]>=0 && p2[i]>=0) {
			newp = (p1[i]==p2[i]?1:0);
			informative_cnt[nbrks]++;
			if (oldp>=0 && newp!=oldp)
				brk_pos[nbrks++] = i;
			oldp = newp;
		}
	
	// remove crossovers with flanking informative_cnt[] < min_informative_flank_mrks
	j = 0;
	for (i=0; i<nbrks; i++)
		if (informative_cnt[i]>=min_informative_flank_mrks && informative_cnt[i+1]>=min_informative_flank_mrks)
			brk_pos[j++] = brk_pos[i];
	
	delete[] informative_cnt;
	
	return j;

}

// simulation offspring genotype with fixed parental haplotypes and randomized crossovers
void breakpoint::simulate(NUCFAMILY *fam, int nc_min, int nc_max) {
	
	int bpos[2][kMaxSibs][kMaxBreakPoints];

	Genotype *pg[2], *og;
	pg[0] = fam->fa()->gt;
	pg[1] = fam->mo()->gt;
	
	int i, j, k, l, m, n, ph;
	PERSONLIST *sl;
	for (n=0,sl=fam->sibs(); sl!=NULL; sl=sl->next,n++) {
		if (sl->node!=NULL && (og=sl->node->gt)!=NULL) {
			for (i=0; i<2; i++) {
				// generate random crossover pos, which is the bottom nb_th a_idx[] with .pos
				n_crossover_sib_sim[i][n] = nc_min + rand_uniform()*(nc_max-nc_min+1);
				bpos[i][n][0] = 0;
				bpos[i][n][n_crossover_sib_sim[i][n]+1] = n_ord-1;
				for (j=1; j<=n_crossover_sib_sim[i][n];) {
					bpos[i][n][j] = 1 + rand_uniform()*(n_ord-2);
					for (k=1; k<j && bpos[i][n][j]>bpos[i][n][k]; k++)
						;
					if (k==j)
						j++;
					else if (bpos[i][n][j]<bpos[i][n][k]) {
						m = bpos[i][n][j];
						for (l=j; l>k; l--)
							bpos[i][n][l] = bpos[i][n][l-1];
						bpos[i][n][k] = m;
						j++;
					}
				}
				// generate og genotypes
				ph = (rand_uniform()>0.5?1:0);
				for (j=0; j<n_crossover_sib_sim[i][n]+1; j++) {
					for (k=bpos[i][n][j]; k<bpos[i][n][j+1]; k++)
						if (og[ord[k]].a[i]>0)
							og[ord[k]].a[i] = pg[i][ord[k]].a[ph];
					ph = 1 - ph;
				}
			}
		}
	}
	
	for (k=0; k<2; k++) {
		for (i=0; i<n; i++) {
			cout << "parent " << k << ", sib " << i << ": ";
			for (j=0; j<n_crossover_sib_sim[k][i]; j++)
				cout << bpos[k][i][j+1] << " ";
			cout << "\n";
		}
	}
	
}

void breakpoint::get_simulate(mystream &os, NUCFAMILY *fam, int nc_min, int nc_max, int ncycles, int min_informative_flank_mrks, float max_drop_rate) {
	
	if (fam==NULL || fam->fa()==NULL || fam->fa()->gt==NULL || fam->mo()==NULL || fam->mo()->gt==NULL)
		return;
	
	if (geno_drop_rate(fam->fa()->gt) >= max_drop_rate || geno_drop_rate(fam->mo()->gt) >= max_drop_rate)
		return;

	int i, j, k, n;
	bool good_sib[kMaxSibs];
	PERSONLIST *sl;
	Genotype *og;
	for (i=0,n=0,sl=fam->sibs(); sl!=NULL; sl=sl->next,i++)
		if (sl->node!=NULL && (og=sl->node->gt)!=NULL && geno_drop_rate(og)<max_drop_rate) {
			n++;
			good_sib[i] = true;
		} else
			good_sib[i] = false;
	if (n<2)
		return;
	
	int nsibs = fam->sib_count();
	for (int i=0; i<ncycles; i++) {
		simulate(fam, nc_min, nc_max);
		for (j=0; j<2; j++) {
			for(avg_crossover[j]=0,k=0; k<nsibs; k++)
				if (good_sib[k])
					avg_crossover[j] += n_crossover_sib_sim[j][k];
			avg_crossover[j] /= n;
		}
		MYOUTPUT(os, i << "\t" << avg_crossover[0] << "\t" << avg_crossover[1] << "\t")
		get(fam, min_informative_flank_mrks, max_drop_rate);
		MYOUTPUT(os, avg_crossover[0] << "\t" << avg_crossover[1] << "\n")
	}
}

bool breakpoint::trio_phase(PERSON *per, int min_informative_flank_mrks) {

	if (per==NULL || per->gt==NULL || per->ascend==NULL)
		return false;
	
	int i, j;
	bool p_missing[2], p_fail[2];
	NUCFAMILY *fam = per->ascend;
	PERSON *p[2];
	p[0] = fam->fa();
	p[1] = fam->mo();
	for (i=0; i<2; i++) {
		p_missing[i] = (p[i]==NULL || p[i]->gt==NULL);
		p_fail[i] = (p_missing[i] || geno_drop_rate(p[i]->gt)>=min_informative_flank_mrks);
	}
	if (p_fail[0] && p_fail[1])
		return false;
	
	Genotype *pg[2];
	for (i=0; i<2; i++)
		pg[i] = (p_missing[i]? new Genotype[nloc] : p[i]->gt);
	
	for (i=0; i<n_ord; i++) {
		j = ord[i];
		if (!og_origin(pg[0][j], pg[1][j], per->gt[j])) 
			per->gt[j].set(0,0);	// !!!!!!!!!!!!! destructive here !!!!!!!!!!!!!!!!
	}
	
	for (i=0; i<2; i++)
		if (p_missing[i])
			delete[] pg[i];
	
	return true;
}

// 1+ parent phased, 1 informative sib
bool breakpoint::get_phased(NUCFAMILY *fam, int min_informative_flank_mrks, float max_drop_rate) {

	nbrks[0]=0; 
	nbrks[1]=0;
	avg_crossover[0] = -1;
	avg_crossover[1] = -1;
	
	if (n_ord < kMinMarkers_for_brk || fam==NULL)
		return false;
	
	int i, j, k, m;
	bool p_missing[2], p_fail[2], p_phased[2];
	PERSON *p[2];
	p[0] = fam->fa();
	p[1] = fam->mo();
	for (i=0; i<2; i++) {
		p_missing[i] = (p[i]==NULL || p[i]->gt==NULL);
		p_fail[i] = (p_missing[i] | geno_drop_rate(p[i]->gt)>=max_drop_rate);
		p_phased[i] = (!p_fail[i] && trio_phase(p[i], min_informative_flank_mrks));
	}
	
	if (!p_phased[0] && !p_phased[1])
		return false;
	
	int nsibs = fam->sib_count();
	if (nsibs<1) 
		return false;
	
	PERSONLIST *sl;
	Genotype *pg[2], *og;
	for (i=0; i<2; i++)
		pg[i] = (p_missing[i]? new Genotype[nloc] : p[i]->gt);

	int ph[2], nb[2], tot[2], new_ph;
	int *info_cnt[2], *bpos[2];
	for (i=0; i<2; i++) {
		bpos[i] = new int[n_ord];
		info_cnt[i] = new int[n_ord];
		tot[i] = 0;
	}
	
	for (k=0,sl=fam->sibs(); sl!=NULL; sl=sl->next) {
		if (sl->node!=NULL && (og=sl->node->gt)!=NULL && geno_drop_rate(og)<max_drop_rate) {
			for (i=0; i<2; i++) {
				ph[i] = -1;
				nb[i] = 0;
				for (j=0; j<n_ord; j++) {
					bpos[i][j] = 0;
					info_cnt[i][j] = 0;
				}
			}
			for (i=0; i<n_ord; i++) {
				j = ord[i];
				if (og_origin(pg[0][j], pg[1][j], og[j])) {
					for (m=0; m<2; m++) {
						if (p_phased[m] && pg[m][j].a[0]>0 && pg[m][j].a[0]!=pg[m][j].a[1]) {
							new_ph = (pg[m][j].a[0]==og[j].a[m]? 0 : 1);
							if (ph[m]>=0 && new_ph!=ph[m])	// crossover
								bpos[m][nb[m]++] = i;
							else
								info_cnt[m][nb[m]]++;
							ph[m] = new_ph;
						}
					}
				}
			}
			// remove breakpoint with <min_informative_flank_mrks at two flanks;
			for (i=0; i<2; i++) {
				if (p_phased[i]) {
					for (j=0,m=0; j<nb[i]; j++)
						if (info_cnt[i][j]>=min_informative_flank_mrks && info_cnt[i][j+1]>=min_informative_flank_mrks)
							bpos[i][m++] = bpos[i][j];
					nb[i] = m;
					tot[i] += m;
				}
			}
			k++;
		}
	}
	
	if (k>0)
		for (i=0; i<2; i++)
			if (p_phased[i])
				avg_crossover[i] = (float)tot[i]/(float)k; // since only one sib is informative;

	for (i=0; i<2; i++) {
		delete[] info_cnt[i];
		delete[] bpos[i];
		if (p_missing[i])
			delete[] pg[i];
	}
	
	return k>0;
}
