/*
 *  sufficient_stat.cpp
 *  fbat
 *
 *  Created by Xin Xu on 12/14/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <cmath>
#include "sufficient_stat.h"
//#include "newmatio.h"
#include "myerror.h"
#include "mystream.h"
#include "rand.h"

static const double very_small_number	= 1.0e-6;
static const int flag_censored			= 1<<18;

Haplotype_List::Haplotype_List(int n) {
	
	for (int i=0; i<2; i++) {
		h[i] = new char[n];
		for (int j=0; j<n; j++)
			h[i][j] = 0;
	}
	nmrk = n;
	next = NULL;
	p = 0.0;
	q = 0.0;
	
}

Haplotype_List::Haplotype_List(int n, char *h1, char *h2, Haplotype_List *nxt) {
	
	int i;
	for (i=0; i<2; i++)
		h[i] = new char[n];
	for (i=0; i<n; i++) {
		h[0][i] = h1[i];
		h[1][i] = h2[i];
	}
	
	nmrk = n;
	next = nxt;
	p = 0.0;
	q = 0.0;

}

bool Haplotype_List::equal(char *h1, char *h2) {
	
	int i;
	
	for (i=0; i<nmrk && h[0][i]==h1[i] && h[1][i]==h2[i]; i++)
		;
	if (i==nmrk)
		return true;
	
	for (i=0; i<nmrk && h[0][i]==h2[i] && h[1][i]==h1[i]; i++)
		;
	
	return (i==nmrk);

}

Haplotype_List* Haplotype_List::find(char *h1, char *h2) {
	
	Haplotype_List *hl=this;
	while (hl!=NULL && !hl->equal(h1,h2))
		hl = hl->next;
		
	return hl;

}
	
Haplotype_List::~Haplotype_List() {
	
	for (int i=0; i<2; i++)
		if (h[i] != NULL)
			delete[] h[i];
	
	if (next)
		delete next;
}

Mating_Haplotype::Mating_Haplotype(int n, Haplotype_List *h1, Haplotype_List *h2, Mating_Haplotype *nxt) {
	
	int i, j;
	
	nmrks = n;
	for (i=0; i<4; i++)
		h[i] = new char[n];
	for (j=0; j<2; j++)
		for (i=0; i<n; i++) {
			h[j][i] = h1!=NULL? h1->h[j][i] : 0;
			h[j+2][i] = h2!=NULL? h2->h[j][i] : 0;
		}
		
	hash_val = (h1!=NULL && h2!=NULL)? hash() : 0;
	for (i=0; i<4; i++)
		pg[i] = 0.0;
	pset = 0.0;
	frq = 0.0;
	
	hl[0] = NULL;
	hl[1] = NULL;

	next = nxt;
	
}

int Mating_Haplotype::hash() {
	
	int i, j, k;
	int a[4];
	
	for (i=0; i<4; i++) {
		a[i] = 0;
		for (j=0; j<nmrks; j++) {
			k = h[i][j];
			a[i] ^= (k<<(j&8));
		}
		a[i] &= 255;
	}
		
	i = a[0]<a[1]? a[0] : a[1];
	j = a[2]<a[3]? a[2] : a[3];
	if (i<j || (i==j && a[0]+a[1]<a[2]+a[3]))
		k = i + ((a[0]+a[1]-i)<<8) + (j<<16) + ((a[2]+a[3]-j)<<24);
	else
		k = j + ((a[2]+a[3]-j)<<8) + (i<<16) + ((a[0]+a[1]-i)<<24);
	
	return k;
}

int Mating_Haplotype::len() {
	
	int i=0;
	
	for (Mating_Haplotype *mh=next; mh!=NULL; mh=mh->next)
		i++;
	
	return i+1;
}

bool Mating_Haplotype::equal(Mating_Haplotype *mh) {
	
	if (mh->hash_val != hash_val)
		return false;
	
	if (((strncmp(h[0],mh->h[0],nmrks)==0 && strncmp(h[1],mh->h[1],nmrks)==0) 
		|| (strncmp(h[0],mh->h[1],nmrks)==0 && strncmp(h[1],mh->h[0],nmrks)==0)) &&
		((strncmp(h[2],mh->h[2],nmrks)==0 && strncmp(h[3],mh->h[3],nmrks)==0) 
		|| (strncmp(h[2],mh->h[3],nmrks)==0 && strncmp(h[3],mh->h[2],nmrks)==0)))
		return true;
	
	if (((strncmp(h[0],mh->h[2],nmrks)==0 && strncmp(h[1],mh->h[3],nmrks)==0) 
		|| (strncmp(h[0],mh->h[3],nmrks)==0 && strncmp(h[1],mh->h[2],nmrks)==0)) &&
		((strncmp(h[2],mh->h[0],nmrks)==0 && strncmp(h[3],mh->h[1],nmrks)==0) 
		|| (strncmp(h[2],mh->h[1],nmrks)==0 && strncmp(h[3],mh->h[0],nmrks)==0)))
		return true;
	
	return false;
}

// return # of homoterozygote parents
int Mating_Haplotype::n_heterozygotes() {
	
	int i, j, n=0;
	for (i=0; i<3; i+=2) {
		for (j=0; j<nmrks && h[i][j]==h[i+1][j];)
			j++;
		if (j==nmrks)
			n++;
	}
	
	return 2-n;

}
		
Mating_Haplotype::~Mating_Haplotype() {
	
	for (int i=0; i<4; i++)
		if (h[i])
			delete[] h[i];
	
	
	if (next)
		delete next;
}

bool Mating_Haplotype::gt_compatible(Genotype *g) {
	
	int i, j, k;
	
	for (i=0; i<2 ; i++)
		for (j=2; j<4; j++) {
			for (k=0; k<nmrks && (g[k].a[0]==0 || g[k].equal(h[i][k],h[j][k])); k++)
				;
			if (k==nmrks)
				return true;
		}
	
	return false;
}

/*
void Mating_Haplotype::update_pg(Genotype *gt[], int ng, int sex[], bool xlink) {
	
	if (ng>4)
		ng = 4;
	
	int i, j, k, l;
	
	for (i=0; i<ng; i++)
		pg[i] = 0.0;
		
	for (i=0; i<(xlink?1:2) ; i++)
		for (j=2; j<4; j++) 
			for (l=0; l<ng; l++) {
				for (k=0; k<nmrks && (gt[l][k].a[0]==0 || gt[l][k].equal((xlink && sex[l]==1?0:h[i][k]),h[j][k])); k++)
					;
				if (k==nmrks) {
					pg[l] += xlink?0.5:0.25;
					break;
				}
			}
	
}
*/

bool Mating_Haplotype::unique(int idx) {
	
	int j;
	
	for (j=0; j<nmrks && h[idx*2][j]!=kotherallele && h[idx*2+1][j]!=kotherallele; j++)
		;
	
	return (j==nmrks);

}

double Mating_Haplotype::get_p_set(int cnt[4], bool is_observed) {
	
	double p = 1.0;
	for (int i=0; i<4; i++)
		if (cnt[i]>0)
			p *= pow(pg[i], cnt[i]);
	
	if (is_observed)	// record to pset
		pset = p;
		
	return p;
}

Haploid::Haploid(int n, char *ht, double f) {
	
	if (n<1) {
		nmrk = 0;
		h = NULL;
	} else {
		h = new char[n];
		for (int i=0; i<n; i++)
			h[i] = ht[i];
		nmrk = n;
	}
	
	p = f;
	q = 0;
}
	
bool Haploid::equal(char *ht) {
	
	int i;
	
	for (i=0; i<nmrk && h[i]==ht[i]; i++)
		;
	
	return i==nmrk;

}

void HaploidTable::maketable(Haplotype_List *hl) {
	
	const double minp = 0.0001;
	
	if (hl==NULL)
		return;
	
	Haplotype_List *hl2;
	Haploid *hld;
	int i;
	size = 0;
	for (hl2=hl; hl2!=NULL; hl2=hl2->next) {
		if (hl2->p>=minp)
			for (i=0; i<2; i++) {
				if (hl2->h[i][0]>0) {	// not missing as in male x-linked markers
					if ((hld=find(hl2->h[i]))==NULL) {
						if (hpl==NULL)
							hpl = new HaploidList(new Haploid(hl2->nmrk, hl2->h[i], hl2->p/2), hpl);
						else
							new HaploidList(new Haploid(hl2->nmrk, hl2->h[i], hl2->p/2), hpl);
						size++;
					} else
						hld->p += hl2->p/2;
				}
			}
	}
	
	// sort hpl by frequency
	if (hpl==NULL || hpl->next==NULL)
		return;
	HaploidList *hpl2, *hpl1, *maxh;
	double maxp;
	for (hpl1=hpl; hpl1->next!=NULL; hpl1=hpl1->next) {
		maxp = hpl1->node->p;
		maxh = hpl1;
		for (hpl2=hpl1->next; hpl2!=NULL; hpl2=hpl2->next)
			if (hpl2->node->p>maxp) {
				maxp = hpl2->node->p;
				maxh = hpl2;
			}
		if (maxh!=hpl1) {	// swap node
			hld = hpl1->node;
			hpl1->node = maxh->node;
			maxh->node = hld;
		}
	}
	
	// normlize the p so it add to 1
	double sum = 0;
	for (hpl1=hpl; hpl1!=NULL; hpl1=hpl1->next)
		sum += hpl1->node->p;
	for (hpl1=hpl; hpl1!=NULL; hpl1=hpl1->next)
		hpl1->node->p /= sum;
}

Haploid* HaploidTable::find(char *h) {
	
	HaploidList *hpl2;
	
	for (hpl2=hpl; hpl2!=NULL; hpl2=hpl2->next)
		if (hpl2->node!=NULL && hpl2->node->equal(h))
			break;
	
	return hpl2==NULL? NULL : hpl2->node;

}

int HaploidTable::allele(char *h) {
	
	HaploidList *hpl2;
	int i = 1;
	for (hpl2=hpl; hpl2!=NULL; hpl2=hpl2->next,i++)
		if (hpl2->node!=NULL && hpl2->node->equal(h))
			break;
	
	return hpl2==NULL? 0 : i;

}

int HaploidTable::gid(Haplotype_List *hl) {
	
	if (hl==NULL || hl->h[0]==NULL || hl->h[1]==NULL)
		return -1;

	int i, a[2];
	for (i=0; i<2; i++)
		a[i] = allele(hl->h[i]);
	
	return a[0]==0 || a[1]==0? -2 : (a[0]<=a[1]? ((a[1]-1)*a[1])/2+a[0]-1 : ((a[0]-1)*a[0])/2+a[1]-1);

}
	
void HaploidTable::dprime() {
	
	if (hpl==NULL || hpl->node==NULL)
		return;
		
	const int kmaxallele = 40;
	
	int i, k, j=hpl->node->nmrk;
	double **f = new double*[j];
	for (i=0; i<j; i++) {
		f[i] = new double[kmaxallele];
		for (k=0; k<kmaxallele; k++)
			f[i][k] = 0.0;
	}
	
	HaploidList *hpl2;
	
	for (hpl2=hpl; hpl2!=NULL; hpl2=hpl2->next)
		if (hpl2->node!=NULL)
			for (i=0; i<j; i++)
				f[i][hpl2->node->h[i]] += hpl2->node->p;
	
	double minf, a, b;
	for (hpl2=hpl; hpl2!=NULL; hpl2=hpl2->next)
		if (hpl2->node!=NULL) {
			a = 1.0;
			minf = 1.0;
			for (i=0; i<j; i++) {
				b = f[i][hpl2->node->h[i]];
				a *= b;
				if (b<minf)
					minf = b;
			}
			b = fabs(minf-a);
			if (b<a)
				b = a;
			hpl2->node->q = (hpl2->node->p-a)/b;
		}
	
	for (i=0; i<j; i++)
		if (f[i]!=NULL)
			delete[] f[i];
	
	delete[] f;

}

double HaploidTable::hapfreq(int a) {
	
	if (hpl==NULL || hpl->node==NULL || a<1 || a>size)
		return 0.0;
	
	int i=1;
	HaploidList *hpl2;
	
	for (hpl2=hpl; hpl2!=NULL && i<a; hpl2=hpl2->next, i++)
		;
	
	return hpl2!=NULL? hpl2->node->p : 0.0;

}

double Off_Genotype_Pattern::sum_log_tao(int n) {
	
	int i;
	double t = 0.0;
	
	for (i=2; i<=n; i++)
		t += log((double)i);
	
	return t;
}

// # of permutation of cnt[] = sum_cnt!/(cnt[0]! * cnt[1]! * cnt[2]! * cnt[3]!)
double Off_Genotype_Pattern::permut_n() {
	
	int i, k;
	double t = 0.0;
	
	for (k=0,i=0; i<4; i++) {
		t += sum_log_tao(cnt[i]);
		k += cnt[i];
	}
	return exp(sum_log_tao(k)-t); 
	
}

// # of permutation of cnt[] conditioning on sex[] = product of permute_n in male and permute_n in female
double Off_Genotype_Pattern::permut_n_x(int sex[]) {
	
	int i, k[2];
	double t = 0.0;
	
	k[0] = 0;
	k[1] = 0;
	for (i=0; i<4; i++) 
		if (cnt[i]>0) {
			t += sum_log_tao(cnt[i]);
			k[sex[i]==1?0:1] += cnt[i];
		}
	return exp(sum_log_tao(k[0])+sum_log_tao(k[1])-t); 
	
}

Off_Genotype_Pattern::Off_Genotype_Pattern(int count[4], double p, Off_Genotype_Pattern *nxt) {
	
	for (int i=0; i<4; i++)
		cnt[i] = count[i];
	pg = p;
	next = nxt;
	
}

Mating_mm_genotype::Mating_mm_genotype(int nmrk, Genotype *mg[2], bool sexlink) {
	
	nMarkers = nmrk;
	xlink = sexlink;
	
	int i, j;
	for (i=0; i<2; i++) {
		gt[i] = mg[i];
		seg[i] = new int[nMarkers];
		for (j=0; j<nMarkers; j++)
			seg[i][j] = gt[i][j].homozygote()? 1 : 0;	// 1 = homozygote parent
	}
	
	next = NULL;
}

// order the two alleles of og by parental orgin such that og.a[0] from father (pg1), og.a[1] from mother
bool Mating_mm_genotype::og_by_origin_order(Genotype &pg1, Genotype &pg2, Genotype &og) {
	
	if (og.a[0]==0) return false;
	
	if (og.homozygote()) return true;
	
	if (og.equal(pg1) && og.equal(pg2)) return false;
	
	int i, j;
	
	j = -1;	// og.a[j] from pg1
	for (i=0; i<2; i++) {
		if (pg1.a[0]!=og.a[i] && pg1.a[1]!=og.a[i]) {
			j = 1-i;
			break;
		}
		else if (pg2.a[0]!=og.a[i] && pg2.a[1]!=og.a[i]) {
			j = i;
			break;
		}
	}
	
	if (j==1)
		og.swap();
	
	return j>=0;

}

// resole the phase using ogt for non-x-linked markers
bool Mating_mm_genotype::resolve(int noffs, Genotype *ogt[]) {
	
	if (nMarkers<2)
		return true;
		
	int i, j, k, l, p, phase;
	bool *oried = new bool[nMarkers];
	for (i=0; i<noffs; i++) {
		for (j=0; j<nMarkers; j++)
			oried[j] = og_by_origin_order(gt[0][j], gt[1][j], ogt[i][j]);
		for (k=0; k<2; k++) {
			for (j=0; j<nMarkers; j++)
				if (oried[j] && seg[k][j]!=1) {	// informative heterzygotes
					phase = gt[k][j].a[0]==ogt[i][j].a[k]? 0 : 1;
					if (seg[k][j]>1 && seg[k][j]!=i+2) {	// has overlap of heterozygote hap fragment
						p = seg[k][j];
						for (l=0; l<nMarkers; l++) {	// combinate overlapped hap
							if (seg[k][l]==p) {
								if (oried[l] && ogt[i][l].a[k]!=gt[k][l].a[phase]) {	// crossover
									if (oried!=NULL)
										delete[] oried;
									return false;
								}
								if (phase==1)
									gt[k][l].swap();
								seg[k][l] = i+2;
							}
						}
					} else {
						if (phase==1)
							gt[k][j].swap();
						seg[k][j] = i+2;
					}
				}
		}
	}
	
	//sorting seg numbers, start from 1001 - an arbitary big number
	for (k=0; k<2; k++)
		for (p=1001,i=0; i<nMarkers; i++)
			if (seg[k][i]>1 && seg[k][i]<1000) {
				for (l=seg[k][i],j=i; j<nMarkers; j++)
					if (seg[k][j]==l)
						seg[k][j] = p;
				p++;
			}
	
	// bring seg number back starting at 2 for #>1000 (2=1001, 3=1002, etc)
	for (k=0; k<2; k++)
		for (i=0; i<nMarkers; i++)
			if (seg[k][i]>1000)
				seg[k][i] -= 999;
	
	if (oried!=NULL)		
		delete[] oried;
		
	return true;
}

// resole the phase of gt[1] (mother) using ogt for x-linked markers
bool Mating_mm_genotype::resolve_x(int noffs, Genotype *ogt[], int *sex) {
	
	if (nMarkers<2)
		return true;
		
	int i, j, k, l, p;
	int phase;
	for (k=0; k<noffs; k++) {
		if (sex[k]==2) {	// phase ogt[] first
			for (i=0; i<nMarkers; i++)
				if (gt[0][i].a[0]==ogt[k][i].a[0])	// here ogt[k][i].a[0] is maternal alleles
					ogt[k][i].swap();
		}
		for (i=0; i<nMarkers; i++) {
			if (seg[1][i]!=1 && seg[1][i]!=k+2 && ogt[k][i].a[0]>0) {	// heterozygotes
				phase = gt[1][i].a[0]==ogt[k][i].a[0]? 0 : 1;
				if (seg[1][i]>1 && seg[1][i]!=k+2) {	// has overlap of heterozygote hap fragment
					j = seg[1][i];
					for (l=0; l<nMarkers; l++) {	// check consistence of overlaped seg[1] and merged into one seg[1]ment
						if (seg[1][l]==j) {
							if (ogt[k][l].a[0]>0 && gt[1][l].a[phase]!=ogt[k][l].a[0])	// inconsistence overlap or crossover happened
								return false;
							if (phase == 1)
								gt[1][l].swap();
							seg[1][l] = k+2;
						}
					}
				} else {
					if (phase == 1)
						gt[1][i].swap();
					seg[1][i] = k+2;
				}
			}
		}
	}
	
	//sorting seg[1] numbers, start from 1001 - an arbitary big number
	for (p=1001,i=0; i<nMarkers; i++)
		if (seg[1][i]>1 && seg[1][i]<1000) {
			for (l=seg[1][i],j=i; j<nMarkers; j++)
				if (seg[1][j]==l) seg[1][j] = p;
			p++;
		}
	
	// bring seg[1] number back starting at 2 for #>1000 (2=1001, 3=1002, etc)
	for (i=0; i<nMarkers; i++)
		if (seg[1][i]>1000)
			seg[1][i] -= 999;
	
	return true;
}

Haplotype_List* Mating_mm_genotype::enumberate_parent_haplotype(bool father) {

	const int kMaxSeg = 15;	
	
	Haplotype_List *hl, *old_hl=NULL;
	int i, j, k, l, maxseg_id, nseg;
	
	if (father && xlink) { // gt[0] already phased
		hl = new Haplotype_List(nMarkers);
		for (i=0; i<nMarkers; i++) {
			hl->h[0][i] = gt[0][i].a[0];
			hl->h[1][i] = 0;
		}
		return hl;
	}
	
	if (nMarkers==1) {	// no need to phase
		hl = new Haplotype_List(nMarkers);
		hl->h[0][0] = gt[father?0:1][0].a[0];
		hl->h[1][0] = gt[father?0:1][0].a[1];
		return hl;
	}
	
	Genotype *genotype = father?gt[0]:gt[1];
	int *segment = father?seg[0]:seg[1];
	
	for (nseg=0,maxseg_id=1,i=0; i<nMarkers; i++)
		if (segment[i]==0)
			nseg++;
		else if (segment[i]>maxseg_id)
			maxseg_id = segment[i];		// max(seg_id,1)
	nseg += maxseg_id - 1;
	
	if (nseg>kMaxSeg)
		return NULL;
	
	for (i=0; i<nMarkers && segment[i]==1; i++)
		;
	int first_het_seg = i;	// fix the phase of first heterzygotic segment
	int segid_off = (i<nMarkers && segment[i]>0)? 3 : 2;
	
	// phase for nseg-1 bit:
	// 0 to maxseq_id-segid_off bit: segment[i-segid_off]>1
	// maxseq_id-segid_off+1 to nseg-1 bit: segment[i>first_het_seq]=0
	if (nseg>1) {
		for (k=0; k<(1<<(nseg-1)); k++) {
			hl = new Haplotype_List(nMarkers);
			for (i=0; i<first_het_seg; i++) {
				hl->h[0][i] = genotype[i].a[0];
				hl->h[1][i] = genotype[i].a[1];
			}
			for (j=maxseg_id-segid_off+1; i<nMarkers; i++) {
				if (segment[i]==1 || (segment[i]>0 && segment[first_het_seg]==segment[i]))
					l = 0;
				else if (segment[i]==0)
					l = (k&(1<<(j++)))? 1 : 0;
				else
					l = (k&(1<<(segment[i]-segid_off)))? 1 : 0;
				hl->h[0][i] = genotype[i].a[l];
				hl->h[1][i] = genotype[i].a[1-l];
			}
			hl->next = old_hl;
			old_hl = hl;
		}
	} else {	// just copy genotype
		old_hl = new Haplotype_List(nMarkers);
		for (i=0; i<nMarkers; i++)
			for (j=0; j<2; j++)
				old_hl->h[j][i] = genotype[i].a[j];
	}
	
	return old_hl;

}

bool Mating_mm_genotype::compatible_mating_haplotype(Haplotype_List *h1, Haplotype_List *h2, Genotype *ogt, int sex) {
	
	if (h1==NULL || h2==NULL || ogt==NULL)
		return false;
	
	int i, j, k;
	
	if (xlink && sex==1) {
		for (j=0; j<2; j++) {
			for (k=0; k<nMarkers && (ogt[k].a[0]==0 || ogt[k].a[0]==h2->h[j][k]); k++)
				;
			if (k==nMarkers) return true;
		}
	} else {
		for (i=0; i<(xlink?1:2); i++)
			for (j=0; j<2; j++) {
				for (k=0; k<nMarkers && (ogt[k].a[0]==0 || ogt[k].equal(h1->h[i][k],h2->h[j][k])); k++)
					;
				if (k==nMarkers)
					return true;
			}
	}
	
	return false;

}

// derive a list of compatible mating haplotypes from the compatible phase genotypes
Mating_Haplotype* Mating_mm_genotype::enumerate_mating_haplo(int noffs, Genotype *ogt[], int *sex, int maxmh) {
	
	Haplotype_List *hl[2], *h1, *h2;
	
	hl[0] = enumberate_parent_haplotype(true);
	
	int i;
	bool identical_pgt = false;
	if (!xlink) {
		for (i=0; i<nMarkers && gt[0][i].equal(gt[1][i]); i++)
			;
		if (i==nMarkers)
			identical_pgt = true;
	}
	// when parental genotypes are identical, their haplotype segments and permutation should be the same
	// since the algorithm can't distinguish which parent. So the mating type permutation can be permuted
	// from h1
	hl[1] = identical_pgt? hl[0] : enumberate_parent_haplotype(false);
	
	Mating_Haplotype *mh = NULL;
	int m = 0;	
	for (h1=hl[0]; h1!=NULL && (maxmh==0 || m<maxmh); h1=h1->next) {
		for (h2=(identical_pgt?h1:hl[1]); h2!=NULL && (maxmh==0 || m<maxmh); h2=h2->next) {
			for (i=0; i<noffs && compatible_mating_haplotype(h1,h2,ogt[i],sex[i]); i++)
				;
			if (i==noffs) {	// compatible mating type
				mh = new Mating_Haplotype(nMarkers,h1,h2,mh);
				m++;
			}
		}
	}

	for (i=0; i<(identical_pgt?1:2); i++)
		if (hl[i])
			delete hl[i];
			
	return mh;

}

Mating_mm_genotype::~Mating_mm_genotype() {
	
	for (int i=0; i<2; i++) {
		if (seg[i])
			delete[] seg[i];
		if (gt[i])
			delete[] gt[i];
	}
	
	if (next) 
		delete next;
}

Sufficient_Stat::Sufficient_Stat(NUCFAMILY *fam, int nmrk, int idx[], bool sexlink, int maxcmht) {
	
	xlink = sexlink;
	mhap = NULL;

	int i;
	n_common_ogt = -1;
	for (i=0; i<4; i++) {
		common_ogt[i] = NULL;
		common_oht[i] = NULL;
	}
		
	ogpat = NULL;
	
	nMarkers = (nmrk>0 && idx!=NULL)? nmrk: 0;
	if (nMarkers==0 || fam==NULL)
		return;
	
	for (i=0; i<kMaxMembers; i++)
		otr[i] = NULL;
		
	for (i=0; i<kMaxMembers; i++)
		ogt[i] = new Genotype[nMarkers];
	for (i=0; i<2; i++)
		pgt[i] = new Genotype[nMarkers];
	
	//!/	
	fa_bool=0; //!/
	mo_bool=0; //!/
	if (fam->fa() && fam->fa()->gt) //!/
	{ //!/
		for (i=0; i<nMarkers; i++) //!/
		{ //!/
			pgt[0][i].set(fam->fa()->gt[idx[i]]); //!/
			if(pgt[0][i].defined(xlink)){fa_bool=1;} //!/
		} //!/
	} //!/
	//!/
	if (fam->mo() && fam->mo()->gt) //!/
	{ //!/
		for (i=0; i<nMarkers; i++) //!/
		{ //!/
			pgt[1][i].set(fam->mo()->gt[idx[i]]); //!/
			if(pgt[1][i].defined(xlink)){mo_bool=1;} //!/
		} //!/
	} //!/
    //!/
	
	n_offspring = 0;
	PERSONLIST *per;
	bool good;
	for (per=fam->sibs(); per!=NULL; per=per->next) {
		good = false;
		for (i=0; i<nMarkers; i++) {
			if (per->node->gt[idx[i]].defined(xlink && per->node->sex==1)) {
				ogt[n_offspring][i].set(per->node->gt[idx[i]]);
				good = true;
			} else
				ogt[n_offspring][i].set(0,0);
		}
		if (good) {	// at least one genotype observed
			otr[n_offspring] = per->node->pt;
			sex[n_offspring++] = per->node->sex;
		}
	}
			
	maxcmh = maxcmht;
	mhap = compatible_mating_haplotypes();
	
	EGX = NULL;
	GX = NULL;
	
}

Sufficient_Stat::Sufficient_Stat(Genotype *fa_gt, Genotype *mo_gt,  Genotype *off_gt[], int sx[], int nmrk, int n_off, bool sexlink, int maxcmht) {

	xlink = sexlink;
	nMarkers = nmrk;
	mhap = NULL;

	int i, k;

	n_common_ogt = 0;
	for (i=0; i<4; i++) {
		common_ogt[i] = NULL;
		common_oht[i] = NULL;
	}
		
	ogpat = NULL;
	
	if (nMarkers==0 || n_off==0)
		return;
		
	for (i=0; i<kMaxMembers; i++)
		otr[i] = NULL;
		
	for (i=0; i<kMaxMembers; i++)
		ogt[i] = new Genotype[nMarkers];
	for (i=0; i<2; i++)
		pgt[i] = new Genotype[nMarkers];
    //!/	
	fa_bool=0; //!/	
	mo_bool=0; //!/	
	if (fa_gt!=NULL) //!/	
	{ //!/	
		for (i=0; i<nMarkers; i++) //!/	
		{ //!/	
			pgt[0][i].set(fa_gt[i]); //!/	
			if(pgt[0][i].defined(xlink)){fa_bool=1;} //!/	
		} //!/	
	} //!/	
	
	if (mo_gt!=NULL) //!/	
	{ //!/	
		for (i=0; i<nMarkers; i++) //!/	
		{ //!/	
			pgt[1][i].set(mo_gt[i]); //!/	
			if(pgt[1][i].defined(xlink)){mo_bool=1;} //!/	
		} //!/	
	} //!/	
    //!/

	bool good;
	n_offspring = 0;
	for (k=0; k<n_off; k++) {
		good = false;
		for (i=0; i<nMarkers && off_gt[k][i].defined((sexlink && sx[k]==1)); i++) {
			ogt[n_offspring][i].set(off_gt[k][i]);
			good = true;
		}
		if (good)	// at least one genotype observed
			sex[n_offspring++] = sx[k];
	}
	
	maxcmh = maxcmht;
	mhap = compatible_mating_haplotypes();
	
	EGX = NULL;
	GX = NULL;

}

Mating_Haplotype* Sufficient_Stat::compatible_mating_haplotypes() {
	
	if (nMarkers==0) return NULL;
		
	Compatible_parent_genotype **cpg = new Compatible_parent_genotype*[nMarkers];
	int i, j;
	for (i=0; i<nMarkers; i++)
		cpg[i] = enumerate_compatible_mating_genotype(i);
	
	if (nMarkers==1) {	// no need to phase, add cpg directly to mhap
		Mating_Haplotype *cmh=NULL;
		for (i=0; i<cpg[0]->ngt; i++) {
			cmh = new Mating_Haplotype(nMarkers, NULL, NULL, cmh);
			cmh->h[0][0] = cpg[0]->gt[i][0].a[0];
			cmh->h[1][0] = xlink?0:cpg[0]->gt[i][0].a[1];
			cmh->h[2][0] = cpg[0]->gt[i][1].a[0];
			cmh->h[3][0] = cpg[0]->gt[i][1].a[1];
			cmh->hash_val = cmh->hash();
		}
		for (i=0; i<nMarkers; i++)
			if (cpg[i])
				delete cpg[i];
		if (cpg)
			delete[] cpg;
			
		return cmh;
	}
	unsigned int iteration_counter=0; //!/
	// enumerate all possible cmg
	Mating_Haplotype *cmh=NULL, *cmh_last_node;
	int ncmh = 0; // added 10/20/04 to discard family with ncmh>maxcmh
	for (i=0; i<nMarkers && cpg[i]->ngt>0; i++)
		;
	if (i==nMarkers) {
		// iterate all recombinations of cmg[i] at all i, and extract haplotype information
		Mating_mm_genotype *mh;
		Genotype *pg[2];
		int *idx = new int[nMarkers];
		for (i=0; i<nMarkers; i++)
			idx[i] = 0;
		while (idx[0]<cpg[0]->ngt) {
			for (i=0; i<2; i++) {
				pg[i] = new Genotype[nMarkers];
				for (j=0; j<nMarkers; j++)
					pg[i][j].set(cpg[j]->gt[idx[j]][i]);
			}
			mh = new Mating_mm_genotype(nMarkers, pg, xlink);
			if (xlink?mh->resolve_x(n_offspring, ogt, sex):mh->resolve(n_offspring, ogt)) {
				if (cmh==NULL) {
					cmh = mh->enumerate_mating_haplo(n_offspring, ogt, sex, 2*maxcmh);
					if (cmh!=NULL)
						for (ncmh=1,cmh_last_node=cmh; cmh_last_node->next!=NULL; cmh_last_node=cmh_last_node->next)
							ncmh++;
				} else {
					cmh_last_node->next = mh->enumerate_mating_haplo(n_offspring, ogt, sex, maxcmh);
					while (cmh_last_node->next!=NULL) {
						cmh_last_node = cmh_last_node->next;
						ncmh++;
					}
				}
			}
			if (mh!=NULL)
				delete mh;
			if (maxcmh!=0 && ncmh>2*maxcmh)
				break;
			for(j=nMarkers-1; j>=0 && ++idx[j]==cpg[j]->ngt; j--)
				;
			if (j<0)
				break;
			while (++j<nMarkers)
				idx[j] = 0;
			iteration_counter++; //!/
			if(iteration_counter==1000000 && ncmh==0) break; //!/
		}
		
		if (idx)
			delete idx;
	

	}
	
	// remove duplicates:
	if (!xlink) {
		Mating_Haplotype *cmh1, *cmh2;
		for (cmh1=cmh; cmh1!=NULL && cmh1->next!=NULL; cmh1=cmh1->next)
			for (cmh2=cmh1; cmh2!=NULL && cmh2->next!=NULL; cmh2=cmh2->next)
				if (cmh2->next->equal(cmh1)) {
					cmh_last_node = cmh2->next;
					cmh2->next = cmh_last_node->next;
					cmh_last_node->next = NULL;
					delete cmh_last_node;
				}
	}
	
	for (i=0; i<nMarkers; i++)
		if (cpg[i])
			delete cpg[i];
	if (cpg)
		delete[] cpg;
	
	return cmh;

}

int Sufficient_Stat::get_common_ogt() {
	
	bool good[4];
	Mating_Haplotype *cmh;
	
	if (mhap==NULL)
		return 0;
	
	int i, j, k, l;
	
	for (i=0; i<4; i++)
		common_ogt[i] = new Genotype[nMarkers];
		
	if (!xlink) {
		for (l=0,i=0; i<2; i++)
			for (j=2; j<4; j++) {
				for (k=0; k<nMarkers; k++)
					common_ogt[l][k].set(mhap->h[i][k],mhap->h[j][k]);
				common_ogt_sex[l] = 2; // no matter
				l++;
			}
	} else {
		for (i=0; i<2; i++) {	// two sons
			for (k=0; k<nMarkers; k++)
				common_ogt[i][k].set(mhap->h[2+i][k],0);
			common_ogt_sex[i] = 1;
		}
		for (i=2; i<4; i++) {	// two daughters
			for (k=0; k<nMarkers; k++)
				common_ogt[i][k].set(mhap->h[0][k],mhap->h[i][k]);
			common_ogt_sex[i] = 2;
		}
	}
	
	for (i=0; i<4; i++) 
		good[i] = (common_ogt[i][0].a[0]>0 && common_ogt[i][0].a[0]!=kotherallele && ((xlink && common_ogt_sex[i]==1) || common_ogt[i][0].a[1]!=kotherallele));
	for (i=0; i<3; i++)	{	// check for duplicate
		for (j=i+1; j<4; j++) {
			if (common_ogt_sex[i]==common_ogt_sex[j] && good[i] && good[j]) {
				for (k=0; k<nMarkers && common_ogt[i][k].equal(common_ogt[j][k]); k++)
					;
				if (k==nMarkers)
					good[j] = false;
			}
		}
	}
	
	for (cmh=mhap->next; cmh!=NULL; cmh=cmh->next)
		for (i=0; i<4; i++)
			if (good[i] && !cmh->gt_compatible(common_ogt[i], common_ogt_sex[i], xlink))
				good[i] = false;
	
	// delete bad (duplicated)
	for (i=0; i<4; i++)
		if (!good[i]) {
			delete[] common_ogt[i];
			common_ogt[i] = NULL;
		}
	n_common_ogt = 0;
	for (i=0; i<4; i++)
		if (good[i]) {
			common_ogt_sex[n_common_ogt] = common_ogt_sex[i];
			common_ogt[n_common_ogt] = common_ogt[i];
			n_common_ogt++;
		}
	for (i=n_common_ogt; i<4; i++)
		common_ogt[i] = NULL;
	
	// update mhap->pg[] & mhap->pset
	if (n_common_ogt>0) {
		int ncnt[4] = {0,0,0,0};	// # of sibs with common_ogt[]
		for (i=0; i<n_offspring; i++)
			if ((j=which_comm_ogt(ogt[i],sex[i]))>=0)
				ncnt[j]++;
		for (cmh=mhap; cmh!=NULL; cmh=cmh->next) {
			cmh->update_pg(common_ogt, common_ogt_sex, n_common_ogt, xlink);
			cmh->get_p_set(ncnt, true);
		}
	}
	
	return n_common_ogt;

}

int Sufficient_Stat::enumerate_genotype_from_allele(int idx, Genotype *ppgt) {
	
	int i, j, k;
	int npa = 0;
	char pa[4], a;
	int missing_indicator=0; //!/
	
	//parents alleles
	for (k=0; k<2; k++) { //!/ modification 8/7/2024, added curly brackets around if-statement
		if (pgt[k][idx].defined(xlink && k==0)){
			for (i=0; i<(xlink && k==0?1:2); i++) {
				a = pgt[k][idx].a[i];
				for(j=0; j<npa && pa[j]!=a; j++)
					;
				if (j==npa) // new allele
					pa[npa++] = a;
			}
		}
		else if ((k==0 && fa_bool==1) || (k==1 && mo_bool==1)) //!/
		{ //!/
				missing_indicator=1; //!/
		} //!/
	}
	
	// offs alleles
	for (k=0; k<n_offspring; k++) { //!/ modification 8/7/2024, added curly brackets around if-statement
		if (ogt[k][idx].defined(xlink && sex[k]==1)){
			
			for (i=0; i<(xlink && sex[k]==1?1:2); i++) {
				a = ogt[k][idx].a[i];
				for(j=0; j<npa && pa[j]!=a; j++)
					;
				if (j==npa) // new allele
					pa[npa++] = a;
			}
		}
		else{ //!/
				missing_indicator=1; //!/
		} //!/
	}
	
	if (npa == 0) return 0;
	
	if (npa==1 && missing_indicator==1)	// added an other_allele for monomorphic marker //!/
		pa[npa++] = kotherallele;
		
	int npg = 0;
	for (i=0; i<npa; i++)
		for (j=i; j<npa; j++)
			ppgt[npg++].set(pa[i], pa[j]);

	return npg;

}

Compatible_parent_genotype* Sufficient_Stat::enumerate_compatible_mating_genotype(int idx) {
	
	int npar = (pgt[0][idx].defined(xlink)?1:0) + (pgt[1][idx].defined(false)?1:0);
	
	Compatible_parent_genotype *cpg = new Compatible_parent_genotype;
	if (npar == 2) {
		cpg->gt[0][0].set(pgt[0][idx]);
		cpg->gt[0][1].set(pgt[1][idx]);
		cpg->ngt = 1;
		return cpg;
	}

	cpg->ngt = 0;
	
	int npa; // possible # of allele for each parent
	char pa[4];
	int npg[2];	// possible # of genotype for each parent, 0=fa, 1=mo
	Genotype ppgt[2][10];

	int i, j, k;
	if (xlink) {
		// determine possible fa genotype
		if (pgt[0][idx].defined(true)) {
			npg[0] = 1;
			ppgt[0][0].set(pgt[0][idx].a[0],0);
		} else {
			// fa allele space will be the allele shared by daughters
			npa = 0;
			for (k=0; k<n_offspring; k++)
				if (sex[k]==2 && ogt[k][idx].defined(false)) {
					if (npa==0) {
						for (i=0; i<2; i++)
							pa[i] = ogt[k][idx].a[i];
						npa = (pa[0]==pa[1]?1:2);
					} else {
						for (i=0; i<npa; i++)
							if (pa[i]>0 && pa[i]!=ogt[k][idx].a[0] && pa[i]!=ogt[k][idx].a[1])	// not the shared allele
								pa[i] = 0;
					}
				}
			npg[0] = 0;
			if (npa == 0)	// no daughter genotypes observed
				ppgt[0][npg[0]++].set(kotherallele, 0);
			else {
				for(i=0; i<npa; i++)
					if (pa[i]>0)
						ppgt[0][npg[0]++].set(pa[i], 0);
				if (npg[0]==0)
					throw (ERROR("mendelian error, this shouldn't happen"));
			}
		}
		
		// determine possible mo genotype
		if (pgt[1][idx].defined(false)) {
			npg[1] = 1;
			ppgt[1][0].set(pgt[1][idx]);
		} else 
			npg[1] = enumerate_genotype_from_allele(idx, &ppgt[1][0]);
		
		cpg->ngt = 0;
		for (i=0; i<npg[0]; i++)
			for (j=0; j<npg[1]; j++)
				if (is_compatible_mating(ppgt[0][i], ppgt[1][j], idx)) {
					cpg->gt[cpg->ngt][0].set(ppgt[0][i]);
					cpg->gt[cpg->ngt][1].set(ppgt[1][j]);
					cpg->ngt++;
				}
	} else { // for non-X-linked markers, do not distinguish fa and mo, use mo instead
		cpg->ngt = 0;
		if (npar == 1) {
			j = (pgt[0][idx].defined(false)? 0 : 1);	// pgt[j] genotyped
			npg[j] = 1;
			ppgt[j][0].set(pgt[j][idx]);
			npg[1-j] = enumerate_genotype_from_allele(idx, &ppgt[1-j][0]);
			for (i=0; i<npg[1-j]; i++)
				if (is_compatible_mating(ppgt[j][0], ppgt[1-j][i], idx)) {
					cpg->gt[cpg->ngt][j].set(ppgt[j][0]);
					cpg->gt[cpg->ngt][1-j].set(ppgt[1-j][i]);
					cpg->ngt++;
				}
		} else {
			npg[1] = enumerate_genotype_from_allele(idx, &ppgt[1][0]);
			for (i=0; i<npg[1]; i++)
				for (j=i; j<npg[1]; j++)
					if (is_compatible_mating(ppgt[1][i], ppgt[1][j], idx)) {
						cpg->gt[cpg->ngt][0].set(ppgt[1][i]);
						cpg->gt[cpg->ngt][1].set(ppgt[1][j]);
						cpg->ngt++;
						// added to distinguish fa from mo
						if (i!=j) {
							cpg->gt[cpg->ngt][1].set(ppgt[1][i]);
							cpg->gt[cpg->ngt][0].set(ppgt[1][j]);
							cpg->ngt++;
						}
					}
		}
	}
	
	return cpg;
	
}

bool Sufficient_Stat::is_compatible_mating(Genotype &fg, Genotype &mg, int idx) {
	
	int k;
	for (k=0; k<n_offspring; k++)
		if (ogt[k][idx].defined(xlink && sex[k]==1)) {
			if (xlink) {
				if (sex[k]==1?(ogt[k][idx].a[0]!=mg.a[0] && ogt[k][idx].a[0]!=mg.a[1]):(!ogt[k][idx].equal(fg.a[0], mg.a[0]) && !ogt[k][idx].equal(fg.a[0], mg.a[1])))
					return false;
			} 
			else if (!ogt[k][idx].equal(fg.a[0], mg.a[0]) && !ogt[k][idx].equal(fg.a[1], mg.a[0]) && !ogt[k][idx].equal(fg.a[0], mg.a[1]) && !ogt[k][idx].equal(fg.a[1], mg.a[1]))
				return false;
		}
	
	return true;
}

int Sufficient_Stat::mhap_hash(int &hash_code) {
	
	hash_code = 0;
	int len = 0;
	for (Mating_Haplotype *mh=mhap; mh!=NULL; mh=mh->next) {
		len++;
		hash_code ^= mh->hash_val;
	}
	
	return len;
}

// derive ogpat, a list of offspring genotype pattern that their mhap is the same as the observed offspring genotype and the ratio of frequency in each mhap is a constant
void Sufficient_Stat::compatible_offspring_gt_set() {
	
	// genotype set of n=count_infomative_offspring() offspring can be enumerated by combination:
	// C(n_comm_ogt, 1) + ... + C(n_comm_ogt, minimal(n,n_comm_ogt))
	// maximum will contain 2^4-1 = 15 genotype combinations
	
	ogpat = NULL;
	if (n_common_ogt<2) return;	// no offs variation
	
	int i, j, k;	
	int n = 0;
	for (k=0; k<n_offspring; k++)
		if ((j=which_comm_ogt(ogt[k],sex[k])) >= 0)
			n++;
	if (n>n_common_ogt)
		n = n_common_ogt;
		
	if (n==0)
		return;
	
	int len1, len2, hash1, hash2;
	len1 = mhap_hash(hash1);	
	
	int b[4];	// indicator for presence of common_ogt[b[]]
	int cogt_flag = 0; // indicator for presence of common_ogt[] in the observed data
	int gsex[4];
	Genotype *g[4];
	for (k=0,j=0,i=0; i<n_offspring; i++)
		if ((j=which_comm_ogt(ogt[i],sex[i]))>=0 && (cogt_flag & (1<<j))==0) {
			cogt_flag |= (1<<j);
			b[k++] = j;
		}
	permute_off_set_distrib(b, k);
			
	Sufficient_Stat *ss = NULL;
	Mating_Haplotype *mh1, *mh2;
	
	for (i=1; i<(1<<n_common_ogt); i++)	// i = cgt_flag
		if (i!=cogt_flag) {
			for (k=0,j=0; j<n_common_ogt; j++)
				if ((i>>j)&1) {
					g[k] = common_ogt[j];
					gsex[k] = common_ogt_sex[j];
					b[k] = j;
					k++;
				}
			if (k<=n) {	// valid set
				ss = new Sufficient_Stat(pgt[0], pgt[1], g, gsex, nMarkers, k, xlink, maxcmh);
				if (ss->mhap!=NULL) {
					len2 = ss->mhap_hash(hash2);
					if (len1==len2 && hash1==hash2)	{ // two set may be eaqual?
						// now have to use brutal force to check if they are equal
						for (mh1=mhap; mh1!=NULL; mh1=mh1->next) {
							for (mh2=ss->mhap; mh2!=NULL && (mh2->hash_val!=mh1->hash_val || !mh1->equal(mh2)); mh2=mh2->next)
								;
							if (mh2==NULL)
								break;
						}
						if (mh1==NULL)	// two sets are equal; permute and add this set
							permute_off_set_distrib(b, k);
					}
				}
				if (ss!=NULL)
					delete ss;
			}
		}
	
	// normalized ogpat->pg so it sums to 1
	if (ogpat!=NULL) {
		Off_Genotype_Pattern *pat;
		double sum = 0.0;	
		for (pat=ogpat; pat!=NULL; pat=pat->next) {
			pat->pg *= (xlink?pat->permut_n_x(common_ogt_sex):pat->permut_n());
			sum += pat->pg;
		}
		if (sum>0.0)
			for (pat=ogpat; pat!=NULL; pat=pat->next)
				pat->pg /= sum;
	}

}

int Sufficient_Stat::which_comm_ogt(Genotype *g, int sx) {
	
	int i, j, k;
	
	for (k=-1,i=0; i<n_common_ogt; i++) 
		if (!xlink || sx==common_ogt_sex[i]) {
			for (j=0; j<nMarkers && (!g[j].defined(xlink && sx==1) || common_ogt[i][j].equal(g[j].a[0],(xlink && sx==1?0:g[j].a[1]))); j++)
				;
			if (j==nMarkers) {	// match
				if (k>=0)
					k = -2;
				else if (k==-1)
					k = i;
			}
		}
	
	return k;
}

// permute offspring gt set: list L_2
void Sufficient_Stat::permute_off_set_distrib(int gid[], int ng) {
	
	int i, j, k, n, n_male=0;

	for (n=0,k=0; k<n_offspring; k++)
		if ((j=which_comm_ogt(ogt[k],sex[k])) >= 0) {
			n++;
			if (sex[k]==1)
				n_male++;	// number of sons
		}
	if (n<ng || ng==0) return;
	
	int cnt[4], gcnt[4];
	int totalcnt[4];
	for (i=0; i<4; i++) {
		gcnt[i] = 1;
		totalcnt[i]=i+1;
	}
	
	while (gcnt[0]<=n-ng+1) {
		gcnt[ng-1] = (ng==1?n:n-totalcnt[ng-2]);
		for (i=0; i<4; i++)
			cnt[i] = 0;
		for (i=0; i<ng; i++)
			cnt[gid[i]] = gcnt[i];
			
		// cnt[] is a valid permutation
		// added 1/15/07 to make sure the number of male and female in  offspring is fixed if xlink==true
		if (xlink) {
			for (i=0,j=0; i<n_common_ogt; i++)
				if (common_ogt_sex[i]==1)
					j += cnt[i];
			if (j==n_male)	// condition met
				check_and_add_offpat(cnt);
		} else
			check_and_add_offpat(cnt);
			
		for (j=ng-2; j>=0 && gcnt[j]==n-(j==0?0:totalcnt[j-1])-ng+1+j; j--)
			;
		if (j<0) break;
		gcnt[j]++;
		totalcnt[j]++;
		while (++j<ng-1) {
			gcnt[j] = 1;
			totalcnt[j] = totalcnt[j-1]+1;
		}
	}
	
}

void Sufficient_Stat::check_and_add_offpat(int cnt[]) {
	
	if (mhap==NULL)
		return;

	double p = mhap->get_p_set(cnt);
	double r = p/mhap->pset;	// p ratio between cnt[] and observed and need to be  constant across all mh
	
	// check if A_kl is a constant
	Mating_Haplotype *mh;
	for (mh=mhap->next; mh!=NULL && fabs(r-mh->get_p_set(cnt)/mh->pset)<0.00001; mh=mh->next) 
		;
	
	if (mh==NULL) // a valid pattern
		ogpat = new Off_Genotype_Pattern(cnt, p, ogpat);
}

// possible haplotypes of common_ogt
void Sufficient_Stat::get_common_oht() {

	int i, j, k, l;
	
	Mating_Haplotype *mh;
	
	for (i=0; i<n_common_ogt; i++) {
		if (xlink && common_ogt_sex[i]==1) {	// only 1 possible phase
			common_oht[i] = new Haplotype_List(nMarkers,mhap->h[2],mhap->h[0],NULL);
			for (j=0; j<nMarkers; j++) {
				common_oht[i]->h[0][j] = common_ogt[i][j].a[0];
				common_oht[i]->h[1][j] = 0;
			}
		} else {
			for (mh=mhap; mh!=NULL; mh=mh->next) {
				for (j=0; j<(xlink?1:2); j++)
					for (k=2; k<4; k++) {
						for (l=0; l<nMarkers && common_ogt[i][l].equal(mh->h[j][l],mh->h[k][l]); l++)
							;
						if (l==nMarkers && (common_oht[i]==NULL || common_oht[i]->find(mh->h[j],mh->h[k])==NULL))	// gt match and novel hap in common_oht[i]
							common_oht[i] = new Haplotype_List(nMarkers,mh->h[j],mh->h[k],common_oht[i]);
					}
			}
		}
	}

}

// derive the possibility of each haplotype type of a genotype given mating type expections
void Sufficient_Stat::update_comm_oht_freq(Haplotype_List *emlist) {
	
	Haplotype_List *hl, *hl1, *hl2;	
	Mating_Haplotype *mh;
	
	double fsum = 0;
	// update prob(mh) given emlist
	for (mh=mhap; mh!=NULL; mh=mh->next) {
		hl1 = emlist->find(mh->h[0],mh->h[1]);	// NULL if not in emlist
		hl2 = emlist->find(mh->h[2],mh->h[3]);
		mh->frq = (hl1==hl2?1:2)*(hl1?hl1->p:very_small_number)*(hl2?hl2->p:very_small_number)*mh->pset;	// give a small prob if not in emlist
		fsum += mh->frq;
	}
	for (mh=mhap; mh!=NULL; mh=mh->next)
		mh->frq /= fsum;
	
	for (int i=0; i<n_common_ogt; i++) {
		fsum = 0.0;
		for (hl=common_oht[i]; hl!=NULL; hl=hl->next) {
			hl->p = 0;
			for (mh=mhap; mh!=NULL; mh=mh->next)
				hl->p += mh->offspring_ht_prob(hl)*mh->frq;
			fsum += hl->p;
		}
		if (fsum==0.0)
			throw ERROR("internal error: update_comm_oht_freq() sum = 0");
			
		for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
			hl->p /= fsum;
	}
		
}

// calculate Var(X_h) where X_h=sum_hgt(P_hgt*X_hgt)
// here pg=p[g1=g] and pgg=p[g1=g,g2=g']
// and pg[4] and pgg[4][4] index from 0...n_comm_ogt-1
// pg and pgg will use the sum of og->pg and disregard of the comm_oht[] for each common_ogt
// comm_oht[] will be used later in the calculation of X[][hap] for each haplotype in htab
void Sufficient_Stat::pgpgg2() {
	
	int i, j;
	
	for (i=0; i<4; i++) {
		pg[i] = 0.0;
		for (j=0; j<4; j++)
			pgg[i][j] = 0.0;
	}
	
	if (ogpat==NULL)
		return;
	
	int noff = 0;	// # of informative offs
	for (i=0; i<n_common_ogt; i++)
			noff += ogpat->cnt[i];
	
	Off_Genotype_Pattern *og;
	for (og=ogpat; og!=NULL; og=og->next) {
		for (i=0; i<n_common_ogt; i++) {
			if (og->cnt[i]>0) {
				pg[i] += og->pg*og->cnt[i]/noff;
				if (og->cnt[i]>1)
					pgg[i][i] += og->pg*og->cnt[i]*(og->cnt[i]-1)/(noff*(noff-1));
				for (j=i+1; j<n_common_ogt; j++)
					if (og->cnt[j]>0) 
						pgg[i][j] += og->pg*og->cnt[i]*og->cnt[j]/(noff*(noff-1));
			}
		}	
	}
	
	for (i=0; i<n_common_ogt; i++)
		for (j=i+1; j<n_common_ogt; j++)
			pgg[j][i] = pgg[i][j];
				
}

void Sufficient_Stat::pgpgg_x() {
	
	int i, j, k;
	
	for (k=0; k<2; k++)
		for (i=0; i<4; i++)
			pg_x[k][i] = 0.0;
	for (k=0; k<3; k++)
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
				pgg_x[k][i][j] = 0.0;
	
	if (ogpat==NULL || !xlink)
		return;
	
	int noff[2] = {0,0};	// # of informative male & female offs
	for (i=0; i<n_common_ogt; i++)
			noff[(common_ogt_sex[i]==1?0:1)] += ogpat->cnt[i];
	
	Off_Genotype_Pattern *og;
	for (og=ogpat; og!=NULL; og=og->next) {
		for (i=0; i<n_common_ogt; i++) {
			if (og->cnt[i]>0) {
				k = (common_ogt_sex[i]==1?0:1);
				pg_x[k][i] += og->pg*og->cnt[i]/noff[k];
				if (og->cnt[i]>1)
					pgg_x[k][i][i] += og->pg*og->cnt[i]*(og->cnt[i]-1)/(noff[k]*(noff[k]-1));
				for (j=i+1; j<n_common_ogt; j++)
					if (og->cnt[j]>0) {
						if (common_ogt_sex[i]==common_ogt_sex[j])
							pgg_x[k][i][j] += og->pg*og->cnt[i]*og->cnt[j]/(noff[k]*(noff[k]-1));
						else
							pgg_x[2][i][j] += og->pg*og->cnt[i]*og->cnt[j]/(2*noff[0]*noff[1]);
					}
			}
		}	
	}
	
	for (k=0; k<3; k++)
		for (i=0; i<n_common_ogt; i++)
			for (j=i+1; j<n_common_ogt; j++)
				pgg_x[k][j][i] = pgg_x[k][i][j];
				
}

void Sufficient_Stat::output_info (mystream &os, HaploidTable *htab, char *fam_id, bool informative_only) {
	
	int i, j, k, l, m;
	
	if (informative_only)
		if (!informative())
			return;

	MYOUTPUT(os, fam_id << "\n")
	MYOUTPUT(os, "observed genotype configuration\n")
	for (k=0; k<2; k++)
		for (j=0; j<(xlink && k==0?1:2); j++) {
			if (j==0)
				MYOUTPUT(os, setw(15) << (k==0?"father":"mother") << "= ")
			else
				MYOUTPUT(os, setw(17) << " ")
			for (i=0; i<nMarkers; i++)
				MYOUTPUT(os, setw(4) << (pgt[k][i].a[j]==0?0:htab->loc_ptr[i]->alleleName[pgt[k][i].a[j]]))
			MYOUTPUT(os, (xlink && k==0 && j==0?"\n\n":"\n"))
		}
	
	for (k=0; k<n_offspring; k++)
		for (j=0; j<2; j++) {
			if (j==0)
				MYOUTPUT(os, "offspring " << setw(5) << k+1 << "= ")
			else
				MYOUTPUT(os, setw(17) << " ")
			if (j==0 || !xlink || sex[k]==2) {
				for (i=0; i<nMarkers; i++)
					MYOUTPUT(os, setw(4) << (ogt[k][i].a[j]==0?0:htab->loc_ptr[i]->alleleName[ogt[k][i].a[j]]))
			}
			MYOUTPUT(os, "\n")
		}
	
	Mating_Haplotype *mh;
	Haplotype_List *hl;
	Off_Genotype_Pattern *pat;
	if (mhap!=NULL && mhap->len()>maxcmh)
		MYOUTPUT(os, "too many compatible mating types\n")
	else if (mhap==NULL)
		MYOUTPUT(os, "uninformative family or compatible mating types not found (recombination occured).\n")
	else if (n_common_ogt<2)
		MYOUTPUT(os, "Only 1 compatible offspring genotypes observed, no genetic variation.\n")
	else {
		for (k=0,mh=mhap; mh!=NULL; mh=mh->next) {
			k++;
			MYOUTPUT(os, "\ncompatible mating haplotype " << k << "\n")
			for (i=0; i<4; i++) {
				if (!(xlink && i==1)) {
					l = htab->allele(mh->h[i]);
					MYOUTPUT(os,  "h" << setw(4) << l << setw(4) << ":")
					for (j=0; j<nMarkers; j++)
						if (mh->h[i][j]==kotherallele)
							MYOUTPUT(os, setw(4) << "x")
						else
							MYOUTPUT(os, setw(4) << htab->loc_ptr[j]->alleleName[mh->h[i][j]])
					MYOUTPUT(os, "\n")
					if (i==1)
						MYOUTPUT(os, "\n")
				} else
					MYOUTPUT(os, "\n")
			}
		}
		MYOUTPUT(os, "\nThere are " << n_common_ogt << " compatible offspring genotypes:\n")
		for (i=0; i<n_common_ogt; i++) {
			if (nMarkers>1)
				for (m=0, hl=common_oht[i]; hl!=NULL; hl=hl->next) {
					m++;
					MYOUTPUT(os, "\ncompatible offspring genotype " << i+1 << " (g" << i+1 << "), phase " << m << "; EM P=" << hl->p << "\n")
					for (k=0; k<2; k++) {
						if (k==0 || !xlink || common_ogt_sex[i]==2) {
							l = htab->allele(hl->h[k]);
							MYOUTPUT(os,  "h" << setw(4) << l << setw(4) << ":")
							for (j=0; j<nMarkers; j++)
								MYOUTPUT(os,  setw(4) << htab->loc_ptr[j]->alleleName[hl->h[k][j]])
						}
						MYOUTPUT(os, "\n")
					}
				}
			else {
				MYOUTPUT(os, "\ncompatible offspring genotype " << i+1 << " (g" << i+1 << ") = ")
				MYOUTPUT(os, htab->loc_ptr[0]->alleleName[common_ogt[i][0].a[0]] << "/")
				if (!xlink || common_ogt_sex[i]==2)
					MYOUTPUT(os, htab->loc_ptr[0]->alleleName[common_ogt[i][0].a[1]])
				MYOUTPUT(os, "\n");
			}
		}			
		MYOUTPUT(os, "\ndistribution of compatible offspring genotype configurations:\n\n")
		for (i=0; i<n_common_ogt; i++)
			MYOUTPUT(os, "#g" << i+1 << "\t")
		MYOUTPUT(os, "P[G]\n")
		if (ogpat==NULL) {
			if (n_common_ogt==1) {
				for (j=0,i=0; i<n_offspring; i++)
					if (which_comm_ogt(ogt[i],sex[i])==0)
						j++;
				MYOUTPUT(os, j << "\t1.000000\n")
			} else {
				MYOUTPUT(os, "compatible offspring genotype configuration not found\n")
				//MYOUTPUT(os, "\n--------------------------------\n")
			}
		} else {
			for (pat=ogpat; pat!=NULL; pat=pat->next) {
				for (j=0; j<n_common_ogt; j++)
					MYOUTPUT(os, pat->cnt[j] << "\t")
				MYOUTPUT(os, pat->pg << "\n")
			}
		}		
		//MYOUTPUT(os, "\n--------------------------------\n")
		if (!xlink) {
			MYOUTPUT(os, "\npg= \n")
			for (i=0; i<n_common_ogt; i++)
				MYOUTPUT(os, "g" << i+1 << "\t" << pg[i] << "\n")
			MYOUTPUT(os, "\npgg = \n")
			for (i=0; i<n_common_ogt; i++) {
				MYOUTPUT(os, "g" << i+1)
				for (j=0; j<n_common_ogt; j++)
					MYOUTPUT(os, "\t" << pgg[i][j])
				MYOUTPUT(os, "\n")
			}
		} else {
			for (k=0; k<2; k++) {
				MYOUTPUT(os, (k==0?"\npg:male= \n":"\npg:female= \n"))
				for (i=0; i<n_common_ogt; i++)
					MYOUTPUT(os, "g" << i+1 << "\t" << pg_x[k][i] << "\n")
			}
			for (k=0; k<3; k++) {
				MYOUTPUT(os, (k==0?"\npgg:male_male= \n":(k==1?"\npgg:female_female= \n":"\npgg:male_female= \n")))
				for (i=0; i<n_common_ogt; i++) {
					MYOUTPUT(os, "g" << i+1)
					for (j=0; j<n_common_ogt; j++)
						MYOUTPUT(os, "\t" << pgg_x[k][i][j])
					MYOUTPUT(os, "\n")
				}
			}
		}
		//MYOUTPUT(os, "\n--------------------------------\n")
	}

}

void Sufficient_Stat::analyze(Haplotype_List *emlist) {
	
	if (n_common_ogt == -1) // haven't called get_common_ogt() before
		get_common_ogt();
	
	compatible_offspring_gt_set();	// obtain ogpat
	
	if (xlink)
		pgpgg_x();
	else
		pgpgg2();
	
	if (nMarkers>1) {
		get_common_oht();
		if (emlist!=NULL)
			update_comm_oht_freq(emlist);
	}
	
}

bool Sufficient_Stat::fbat_stat(ColumnVector &EX, Matrix &VX, Matrix &CX, double gx[][3], HaploidTable *htab) {
	
	int na = VX.Nrows();
	
	if (htab && htab->size!=VX.Nrows())
		throw Param_Error("Sufficient_Stat::fbat_stat wrong matrix dimension");
	
	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat htab not specified");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL)
		return false;
	
	int i, j;
	int a[2];
	Haplotype_List *hl;

	Matrix X(na, n_common_ogt);
	X = 0;
	for (i=0; i<n_common_ogt; i++) {
		if (nMarkers>1) {
			for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
				if (hl->p>0) {
					for (j=0; j<2; j++)
						a[j] = htab->allele(hl->h[j]);
					if (a[0]>0) {
						X(a[0],i+1) += hl->p*gx[a[0]-1][(a[0]==a[1] && (!xlink || common_ogt_sex[i]==2)?2:1)];
						if (a[1]>0 && a[0]!=a[1] && (!xlink || common_ogt_sex[i]==2))
							X(a[1],i+1) += hl->p*gx[a[1]-1][1];
					}
				}
		} else {
			for (j=0; j<2; j++)
				a[j] = (htab? htab->allele(&common_ogt[i][0].a[j]) : common_ogt[i][0].a[j]);
			if (a[0]>0) {
				X(a[0],i+1) += gx[a[0]-1][(a[0]==a[1] && (!xlink || common_ogt_sex[i]==2)?2:1)];
				if (a[1]>0 && a[0]!=a[1] && (!xlink || common_ogt_sex[i]==2))
					X(a[1],i+1) += gx[a[1]-1][1];
			}
		}
	}
	
	ColumnVector P(n_common_ogt);
	for (i=0; i<n_common_ogt; i++)
		P(i+1) = pg[i];
		
	DiagonalMatrix Q(n_common_ogt);
	for (i=0; i<n_common_ogt; i++)
		Q(i+1,i+1) = P(i+1);
	
	SymmetricMatrix PP(n_common_ogt);
	for (i=0; i<n_common_ogt; i++)
		for (j=i; j<n_common_ogt; j++)
			PP(i+1,j+1) = pgg[i][j];
	
	EX = X * P;
	VX = X * Q * X.t() - EX * EX.t();
	CX = X * PP * X.t() - EX * EX.t();
	
	return true;
	
}

// return # of informative offs
int Sufficient_Stat::fbat_stat_observed(int trt_id, double offset, double &tsum, double &tssq, ColumnVector &S, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab) {
	
	tsum = 0;
	tssq = 0;
	double t;
	int i, j;
	int n = 0;
	
	S = 0;
	X = 0;
	ColumnVector Y = X;
	for (i=0; i<n_offspring; i++) {
		if (otr[i]!=NULL && otr[i][trt_id]<kVeryLargeNumber && (j=which_comm_ogt(ogt[i], sex[i]))>=0 ) {
			t = otr[i][trt_id] - ( (!(flag&flag_censored) || otr[i][0]>0.5)? offset : 0 );	// for cencored data, t-u in affecteds, t in censored group
			tsum += t;
			tssq += t*t;
			comm_ogt_X(Y, j, gx, htab);
			X += Y;
			S += t*Y;			
			n++;
		}
	}
	
	return n;
}

/*
bool Sufficient_Stat::fbat_stat2(int trt_id, double offset, ColumnVector &ES, ColumnVector &S, Matrix &V, double gx[], int flag, HaploidTable *htab) {
	
	int na = V.Nrows();
	
	if (htab && htab->size!=na)
		throw Param_Error("Sufficient_Stat::fbat_stat wrong matrix dimension");
	
	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat htab not specified");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL || trt_id<0)
		return false;

	double tsum = 0;
	double tssq = 0;
	S = 0;
	
	ColumnVector X(na);
	int n = fbat_stat_observed(trt_id, offset, tsum, tssq, S, X, gx, flag, htab);
	if (n>0 && tssq>0) {
		ColumnVector EX(na);
		Matrix VX(na, na);
		Matrix CX(na, na);
		if (fbat_stat(EX, VX, CX, gx, htab)) {
			ES = tsum*EX;
			V = tssq*VX + (tsum*tsum-tssq)*CX;
			return true;
		}
	}
	
	return false;

}
*/

// added A & B for calculation of nuisance parameter u, where total Var(S') = A*u^2 - 2*B*u + Var(S), u=sum(B)/sum(A) minimize Var(S')
// A = n*[VX+(n-1)*CX], B=tsum*[VX+(n-1)*CX]
// S'-E(S') = S-E(S) + u*(X-EX)
bool Sufficient_Stat::fbat_stat3(int trt_id, double offset, ColumnVector &ES, ColumnVector &S, Matrix &V, Matrix &A, Matrix &B, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab) {
	
	int na = V.Nrows();
	
	if (htab && htab->size!=na)
		throw Param_Error("Sufficient_Stat::fbat_stat wrong matrix dimension");
	
	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat htab not specified");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL || trt_id<0)
		return false;
	
	S = 0;
	ES = 0;
	V = 0;
	A = 0;
	B = 0;
	X = 0;
	
	if (xlink)
		return fbat_stat3_x(trt_id, offset, ES, S, V, A, B, X, gx, flag, htab);

	double tsum = 0;
	double tssq = 0;
	
	int n = fbat_stat_observed(trt_id, offset, tsum, tssq, S, X, gx, flag, htab);
	if (n>0) {	// can't add tssq>0 condition if you want to estimate nuisance parameter
		ColumnVector EX(na);
		Matrix VX(na, na);
		Matrix CX(na, na);
		if (fbat_stat(EX, VX, CX, gx, htab)) {
			ES = tsum*EX;
			V = tssq*VX + (tsum*tsum-tssq)*CX;
			X -= n*EX;
			A = n*(VX + (n-1)*CX);
			B = tsum*(VX + (n-1)*CX);
			return true;
		}
	}
	
	return false;

}

/*
// added to recount informative fam# given nuisance parameter mu
bool Sufficient_Stat::fbat_stat4(int trt_id, Matrix &VX, Matrix &CX, double &tsum, double &tssq, int &noff, double gx[], int flag, HaploidTable *htab) {
	
	if (!informative())
		return false;
		
	tsum = 0;
	tssq = 0;
	noff = 0;
	
	int i, j;
	double t;
	for (i=0; i<n_offspring; i++) {
		if (otr[i]!=NULL && otr[i][trt_id]<kVeryLargeNumber && (j=which_comm_ogt(ogt[i], sex[i]))>=0 ) {
			t = otr[i][trt_id];	// for cencored data, t-u in affecteds, t in censored group
			tsum += t;
			tssq += t*t;
			noff++;
		}
	}	
	
	ColumnVector EX(VX.Nrows());
	if (noff==0 || !fbat_stat(EX, VX, CX, gx, htab))
		return false;
	
	return true;
			
}
*/

bool Mating_Haplotype::gt_compatible(Genotype *g, int sex, bool xlink) {
	
	int i, j, k;
	
	if (!xlink) {
		for (i=0; i<2 ; i++)
			for (j=2; j<4; j++) {
				for (k=0; k<nmrks && (g[k].a[0]==0 || g[k].equal(h[i][k],h[j][k])); k++)
					;
				if (k==nmrks)
					return true;
			}
	}
	else if (sex==1) {
		for (j=2; j<4; j++) {
			for (k=0; k<nmrks && (g[k].a[0]==0 || g[k].a[0]==h[j][k]); k++)
				;
			if (k==nmrks)
				return true;
		}
	} 
	else {
		for (j=2; j<4; j++) {
			for (k=0; k<nmrks && (g[k].a[0]==0 || g[k].equal(h[0][k],h[j][k])); k++)
				;
			if (k==nmrks)
				return true;
		}
	} 

	return false;
}

// presume no missing, so pls use common_gt only
void Mating_Haplotype::update_pg(Genotype *gt[], int sex[], int ng, bool sexlink) {
	
	if (ng>4)
		ng = 4;
	
	int i, j, k, l;
	
	for (i=0; i<ng; i++)
		pg[i] = 0.0;
	
	for (l=0; l<ng; l++) {
		if (sexlink && sex[l]==1) {
			for (j=2; j<4; j++) {
				for (k=0; k<nmrks && gt[l][k].a[0]==h[j][k]; k++)	// presume no missing, so pls use common_gt only
					;
				if (k==nmrks)
					pg[l] += 0.5;
			}
		} else {
			for (i=0; i<(sexlink?1:2); i++)
				for (j=2; j<4; j++) {
					for (k=0; k<nmrks && (gt[l][k].equal(h[i][k],h[j][k])); k++)
						;
					if (k==nmrks)
						pg[l] += (sexlink?0.5:0.25);
				}
		}
	}
	
}

bool Haplotype_List::identical(char *h1, char *h2) {
	
	int i;
	for (i=0; i<nmrk && h1[i]==h2[i]; i++)
		;
	
	return (i==nmrk);

}

bool Haplotype_List::compatible(char *h1, char *h2) {
	
	int i;
	
	for (i=0; i<nmrk && (h[0][i]!=kotherallele?(h[0][i]==h1[i]):(h1[i]!=h[1][i])) && (h[1][i]!=kotherallele?(h[1][i]==h2[i]):(h2[i]!=h[0][i])); i++)
		;
		
	if (i==nmrk)
		return true;
	
	for (i=0; i<nmrk && (h[0][i]!=kotherallele?(h[0][i]==h2[i]):(h2[i]!=h[1][i])) && (h[1][i]!=kotherallele?(h[1][i]==h1[i]):(h1[i]!=h[0][i])); i++)
		;
	
	return (i==nmrk);

}

	
double Mating_Haplotype::offspring_ht_prob(Haplotype_List *hl, int sex, bool sexlink) {
	
	if (hl==NULL)
		return 0.0;
	else if (!sexlink)
		return (hl->equal(h[0],h[2])?0.25:0)+(hl->equal(h[0],h[3])?0.25:0)+(hl->equal(h[1],h[2])?0.25:0)+(hl->equal(h[1],h[3])?0.25:0);
	else if (sex==1)
		return (hl->identical(hl->h[0],h[2])?0.5:0)+(hl->identical(hl->h[0],h[3])?0.5:0);
	else {
		int j;	// hl->h[j] from father
		for (j=0; j<2 && !hl->identical(hl->h[j],h[0]); j++)
			;
		return (j==2? 0 : (hl->identical(hl->h[1-j],h[2])?0.5:0)+(hl->identical(hl->h[1-j],h[3])?0.5:0));
	}

}

HaploidTable::HaploidTable(Locus_Info *loc, int n, int idx[], Haplotype_List *hl) {
	
	size = 0;
	hpl = NULL;
	nmrk = n;
	if (n>0 && loc!=NULL) {
		loc_ptr = new Locus_Info*[n];
		for (int i=0; i<n; i++)
			loc_ptr[i] = loc + idx[i];
	} else
		loc_ptr = NULL;
	
	if (hl!=NULL)
		maketable(hl);

}

void HaploidTable::print(mystream &os) {
	
	int j, k;
	
	MYOUTPUT(os, "\nestimating haplotype frequencies for the selected markers:\n")
	for (k=0; k<nmrk; k++)
		MYOUTPUT(os, loc_ptr[k]->name << "  ")
	MYOUTPUT(os, "\n")
	
	MYOUTPUT(os, "\n" << setw(5) << "a#" << setw(nmrk*4) << "Hapl" << setw(10) << "freq" << "\n")
	MYOUTPUT(os, "---------------------------------------\n")
	
	HaploidList *hpl2;
	for (k=1,hpl2=hpl; hpl2!=NULL; hpl2=hpl2->next,k++) {
		MYOUTPUT(os, "a" << setw(4) << k)
		for (j=0; j<nmrk; j++)
			MYOUTPUT(os, setw(4) << loc_ptr[j]->alleleName[hpl2->node->h[j]])
		MYOUTPUT(os, setw(10) << hpl2->node->p << "\n")
	}


}

Sufficient_Stat::~Sufficient_Stat() {
	
	int i;
	for (i=0; i<2; i++)
		if (pgt[i])
			delete[] pgt[i];
			
	for (i=0; i<kMaxMembers; i++)
		if (ogt[i])
			delete[] ogt[i];
			
	for (i=0; i<4; i++) {
		if (common_ogt[i])
			delete[] common_ogt[i];
		if (common_oht[i])
			delete common_oht[i];
	}
	
	if (mhap)
		delete mhap;
		
		
	if (ogpat)
		delete ogpat;
	
	if (GX)
		delete GX;
	
	if (EGX)
		delete EGX;

}

// calculate expected genotype probability in offspring based on pg[]
// return false if this family is not informative and pg[] is not calculated
bool Sufficient_Stat::marker_prob(int a_idx, double mrk_prob[]) {
	
	//if (nMarkers!=1 || a_idx<1 || a_idx>2)
	//	throw ERROR("Sufficient_Stat::marker_prob only works for bi-allelic single-marker");
	
	if (!informative())
		return false;
	
	int i, k;
	for (i=0; i<2; i++)
		mrk_prob[i] = 0;
	for (i=0; i<n_common_ogt; i++) {
		k = common_ogt[i][0].allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
		if (k>0)
			mrk_prob[k-1] += pg[i];
	}
	
	return true;

}

// return # of informative offs
int Sufficient_Stat::fbat_stat_observed_x(int trt_id, double offset, double tsum[], double tssq[], int noff[], ColumnVector &S, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab) {
	
	double t;
	int i, j;
	
	for (i=0; i<2; i++) {
		tsum[i] = 0;
		tssq[i] = 0;
		noff[i] = 0;
	}
	
	S = 0;
	X = 0;
	ColumnVector Y = X;
	for (i=0; i<n_offspring; i++) {
		if (otr[i]!=NULL && otr[i][trt_id]<kVeryLargeNumber && (j=which_comm_ogt(ogt[i], sex[i]))>=0 ) {
			t = otr[i][trt_id] - ( (!(flag&flag_censored) || otr[i][0]>0.5)? offset : 0 );	// for cencored data, t-u in affecteds, t in censored group
			tsum[sex[i]==1?0:1] += t;
			tssq[sex[i]==1?0:1] += t*t;
			noff[sex[i]==1?0:1]++;
			comm_ogt_X(Y, j, gx, htab);
			X += Y;
			S += t*Y;
		}
	}
	
	return noff[0]+noff[1];
}

bool Sufficient_Stat::fbat_stat3_x(int trt_id, double offset, ColumnVector &ES, ColumnVector &S, Matrix &V, Matrix &A, Matrix &B, ColumnVector &X, double gx[][3], int flag, HaploidTable *htab) {
	
	int na = V.Nrows();
	
	if (htab && htab->size!=na)
		throw Param_Error("Sufficient_Stat::fbat_stat wrong matrix dimension");
	
	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat htab not specified");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL || trt_id<0)
		return false;

	double tsum[2] = {0, 0};
	double tssq[2] = {0, 0};
	int noff[2] = {0,0};
	
	S = 0;
	ES = 0;
	V = 0;
	A = 0;
	B = 0;
	fbat_stat_observed_x(trt_id, offset, tsum, tssq, noff, S, X, gx, flag, htab);
	if (noff[0]+noff[1]==0)
		return false;
	
	int i, j;
	int a[2];
	Haplotype_List *hl;
	Matrix Xm(na, 4);
	Matrix Xf(na, 4);
	Xm = 0;
	Xf = 0;
	if (nMarkers>1) {	// haplotype test
		for (i=0; i<n_common_ogt; i++) {
			for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
				if (hl->p>0) {
					for (j=0; j<2; j++)
						a[j] = htab->allele(hl->h[j]);
					if (a[0]>0) {
						if (common_ogt_sex[i]==1)
							Xm(a[0], i+1) += hl->p;
						else {
							Xf(a[0], i+1) += hl->p*gx[a[0]-1][(a[0]==a[1]?2:1)];
							if (a[1]>0 && a[0]!=a[1])
								Xf(a[1], i+1) += hl->p*gx[a[1]-1][1];
						}
					}
				}
		} 
	} else {
		for (i=0; i<n_common_ogt; i++) {
			for (j=0; j<2; j++)
				a[j] = (htab? htab->allele(&common_ogt[i][0].a[j]) : common_ogt[i][0].a[j]);
			if (a[0]>0) {
				if (common_ogt_sex[i]==1)
					Xm(a[0], i+1) += 1;
				else {
					Xf(a[0], i+1) += gx[a[0]-1][(a[0]==a[1]?2:1)];
					if (a[1]>0 && a[0]!=a[1])
						Xf(a[1], i+1) += gx[a[1]-1][1];
				}
			}
		}
	}
	
	ColumnVector EXm(na);
	ColumnVector EXf(na);
	Matrix EX2m(na, na);
	Matrix EX2f(na, na);
	Matrix EXmXm(na, na);
	Matrix EXfXf(na, na);
	Matrix EXmXf(na, na);
	EXm = 0;
	EXf = 0;
	EX2m = 0;
	EX2f = 0;
	EXmXm = 0;
	EXfXf = 0;
	EXmXf = 0;
	for (i=0; i<n_common_ogt; i++) {
		if (common_ogt_sex[i]==1) {
			EXm += pg_x[0][i]*Xm.Column(i+1);
			EX2m += pg_x[0][i]*(Xm.Column(i+1)*Xm.Column(i+1).t());
		} else {
			EXf += pg_x[1][i]*Xf.Column(i+1);
			EX2f += pg_x[1][i]*(Xf.Column(i+1)*Xf.Column(i+1).t());
		}
		for (j=0; j<n_common_ogt; j++) {
			if (common_ogt_sex[i]==common_ogt_sex[j]) {
				if (common_ogt_sex[i]==1)
					EXmXm += pgg_x[0][i][j]*(Xm.Column(i+1)*Xm.Column(j+1).t());
				else
					EXfXf += pgg_x[1][i][j]*(Xf.Column(i+1)*Xf.Column(j+1).t());
			} else {
				if (common_ogt_sex[i]==1)
					EXmXf += pgg_x[2][i][j]*(Xm.Column(i+1)*Xf.Column(j+1).t());
				else
					EXmXf += pgg_x[2][i][j]*(Xf.Column(i+1)*Xm.Column(j+1).t());
			}
		}
	}
	
	double t0 = tsum[0]*tsum[0];
	double t1 = tsum[1]*tsum[1];
	double tt = (tsum[0]+tsum[1])*(tsum[0]+tsum[1]);
	ES = tsum[0]*EXm + tsum[1]*EXf;
	V = tssq[0]*EX2m + (t0-tssq[0])*EXmXm + tssq[1]*EX2f + (t1-tssq[1])*EXfXf + (tt-t0-t1)*EXmXf - ES*ES.t();
	ColumnVector EX = (noff[0]*EXm + noff[1]*EXf);
	X -= EX;
	B = tsum[0]*(EX2m+(noff[0]-1)*EXmXm+noff[1]*EXmXf) + tsum[1]*(EX2f+(noff[1]-1)*EXfXf+noff[0]*EXmXf) - 0.5*(ES*EX.t()+EX*ES.t());
	A = noff[0]*(EX2m+(noff[0]-1)*EXmXm+noff[1]*EXmXf) + noff[1]*(EX2f+(noff[1]-1)*EXfXf+noff[0]*EXmXf) - EX*EX.t();

	return true;

}

// P(G_i|Y_i,M) = P(Y_i|G_i)*P(G_i|M)/sum_j[P(Y_i|G_j)*P(G_j|M)]
// E(S|Y,M) = sum_j(sum_i[P(G_i|Y_j,M)*X(G_i)]*Y_j)
// P(M_i|Y,Pa,G) = Product_k(P(G_k|M_i))*P(Y|M_i)*P(M_i|Pa)/sum_j[P(Y|M_j)*P(M_j|Pa)]
// P(Y|M_i) = P(Y_fa|M_i)*P(Y_mo|M_i)
// E(S|Y) = sum_i(P(M_i|Y,Pa)*E(S|Y,M_i))

// P(Y|G) = tr_info->prob_y_g()
double Sufficient_Stat::fbat_es_h1(Trait_Info *tr_info, int trt_id, double offset, int a_idx, double gx[], double pa) {
	
	if (nMarkers!=1 || a_idx>1)
		throw Param_Error("Sufficient_Stat::fbat_es_h1 requires bi-allelic single marker");
	
	int i, j, n;
	int sx[kMaxMembers];
	int gid[kMaxMembers];
	
	Genotype *gt[kMaxMembers];
	double y[kMaxMembers];
	double pg1[kMaxMembers][4];
	double x[4];
	
	// set gt score for each common_ogt
	for (i=0; i<n_common_ogt; i++)
		x[i] = (xlink && common_ogt_sex[i]==1?(common_ogt[i][0].a[0]==a_idx?1:0):gx[common_ogt[i][0].allele_cnt(a_idx)]);
		
	n = 0;
	for (i=0; i<n_offspring; i++) {
		if (otr[i]!=NULL && otr[i][trt_id]<kVeryLargeNumber && (j=which_comm_ogt(ogt[i], sex[i]))>=0 ) {
			sx[n] = sex[i];
			gid[n] = j;
			gt[n] = &common_ogt[j][0];
			y[n] = otr[i][trt_id] - offset;
			n++;
		}
	}
	
	double sum, psum = 0;
	Mating_Haplotype *mh;
	double q[3];	// genotype frequency from HWE
	q[0] = (1-pa)*(1-pa);
	q[1] = 2*pa*(1-pa);
	q[2] = pa*pa;
	
	for (mh=mhap; mh!=NULL; mh=mh->next) {
		// calculate P(G_i|Y_i,M) = pg1[person][comm_ogt] in each person
		for (i=0; i<n; i++) {
			sum = 0;
			for (j=0; j<n_common_ogt; j++) {
				pg1[i][j] = tr_info->prob_y_g(y[n], gt[i]->allele_cnt(a_idx, (sx[i]==1 && xlink)))*mh->pg[j];
				sum += pg1[i][j];
			}
			for (j=0; j<n_common_ogt; j++)
				pg1[i][j] /= sum;
		}
		// calculate E(S|Y,M) for allele a_idx
		sum = 0;
		for (i=0; i<n; i++)
			for (j=0; j<n_common_ogt; j++)
				sum += pg1[i][j]*y[i]*x[j];
		mh->es = sum;
		
		// calculate mating type frequency from allele grequency and observed offspring genotype under H1
		int a1, a2;
		a1 = (mh->h[0][0]==a_idx?1:0) + (!xlink && mh->h[1][0]==a_idx?1:0);
		a2 = (mh->h[2][0]==a_idx?1:0) + (mh->h[3][0]==a_idx?1:0);
		if (xlink)
			mh->q = (a1==1?pa:1-pa)*q[a2];
		else
			mh->q = (a1==a2?1:2)*q[a1]*q[a2];
		for (i=0; i<n; i++)
			mh->q *= pg1[i][gid[i]];
		psum += mh->q;
	}
	double es1 = 0;	// E(S) under H1
	for (mh=mhap; mh!=NULL; mh=mh->next) {
		mh->q /= psum;
		es1 += mh->es*mh->q;
	}
	
	// now calculate E(S|Y) under H0
	double tsum[2] = {0,0};
	double ex[2] = {0,0};
	for (i=0; i<n; i++) 
		tsum[(xlink && sx[i]==2? 1:0)] += y[i];
	
	for (i=0; i<n_common_ogt; i++)
		ex[(xlink && common_ogt_sex[i]==2? 1:0)] += pg[i];
	
	double es0 = tsum[0]*ex[0] + (xlink?tsum[1]*ex[1]:0);
	
	return es1 - es0;
	
}

// return var(s) for allele a_idx, s=S-ES for allele idx, used for mm-test
double Sufficient_Stat::fbat_stat1(int trt_id, int na, int a_idx, double offset, double gx[][3], int flag, double &s, HaploidTable *htab) {
	
	if (htab && htab->size!=na)
		throw Param_Error("Sufficient_Stat::fbat_stat1 wrong matrix dimension");
	
	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat1 htab not specified");
	
	if (a_idx<1 || a_idx>na)
		throw Param_Error("Sufficient_Stat::fbat_stat1 invalid a_idx");
		
	s = 0;
	
	ColumnVector FS(na);
	ColumnVector FES(na);
	Matrix FV(na,na);	
	Matrix A(na,na);	// these six are used for nuisance parameter etimation
	Matrix B(na,na);
	ColumnVector X(na);
	
	if (trt_id<0 || !informative() || !fbat_stat3(trt_id, offset, FES, FS, FV, A, B, X, gx, flag, htab))
		return 0;
	
	s = FS(a_idx) - FES(a_idx);
	
	return FV(a_idx,a_idx);
	
}

double Sufficient_Stat::ex(int a_idx, double gx[][3], HaploidTable *htab) {
	
	if (a_idx<1)
		throw Param_Error("invalid a_idx in Sufficient_Stat::ex");
	if (nMarkers>1 && htab==NULL)
		throw Param_Error("htab not defined in Sufficient_Stat::ex");
		
	int i, j, n;
	Genotype gt;
	Haplotype_List *hl;
	
	double ex = 0;
	
	for (i=0; i<n_common_ogt; i++) {
		if (nMarkers>1) {
			for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
				if (hl->p>0) {
					for (j=0; j<2; j++)
						gt.a[j] = htab->allele(hl->h[j]);
					n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
					ex += pg[i]*hl->p*(xlink && common_ogt_sex[i]==1? (n>0?1:0) : gx[a_idx-1][n]);
				}
		} else {
			for (j=0; j<2; j++)
				gt.a[j] = (htab? htab->allele(&common_ogt[i][0].a[j]) : common_ogt[i][0].a[j]);
			n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
			ex += pg[i]*(xlink && common_ogt_sex[i]==1? (n>0?1:0) : gx[a_idx-1][n]);
		}
	}
	
	return ex;

}

// return E(X) for male & female (into mx[]) seperately to accomandate X-markers
void Sufficient_Stat::ex(double mx[2], int a_idx, double gx[][3], HaploidTable *htab) {
	
	if (a_idx<1)
		throw Param_Error("invalid a_idx in Sufficient_Stat::ex");
	if (nMarkers>1 && htab==NULL)
		throw Param_Error("htab not defined in Sufficient_Stat::ex");
		
	int i, j, n;
	Genotype gt;
	Haplotype_List *hl;
	
	mx[0] = mx[1] = 0; 
	for (i=0; i<n_common_ogt; i++) {
		if (nMarkers>1) {
			for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
				if (hl->p>0) {
					for (j=0; j<2; j++)
						gt.a[j] = htab->allele(hl->h[j]);
					n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
					if (xlink) {
						if (common_ogt_sex[i]==1)
							mx[0] += pg_x[0][i]*hl->p*(n>0?1:0);
						else
							mx[1] += pg_x[1][i]*hl->p*gx[a_idx-1][n];
					} else
						mx[0] += pg[i]*hl->p*gx[a_idx-1][n];
				}
		} else {
			for (j=0; j<2; j++)
				gt.a[j] = (htab? htab->allele(&common_ogt[i][0].a[j]) : common_ogt[i][0].a[j]);
			n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
			if (xlink) {
				if (common_ogt_sex[i]==1)
					mx[0] += pg_x[0][i]*(n>0?1:0);
				else
					mx[1] += pg_x[1][i]*gx[a_idx-1][n];
			} else
				mx[0] += pg[i]*gx[a_idx-1][n];
		}
	}
	
	if (!xlink)
		mx[1] = mx[0];

}


//distribution [E(x), Var(x), Cov(x1,x2)] of x for allele a_idx
bool Sufficient_Stat::xdist(int a_idx, double &ex, double &vx, double &cx, double gx[][3], HaploidTable *htab) {
	
	if (a_idx<1)
		throw Param_Error("invalid a_idx in Sufficient_Stat::ex");
	if (nMarkers>1 && htab==NULL)
		throw Param_Error("htab not defined in Sufficient_Stat::ex");
		
	int i, j, n;
	Genotype gt;
	Haplotype_List *hl;
	
	ex = 0;
	vx = 0;
	cx = 0;
	
	if (!informative())
		return false;
	
	ColumnVector X(n_common_ogt);
	X = 0;
	for (i=0; i<n_common_ogt; i++) {
		if (nMarkers>1) {
			for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
				if (hl->p>0) {
					for (j=0; j<2; j++)
						gt.a[j] = htab->allele(hl->h[j]);
					n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
					X(i+1) += hl->p*(xlink && common_ogt_sex[i]==1? (n>0?1:0) : gx[a_idx-1][n]);
				}
		} else {
			for (j=0; j<2; j++)
				gt.a[j] = (htab? htab->allele(&common_ogt[i][0].a[j]) : common_ogt[i][0].a[j]);
			n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
			X(i+1) = (xlink && common_ogt_sex[i]==1? (n>0?1:0) : gx[a_idx-1][n]);
		}
	}
	
	for (i=0; i<n_common_ogt; i++) {
		ex += pg[i]*X(i+1);
		vx += pg[i]*X(i+1)*X(i+1);
		for (j=0; j<n_common_ogt; j++)
			cx += pgg[i][j]*X(i+1)*X(j+1);
	}
	
	vx -= ex*ex;
	cx -= ex*ex;
	
	return true;

}

// S & V for allele a_idx and multiple phenotypes
bool Sufficient_Stat::fbat_stat_mp(int ntr, int trt_id[], int a_idx, ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab) {
	
	double ex, vx, cx;
	
	S = 0;
	V = 0;
	
	if (xlink)
		return fbat_stat_x_mp(ntr, trt_id, a_idx, S, V, gx, flag, htab);
	
	if (V.Nrows()!=ntr)
		throw Param_Error("wrong dimension of V in Sufficient_Stat::fbat_stat_mp");
	
	if (flag&flag_censored)
		throw Param_Error("can't do censored trait in Sufficient_Stat::fbat_stat_mp");
		
	if (!xdist(a_idx, ex, vx, cx, gx, htab))
		return false;
	
	ColumnVector T(ntr), U(ntr), O(ntr);	// t_i, sum of t_i, sum of observed S for each trait
	Matrix Q(ntr, ntr);	// sum of t_i*t_j
	
	int i, j, k, n;
	Genotype gt;
	Haplotype_List *hl;
	ColumnVector X(n_common_ogt);
	X = 0;
	for (i=0; i<n_common_ogt; i++) {
		if (nMarkers>1) {
			for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
				if (hl->p>0) {
					for (j=0; j<2; j++)
						gt.a[j] = htab->allele(hl->h[j]);
					n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
					X(i+1) += hl->p*(xlink && common_ogt_sex[i]==1? (n>0?1:0) : gx[a_idx-1][n]);
				}
		} else {
			for (j=0; j<2; j++)
				gt.a[j] = (htab? htab->allele(&common_ogt[i][0].a[j]) : common_ogt[i][0].a[j]);
			n = gt.allele_cnt(a_idx, (xlink && common_ogt_sex[i]==1));
			X(i+1) = (xlink && common_ogt_sex[i]==1? (n>0?1:0) : gx[a_idx-1][n]);
		}
	}

	O = 0;
	U = 0;
	Q = 0;
	for (i=0; i<n_offspring; i++) {
		if (otr[i]!=NULL && (k=which_comm_ogt(ogt[i], sex[i]))>=0) {
			for (j=0; j<ntr; j++)
				T(j+1) = (otr[i][trt_id[j]]<kVeryLargeNumber?otr[i][trt_id[j]]:0);						
			U += T;
			O += X(k+1)*T;
			Q += T*T.t();
		}
	}
		
	S = O - ex*U;
	V = vx*Q + cx*((U*U.t())-Q);
	
	return true;
}

void Sufficient_Stat::init_permute_X(int na, double gx[][3], HaploidTable *htab) {
	
	if (GX)
		delete GX;		
	GX = new Matrix(na,n_common_ogt);
	if (EGX)
		delete EGX;
	EGX = new ColumnVector(na);
	
	int i, j;
	int a[2];
	Haplotype_List *hl;

	(*GX) = 0;
	for (i=0; i<n_common_ogt; i++) {
		if (nMarkers>1) {
			for (hl=common_oht[i]; hl!=NULL; hl=hl->next)
				if (hl->p>0) {
					for (j=0; j<2; j++)
						a[j] = htab->allele(hl->h[j]);
					if (a[0]>0) {
						(*GX)(a[0],i+1) += hl->p*gx[a[0]-1][(a[0]==a[1] && (!xlink || common_ogt_sex[i]==2)?2:1)];
						if (a[1]>0 && a[0]!=a[1] && (!xlink || common_ogt_sex[i]==2))
							(*GX)(a[1],i+1) += hl->p*gx[a[1]-1][1];
					}
				}
		} else {
			for (j=0; j<2; j++)
				a[j] = (htab? htab->allele(&common_ogt[i][0].a[j]) : common_ogt[i][0].a[j]);
			if (a[0]>0) {
				(*GX)(a[0],i+1) += gx[a[0]-1][(a[0]==a[1] && (!xlink || common_ogt_sex[i]==2)?2:1)];
				if (a[1]>0 && a[0]!=a[1] && (!xlink || common_ogt_sex[i]==2))
					(*GX)(a[1],i+1) += gx[a[1]-1][1];
			}
		}
	}
	
	(*EGX) = 0;
	for (i=1; i<=na; i++)
		for (j=0; j<n_common_ogt; j++)
			(*EGX)(i) += pg[j]*(*GX)(i,j+1);
	
	// move offspring with good genotype to the front of the list 
	int bad_sex[kMaxMembers];
	double *bad_otr[kMaxMembers];
	Genotype *bad_ogt[kMaxMembers];
	int nbad = 0;
	for (i=0,j=0; i<n_offspring; i++) {
		if (which_comm_ogt(ogt[i],sex[i])>=0) {
			if (i!=j) {
				ogt[j] = ogt[i];
				otr[j] = otr[i];
				sex[j] = sex[i];
			}
			j++;
		} else {
			bad_ogt[nbad] = ogt[i];
			bad_otr[nbad] = otr[i];
			bad_sex[nbad] = sex[i];
			nbad++;
		}
	}
	for (i=0; i<nbad; i++) {
		ogt[j] = bad_ogt[i];
		otr[j] = bad_otr[i];
		sex[j] = bad_sex[i];
		j++;
	}

}

void Sufficient_Stat::permute(ColumnVector &S, int trt_id, double offset, int flag) {
	
	S = 0;
	if (!informative())
		return;
	
	Off_Genotype_Pattern *op;
	// randomly select op
	double d = rand_uniform();
	for (op=ogpat; op->next!=NULL && d>op->pg; op=op->next)
		d -= op->pg;
		
	// ranlomly permute genotype of op;
	int i, j;
	int n = 0;
	int gid[kMaxMembers];
	double r[kMaxMembers];
	for (i=0; i<n_common_ogt; i++)
		for (j=0; j<op->cnt[i]; j++)
			gid[n++] = i;
	for (i=0; i<n; i++)
		r[i] = rand_uniform();
	// brute force sort of r[i] since n is small
	int imax;
	double rmax;
	for (i=0; i<n-1; i++) {
		rmax = r[i];
		imax = i;
		for (j=i+1; j<n; j++)
			if (r[j]>rmax) {
				imax = j;
				rmax = r[j];
			}
		if (i!=imax) {	//swap r[i], r[imax]; gid[i], gid[imax]
			j = gid[i];
			gid[i] = gid[imax];
			gid[imax] = j;
			r[imax] = r[i];
			r[i] = rmax;
		}
	}
	
	double t;
	int na = GX->Nrows();
	for (i=0; i<n; i++) {
		t = (otr[i][trt_id]<kVeryLargeNumber? otr[i][trt_id]-((!(flag&flag_censored) || otr[i][0]>0.5)? offset : 0) : 0);
		for (j=1; j<=na; j++)
			S(j) += t*((*GX)(j,gid[i]+1) - (*EGX)(j));
	}
}

// S & V of fbat-gee for multiple phenotypes (S=S_OBS-S_EXP)
bool Sufficient_Stat::fbat_gee_mp(int na, int ntr, int trt_id[], ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab) {
	
	S = 0;
	V = 0;
	
	if (xlink)
		return fbat_gee_x_mp(na, ntr, trt_id, S, V, gx, flag, htab);
		//throw ERROR("fbat-gee not implemented for X-link markers");
		
	if (nMarkers>1 && htab!=NULL) {
		if (na==0)
			na = htab->size;
		else if (na!=htab->size)
			throw ERROR("wrong number of allele argument in fbat-gee");
	}
	
	if (V.Nrows()!=na*ntr)
		throw Param_Error("wrong dimension of V in Sufficient_Stat::fbat_gee_mp");
	
	if (flag&flag_censored)
		throw Param_Error("can't do censored trait in Sufficient_Stat::fbat_gee_mp");
		
	ColumnVector U(ntr), O(ntr);	// sum of t_i
	Matrix Q(ntr, ntr);	// sum of t_i*t_j
	
	ColumnVector X(na);
	ColumnVector EX(na);
	Matrix VX(na,na), CX(na,na);
	
	int i, j, k;
	if (fbat_stat(EX,VX,CX,gx,htab)) {
		U = 0;
		O = 0;
		Q = 0;
		for (i=0; i<n_offspring; i++) {
			if (otr[i]!=NULL && (k=which_comm_ogt(ogt[i], sex[i]))>=0) {
				/*
				for (j=0; j<ntr; j++)
					if ((O(j+1)=otr[i][trt_id[j]])>=kVeryLargeNumber)
						break;
				if (j==ntr) {	// no missing of T and G
					U += O;
					Q += O*O.t();
					comm_ogt_X(X, k, gx, htab);
					S += KP(O,X);
				}
				*/
				// changed 10/19/2007 to allow missing trait in mp (treat as t=0)
				for (j=0; j<ntr; j++)
					if ((O(j+1)=otr[i][trt_id[j]])>=kVeryLargeNumber)
						O(j+1) = 0;
				U += O;
				Q += O*O.t();
				comm_ogt_X(X, k, gx, htab);
				S += KP(O,X);
			}
		}
		S -= KP(U,EX);
		V = KP(Q,VX)+KP(U*U.t()-Q,CX);
		return true;
	}
	
	return false;

}

void Sufficient_Stat::comm_ogt_X(ColumnVector &X, int which, double gx[][3], HaploidTable *htab) {
	
	if (which<0 || which>=n_common_ogt)
		throw Param_Error("which_common_ogt out of range in Sufficient_Stat::comm_ogt_X");
	
	X = 0;
	bool female_or_auto =  (!xlink || common_ogt_sex[which]==2);
	int k, a[2];
	Haplotype_List *hl;
	if (nMarkers>1) {
		for (hl=common_oht[which]; hl!=NULL; hl=hl->next)
			if (hl->p>0) {
				for (k=0; k<2; k++)
					a[k] = (htab? htab->allele(hl->h[k]) : a[k]);
				if (a[0]>0) {
					X(a[0]) += hl->p*(female_or_auto?gx[a[0]-1][a[0]==a[1]?2:1]:1);
					if (a[1]>0 && a[0]!=a[1] && female_or_auto)
						X(a[1]) += hl->p*gx[a[1]-1][1];;
				}
			}
	} else {
		for (k=0; k<2; k++)
			a[k] = (htab? htab->allele(&common_ogt[which][0].a[k]) : common_ogt[which][0].a[k]);
		if (a[0]>0) {
			X(a[0]) += (female_or_auto?gx[a[0]-1][a[0]==a[1]?2:1]:1);
			if (a[1]>0 && a[0]!=a[1] && female_or_auto) {
				X(a[1]) += gx[a[1]-1][1];
			}
		}
	}

}

double Sufficient_Stat::comm_ogt_X(int a_idx, int which, double gx[3], HaploidTable *htab) {
	
	if (which<0 || which>=n_common_ogt)
		throw Param_Error("which_common_ogt out of range in Sufficient_Stat::comm_ogt_X");
	
	if (nMarkers>1 && htab==NULL)
		throw Param_Error("Haploidtable htab undefined in Sufficient_Stat::comm_ogt_X");
	
	if (a_idx<1)
		throw Param_Error("a_idx out of range in Sufficient_Stat::comm_ogt_X");
		
	double x = 0;
	bool female_or_auto =  (!xlink || common_ogt_sex[which]==2);
	int n;
	Haplotype_List *hl;
	if (nMarkers>1) {
		for (hl=common_oht[which]; hl!=NULL; hl=hl->next)
			if (hl->p>0) {
				n = (htab->allele(hl->h[0])==a_idx?1:0);
				x += hl->p * (female_or_auto? gx[n+(htab->allele(hl->h[1])==a_idx?1:0)] : n);
			}
	} 
	else if (htab) {
		n = (htab->allele(&common_ogt[which][0].a[0])==a_idx?1:0);
		x = (female_or_auto? gx[n+(htab->allele(&common_ogt[which][0].a[1])==a_idx?1:0)] : n);
	}
	else {
		n = (common_ogt[which][0].a[0]==a_idx?1:0);
		x = (female_or_auto? gx[n+(common_ogt[which][0].a[1]==a_idx?1:0)] : n);
	}
	
	return x;
}



bool Sufficient_Stat::fbat_ex(ColumnVector &EX, double gx[][3], HaploidTable *htab) {
	
	if (htab && htab->size!=EX.Nrows())
		throw Param_Error("Sufficient_Stat::fbat_stat wrong matrix dimension");
	
	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat htab not specified");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL)
		return false;
	
	EX = 0;
	ColumnVector X = EX;
	int i;
	for (i=0; i<n_common_ogt; i++) {
		comm_ogt_X(X, i, gx, htab);
		EX += pg[i]*X;
	}
	
	return true;

}

bool Sufficient_Stat::fbat_stat_x(ColumnVector &EXm, ColumnVector &EXf, Matrix &VXm, Matrix &VXf, Matrix &CXmm, Matrix &CXff, Matrix &CXmf, double gx[][3], HaploidTable *htab) {

	if (htab && htab->size!=VXm.Nrows())
		throw Param_Error("Sufficient_Stat::fbat_stat wrong matrix dimension");
	
	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat htab not specified");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL)
		return false;

	EXm = 0;
	EXf = 0;
	VXm = 0;
	VXf = 0;
	CXmm = 0;
	CXff = 0;
	CXmf = 0;
	
	ColumnVector Y = EXm;
	Matrix X(Y.Nrows(),n_common_ogt);
	Matrix EX2m = VXm;
	Matrix EX2f = VXm;
	Matrix EXmXm = VXm;
	Matrix EXfXf = VXm;
	Matrix EXmXf = VXm;
	
	int i, j;
	for (i=0; i<n_common_ogt; i++) {
		comm_ogt_X(Y, i, gx, htab);
		X.Column(i+1) = Y;
		if (common_ogt_sex[i]==1) {
			EXm += pg_x[0][i]*Y;
			EX2m += pg_x[0][i]*(Y*Y.t());
		} else {
			EXf += pg_x[1][i]*Y;
			EX2f += pg_x[1][i]*(Y*Y.t());
		}
	}
	
	for (i=0; i<n_common_ogt; i++)
		for (j=0; j<n_common_ogt; j++) {
			if (common_ogt_sex[i]==common_ogt_sex[j]) {
				if (common_ogt_sex[i]==1)
					EXmXm += pgg_x[0][i][j]*(X.Column(i+1)*X.Column(j+1).t());
				else
					EXfXf += pgg_x[1][i][j]*(X.Column(i+1)*X.Column(j+1).t());
			} else 
				EXmXf += pgg_x[2][i][j]*(X.Column(i+1)*X.Column(j+1).t());
		}
	
	VXm = EX2m - EXm*EXm.t();
	VXf = EX2f - EXf*EXf.t();
	CXmm = EXmXm - EXm*EXm.t();
	CXff = EXfXf - EXf*EXf.t();
	CXmf = EXmXf - 0.5*(EXm*EXf.t()+EXf*EXm.t());
	
	return true;
}

bool Sufficient_Stat::fbat_stat_x(int a_idx, double &exm, double &exf, double &vxm, double &vxf, double &cxmm, double &cxff, double &cxmf, double gx[3], HaploidTable *htab) {

	if (htab==NULL && nMarkers>1)
		throw Param_Error("Sufficient_Stat::fbat_stat_x htab not specified");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL)
		return false;

	exm = 0;
	exf = 0;
	vxm = 0;
	vxf = 0;
	cxmm = 0;
	cxff = 0;
	cxmf = 0;
	
	double ex2m = 0;
	double ex2f = 0;
	double exmxm = 0;
	double exfxf = 0;
	double exmxf = 0;
	
	int i, j;
	double x[4];
	for (i=0; i<n_common_ogt; i++) {
		x[i] = comm_ogt_X(a_idx, i, gx, htab);
		if (common_ogt_sex[i]==1) {
			exm += pg_x[0][i]*x[i];
			ex2m += pg_x[0][i]*x[i]*x[i];
		} else {
			exf += pg_x[1][i]*x[i];
			ex2f += pg_x[1][i]*x[i]*x[i];
		}
	}
	
	for (i=0; i<n_common_ogt; i++)
		for (j=0; j<n_common_ogt; j++) {
			if (common_ogt_sex[i]==common_ogt_sex[j]) {
				if (common_ogt_sex[i]==1)
					exmxm += pgg_x[0][i][j]*x[i]*x[j];
				else
					exfxf += pgg_x[1][i][j]*x[i]*x[j];
			} else 
				exmxf += pgg_x[2][i][j]*x[i]*x[j];
		}
	
	vxm = ex2m - exm*exm;
	vxf = ex2f - exf*exf;
	cxmm = exmxm - exm*exm;
	cxff = exfxf - exf*exf;
	cxmf = exmxf - exm*exf;
	
	return true;
}

bool Sufficient_Stat::fbat_gee_x_mp(int na, int ntr, int trt_id[], ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab) {
	
	S = 0;
	V = 0;
	
	if (nMarkers>1 && htab!=NULL) {
		if (na==0)
			na = htab->size;
		else if (na!=htab->size)
			throw ERROR("wrong number of allele argument in fbat-gee");
	}
	
	if (V.Nrows()!=ntr*na)
		throw Param_Error("wrong dimension of V in Sufficient_Stat::fbat_gee_mp");
	
	if (flag&flag_censored)
		throw Param_Error("can't do censored trait in Sufficient_Stat::fbat_gee_mp");
	
	ColumnVector Um(ntr), Uf(ntr), O(ntr);	// sum of t_i
	Matrix Qmm(ntr, ntr), Qff(ntr, ntr), Qmf(ntr,ntr);	// sum of t_i*t_j
	
	ColumnVector X(na), EXm(na), EXf(na);
	Matrix VXm(na,na), VXf(na,na), CXmm(na,na), CXff(na,na), CXmf(na,na);
	
	int i, j, k;
	if (fbat_stat_x(EXm, EXf, VXm, VXf, CXmm, CXff, CXmf, gx, htab)) {
		Um = 0;
		Uf = 0;
		Qmm = 0;
		Qff = 0;
		Qmf = 0;
		for (i=0; i<n_offspring; i++) {
			if (otr[i]!=NULL && (k=which_comm_ogt(ogt[i], sex[i]))>=0) {
				for (j=0; j<ntr; j++)
					if ((O(j+1)=otr[i][trt_id[j]])>=kVeryLargeNumber)
						O(j+1) = 0;
				comm_ogt_X(X, k, gx, htab);
				if (sex[i]==1) {
					Um += O;
					Qmm += O*O.t();
					S += KP(O,(X-EXm));
				} else {
					Uf += O;
					Qff +=  O*O.t();
					S += KP(O,(X-EXf));
				}
			}
		}
		
		
		V = KP(Qmm,VXm) + KP(Qff, VXf) + KP(Um*Um.t()-Qmm,CXmm) + KP(Uf*Uf.t()-Qff,CXff) + KP(0.5*(Um*Uf.t()+Uf*Um.t()),CXmf);
		
		return true;
	}
	
	return false;

}

// S & V for allele a_idx and multiple phenotypes for x-linked markers
bool Sufficient_Stat::fbat_stat_x_mp(int ntr, int trt_id[], int a_idx, ColumnVector &S, Matrix &V, double gx[][3], int flag, HaploidTable *htab) {
	
	S = 0;
	V = 0;
	
	int na;
	
	if (nMarkers>1 && htab!=NULL) {
		if (na==0)
			na = htab->size;
		else if (na!=htab->size || a_idx<1 || a_idx>na)
			throw ERROR("wrong number of allele argument in FBAT_LC_MP_X");
	}
	
	if (V.Nrows()!=ntr)
		throw Param_Error("wrong dimension of V in Sufficient_Stat::fbat_stat_x_mp");
	
	if (flag&flag_censored)
		throw Param_Error("can't do censored trait in Sufficient_Stat::fbat_stat_x_mp");
	
	ColumnVector Um(ntr), Uf(ntr), O(ntr);	// sum of t_i
	Matrix Qmm(ntr, ntr), Qff(ntr, ntr), Qmf(ntr,ntr);	// sum of t_i*t_j
	
	int i, j, k;
	double exm;
	double exf;
	double vxm;
	double vxf;
	double cxmm;
	double cxff;
	double cxmf;
	double x;
	if (fbat_stat_x(a_idx, exm, exf, vxm, vxf, cxmm, cxff, cxmf, gx[a_idx-1], htab)) {
		Um = 0;
		Uf = 0;
		Qmm = 0;
		Qff = 0;
		Qmf = 0;
		for (i=0; i<n_offspring; i++) {
			if (otr[i]!=NULL && (k=which_comm_ogt(ogt[i], sex[i]))>=0) {
				for (j=0; j<ntr; j++)
					if ((O(j+1)=otr[i][trt_id[j]])>=kVeryLargeNumber)
						O(j+1) = 0;
				x = comm_ogt_X(a_idx, k, gx[a_idx-1], htab);
				if (sex[i]==1) {
					Um += O;
					Qmm += O*O.t();
					S += O*(x-exm);
				} else {
					Uf += O;
					Qff +=  O*O.t();
					S += O*(x-exf);
				}
			}
		}
		
		V = vxm*Qmm + cxmm*((Um*Um.t())-Qmm) + vxf*Qff + cxff*((Uf*Uf.t())-Qff) + cxmf*(0.5*(Um*Uf.t()+Uf*Um.t()));
		
		return true;
	}
	
	return false;

}


// return true if vx>0, ex[2] for male and female seperately if xlink
/*
bool Sufficient_Stat::fbat_ex(int a_idx, double ex[], double gx[]) {
	
	if (nMarkers!=1)
		throw Param_Error("Sufficient_Stat::fbat_ex works for single markers only");
	
	if (mhap==NULL || n_common_ogt<2 || ogpat==NULL)
		return false;
		
		
		     
			 
			 
			 
	
	int i, j, k;
	
	double vx[2] = {0,0};
	double cx[2] = {0,0};
	double x[4];
	if (!xlink) {
		for (i=0; i<n_common_ogt; i++) {
			x[i] = gx[common_ogt[i][0].allele_cnt(a_idx, false)];
			ex[0] += pg[i]*x[i];
			vx[0] += pg[i]*x[i]*x[i];
		}
		for (i=0; i<n_common_ogt; i++)
			for (j=0; j<n_common_ogt; j++)
				cx[0] += pgg[i][j]*x[i]*x[j];		
		vx[0] -= ex[0]*ex[0];
		cx[0] -= ex[0]*ex[0];
	} else {
		for (i=0; i<n_common_ogt; i++) {
			k = common_ogt_sex[i];
			x[i] = (k==1?common_ogt[i][0].allele_cnt(a_idx, true):gx[common_ogt[i][0].allele_cnt(a_idx, false)]);
			ex[k-1] += pg_x[k-1][i]*x[i];
			vx[k-1] += pg[k-1][i]*x[i]*x[i];
		}
		for (k=0; k<2; k++) {
			for (i=0; i<n_common_ogt; i++)
				for (j=0; j<n_common_ogt; j++)
					cx[k] += pgg_x[k][i][j]*x[i]*x[j];
			vx[k] -= ex[k]*ex[k];
			cx[k] -= ex[k]*ex[k];
		}
	}
	
	int n[2] = {0,0};
	for (i=0; i<n_offspring; i++)
		if (ogt[i][0].a[0]>0 && otr[i]
	
	return true;
	
}
*/