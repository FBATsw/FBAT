#include <cstdio>
#include "pedigree.h"
#include "myerror.h"
#include "mystream.h"
#include "sufficient_stat.h"

char Locus_Info::allele_id(int aName) {
	
	char i;
	
	if (aName==0)
		return 0;
	
	for (i=1; i<=nAlleles && alleleName[i]!=aName; i++)
		;
	if (i>nAlleles) {
		if (nAlleles>=kMaxAlleles) {
			char s[255];
			sprintf(s, "marker %s: too many alleles (>%d)", name, kMaxAlleles);
			throw ERROR( s);
		}
		alleleName[++nAlleles] = aName;
	}
	
	return i;

}

double Trait_Info::prob_y_g(double y, int g) {
	
	if (g<0 || g>2)
		throw Param_Error("invalid genotype (allele acount) for P(Y|G)");
	
	if (nlevels!=2 && nlevels!=0)
		throw Param_Error("P(Y|G) only good for dichotamous trait or normally distributed quantitative trait");
		
	
	return (nlevels==2?penetrance[g]:pdf_norm(y, mean[g], sd));

}
	
PERSON::PERSON(char *s) { 
	
	if (s)
		strncpy(id,s,kMaxIDLen-1); 
	sex=0; 
	affection=0; 
	flag=0; 
	gt = NULL; 
	pt = new double[1]; // pt[0] will always be affection trait
	*pt = kEvenLargerNumber;
	ascend=NULL;
	for (int i=0; i<kMaxMarriages; i++)
		descend[i] = NULL; 

}

void PERSON::sex_status(char status) {
	
	if (status<0 || status>2) {
		char s[80];
		sprintf(s, "id %s : unknown sex code", id);
		throw ERROR(s);
	}
	if (sex==0)
		sex = status;
	else if (sex!=status && status) {
		char s[80];
		sprintf(s, "id %s : inconsistent sex code", id);
		throw ERROR(s);
	}
}

void PERSON::affection_status(char status) {
	
	if (status<0 || status>2) {
		char s[80];
		sprintf(s, "id %s : unknown affection code", id);
		throw ERROR(s);
	}
	if (affection==0) {
		affection = status;
		*pt = (affection==0?kVeryLargeNumber:affection-1.0);
	}
	else if (affection!=status && status) {
		char s[80];
		sprintf(s, "id %s : inconsistent affection code", id);
		throw ERROR(s);
	}
}

void PERSON::descends(NUCFAMILY *fam) {
	
	if (fam==NULL)
		return;
	
	int i;
	for (i=0; i<kMaxMarriages && descend[i]!=NULL; i++)
			;
	if (i>=kMaxMarriages) {
		char s[80];
		sprintf(s, "too many marriages, id %s; a maximum of %d allowed", id, kMaxMarriages);
		throw Data_Error(s);
	} else
		descend[i] = fam;

}

void  PERSON::ascends(NUCFAMILY *fam) {
	
	if (ascend!=NULL && ascend!=fam) {
		char s[80];
		sprintf(s, "more than one pair of parents, id %s", id);
		throw Data_Error(s);
	}
	ascend = fam;

}

NUCFAMILY::NUCFAMILY(PERSON *f, PERSON *m) { 
	
	members = new PERSONLIST(f, NULL);
	new PERSONLIST(m, members);
	
	if (f!=NULL)
		f->descends(this);
	
	if (m!=NULL)
		m->descends(this);
	
	stat = NULL;	
	flags = 0;
}

PERSON* NUCFAMILY::id_lookup(char *id) {
	
	PERSONLIST *plist;
	for (plist=members; plist!=NULL; plist=plist->next)
		if (plist->node && plist->node->match_id(id))
			return plist->node;
	
	return NULL;

}
	
bool NUCFAMILY::mendelerr(int loc_index, bool xlink) {
	
	int i, j;
	int nf, nm, na;
	char aa[4], bb[4], *a;	// bb[] is mark for presence of female aa[]
	Genotype fg[10], mg[10];
	Genotype *gt;
	PERSONLIST *plist;
	

	if ((fa() && (gt=fa()->gt) && gt[loc_index].defined(xlink))) {
		nf = 1;
		fg[0].set(gt[loc_index]);
	} else
		nf = 0;
	if ((mo() && (gt=mo()->gt) && gt[loc_index].defined())) {
		nm = 1;
		mg[0].set(gt[loc_index]);;
	} else
		nm = 0;
	if (nf==0 || nm==0) {	// at least 1 parent unknown
		for (i=0; i<4; i++)
			bb[i] = 0; 
		for (na=0,plist=sibs(); plist!=NULL; plist=plist->next)
			if (plist->node && plist->node->good_gt(loc_index,xlink)) {
				a = &plist->node->gt[loc_index].a[0];
				for (j=0; j<(xlink && plist->node->sex==1?1:2); j++) {
					for (i=0; i<na && a[j]!=aa[i]; i++)
						;
					if (i==na) // new allele
						if (na==(xlink?3:4)) 
							return true;	// more than 4 alleles, or 3 alleles for xlinked
						else
							aa[na++] = a[j];
					if (xlink && plist->node->sex==2)	// 
						bb[i] = 1;
				}
			}
		
		// fill permutation of unknow parent's genotypes with b[]			
		if (nf==0) {
			if (xlink) {
				for (i=0; i<na; i++) {
					if (bb[i])
						fg[nf++].set(aa[i],aa[i]);
				}
				if (nf==0) 	// no daughter observed, set fg to any gt since its irrelevant
					fg[nf++].set(1,1);
			} else {
				for (i=0; i<na; i++)
					for (j=i; j<na; j++)
						fg[nf++].set(aa[i],aa[j]);
			}
		}
		if (nm==0) {
			for (i=0; i<na; i++)
				for (j=i; j<na; j++)
					mg[nm++].set(aa[i],aa[j]);
		}
	}
	
	if (nf==0 || nm==0)
		return false;
	
	for (i=0; i<nf; i++)
		for (j=0; j<nm; j++)
			if (compatible_mating(fg[i], mg[j], loc_index, xlink))
				return false;
				
	return true;

}

void NUCFAMILY::reset_genotype(int loc_index) {
	
	PERSONLIST *plist;
	
	for (plist=members; plist!=NULL; plist=plist->next)
		if (plist->node && plist->node->gt)
			plist->node->gt[loc_index].set(0,0);		
	
}

// return number of sibs genotypes at loc_index
int NUCFAMILY::sibsGenotyped(int loc_index, bool xlink) {
	
	int n;
	PERSONLIST *sl;
	
	for (n=0,sl=sibs(); sl!=NULL; sl=sl->next)
		if (sl->node->good_gt(loc_index,xlink))
			n++;
	
	return n;

}
		
int NUCFAMILY::sib_count() {
	
	int n;
	PERSONLIST *sl;
	
	for (n=0,sl=sibs(); sl!=NULL; sl=sl->next)
		n++;
	
	return n;
	
}

NUCFAMILY::~NUCFAMILY() {
	
	int i;
	if (fa()) {
		for (i=0; i<kMaxMarriages && fa()->descend[i]!=this; i++)
			;
		for (; i<kMaxMarriages-1; i++)
			fa()->descend[i] = fa()->descend[i+1];
		fa()->descend[kMaxMarriages-1] = NULL;
	}
	if (mo()){
		for (i=0; i<kMaxMarriages && mo()->descend[i]!=this; i++)
			;
		for (; i<kMaxMarriages-1; i++)
			mo()->descend[i] = mo()->descend[i+1];
		mo()->descend[kMaxMarriages-1] = NULL;
	}
	for (PERSONLIST *slist=sibs(); slist!=NULL; slist=slist->next)
		slist->node->ascend = NULL;
	
	members->detachnode();
	delete members;
	members = NULL;
	
	if (stat!=NULL)
		delete stat;
	stat = NULL;

}
	
NUCFAMILY* PEDIGREE::fam_lookup(char *fid, char *mid) {
	
	NUCFAMILYLIST *flist;
	NUCFAMILY *fam;
	for (flist=nucfams; flist!=NULL; flist=flist->next) {
		fam = flist->node;
		if (fam && fam->fa() && fam->fa()->match_id(fid) && fam->mo() && fam->mo()->match_id(mid))
			return fam;
	}
	
	return NULL;

}
	
PERSON* PEDIGREE::id_lookup(char *id) {
	
	PERSONLIST *plist;
	for (plist=members; plist!=NULL; plist=plist->next)
		if (plist->node && plist->node->match_id(id))
			return plist->node;
	
	return NULL;

}

void PEDIGREE::removeTraits() {
	
	PERSONLIST *plist;
	for (plist=members; plist!=NULL; plist=plist->next)
		if (plist->node && plist->node->pt) {
			delete[] plist->node->pt;
			plist->node->pt = NULL;
		}
	
}

void PEDIGREE::print_marker(int idx, mystream &os, Locus_Info *ti, bool xlink) {
	
	for (PERSONLIST *pl=members; pl!=NULL; pl=pl->next)
		if (pl->node && pl->node->good_gt(idx,xlink)) {
			MYOUTPUT(os, setw(16) << name << setw(16) << pl->node->id << setw(16) << (pl->node->ascend?pl->node->ascend->fa()->id:0) << setw(16) << (pl->node->ascend?pl->node->ascend->mo()->id:0) << setw(4) << pl->node->sex)
			MYOUTPUT(os, ti[idx].alleleName[pl->node->gt[idx].a[0]] << "/" << ti[idx].alleleName[pl->node->gt[idx].a[0]])
		}

}

bool NUCFAMILY::compatible_mating(Genotype &fg, Genotype &mg, int locidx, bool xlink) {
	
	PERSONLIST *per;
	Genotype *g;
	for (per=sibs(); per!=NULL; per=per->next)
		if (per->node->good_gt(locidx,xlink)) {
			g = per->node->gt + locidx;
			if (xlink) {
				if (per->node->sex==1?(g->a[0]!=mg.a[0] && g->a[0]!=mg.a[1]):(!g->equal(fg.a[0], mg.a[0]) && !g->equal(fg.a[0], mg.a[1])))
					return false;
			} 
			else if (!g->equal(fg.a[0], mg.a[0]) && !g->equal(fg.a[1], mg.a[0]) && !g->equal(fg.a[0], mg.a[1]) && !g->equal(fg.a[1], mg.a[1]))
				return false;
		}
	
	return true;

}

