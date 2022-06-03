#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ostream>
#include <iomanip>
#include "myerror.h"
#include "myutil.h"
#include "data.h"
#include "breakpoint.h"

const int kMaxLineSize = 1000000;

void DATA::init() {
	
	static COMMAND cmds[] = {
		{	"load", 
 			"load [{[-x] ped, phe, map}] file_name  -- load formated data file; -x: sex-linked markers",
  			static_cast<funcp> (&DATA::read_data_file)
 		},
 		{	"afreq",
 			"afreq [marker...]  --  print out markers' allele frequency",
 			static_cast<funcp> (&DATA::allelefrequency)
 		},
 		{	"genotype",
 			"genotype pedigree mrk1 [mrk2..]  --  show markers' genotypes in pedigree",
 			static_cast<funcp> (&DATA::outgenotype)
 		},
		{	"trait", 
 			"trait [trait_name...] -- select trait or traits for analysis",
  			static_cast<funcp> (&DATA::select_trait)
 		},
 		{	"viewmarker", 
 			"viewmarker marker [pedigree] -- display detailed maker genotype distribution",
  			static_cast<funcp> (&DATA::view_marker)
 		},
		{	"heterozygosity", 
 			"heterozygosity [marker...]  --  print out markers' heterozygosity",
  			static_cast<funcp> (&DATA::do_hetoerozygosity)
 		},
		{	"structure", 
 			"structure  --  print out count of family types by parent# and sib#",
  			static_cast<funcp> (&DATA::structure_summary)
 		},
		{	"breakpoint", 
 			"breakpoint [-i# -d# -s#,#,#]  --  breakpoint analysis; -i#: # of support markers, -d#: geno drop rate",
  			static_cast<funcp> (&DATA::do_breakpoint)
 		},
	};

	addcmd(cmds, 8);
	
}

PEDIGREE* DATA::ped_lookup(char *name) {
	
	PEDIGREELIST* plist;
	PEDIGREE *p;
	for (plist=peds; plist!=NULL; plist=plist->next) {
		p = plist->node;
		if (p && p->match_name(name))
			return p;
	}
	return NULL;

}

PERSON* DATA::id_lookup(char *fname, char *id) {
	
	PEDIGREE* ped=ped_lookup(fname);
	
	return ped?ped->id_lookup(id):NULL;

}

void DATA::read_data_file(char *s) {
	
	char name[640], format[640], *c;
	
	int wc = WordCount(s);
	int i;
	for (i=2; i<=wc && GetWord(s, i, name) && *name=='-'; i++) {
		switch (name[1]) {
			case 'x': case 'X':
				xlink = true;
				break;
			default:
				MYOUTPUT(os, "unknown option\n")
				usage(s);
				return;
		}
	}
	
	if (i!=wc && i!=wc-1) {
		usage(s);
		return;
	}
	
	if (wc==i) {	// file name with extension
		GetWord(s, wc, name);
		for (c=name+strlen(name)-1; *c!='.' && c>name; c--)
			;
		if (c<=name) {
			MYOUTPUT(os, "error: file extension expected" << endl)
			return;
		}
		strncpy(format, c+1, 639);
	} else {
		GetWord(s, wc-1, format);
		GetWord(s, wc, name);
	}
	if (strcmp(format, "ped")==0)
		ReadPedFile(name);
	else if (strcmp(format, "phe")==0)
		ReadPhenoFile(name);
	else if (strcmp(format, "map")==0)
		ReadMapFile(name);
	else
		MYOUTPUT(os, "error: unknown format/extension " << format << endl);

}

void DATA::ReadPedFile(char *fname) {
	
	fstream fs;

	if (peds)
		throw ERROR("pedigree data already exists");
/*
	if (loc)
		throw ERROR("marker info already exists");
*/	
	fs.open(fname, ios::in);
	if (fs.fail())
		throw File_Error("error: openning ped file failed");
	
	try {
		ReadPedFile(fs);
	}
	catch (ERROR &err) {
		fs.close();
/*
		if (loc) {
			delete[] loc;
			loc = NULL;
			nLoci = 0;
		}
*/
		if (peds) {
			delete peds;
			peds = NULL;
			nPersons = 0;
			nPeds = 0;
			nNucfams = 0;
			nTraits = 0;
		}
		err.rethrow();
	}
	catch (...) {
		fs.close();
/*
		if (loc) {
			delete[] loc;
			loc = NULL;
			nLoci = 0;
		}
*/
		if (peds) {
			delete peds;
			peds = NULL;
			nPersons = 0;
			nPeds = 0;
			nNucfams = 0;
			nTraits = 0;
		}
		throw ERROR("a library run time error occured");
	}
	
	fs.close();
	
	for (int i=0; i<nLoci; i++)
		sort_alleles(i);
	
	mendelerr("mendelerr -a");
	
	trt = new Trait_Info[1];
	trt[0].set("affection");
	nTraits = 1;
	
	update_affection_trait();
	
}

// if nLoci==0, need the first line defination
void DATA::ReadPedFile(fstream &fs) {
	
	char *s, name[64], id[64], fid[64], mid[64], *c1, *c2;
	char a[2], b;
	int i, j;	
	PEDIGREE *ped;
	NUCFAMILY *fam;
	PERSON *per, *fa, *mo;
	
	
	long linesize;
	char ter = line_terminator(fs, linesize);
	// detect if the ped file has the locus name header and determine if has liability class
	s = new char[linesize];
	if (s==NULL)
		throw("Out of Memory");
	
	if (nLoci==0 && loc==NULL) {	// read in header
		fs.seekg(0);
		fs.getline(s, linesize, ter);	// 1st line
		if (fs.gcount()>=linesize-1)
			throw Data_Error("line size exceeding the limit");
		j = WordCount(s);			
		fs.getline(s, linesize, ter);
		i = WordCount(s);
		if (i==j)
			throw Data_Error("no header in the ped file, please load the map file first");
		if (i!=2*j+6)
			throw Data_Error("number of fields inconsistent with number of markers");		
		if (j<1)
			throw Data_Error("no genotype data defined");
		if ((loc=new Locus_Info[j])==NULL)
			throw Mem_Error("out of memory");
		nLoci = j;		
		fs.seekg(0);
		fs.getline(s, linesize, ter);
		char *c = s;
		for (i=0; i<nLoci; i++) {
			c = GetWord(c, 1, loc[i].name);
			if (loc[i].name[0]==0) {
				sprintf(name, "error: null string name for marker %d", i+1);
				throw Data_Error(name);
			}
			loc[i].sexlinked(xlink);
			c += strlen(loc[i].name);
			while (*c==' ')
				c++;
			if (*c=='\t')
				c++;
		}
	}
	
	nPersons = 0;
	nPeds = 0;
	nNucfams = 0;
	if (peds!=NULL)
		throw ERROR("pedigree stack not cleaned up correctly, restart the program and try again");
		
	while (!fs.eof() && fs.getline(s, linesize, ter)) {
		if ((i=WordCount(s))==0)
			continue;
		else if (2*nLoci!=(i-6)) {	// modified 12/07/04 for strict consistency of column numbers
			sprintf(s+(i<30?i:30), " ... ==> inconsistent number of fields");
			throw Data_Error(s);
		}
		try {
			GetWord(s, 1, name);
			ped = ped_lookup(name);
			if (ped==NULL) {	// new pedigree
				if (peds)
					new PEDIGREELIST(ped=new PEDIGREE(name), peds);
				else
					peds = new PEDIGREELIST(ped=new PEDIGREE(name), peds);
				if (ped==NULL || peds==NULL)
					throw Mem_Error("out of memory");
				nPeds++;
			}
			GetWord(s, 2, id);
			if ((per=ped->id_lookup(id)) == NULL) {
				if  ((per=new PERSON(id))==NULL) 
					throw Mem_Error("out of memory");
				ped->addperson(per);
			}
			if (per->gt==NULL) {
				if ((per->gt=new Genotype[nLoci])==NULL)
					throw Mem_Error("out of memory");
				nPersons++;
			}				
			GetWord(s, 3, fid);
			GetWord(s, 4, mid);
			if (!test_root(fid)) {
				if ((fa=ped->id_lookup(fid)) == NULL) {
					if ((fa=new PERSON(fid))==NULL)
						throw Mem_Error("out of memory");
					ped->addperson(fa);
				}
				fa->sex_status(1);
			}		
			if (!test_root(mid)) {
				if ((mo=ped->id_lookup(mid)) == NULL) {
					if ((mo=new PERSON(mid))==NULL)
						throw Mem_Error("out of memory");
					ped->addperson(mo);
				}
				mo->sex_status(2);
			}		
			if (!(test_root(fid) && test_root(mid))) {
				if ((fam=ped->fam_lookup(fid,mid))==NULL) {	// new NUCFAMILY
					fam = new NUCFAMILY(ped->id_lookup(fid), ped->id_lookup(mid));
					if (fam==NULL)
						throw Mem_Error("out of memory");
					ped->addfamily(fam);
					nNucfams++;
				}
				if (fam->id_lookup(id)) {
					sprintf(s, "pedigree %s, parental id (%s,%s), id %s redefined", ped->name, fid, mid, id);
					throw Data_Error(s);
				}
				fam->sib(per);
			}
			switch(*(GetWord(s, 5, NULL))) {	// sex
				case '1' : case 'M' : case 'm':
					per->sex_status(1); break;
				case '2' : case 'F' : case 'f':
					per->sex_status(2); break;
				case '0' :
					per->sex_status(0); break;
				default:
					throw Data_Error("unknown sex code");
			}
			c1 = GetWord(s, 6, NULL);	// affection status
			if (*c1>='0' && *c1<='2')
				per->affection_status(*c1-'0');
			else 
				throw Data_Error("unknown affection code");
			c1 = GetWord(s, 7, NULL);	// genotypes
			for (i=0; i<nLoci; i++) {
				a[0] = loc[i].allele_id((int)strtol(c1, &c2, 10));
				a[1] = loc[i].allele_id((int)strtol(c2, &c1, 10));
				if (a[0]<a[1]) {	// swap so a[0]>=a[1]
					b = a[1];
					a[1] = a[0];
					a[0] = b;
				}
				if (loc[i].sexlinked() && per->sex==1 && a[1]>0 && a[1]!=a[0]) {
	//					sprintf(s, "invalid sex-linked genotype: marker %s", loc[i].name);
	//					throw Data_Error(s);
					MYOUTPUT(os, "heterozygotic sex-linked genotype in men (ped " << ped->name << ",id " << per->id << ") at marker " << loc[i].name << ", reset to 0\n")
					a[1] = 0;
					a[0] = 0;
				}
				if (per->gt[i].defined() && !per->gt[i].equal(a[0], a[1])) {
					sprintf(s, "inconsistent genotype: marker %s", loc[i].name);
					throw Data_Error(s);
				}
				per->gt[i].set(a[0],a[1]);
				if (c1==c2) {
					sprintf(s, "numbered alleles expected: marker %s", loc[i].name);
					throw Data_Error(s);
				}
				// added 12/31/03 to disallow only one allele defined
				if ((!loc[i].sexlinked() || per->sex==2) && ((a[0]>0 && a[1]==0) || (a[0]==0 && a[1]>1))) {
					a[1] = 0;
					a[0] = 0;
					//sprintf(s, "ilegal genotype (only one allele defined): marker %s", loc[i].name);
					//throw Data_Error(s);
				}
			}			
		}
		catch (ERROR &err) {
			if (s)
				delete[] s;
			char msg[255];
			sprintf(msg, ": pedigree %s, id %s", ped->name, id);
			err.append(msg);
			err.rethrow();
		}
	}
	

	/* deleted 2/8/2004 to make outgenotype() output correct fid and mid (Kristel) */
	// reactivated the changed ped_cleanup with flag_uninformative 
	i = ped_cleanup();
	if (i>0) {
		nNucfams -= i;
		MYOUTPUT(os, "pedigree cleanup: " << i << " singleton nuclear families are removed\n")
	}
	
	MYOUTPUT(os, "read in: " << nLoci << " markers from " << nPeds)
	MYOUTPUT(os, " pedigrees (" << nNucfams << " nuclear families," << nPersons << " persons)\n")
	
	if (s)
		delete[] s;
}

void DATA::ReadPhenoFile(char *fname) {

	char *s, *c;
	char fno[256], id[256], name[256];
	int i;
	PEDIGREELIST *pl;
	PERSON *per;
	fstream fs;
	
	fs.open(fname, ios::in);
	if (fs.fail())
		throw File_Error("error: openning phenotype file failed");
	
	if (nPeds<1 || peds==NULL)
		throw ERROR("pedigree file must be loaded first");
		
		
	// clear previous phenotype
	if (nTraits>0) {
		for (pl=peds; pl!=NULL; pl=pl->next)
			pl->node->removeTraits();
		if (trt) {
			delete[] trt;
			trt = NULL;
		}
		nTraits = 0;
	}
	
	long linesize;
	char ter = line_terminator(fs, linesize);
		
	s = new char[linesize];
	if (s==NULL)
		throw("Out of Memory");
	
	int nHasTrait = 0;
	int nLines = 1;
	int nMissing = 0;

	try {
		// read nTraits
		fs.getline(s,  linesize, ter);
		nTraits = WordCount(s);
		if (nTraits==0)
			throw Data_Error("error: no trait defined at line 1");
		
		nTraits++;
		trt = new Trait_Info[nTraits];
		trt[0].set("affection", 2);
		for (i=1; i<nTraits; i++) {
			GetWord(s, i, name);
			trt[i].set(name);
		}
		
		while (!fs.eof() && fs.getline(s, linesize, ter)) {
			nLines++;
			if ((i=WordCount(s))<nTraits+1 && i>0) {
				MYOUTPUT(os, "warning at line " << nLines << ": incomplete traits\n")
				continue;
			}
			GetWord(s, 1, fno);
			GetWord(s, 2, id);
			if ((per=id_lookup(fno,id))==NULL) {
				nMissing++;
				continue;
			}
			
			if ((per->pt=new double[nTraits])==NULL)
				throw Mem_Error("error: out of memory");
			
			nHasTrait++;
			for (i=1; i<nTraits; i++) {
				GetWord(s, i+2, id);
				if (*id==0 || (*id=='-' && *(id+1)==0))
					per->notrait(i);
				else {
					per->pt[i] = strtod(id, &c);
					if (c==id) {
						per->notrait(i);
						sprintf(s, "invalid trait: pedigree %s, id %s, trait %d", fno, per->id, i+1);
						throw Data_Error(s);
					}
				}
			}
		}
	}
	catch (ERROR &err) {
		if (nTraits>0) {
			for (pl=peds; pl!=NULL; pl=pl->next)
				pl->node->removeTraits();
			if (trt) {
				delete[] trt;
				trt = NULL;
			}
			nTraits = 0;
		}
		if (s)
			delete[] s;
		s = NULL;
		
		err.rethrow();
	}

	MYOUTPUT(os, nTraits-1 << " quantitative traits have been successfully read\n")
	MYOUTPUT(os, nHasTrait << " persons have been phenotyped\n")
	if (nMissing)
		MYOUTPUT(os, "warning: " << nMissing << " persons are not in any pedigrees\n")
	
	if (s)
		delete[] s;
	
	update_affection_trait();	// set affection trait
}

// sort the alleles of a locus according to LOCUS::alleleName, update the genotype data as well
void DATA::sort_alleles(int loc_index) {
	
	char i, j, m, n;
	char id[kMaxAlleles];
	int name[kMaxAlleles];
	
	n = loc[loc_index].nAlleles;
	for (i=1; i<=n; i++) {
		id[i] = i;
		name[i] = loc[loc_index].alleleName[i];
	}
	
	// since nAllels typically < 20, a inefficient sort like this will be fine
	for (i=1; i<n; i++) {
		for (m=i,j=i+1; j<=n; j++)
			if (name[j]<name[m])
				m = j;
		if (m!=i) {	// swap m and i
			id[0] = id[i]; id[i] = id[m]; id[m] = id[0];
			name[0] = name[i]; name[i] = name[m]; name[m] = name[0];
		}
	}
	name[0] = 0;
	for (i=0; i<=n; i++)
		loc[loc_index].alleleName[i] = name[i];
	
	char order[kMaxAlleles];
	order[0] = 0;
	for (i=1; i<=n; i++)
		order[id[i]] = i;		
	// update genotype datas using new allele order
	PEDIGREELIST* dlist;
	PERSONLIST *rlist;
	Genotype *gt;
	for (dlist=peds; dlist!=NULL; dlist=dlist->next)
		for (rlist=dlist->node->members; rlist!=NULL; rlist=rlist->next)
			if (rlist->node && rlist->node->gt) {
				gt = rlist->node->gt+loc_index;
				gt->set(order[gt->a[0]], order[gt->a[1]]);
			}

}

/* 	mendelerr -a  : reset genotypes of all members of family with mendelian error to 0
	mendelerr -s  : if parental genotypes are available, reset only the errored sib's genotype to 0;
					else reset all members of family with mendelian error to 0
*/
void DATA::mendelerr(char *s) {
	
	int	i, n=0;
	
	PEDIGREELIST *plist;
	NUCFAMILYLIST *flist;
	
	if ((i=WordCount(s)) > 2)
		throw Param_Error("wrong number of command arguments");
	
	int cleanit = 0;
	if (i==2) {
		if (strcmp(GetWord(s, 2, NULL), "-s")==0)
			cleanit = 1;
		else if (strcmp(GetWord(s, 2, NULL), "-a")==0)
			cleanit = 2;
		else
			throw Param_Error("unknown command argument");
	}
	
	if (nNucfams<1 || peds==NULL || nLoci<1 || loc==NULL) 
		return;
	
	// added 1/14/09 to report excessive family/marker memdelian errors
	const float k_max_fam_err = 0.05;
	const float k_max_mrk_err = 0.05;
	const int min_fam_cnt = 100;
	const int min_mrk_cnt = 100;
	int *mrk_err_cnt = new int[nLoci];
	for (i=0; i<nLoci; i++)
		mrk_err_cnt[i] = 0;
	int fam_err_cnt;
	for (plist=peds; plist!=NULL; plist=plist->next) {
		for (flist=plist->node->nucfams; flist!=NULL; flist=flist->next) {
			fam_err_cnt = 0;
			for (i=0; i<nLoci; i++)
				if (loc[i].chrNum<24 && flist->node->mendelerr(i, loc[i].sexlinked())) {
					MYOUTPUT(os, "mendelian error: locus " << loc[i].name << ", pedigree " << plist->node->name << " [")
					if (flist->node->fa())
						MYOUTPUT(os, flist->node->fa()->id)
					if (flist->node->mo())
						MYOUTPUT(os, "," << flist->node->mo()->id << "]\n")
					else
						MYOUTPUT(os, ",]\n")
					n++;
					if (cleanit) 
						flist->node->reset_genotype(i);
					fam_err_cnt++;
					mrk_err_cnt[i]++;
				}
			if (nLoci >= min_mrk_cnt && fam_err_cnt >= k_max_fam_err*nLoci)
				MYOUTPUT(os, "Warning: the above family has excessive mendelian errors (" << 100.0*(float)fam_err_cnt/(float)nLoci << "%)\n")
		}
	}
	if (nNucfams >= min_fam_cnt) {
		for (i=0; i<nLoci; i++)
			if (mrk_err_cnt[i] >= k_max_mrk_err*nNucfams)
				MYOUTPUT(os, "Warning: " << loc[i].name << " has excessive mendelian errors (" << 100.0*(float)mrk_err_cnt[i]/(float)nNucfams << "%)\n")
	}
	
	delete[] mrk_err_cnt;
	
	if (n>0) {
		MYOUTPUT(os, "A total of " << n << " mendelian errors have been found\n")
		if (cleanit)
			MYOUTPUT(os, "genotypes of families with mendelian error have been reset to 0\n");
	}
}

// Marker nAlleles most_freq_allele median_allele allele_var
void DATA::draw_allelefrq(int loc_index) {
	
	double *frq = allelefrq(loc_index);
	
	int i;
	
	MYOUTPUT(os, loc[loc_index].name << "\n")
	for (i=1; i<=loc[loc_index].nAlleles; i++) {
		MYOUTPUT(os, setw(6) << loc[loc_index].alleleName[i] << setw(16) << frq[i])
		for (float f=frq[i]; f>0.0; f-=0.02)
			MYOUTPUT(os, "*");
		MYOUTPUT(os, "\n");
	}
	MYOUTPUT(os, "\n");
}

void DATA::allelefrequency(char *s) {
	
	int i=WordCount(s);
	
	if (i==1)	// do all loci
		for (i=0; i<nLoci; i++)
			draw_allelefrq(i);
	else {
		char name[256];
		int j, k;
		for (j=2; j<=i; j++) {
			GetWord(s, j, name);
			if ((k=find_locus(name))>=0)
				draw_allelefrq(k);
			else
				MYOUTPUT(os, "locus " << name << " not found\n");
		}
	}

}

int DATA::find_locus(char *name) {
	
	int i;
	for (i=0; i<nLoci && strcmp(loc[i].name, name); i++)
		;
	
	return i>=nLoci?-1:i;

}

/*
int DATA::find_locus(int chrNum, int order) {
	
	if (chrNum<1 || chrNum>23)
		return -1;
	
	int i;
	for (i=0; i<nLoci && (loc[i].chrNum!=chrNum || loc[i].order!=order); i++)
		;
	
	return i>=nLoci?-1:i;

}
*/

int DATA::find_trait(char *name) {
	
	int i;
	for (i=0; i<nTraits && !trt[i].match_name(name); i++)
		;
	
	return i==nTraits?-1:i;

}

// printout genotype of specified markers in a pedigree
// genotype ped marker1 ...
void DATA::outgenotype(char *s) {
	
	int i = WordCount(s);
	
	if (i<3)
		throw Param_Error("wrong number of command arguments");
	
	char ss[255];
	GetWord(s, 2, ss);
	PEDIGREE *ped = ped_lookup(ss);
	if (ped==NULL)
		throw Param_Error("specified pedigree not found");
		
	int *index = new int[i-2];
	if (index==NULL)
		throw Mem_Error("out of memory");
		
	int j, n;
	for (n=0,j=3; j<=i; j++) {
		GetWord(s, j, ss);
		index[n] = find_locus(ss);
		if (index[n]>=0)
			n++;
		else
			MYOUTPUT(os, "warning: marker " << ss << " not found\n");
	}
	
	if (n>0) {	
		for (i=0; i<n; i++)
			MYOUTPUT(os, loc[index[i]].name << " ");
		MYOUTPUT(os, "\n")
		
		PERSON *per;
		for (PERSONLIST *plist=ped->members; plist!=NULL; plist=plist->next) {
			per = plist->node;
			MYOUTPUT(os, per->id << " ")
			if (per->ascend && per->ascend->fa())
				MYOUTPUT(os, per->ascend->fa()->id << " ")
			else 
				MYOUTPUT(os, "0 ")
			if (per->ascend && per->ascend->mo())
				MYOUTPUT(os, per->ascend->mo()->id << " ")
			else 
				MYOUTPUT(os, "0 ")
			for (i=0; i<n; i++) {
				for (j=0; j<2; j++)
					if (per->gt==NULL)
						MYOUTPUT(os, "0 ")
					else
						MYOUTPUT(os, loc[index[i]].alleleName[per->gt[index[i]].a[j]] << " ")
				MYOUTPUT(os, "\n")
			}
		}
	}
	
	delete[] index;
	
}

// EM calculation for genotype frequency and then allele frequency from pedigree
double* DATA::allelefrq_em(int loc_index, int maxCycles) {

	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	PERSONLIST *rlist;
	Genotype *gt;
	int a[2], b[2], c[2];
	int i, j, k, m, n;
	int na = loc[loc_index].nAlleles;
	int ng = na*(na+1)/2;
	double *pg = new double[ng];
	double *qg = new double[ng];
	double *tq = new double[ng];
	double tqq[kMaxAlleles][kMaxAlleles];
	double tmp, tqsum;
	
	if (  loc_index<0 || loc_index>=nLoci)
		throw ERROR("ilegal locus #");
	
	if (loc[loc_index].frq)
		return loc[loc_index].frq;
		
	double *frq = new double[loc[loc_index].nAlleles+1];
	if (frq==NULL)
		throw Mem_Error("out of memory");
	
	loc[loc_index].frq = frq;
	for (i=0; i<=na; i++)
		frq[i] = 0.0;
	
	for (tmp=1.0/ng,i=0; i<ng; i++)
		pg[i] = tmp;
	
	for (n=0; n<maxCycles; n++) {
		for (i=0; i<ng; i++)
			qg[i] = 0.0;
	
		for (ped=peds; ped!=NULL; ped=ped->next) {
			flist=ped->node->nucfams;
			if (flist==NULL || flist->next!=NULL) {	// no families or complex pedigree, use roots only
				for (rlist=ped->node->members; rlist!=NULL; rlist=rlist->next)
					if (rlist->node && rlist->node->ascend==NULL && rlist->node->good_gt(loc_index))
						qg[rlist->node->gt[loc_index].index()]++; //qg[gtindex(rlist->node->gt[loc_index])]++;
			}
			else if ((i=flist->node->parentsGenotyped(loc_index))==2) {	// nuclear family with both parents
				qg[flist->node->fa()->gt[loc_index].index()]++;	// qg[gtindex(flist->node->fa()->gt[loc_index])]++;
				qg[flist->node->mo()->gt[loc_index].index()]++;	//qg[gtindex(flist->node->mo()->gt[loc_index])]++;
			}
			else if (i==1) {
				if (flist->node->mo() && flist->node->mo()->good_gt(loc_index)) {
					a[0] = flist->node->mo()->gt[loc_index].a[0];
					a[1] = flist->node->mo()->gt[loc_index].a[1];
					qg[flist->node->mo()->gt[loc_index].index()]++;	//qg[gtindex(flist->node->mo()->gt[loc_index])]++;
				} else {
					a[0] = flist->node->fa()->gt[loc_index].a[0];
					a[1] = flist->node->fa()->gt[loc_index].a[1];
					qg[flist->node->fa()->gt[loc_index].index()]++;	//qg[gtindex(flist->node->fa()->gt[loc_index])]++;
				}
				if (flist->node->sibsGenotyped(loc_index, loc[loc_index].sexlinked())==0)
					continue;
				for (i=0; i<ng; i++)
					tq[i] = pg[i];
				for (rlist=flist->node->sibs(); rlist!=NULL; rlist=rlist->next)
					if (rlist->node && rlist->node->good_gt(loc_index, loc[loc_index].sexlinked())) {
						gt = rlist->node->gt+loc_index;
						c[0] = gt->a[0];
						c[1] = gt->a[1];
						for (b[0]=1; b[0]<=na; b[0]++)
							for (b[1]=b[0]; b[1]<=na; b[1]++) {
								i = gtindex(b[0],b[1]);
								if (tq[i]==0.0)
									continue;
								for (j=0,k=0; k<2; k++)
									for (m=0; m<2; m++)
										if ((a[k]==c[0] && b[m]==c[1]) || (a[k]==c[1] && b[m]==c[0]))
											j++;
								tq[i] *= 0.25*j;
							}
					}
				for (tqsum=0.0,i=0; i<ng; i++)
					tqsum += tq[i];
				for (i=0; i<ng; i++)
					qg[i] += tq[i]/tqsum;
			}
			else if (flist->node->sibsGenotyped(loc_index, loc[loc_index].sexlinked())) {
				for (rlist=flist->node->sibs(); rlist!=NULL; rlist=rlist->next)
					if (rlist->node && rlist->node->good_gt(loc_index, loc[loc_index].sexlinked()))
						break;
				gt = rlist->node->gt+loc_index;
				a[0] = gt->a[0];
				b[0] = gt->a[1];
				for (i=1; i<=na; i++) {
					k = gtindex(a[0],i);
					tmp = pg[k];
					for (j=1; j<=na; j++) {
						m = gtindex(b[0],j);
						tqq[i-1][j-1] = tmp*pg[m]*(k==m?1:2);
					}
				}
				for (; rlist!=NULL; rlist=rlist->next)
					if (rlist->node && rlist->node->good_gt(loc_index, loc[loc_index].sexlinked())) {
						gt = rlist->node->gt+loc_index;
						c[0] = gt->a[0];
						c[1] = gt->a[1];
						for (a[1]=1; a[1]<=na; a[1]++)
							for (b[1]=1; b[1]<=na; b[1]++) {
								if (tqq[a[1]-1][b[1]-1]==0.0)
									continue;
								for (j=0,k=0; k<2; k++)
									for (m=0; m<2; m++)
										if ((a[k]==c[0] && b[m]==c[1]) || (a[k]==c[1] && b[m]==c[0]))
											j++;
								tqq[a[1]-1][b[1]-1] *= 0.25*j;
							}
					}
			
				for (tqsum=0.0,i=0; i<na; i++)
					for (j=0; j<na; j++)
						tqsum += tqq[i][j];
				for (i=0; i<na; i++)
					for (j=0; j<na; j++) {
						tmp = tqq[i][j]/tqsum;
						qg[gtindex(a[0],i+1)] += tmp;
						qg[gtindex(b[0],j+1)] += tmp;
					}
			}
		}
		
		for (tqsum=0.0,i=0; i<ng; i++)
			tqsum += qg[i];
			
		if (n==0)
			frq[0] = 2*tqsum;
			
		for (i=0; i<ng; i++)
			qg[i] /= tqsum;
			
		for (i=0; i<ng && fabs(pg[i]-qg[i])<kConvFreq; i++)
			;
		if (i==ng)
			break;
		
		for (i=0; i<ng; i++)
			pg[i] = qg[i];
	}
	
	for (i=1; i<=na; i++)
		for (j=i; j<=na; j++) {
			k = gtindex(i, j);
			frq[i] += pg[k];
			frq[j] += pg[k];
		}
	
	for (i=1; i<=na; i++)
		frq[i] *= 0.5;
	
	// delete[] pg;
	delete[] qg;
	delete[] tq;
	
	// added 10/14/99 for genotype frequency
	if (loc[loc_index].frq_gt)
		delete[] loc[loc_index].frq_gt;
	loc[loc_index].frq_gt = pg;
	
	return frq;

}

// remove nuclear families that 
// both parents are root and not genotyped and have only 1 offspring
// return number of nuclear families removed
int DATA::ped_cleanup() {
	
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	//PERSONLIST *rlist;
	int n = 0;
	/* changed 3/8/2004 to setflag(flag_uninformative) insted of delete uninformative nucfams
	for (ped=peds; ped!=NULL; ped=ped->next) {
		for (flist=ped->node->nucfams; flist!=NULL;) {	// modified 9/30/03 to add ->fa()==NULL, and ->mo()==NULL
			if ((flist->node->fa()==NULL || (flist->node->fa()->ascend==NULL && flist->node->fa()->gt==NULL)) && 
					(flist->node->mo()==NULL || (flist->node->mo()->ascend==NULL && flist->node->mo()->gt==NULL)) && 
					flist->node->sibs()->next==NULL) {
				alist = flist->next;
				ped->node->nucfams = ped->node->nucfams->detach(flist);
				delete flist;
				flist = alist;
				n++;
			} else
				flist=flist->next;
		}
	}
	*/
	
	for (ped=peds; ped!=NULL; ped=ped->next)
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next){
			if (flist->node->sibs()==NULL ||
					(flist->node->sibs()->next==NULL && 
						((flist->node->fa()==NULL || (flist->node->fa()->ascend==NULL && flist->node->fa()->gt==NULL)) || 
							(flist->node->mo()==NULL || (flist->node->mo()->ascend==NULL && flist->node->mo()->gt==NULL))))) {
				flist->node->setflag(flag_uninformative);
				n++;
				// MYOUTPUT(os, "[" << ped << "-" << flist << "](U)\n")
			}
		     else
		     {
			// MYOUTPUT(os, "[" << ped << "-" << flist << "]\n")
		     }
		}
	return n;
}
/*
void DATA::haplotypefrequency(char *s) {
	
	int i = WordCount(s);
	
	char name[256];
	
	int j, k;
	for (j=2; j<=i && GetWord(s,j,name) && *name=='-'; j++)
		switch (name[1]) {
			default:
				throw Param_Error("unknown option");
		}
	
	int nsel = j>=i? nLoci : i-j+1;
	int *idx = new int[nsel];
	if (j>i)
		for (nsel=0; nsel<nLoci; nsel++)
			idx[nsel] = nsel;
	else {
		for (nsel=0; j<=i; j++) {
			GetWord(s, j, name);
			if ((k=find_locus(name))>=0)
				idx[nsel++] = k;
		}
	}
	
	MYOUTPUT(os, "estimating phase genotype frequencies for the selected markers:\n")
	
	for (k=0; k<nsel; k++)
		MYOUTPUT(os, loc[idx[k]].name << "  ")
	MYOUTPUT(os, "\n")
	
	MYOUTPUT(os, "\n" << setw(nsel*4) << "Hapl" << setw(10) << "freq" << "\n")
	MYOUTPUT(os, "---------------------------------------\n")
	
	Haplotype_List *hl = hapfreq_ss(nsel, idx);
	
	if (hl!=NULL) {
		Haplotype_List *hl2;
		for (hl2=hl; hl2!=NULL; hl2=hl2->next) {
			if (hl2->p>=0.0001) {
				MYOUTPUT(os, "\n" << hl2->p << "\n")
				for (i=0; i<2; i++) {
					for (j=0; j<nsel; j++)
						MYOUTPUT(os, setw(4) << loc[idx[j]].alleleName[hl2->h[i][j]])
					MYOUTPUT(os, "\n")
				}
			}
		}
		
		while (hl!=NULL) {
			hl2 = hl->next;
			delete hl;
			hl = hl2;
		}
	}
	
	if (idx)
		delete[] idx;

}
*/

// to detect which line terminator is used: "\n" or "\r"
char DATA::line_terminator(fstream &fs, long &linesize) {
	
	streampos fp = fs.tellg();
	char *s;
	long sn, sr;

	linesize = kMaxLineSize;	
	while (true) {
		s = new char[linesize];
		if (s==NULL)
			throw("Out of Memory");
		fs.seekg(0);
		*s = 0;
		fs.getline(s, linesize, '\n');
		fs.clear();
		sn = fs.gcount();
		fs.seekg(0);
		*s = 0;
		fs.getline(s, linesize, '\r');	
		fs.clear();
		sr = fs.gcount();
		delete[] s;
		
		if (sn>=linesize-1 && sr>=linesize-1)
			linesize *= 10;
		else
			break;
	}
	
	fs.seekg(fp);
	
	if (sn==sr-1)
		return '\r';
	else if (sr==sn-1)
		return '\n';
	else
		return (sn<sr)? '\n' : '\r';

}

void DATA::select_trait(char *s) {
	
	int wcnt = WordCount(s);
	
	char word[255];
	
	if (wcnt==1)
		display_trait();
	else if (wcnt-1>kMaxSelectTraits) {
		sprintf(word, "over the limit of trait number (%d)",  kMaxSelectTraits);
		throw Param_Error(word);
	} else {
		int idx[kMaxSelectTraits];
		for (int i=0; i<wcnt-1; i++) {
			GetWord(s, i+2, word);
			idx[i] = find_trait(word);
			if (idx[i]<0) {
				MYOUTPUT(os, "error: trait \"" << word << "\" is unknown\n");
				return;
			}
		}
		// reset previous selection
		for (int i=0; i<n_sel_trt; i++)
			if (trait_id[i]>=0 && trait_id[i]<nTraits)
				trt[trait_id[i]].resetflag(flag_selected);
		// update new selection
		n_sel_trt = wcnt-1;
		for (int i=0; i<wcnt-1; i++) {
			trait_id[i] = idx[i];
			trt[trait_id[i]].setflag(flag_selected);
		}
		display_trait();
	}

}

void DATA::display_trait() {
	
	int i;
	
	for (i=0; i<nTraits; i++)
		MYOUTPUT(os, trt[i].name << (trt[i].testflag(flag_selected)?"** ":" "))
	
	MYOUTPUT(os, "\n")
}

void DATA::view_marker(char *s) {
	
	int i = WordCount(s);
	
	if (i!=3 && i!=2)
		throw Param_Error("wrong number of arguments");
	
	char name[255];
	GetWord(s, 2, name);
	int j = find_locus(name);
	if (j<0)
		throw ERROR("specified marker not found");
	
	if (i==2)	// do all families
		for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next)
			if (ped->node)
				ped->node->print_marker(j, os, loc, loc[j].sexlinked());
	else {
		GetWord(s, 3, name);
		PEDIGREE *ped = ped_lookup(name);
		if (ped==NULL)
			throw ERROR("specified pedigree not found");
		else
			ped->print_marker(j, os, loc, loc[j].sexlinked());
	}
}

void DATA::update_affection_trait() {
	
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next)
		for (PERSONLIST *pl=ped->node->members; pl!=NULL; pl=pl->next)
			if (pl->node && pl->node->pt)
				pl->node->pt[0] = affection_trait[pl->node->affection];
	
}

double* DATA::allelefrq(int loc_index, bool weight_by_fam) {
	
	Genotype *gt;

	int na = loc[loc_index].nAlleles;

	if (loc_index<0 || loc_index>=nLoci)
		throw ERROR("ilegal locus #");
	
	if (loc[loc_index].frq)
		return loc[loc_index].frq;

	int *pcnt = new int[na+1];
	int n_roots;
	
	double *frq = new double[na+1];
	if (frq==NULL)
		throw Mem_Error("out of memory");
	
	loc[loc_index].frq = frq;
	
	int i;
	for (i=0; i<=na; i++)
		frq[i] = 0.0;
	
	// changed 6/25/08 to weight offspring genotype on # of missing parents
	double miss_founder_allele;
	double w=0.;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		n_roots = 0;
		miss_founder_allele = 0;
		for (i=0; i<=na; i++)
			pcnt[i] = 0;	// allele cnt in founders
		for (PERSONLIST *rlist=ped->node->members; rlist!=NULL; rlist=rlist->next)
			if (rlist->node) {
				if (rlist->node->good_gt(loc_index,loc[loc_index].sexlinked())) {
					gt = rlist->node->gt + loc_index;
					if (rlist->node->testroot()) {
						n_roots++;
						frq[gt->a[0]]++;
						frq[0]++;
						if (!loc[loc_index].sexlinked() || rlist->node->sex==2)
							frq[gt->a[1]]++;
					} else {
						pcnt[gt->a[0]]++;
						pcnt[0]++;
						if (!loc[loc_index].sexlinked() || rlist->node->sex==2) {
							pcnt[gt->a[1]]++;
							pcnt[0]++;
						}
					}
				}
				else if (rlist->node->testroot()) // missing founder
					miss_founder_allele += ((!loc[loc_index].sexlinked() || rlist->node->sex==2)?2:1);
			}
			
		//Weight pcnt on minimal(miss_founder_allele, pcnt[0])
		if (pcnt[0]>0 && miss_founder_allele>0)
			w = (miss_founder_allele > pcnt[0]? 1 : miss_founder_allele/pcnt[0]);
		for (i=1; i<=na; i++)
			frq[i] += pcnt[i]*w;
	}
	
	double sum = 0;
	for (i=1; i<=na; i++)
		sum += frq[i];
	
	for (i=1; i<=na; i++)
		frq[i] /= sum;
	
	if (pcnt)
		delete[] pcnt;
	
	return frq;
	
}

// marker chr cm base sex_link
// mapfile should be read first if used
void DATA::ReadMapFile(char *fname) {
	
	fstream fs;
		
	fs.open(fname, ios::in);
		if (fs.fail())
			throw File_Error("error: opening map file failed");

	long linesize;
	char ter = line_terminator(fs, linesize);

	if (nLoci>0 && loc!=NULL) {
		// remove Locus_Info
		delete[] loc;
		loc = NULL;
		nLoci = 0;
	}
	
	try {		
		char s[256], *c1, *c2;	
		// Get total # of markers
		nLoci = 0;
		while (!fs.eof() && fs.getline(s, 256, ter))
			if (WordCount(s)==5)
				nLoci++;
		
		fs.clear();
		fs.seekg(0);
		
		loc = new Locus_Info[nLoci];
		if (loc==NULL)
			throw ERROR("out of memory");
			
		int	i=0, j, n=0;
		while (!fs.eof() && fs.getline(s, 256, ter)) {
			n++;
			if ((j=WordCount(s))==0)
				continue;
			else if (j!=5) {
				sprintf(s, "line %d: wrong field #", n);
				throw ERROR(s);
			}
			else {
				GetWord(s,1,loc[i].name);
				if (loc[i].name==NULL) {
					sprintf(s, "line %d: null marker nmae", n);
					throw ERROR(s);
				}
				loc[i].chrNum = strtol(c1=GetWord(s,2,NULL), &c2, 10);
				if (c1==c2 || loc[i].chrNum<0 || loc[i].chrNum>24) {
					sprintf(s, "line %d: invalid chr#", n);
					throw ERROR(s);
				}
				loc[i].cm = strtod(c2, &c1);
				if (c1==c2) {
					sprintf(s, "line %d: invalid genetic map position", n);
					throw ERROR(s);
				}
				loc[i].pos = strtol(c1, &c2, 10);
				if (c1==c2 || loc[i].pos<0) {
					sprintf(s, "line %d: invalid physical map position", n);
					throw ERROR(s);
				}
				j = strtol(c2, &c1, 10);
				if (c1==c2 || j<0 || j>1) {
					sprintf(s, "line %d: invalid sexlink code", n);
					throw ERROR(s);
				}
				loc[i].sexlinked(j==1);
				if (loc[i].chrNum>0 && loc[i].chrNum<23 && loc[i].sexlinked()) {
					sprintf(s, "line %d: inconsistent chr# & sexlink code", n);
					throw ERROR(s);
				}
				i++;
			}
		}
		
	}
	catch (ERROR &err) {
		fs.close();
		if (loc) {
			delete[] loc;
			loc = NULL;
			nLoci = 0;
		}
		err.rethrow();
	}
	
	fs.close();	
	MYOUTPUT(os, "read in " << nLoci << " markers' info\n");
	
}

// return proposion of heterozygotes in founders
double DATA::hetoerozygosity(int loc_index) {
	
	if (loc_index<0 || loc_index>=nLoci)
		throw ERROR("ilegal locus #");
	
	int cnt1 = 0, cnt2 = 0;
	
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (PERSONLIST *rlist=ped->node->members; rlist!=NULL; rlist=rlist->next)
			if (rlist->node && rlist->node->testroot() && rlist->node->good_gt(loc_index,loc[loc_index].sexlinked())) {
				cnt2++;
				if (rlist->node->gt[loc_index].homozygote())
					cnt1++;
			}
	}
	
	return 1.0 - (double)cnt1/(double)cnt2;
	
}

void DATA::do_hetoerozygosity(char *s) {
	
	int i=WordCount(s);
	
	if (i==1)	// do all loci
		for (i=0; i<nLoci; i++) {
			MYOUTPUT(os, loc[i].name << "\t" << loc[i].nAlleles << "\t" << hetoerozygosity(i) << "\n")
	} else {
		char name[256];
		int j, k;
		for (j=2; j<=i; j++) {
			GetWord(s, j, name);
			if ((k=find_locus(name))>=0)
				MYOUTPUT(os, loc[k].name << "\t" << loc[k].nAlleles << "\t" << hetoerozygosity(k) << "\n")
			else
				MYOUTPUT(os, "locus " << name << " not found\n");
		}
	}

}

void DATA::do_breakpoint(char *s) {
	
	int min_informative_flank_mrks = 5;
	float g_drop_rate = 0.1;
	bool simulate = false;
	int min_nc = 0;
	int max_nc = 6;
	int ncycle = 100;
	
	int i=WordCount(s);
	char *c, *cc;
	int j;
	for (j=2; j<=i && (c=GetWord(s,j,NULL)) && *c=='-'; j++) {
		switch (c[1]) {
			case 'i':
				min_informative_flank_mrks = strtol(c+2, &cc, 10);
				break;
			case 'd':
				g_drop_rate = strtod(c+2, &cc);
				break;
			case 's':
				if (sscanf(c+2,"%d,%d,%d", &min_nc, &max_nc, &ncycle)!=3)
					throw Param_Error("unknown simultion options");
				simulate = true;
				break;
			default:
				throw Param_Error("unknown option");
		}
	}
	
	if (j<=i)
		throw Param_Error("unknown arguments");
		
	int k;
	breakpoint brk(loc, nLoci);
	for (i=1; i<23; i++) {
		brk.order_marker(i);
		for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
			// do nuc_fam with multiple sibs
			for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
				if (flist->node) {
					flist->node->resetflag(flag_processed);
					if (!simulate) {
						if (brk.get(flist->node, min_informative_flank_mrks, g_drop_rate)) {
							j = (brk.avg_crossover[0]>=0?1:0) + (brk.avg_crossover[1]>=0?1:0);
							k = flist->node->sib_count();
							if (brk.avg_crossover[0]>=0)
								MYOUTPUT(os, ped->node->name << "\t" << flist->node->fa()->id << "\t" << j << "\t1\t" << i << "\t" << k << "\t" <<  brk.avg_crossover[0] << "\n")
							if (brk.avg_crossover[1]>=0)
								MYOUTPUT(os, ped->node->name << "\t" << flist->node->mo()->id << "\t" << j << "\t2\t" << i << "\t" << k << "\t" <<  brk.avg_crossover[1] << "\n")
							flist->node->setflag(flag_processed);
						}
					} else
						brk.get_simulate(os, flist->node, min_nc, max_nc, ncycle, min_informative_flank_mrks, g_drop_rate);
				}
			// do nuc_fam with 1+ phased parent and 1 offspring
			if (!simulate)
				for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
					if (flist->node && !flist->node->testflag(flag_processed) && brk.get_phased(flist->node, min_informative_flank_mrks, g_drop_rate)) {
						if (brk.avg_crossover[0]>=0)
							MYOUTPUT(os, ped->node->name << "\t" << flist->node->fa()->id << "\t3\t1\t" << i << "\t" << k << "\t" <<  brk.avg_crossover[0] << "\n")
						if (brk.avg_crossover[1]>=0)
							MYOUTPUT(os, ped->node->name << "\t" << flist->node->mo()->id << "\t3\t2\t" << i << "\t" << k << "\t" <<  brk.avg_crossover[1] << "\n")
					}
		}
	}
}

// report count of different type of nuc_families
void DATA::structure_summary(char *s) {
	
	int i, j;
	int cnt[3][6]; // cnt[parents#][sib#];
	
	for (i=0; i<3; i++)
		for (j=0; j<6; j++)
			cnt[i][j] = 0;
	
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			i = (flist->node->fa() && flist->node->fa()->gt? 1 : 0);
			if (flist->node->mo() && flist->node->mo()->gt)
				i++;
			j = 0;
			for (PERSONLIST *rlist=flist->node->sibs(); rlist!=NULL; rlist=rlist->next)
				if (rlist->node && rlist->node->gt)
					j++;
			if (j>5)
				j = 5;
			cnt[i][j]++;
		}
	}
	MYOUTPUT(os, "P#/S#   1       2       3       4       5\n")
	for (i=0; i<3; i++) {
		MYOUTPUT(os, setw(8) << i);
		for (j=1; j<6; j++)
			MYOUTPUT(os, setw(8) << cnt[i][j])
		MYOUTPUT(os, "\n")
	}
	
}

			
