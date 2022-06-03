#pragma once
#include <cstring>
#include "mystream.h"
#include "pedigree.h"
#include "shell.h"
#include "mysimplepointerlist.h"
#include "xfunc.h"

const int kMaxSelectTraits = 20;
const int kMaxSelectMarkers = 1000;

const double kConvFreq = 0.000001;

const char kRoot_Id[] = "0";

class DATA : public CMD_SHELL {
	public:
	PEDIGREELIST *peds;
	Locus_Info *loc;
	Trait_Info *trt;
	int nLoci;
	int nTraits;	
	bool xlink;
	int nPeds;
	int nNucfams;
	int nPersons;
	
	int n_sel_trt;
	int	trait_id[kMaxSelectTraits];
	int n_sel_mrk;
	int marker_id[kMaxSelectMarkers];

	double affection_trait[3];

	DATA() { peds=NULL;loc=NULL,trt=NULL;nPeds=0;nNucfams=0;nPersons=0;nLoci=0;nTraits=0;xlink=false;n_sel_trt=0;n_sel_mrk=0;affection_trait[0]=kEvenLargerNumber;affection_trait[1]=0;affection_trait[2]=1;}
	void init();
	
	virtual void read_data_file(char *s);
	
	void ReadPedFile(char *fname);
	void ReadPedFile(fstream &fs);
	void ReadPhenoFile(char *fname);
	void ReadMapFile(char *fname);
	void select_trait(char *s);	
	void view_marker(char *s);
	
	void sort_alleles(int loc_index);
	void mendelerr(char *s);
	void allelefrequency(char *s);
	void outgenotype(char *s);
	void haplotypefrequency(char *s);
	
	void do_breakpoint(char *s);
	void structure_summary(char *s);
	void display_trait();
	
	bool test_root(char *id) { return strcmp(id, kRoot_Id)==0; }
	PEDIGREE* ped_lookup(char *name);
	PERSON* id_lookup(char *fname, char *id);
	double* allelefrq_em(int loc_index, int maxCycles=50);	// EM based, not for X-linked
	double* allelefrq(int loc_index, bool weight_by_fam=true);	// simple tabulation of frq, optionally weight by root# of a pedigree
	void draw_allelefrq(int loc_index);
	int find_locus(char *name);
	int find_locus(int chrNum, int order);
	int find_trait(char *name);
	int ped_cleanup();
	void gtindex2allele(int gtidx, char a[]);
	void update_affection_trait();
	
	//void ld_matrix_print(int nmrk, int idx[], HaploidTable *htab);
	
	char line_terminator(fstream &fs, long &linesize);
	//double rsquare(int loc1, int loc2);
	
	double hetoerozygosity(int loc_index);
	void do_hetoerozygosity(char *s);
	
};
