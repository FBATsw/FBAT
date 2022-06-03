#pragma once
#include "data.h"
#include "newmat.h"
#include "sufficient_stat.h"
#include "ires.h"

enum TEST_MODE { by_allele, by_marker, by_both };
//enum GEN_MODEL { model_additive, model_dominant, model_recessive, model_genotype, model_user_defined };
//enum TRAIT_MODE { single_trait, multi_trait };

static const int kmaxcmh = 1000;
static const int flag_full_phase		= 1<<17;
static const int flag_censored			= 1<<18;
static const int flag_empirical_var		= 1<<19;
static const int flag_optimize_offset	= 1<<20;

class IRES;

class FBAT : public DATA {
public:
	int min_size;
	int rv_min_size;	// min_size for rare variants
	
	GEN_MODEL gen_model;
	TEST_MODE mode;
	int maxcmh;
	double alpha;
	double tr_offset;
	double min_freq;
	unsigned long num_sim_region; //!/
	unsigned long max_sim_region; //!/
	
    FBAT() { min_size=10; rv_min_size=1; gen_model=model_additive; mode=by_allele; alpha=1.0; tr_offset=0; min_freq=0; maxcmh=kmaxcmh;  num_sim_region=10000; max_sim_region=10000;} //!/ added num_sim_region and max_sim_region
	void init();

	void select_model(char *s);
	void select_mode(char *s);
	void select_alpha(char *s);
	void select_minsize(char *s);
	void select_minfreq(char *s);
	void select_offset(char *s);
	void select_maxcmh(char *s);
	void select_rvminsize(char *s);
	void set_affection_trait(char *s);
	
	void select_nsim_region(char *s); //!/
	void select_maxsim_region(char *s); //!/
	void gen_rvtest(char *s); //!/
	void gen_cond_rvtest(char *s);//!/
	int* do_rand_sim_fam(Sufficient_Stat* ss); //!/
	int draw_index_pat(double* probs,int len); //!/
	bool check_consistent_cond(Sufficient_Stat* ss, int index, int* random_index); //!/

	void fbat(char *s);
	void hbat(char *s);
	void hapview_ss(char *s);
	void pbat(char *s) ;
	void haplotypefrequency(char *s);
	
	void ld_matrix_print(int nmrk, int idx[], HaploidTable *htab, bool print_d=true);

	void print_settings();
	
	void gcode(double *gx, GEN_MODEL model);
	void hdt(int nmrk, int idx[], double gx[][3], int flag);
	void fdt(int locidx, double gx[][3], int flag, int width, IRES* ires, bool rvtest, int rvoption, bool print_header=true);	//print header with marker column width
	void x_fbat_stat(int loc_idx, int trt_id, ColumnVector &S, ColumnVector &X, Matrix &V, Matrix &A, Matrix &B, int *fcnt, double gx[][3], double offset, int flag, IRES* ires=NULL);
	void count_infomative_fam(double mu[], int na, int *cnt, double gx[][3], int flag, HaploidTable *htab=NULL); // count informative fam# given offset mu[]
	
	double uncensored_trait_mean(int tr_idx);
	double calc_chisq(ColumnVector &U, Matrix &V, double &chisq, double &df, int idx[]);
	double estimate_offset(int loc_index, int a_idx, double gscore[][3], double &vsum);
	double estimate_offset_mm(int nmrk, int loc_index[], int a_idx[], double gscore[][3]);
	double estimate_dominance_cm(int loc_index, int a_idx, int trt_id, int flag, double offset, double &tval);	// estimate dominance parameter using the conditional mean method
	//double beta_conditional_mean(int loc_index, int a_idx, int trt_id, int flag, double offset, double gscore[][3]);
	//void stat_conditional_mean(int loc_index, int a_idx, int trt_id, int flag, double offset, double gscore[], double &zs, double &zb);
	double pop_regression(int loc_index, int a_idx, int trt_id, int flag, double offset, double gscore[][3], bool use_cmm, double *b=NULL, double *v=NULL);
	double fbat_regression(int loc_index, int a_idx, int trt_id, int flag, double offset, double gscore[][3], double *b=NULL, double *v=NULL);
	
	void gt2x(Genotype &gt, ColumnVector &X, bool female, bool sexlink, double gx[][3]);
	void x_pop_mean(int loc_idx, int trt_id, ColumnVector &X, double gx[][3], int flag, bool use_cmm);
	void x_pop_stat(int loc_idx, int trt_id, ColumnVector &S, Matrix &V, Matrix &A, Matrix &B, int *fcnt, double gx[][3], double offset, int flag, bool use_cmm);
	void pop_asso(int locidx, double gx[][3], int flag, bool use_cmm, int width, bool print_header);
	void pop_count_infomative_fam(int loc_idx, double mu[], int *cnt, double gx[][3], int flag, bool use_cmm); // count informative fam# given offset mu[]
	
	void joint_asso(int locidx, double gx[][3], int flag, int width, bool print_header);
	void fdt_reg(int locidx, double gx[][3], int flag, int width, bool print_header);
	
	void var2cor(SymmetricMatrix &V);
	void var2cor(Matrix &V);
	
	double tr_anova(int trt_id, double &rho);	// fitting Y_ij = u_i + e_ij, u_i=N(0,Vc), e_ij=N(0,Ve), return Vt=Vc+Ve, rho=Vc/Vt
	
	Haplotype_List* hapfreq_ss(int nmrk, int idx[]);

	int mc_hbat(int nmrk, int idx[], double gx[][3], int flag, int max_runs=100000);
	
	int fbat_mm_test(int nmrk, int loc_idx[], double gscore[][3], int flag);
	int fbat_lc_test(int nmrk, int loc_idx[], double gscore[][3], int flag);
	int fbat_minp_test(int nmrk, int loc_idx[], double gscore[][3], int flag, int ncycles);
	//int fbat_lc_test_mp(int loc_idx, double gscore[][3], int flag);
	int fbat_lc_test_mp(int loc_idx, int a_idx, double gscore[][3], int flag);
	int fbat_gee_mp(int loc_index, double gscore[][3], int flag, bool print_header);
	
	double empirical_min_p(double maxz, SymmetricMatrix &R, int maxCycle);

/*
	void estimate_nuisance_param(int loc_index, double mu[]);
	void estimate_nuisance_param_censored(int loc_index, double mu[]);
	//void estimate_dominance_param(int loc_index, double dom[]=NULL);
	void dump_stat(int loc_index, bool by_family);
	void dumpstat_nu(int loc_index, bool by_family);
	void dumpstat_censored(int loc_index, bool by_family);
	void dump_empiric_stat(int loc_index, bool by_family);
	void dumpstat_mp(int loc_index, int ai, int tr_index[], int ntr, bool by_family);
	void dump_stat3(int loc_index, bool by_family);
	//void estimate_nuisance_param_hap(int idx[], double nu[], HaploidTable *htab, Haplotype_List *emhl);
	
	void count_informative(Haplotype_List *emhl, HaploidTable *htab, int idx[], ColumnVector &Mu, ColumnVector &W, int flag);
	
	//void hap(NUCFAMILY *fam, int nmrk, int *idx);
	int fdt1c(int loc_index);
	int fdt1(int locus_index, bool optimize_offset=false);	// allelic test
	int fdt2(int locus_index, bool optimize_offset=false);	// multi-allelic test
	int fdt3(int loc_index);								// test under genotype model
	int fdt4(int loc_index, bool optimize_offset);
	int fdt1mp(int loc_index, int ai, int tr_index[], int ntr);
	int fdt2mp(int loc_index, int tr_index[], int ntr);
	int empiric_fdt1(int loc_index);
	int empiric_fdt2(int loc_index);
	int sdt1(int locus_index);
	int sdt2(int locus_index);
	int empiric_hdt(int nmrk, int idx[], bool full_phased_only=false);
	int hdt(int nmrk, int idx[], bool full_phased_only=false);
	
	int TDTtransmit(int loc_index);
		
	//int nucfam_hbat_stat(HaploidTable *htab, Haplotype_List *emhl, NUCFAMILY *fam, double &tsum, double &tssq, ColumnVector &S, Matrix &X, ColumnVector &EX, ColumnVector &W, SymmetricMatrix &A, SymmetricMatrix &B, bool full_phased_only);
	//inline double xfunc(int a1, int a2, int a, double dom[]=NULL);
	//inline double xfunc(int g1, int a);

	double pexact(int a, int n);
	void fdt(char *s);
	void sdt(char *s);
	void view_stat(char *s);
	void fbat_all(char *s);
	void hapview(char *s);
	void small(char *s);
	void tdt(char *s);
	
	
	double censored_offset(int trait_id);	// return the offset for censored trait
	double stat_conditional_mean(int loc_index, int a_idx, int trait_id, int flag, double offset, double gscore[]);
	double estimate_nuisance_param_mm(int nmrk, int loc_index[], int a_idx[], double gscore[]);
	double stat_conditional_mean_offset(int loc_index, int a_idx, int trait_id, int flag, double gscore[], double &mu);
	double stat_conditional_mean_d(int loc_index, int a_idx, int trait_id, int flag, double &d);
*/
};

