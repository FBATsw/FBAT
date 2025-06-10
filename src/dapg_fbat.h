#pragma once
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "mlr_fbat.h"
#include "mystream.h" //!/

//!/ new
class finemap_results{
	public:
		vector<string> variant_names;
		vector<double> incl_probs;
		vector<double> log_bfs;
		vector<double> z_scores;
		vector<int> cluster_ids;
		
		vector<double> cluster_pips;
		vector<double> cluster_r2s;
		vector<int> cluster_counts;
	
};
//!/
class NSNP {

    public:

        string name;
        double incl_prob;
        int cluster; //signal cluster membership
};

bool sort_nsnp_dec_by_ip(const NSNP &lhs, const NSNP &rhs);



class Nmodel {

    public: 
        string id;
        double prob;
        double post_score;
        int size;
};
bool sort_nmodel_dec(const Nmodel &lhs, const Nmodel &rhs);



// set of models with same size
class size_model {

    public:
        int size;
        double log10_sum_post;
        map<string, double> post_map;
        vector<vector<int> > mvec; // models to be expanded for the next round
        vector<int> snp_cluster;
};

class controller {

    public:
		
        MLR mlr;

		map<int, string> geno_map; //needs to be initialized
		string locus; //init
        int n_var; //init
        
		
		vector<double> kv_vec; //init
		int max_size; //init
		
		double size_select_thresh; //init
		double cluster_pip_thresh; //init
		int    priority_msize; //init
		double log10_snp_thresh; //init
        double ld_control_thresh; //init
		int    size_limit; //init
		
        map<int, int> snp2cluster_map;
		vector<double>  pi_vec; //init
		vector<int>  null_config; //init
		double prior_ratio; //init
		vector<size_model> szm_vec;
		vector<int> cand_set;  
        map<int, int> cand_map;
		map<string,double> single_log10_abfv;
		map<string,double> single_z_score;
		
		// log10 normalizing constant
        double log10_pnorm;
		
		// for reporting
        map<string,int> nsnp_map;
        vector<NSNP> nsnp_vec;
        vector<Nmodel> nmodel_vec;
		
		
		ColumnVector *Z=NULL; //init
		Matrix *R=NULL; //init

		void initialize(ColumnVector &Zp, Matrix &Rp, int n_var_p, map<int, string> geno_map_p, double size_select_thresh_p, double log10_snp_thresh_p);
		void set_default_grid();
        void set_default_options(double correlation_control_thresh_p, double log10_snp_thresh_p);
		
		void set_prior(double pi1);
		void init();
		
		double compute_log10_prior(vector<int> &mcfg);
		size_model compute_post_model_single(vector<int>& bm);
		size_model compute_post_model(int size, int use_abs_cutoff);
		void fine_map();
		int backward_checking(vector<int>& bm, double log10_post);
		double conditional_est(vector<int>& bm);
		size_model append_post_model(int size, map<int, int> &black_list);
		finemap_results summarize_approx_posterior(bool verbose);
		void parse_nmodel(Nmodel nmod);
		
		double compute_r2(int i, int j);
		double compute_average_r2(const vector<int> & vec1, const vector<int> & vec2);
		double compute_average_r2(const vector<int> & vec);


};  




