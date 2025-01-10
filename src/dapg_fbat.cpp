#include "dapg_fbat.h"
using namespace std;
// independent function, does not call other function, okay
double log10_weighted_sum(vector<double> &vec, vector<double> &wts){

	size_t i;
    double max = vec[0];
    for(i=0;i<vec.size();i++){
        if(vec[i]>max)
            max = vec[i];
    }
    double sum = 0;
    for(i=0;i<vec.size();i++){
        sum += wts[i]*pow(10, (vec[i]-max));
    }

    return (max+log10(sum));
}

void controller::initialize(ColumnVector &Zp, Matrix &Rp, int n_var_p, map<int, string> geno_map_p){

    
    set_default_grid();  
	n_var=n_var_p;
	
	set_default_options();
    
	int i,j;
    Z=new ColumnVector(n_var); //non-standard
	R=new Matrix(n_var, n_var); //non-standard
	
	for (i=1; i<=n_var; i++){
		(*Z)(i) = Zp(i);
		for (j=1; j<=n_var; j++){
			(*R)(i,j) = Rp(i,j);
		}
	}
	double prior=1.0/(double)n_var;
	set_prior(prior);
	geno_map=geno_map_p;
	init();
}
// okay
///////////////////////////////////////////////////////////////////////////////////////

/*
void controller::set_default_grid(){
    if(use_ss!=1){
        omg2_vec.push_back(0.04);
        omg2_vec.push_back(0.16);
        omg2_vec.push_back(0.64);
    }else{
        kv_vec.push_back(4);
        kv_vec.push_back(16);
        kv_vec.push_back(25);
    }


}*/
void controller::set_default_grid(){
    
    kv_vec.push_back(4);
    kv_vec.push_back(16);
    kv_vec.push_back(25);
   
}
//mainly unchanged
///////////////////////////////////////////////////////////////////////////////////////

/*
void controller::set_default_options(){

    max_size = p;
    // default thread 
    nthread = 1;

    use_dap1 = 0;
    output_all = 0;
    size_select_thresh = 0.01;  // set to big positive numbers to enforce the adaptive stopping rule
    //snp_select_thresh = 0.01;


    cluster_pip_thresh = 0.25; // for output cluster purpose
    priority_msize = 100;
    log10_snp_thresh = 2;
    
    ld_control_thresh = 0.25; // by default  ld control with r^2 = 0.25

}*/
void controller::set_default_options(){

    max_size = n_var;
    size_limit = n_var;
    size_select_thresh = 0.01;  

    cluster_pip_thresh = 0.25; // for output cluster purpose
    priority_msize = 100;
    log10_snp_thresh = 2;
    
    ld_control_thresh = 0.25; // by default  ld control with r^2 = 0.25

}
//okay

///////////////////////////////////////////////////////////////////////////////////////

/*
void controller::set_prior(double pi1){

    if(pi1 > 1 - 1e-3)
        pi1 = 1-1e-3;

    for(int i=0;i<p;i++){
        pi_vec.push_back(pi1);
    }
    return;
}*/
void controller::set_prior(double pi1){

	int i;
    if(pi1 > 1 - 1e-3)
        pi1 = 1-1e-3;

    for(i=0;i<n_var;i++){
        pi_vec.push_back(pi1);
    }
    return;
}
// mainly unchanged
///////////////////////////////////////////////////////////////////////////////////////
/*
double controller::compute_log10_prior(vector<int> &mcfg){

    double lp=0;
    for(int i=0;i<p;i++){
        if(mcfg[i]==0) 	 
            lp += log(1-pi_vec[i]);
        else
            lp += log(pi_vec[i]);
    }

    return lp/log(10);

}*/
double controller::compute_log10_prior(vector<int> &mcfg){

	int i;
    double lp=0;
    for(i=0;i<n_var;i++){
        if(mcfg[i]==0) 	 
            lp += log(1-pi_vec[i]);
        else
            lp += log(pi_vec[i]);
    }

    return lp/log(10);
}
//mainly unchanged
// compared
///////////////////////////////////////////////////////////////////////////////////////

void controller::init(){  
    
	int i;
    mlr.init(Z,R);
	mlr.set_effect_vec(kv_vec);//!/ call by reference not necessary?
    
    null_config = vector<int>(n_var,0);  
    
    double sum = 0;
    for(i=0;i<n_var;i++){
        sum += (pi_vec[i])/(1-pi_vec[i]);
    }

    prior_ratio = sum/n_var;

}
void controller::fine_map(){

	int i;
    init();
	
    // initialize snp2cluster_map
    for(i=0;i<n_var;i++){
        snp2cluster_map[i] = -1;
    }

    vector<double> log10_pmass_vec;

    double curr_val = compute_log10_prior(null_config);
    log10_pmass_vec.push_back(curr_val);

	
	// record best model
    vector<int> bm;
    // single SNP scan, label candidate SNPs for higher order models
    szm_vec.push_back(compute_post_model_single(bm));
    log10_pmass_vec.push_back(szm_vec[0].log10_sum_post);
	
    vector<double> wv(log10_pmass_vec.size(),1.0);
    double val = log10_weighted_sum(log10_pmass_vec,wv);
    int total_snp = szm_vec[0].mvec.size();
    
   
    double prev_val = val;
    int size = 1;
    int no_bc = 0;
    double increment;
	
	
	while(1){

        // if already at allowed max size, quit search
        if(bm.size() == max_size)
            break;
    

        // attempt to add one more variable
        double log10_post = conditional_est(bm);
		

        // if no more variable to add, then quit search
        if(size == bm.size())
            break;
        else
            size++;

        
        // perform backwards checking, see if can eliminate one variable 
        int elim_index = -1;
        if(size >= 3 && !no_bc){
            elim_index = backward_checking(bm, log10_post); 
            if(elim_index !=-1){
                size--;
            }
        }

        size_model szm;

        double cps = szm_vec[szm_vec.size()-1].log10_sum_post;
        // append an equal size model to correct for necessary backward checking
        if(elim_index != -1){
            
            map<int, int> black_list;
            for(i=0;i<szm_vec[elim_index].snp_cluster.size();i++){
                int snp = szm_vec[elim_index].snp_cluster[i];
                black_list[snp] = 100;
            }
            szm = append_post_model(size, black_list); 
            
            szm.snp_cluster = cand_set;
            szm_vec.push_back(szm);
            for(i=0;i<cand_set.size();i++){
                snp2cluster_map[cand_set[i]] = szm_vec.size();
            }
            log10_pmass_vec.push_back(szm.log10_sum_post);      
            // update log10(NC)
            vector<double> nwv(log10_pmass_vec.size(),1.0);
            val = log10_weighted_sum(log10_pmass_vec,nwv);
            total_snp += cand_set.size();
            
            //  stop continuous backward checking for tiny improvement
            if(val - prev_val < size_select_thresh)
                no_bc = 1;  

        }else{
            //allow backward checking if adding a new variable    
            if(no_bc)
                no_bc = 0;

            // append a size 
            int use_abs_cutoff = 0;
            if(increment > 1 || size == 2)
                use_abs_cutoff = 1;


            szm = compute_post_model(size, use_abs_cutoff); 
            szm.snp_cluster = cand_set;
            szm_vec.push_back(szm);
            for(i=0;i<cand_set.size();i++){
                snp2cluster_map[cand_set[i]] = szm_vec.size();
            }
            log10_pmass_vec.push_back(szm.log10_sum_post);      
            // update log10(NC)
            vector<double> nwv(log10_pmass_vec.size(),1.0);
            val = log10_weighted_sum(log10_pmass_vec,nwv);
            total_snp += cand_set.size();
            
		}


		// stopping rule

        int cs = size;
        // if the prior expectation is high for number of signals, carry on

        // expected contribution of next size partition, assuming model is saturated 
        double project_ratio = (n_var-cs+1)*prior_ratio/cs;
        if(project_ratio > 1){
            increment = val - prev_val;
            prev_val = val;
            continue;
        }
        // else check stopping point
        double ncps = szm.log10_sum_post;
        double rb = log10(double(n_var)-cs+1)+log10(prior_ratio) + cps;
        double lb = log10(double(n_var-2*cs+2)/cs) + log10(prior_ratio) + cps;
        if( ncps <= rb && val - prev_val <= size_select_thresh){
            break;
        }

        increment = val - prev_val;
        prev_val = val;
        
    }
    
    
	
}
// mainly unchanged: loop indices, p -> n_var
// compared
size_model controller::compute_post_model(int size, int use_abs_cutoff){


    size_model smod;
    smod.size = size;
    smod.log10_sum_post = 0; 


    vector<vector<int> > & model_vec = szm_vec[szm_vec.size()-1].mvec;
    vector<vector<int> > mc_vec;
	
	int i,j;
    for(i=0; i<model_vec.size(); i++){
        for(j=0;j<cand_set.size();j++){
            vector<int> cm = model_vec[i];
            cm.push_back(cand_set[j]);
            mc_vec.push_back(cm);

        }
    }

    vector<double> post_vec(mc_vec.size());
    int ms = mc_vec.size();


    map<string, double> post_map;
    vector<string> name_vec(ms);

    for(i=0;i<ms;i++){

        // empty sets to start
        vector<int>  mcfg = null_config;
        vector<int> cm = mc_vec[i];

        string name ="";
        for(j=0;j<size;j++){
            int index = cm[j];
            name += geno_map[index]; 
            if(j!=size-1){
                name += "&";
            }

            mcfg[index]=1;
        }



        MLR local_mlr; 

        local_mlr.copy(mlr,mcfg); 
        double log10_abf  = local_mlr.compute_log10_ABF();
        double log10_post =  log10_abf + compute_log10_prior(mcfg);
        
        post_vec[i] = log10_post;
        name_vec[i] = name;
        post_map[name] = log10_post;
		
        
    }


    vector<double> rst_post_vec;
    map<string, double> rst_post_map;
    vector<vector<int> > rst_mc_vec;

    vector<double> post_vec_sort = post_vec;
    sort(post_vec_sort.rbegin(), post_vec_sort.rend());

    double cutoff = post_vec_sort[0] - 2;
    if(!use_abs_cutoff){
        if(post_vec.size()>priority_msize){
            cutoff = post_vec_sort[priority_msize-1];
        }else{
            cutoff = post_vec_sort[post_vec_sort.size()-1];
        }
    }

    for(i=0;i<post_vec.size();i++){
        if(post_vec[i]>= cutoff){
            rst_post_vec.push_back(post_vec[i]);
            rst_mc_vec.push_back(mc_vec[i]);
            rst_post_map[name_vec[i]] = post_vec[i];
        }
    }

    vector<double> wv(rst_post_vec.size(),1.0);
    smod.log10_sum_post = log10_weighted_sum(rst_post_vec,wv); 
    smod.post_map = rst_post_map;
    smod.mvec = rst_mc_vec;
    return smod;

}
// mainly unchanged: loop indices and omp code
// compared



size_model controller::append_post_model(int size, map<int, int> &black_list){

    size_model smod;
    smod.size = size;
    smod.log10_sum_post = 0;

    vector<vector<int> > & model_vec = szm_vec[szm_vec.size()-1].mvec;
    vector<vector<int> > mc_vec;

    map<string, int> saved_model_string_map;
	int i,j,k;
    for(i=0; i<model_vec.size(); i++){
        for(j=0;j<cand_set.size();j++){
            vector<int> cm = model_vec[i];

            for(k=0;k<cm.size();k++){
                if(black_list[cm[k]] == 100){
                    cm.erase(cm.begin()+k);
                    break;
                }
            }
            if(cm.size()==size-1){

                cm.push_back(cand_set[j]);

                stringstream mstream;
                
                for(k=0;k<cm.size();k++){
                    mstream<< cm[k]<<":";
                }

                string model_string = mstream.str(); 
                if(saved_model_string_map[model_string] != 100){;
                    mc_vec.push_back(cm);
                    saved_model_string_map[model_string] = 100;
                }

            }
        }
    }

    vector<double> post_vec(mc_vec.size());
    int ms = mc_vec.size();

    map<string, double> post_map;
    vector<string> name_vec(ms);

    for(i=0;i<ms;i++){

        // empty sets to start
        vector<int>  mcfg = null_config;
        vector<int> cm = mc_vec[i];

        string name ="";
        for(j=0;j<size;j++){
            int index = cm[j];
            name += geno_map[index]; 
            if(j!=size-1){
                name += "&";
            }

            mcfg[index]=1;
        }



        MLR local_mlr; 

        local_mlr.copy(mlr,mcfg);
        double log10_abf  = local_mlr.compute_log10_ABF(); 
        double log10_post =  log10_abf + compute_log10_prior(mcfg);
        
        post_vec[i] = log10_post;
        name_vec[i] = name;
        post_map[name] = log10_post;
		
        
    }


    vector<double> rst_post_vec;
    map<string, double> rst_post_map;
    vector<vector<int> > rst_mc_vec;

    vector<double> post_vec_sort = post_vec;
    sort(post_vec_sort.rbegin(), post_vec_sort.rend());

    double cutoff = post_vec_sort[0] - 2;

    for(i=0;i<post_vec.size();i++){
        if(post_vec[i]>= cutoff){
            rst_post_vec.push_back(post_vec[i]);
            rst_mc_vec.push_back(mc_vec[i]);
            rst_post_map[name_vec[i]] = post_vec[i];
        }
    }

    vector<double> wv(rst_post_vec.size(),1.0);
    smod.log10_sum_post = log10_weighted_sum(rst_post_vec,wv);//!/
    smod.post_map = rst_post_map;
    smod.mvec = rst_mc_vec;

    return smod;

}
// mainly unchanged: loop indices and omp code
//compared

double controller::conditional_est(vector<int>& bm){

    cand_set.clear();

    double max = -9999;

    vector<int>  mcfg = null_config;

	int i;
    for(i=0;i<bm.size();i++){
        int index = bm[i];
        mcfg[index] = 1;
    }
    

    vector<double> post_vec;

    int max_id;
    for(i=0;i<n_var;i++){

        if(cand_map[i]==1){
            post_vec.push_back(-99999);
			
            continue;
        }
        vector<int>  cmcfg = mcfg;
        cmcfg[i] = 1;


        double rst = mlr.compute_log10_ABF(cmcfg) + compute_log10_prior(cmcfg); 
        post_vec.push_back(rst);

        if(rst > max){
            max= rst;
            max_id = i;
        }

    }

    double thresh = -99999;


    if(post_vec.size()>size_limit){
        vector<double> post_vec_sort = post_vec;
        sort(post_vec_sort.rbegin(), post_vec_sort.rend());
        thresh = post_vec_sort[size_limit-1];
    }

    int flag = 0;
    if(max > -9999){
        for(i=0;i<n_var;i++){
            if( max - post_vec[i]  <= log10_snp_thresh  && post_vec[i]>=thresh){
                if(compute_r2(max_id,i)< ld_control_thresh)
                    continue;
                cand_set.push_back(i);
                cand_map[i]=1;
                flag =1;
				
            }

        }
    }

    if(flag == 1)
        bm.push_back(max_id);

    return max;
}
// mainly unchanged: loop indices, p->n_var
// compared

int controller::backward_checking(vector<int>& bm, double log10_post){


    double max = -99999;
    int max_snp = -1;
    int max_id = -1;

	int i,j;
    for(i=0; i<bm.size()-2; i++){
        vector<int> tm = bm;
        tm.erase(tm.begin()+i);

        vector<int> mcfg = null_config;

        for(j=0;j<tm.size();j++){
            int index = tm[j];
            mcfg[index] = 1;
        }

        double rst =  mlr.compute_log10_ABF(mcfg) + compute_log10_prior(mcfg);
        if(rst > max){
            max= rst;
            max_snp = bm[i];
            max_id=i;
        }

    }

    if(max-log10_post>0){
        bm.erase(bm.begin()+max_id);
        return snp2cluster_map[max_snp];
    }


    return -1;

}
// mainly unchanged: loop indices
// compared 

size_model controller::compute_post_model_single(vector<int>& bm){


    size_model smod;
    smod.size = 1;
    smod.log10_sum_post = 0;

    vector<double> post_vec;
    vector<double> abf_vec;

    double max_log10_abf = -9999;
    double max_log10_post = -9999;
    int max_id = -1;
    int i,index;
	
    for(index=0;index<n_var;index++){

        string name = geno_map[index]; 
        vector<int>  mcfg = null_config;
        mcfg[index]=1;
        double log10_abf = mlr.compute_log10_ABF(mcfg); 
        single_log10_abfv[name] = log10_abf;
		single_z_score[name] = mlr.get_z_score(index+1);
       
        
        abf_vec.push_back(log10_abf);
        double log10_post =  log10_abf + compute_log10_prior(mcfg); 
        post_vec.push_back(log10_post);

        if(log10_post>max_log10_post){
            max_log10_post = log10_post;
            max_id = index;
        }
        smod.post_map[name] = log10_post;
    }


    vector<double> wv(post_vec.size(),1.0);
    smod.log10_sum_post = log10_weighted_sum(post_vec,wv); 

    double thresh = -99999;

    if(post_vec.size()>size_limit){
        vector<double> post_vec_sort = post_vec;
        sort(post_vec_sort.rbegin(), post_vec_sort.rend());
        thresh= post_vec_sort[size_limit-1];
    }


    for(i=0;i<n_var;i++){
        if(max_log10_post - post_vec[i]  <= log10_snp_thresh && post_vec[i]>= thresh) { 
            if(compute_r2(max_id,i)< ld_control_thresh)
                continue;
            vector<int> mv;
            mv.push_back(i);
            smod.mvec.push_back(mv);
            cand_set.push_back(i);
            cand_map[i] = 1;
            snp2cluster_map[i]=0;
        }
    }
    

    smod.snp_cluster = cand_set;
    bm.push_back(max_id);
    return smod;

}
// mainly unchanged: loop indices, single_bhat, single_se, p->n_var
// compared
///////////////////////////////////////////////////////////////////////////////////////

/*
double controller::compute_average_r2(const vector<int> & vec){
    if(vec.size()==1)
        return 1.0;
    double sum_r2 = 0;
    int count = 0;
    for(int i=0;i<vec.size();i++){
        for(int j=i+1;j<vec.size();j++){
            sum_r2 += compute_r2(vec[i],vec[j]);
            count++;
        }
    }
    return sum_r2/count;
}*/
double controller::compute_average_r2(const vector<int> & vec){
    if(vec.size()==1)
        return 1.0;
    double sum_r2 = 0;
    int count = 0;
	int i,j;
    for(i=0;i<vec.size();i++){
        for(j=i+1;j<vec.size();j++){
            sum_r2 += compute_r2(vec[i],vec[j]);
            count++;
        }
    }
    return sum_r2/count;
}
//slight modificiation, mainly unchanged

///////////////////////////////////////////////////////////////////////////////////////

/*
double controller::compute_average_r2(const vector<int> & vec1, const vector<int> & vec2){

    double sum_r2 = 0;
    int count = 0;
    for(int i=0;i<vec1.size();i++){
        for(int j=0;j<vec2.size();j++){
            sum_r2 += compute_r2(vec1[i],vec2[j]);
            count++;
        }
    }
    return sum_r2/count;
}*/
double controller::compute_average_r2(const vector<int> & vec1, const vector<int> & vec2){

    double sum_r2 = 0;
    int count = 0;
	int i,j;
    for(i=0;i<vec1.size();i++){
        for(j=0;j<vec2.size();j++){
            sum_r2 += compute_r2(vec1[i],vec2[j]);
            count++;
        }
    }
    return sum_r2/count;
}
//slight modificiation, mainly unchanged

///////////////////////////////////////////////////////////////////////////////////////

/*double controller::compute_r2(int i, int j){
    double r2;
    if(use_ss==0){
        double *gi = &pars.geno_vec[0][i][0];
        double *gj = &pars.geno_vec[0][j][0];
        r2 = pow(gsl_stats_correlation(gi,1,gj, 1, pars.geno_vec[0][i].size()),2);
    }else{
        r2 = pow(gsl_matrix_get(pars.ld_matrix,i,j),2);
    }
    return r2;
}*/
double controller::compute_r2(int i, int j){
    double r2;
    r2 = pow((*R)(i+1,j+1),2.0);
    return r2;
}
//extracted use_ss=1 part and replaced gsl operation
///////////////////////////////////////////////////////////////////////////////////////


finemap_results controller::summarize_approx_posterior(bool verbose=false){

    vector<double>  log10_pmass_vec;
    log10_pmass_vec.push_back(compute_log10_prior(null_config));

	int i,j,k;
    for(i=0;i<szm_vec.size();i++){
        log10_pmass_vec.push_back(szm_vec[i].log10_sum_post);
    }

    vector<double> wv1(log10_pmass_vec.size(),1.0);
    double log10_pip_NC = log10_weighted_sum(log10_pmass_vec,wv1);

    double val = log10_pmass_vec[log10_pmass_vec.size()-1];

    double sum = 0;
    double ratio = 1;
    int max_model_size = szm_vec[szm_vec.size()-1].size;
    if(max_size == n_var){
        for(k=max_model_size+1;k<=n_var;k++){
            ratio *= (n_var-k+1)*prior_ratio/k;
            sum += ratio;
        }

        log10_pmass_vec.push_back(val+log10(sum));
    }

    vector<double> wv(log10_pmass_vec.size(),1.0);
    log10_pnorm  = log10_weighted_sum(log10_pmass_vec,wv);


    Nmodel nm;
    nm.id = "NULL";
    nm.prob = pow(10, log10_pmass_vec[0]-log10_pnorm);
    double null_prob = nm.prob;
    nm.size = 0;
    nm.post_score = log10_pmass_vec[0];
    nmodel_vec.push_back(nm);

    for(i=0;i<szm_vec.size();i++){
        
        map<string, double>::iterator iter;
        for(iter = szm_vec[i].post_map.begin(); iter != szm_vec[i].post_map.end(); iter++){
            Nmodel nm;
            nm.id = iter->first;
            nm.post_score = iter->second;
            nm.prob = pow(10, iter->second-log10_pnorm);
            nm.size = szm_vec[i].size;
            nmodel_vec.push_back(nm);
            parse_nmodel(nm);
        }
    }

    std::sort(nmodel_vec.begin(),nmodel_vec.end(),sort_nmodel_dec);
    double cump = 0;

    double msize_mean = 0;
    double msize_var = 0;


    for(i=0;i<nmodel_vec.size();i++){

        string name = nmodel_vec[i].id;

        size_t start_pos = 0;
        while((start_pos = name.find("&", start_pos)) != std::string::npos) {
            name.replace(start_pos, 1 , "] [");
            start_pos += 3; 
        }

        //if(nmodel_vec[i].prob >= 1e-5 || name == "NULL"){ //!/
            //std::fprintf(outfd, "%5d   %7.4e    %d    %7.3f   [%s]\n",i+1, nmodel_vec[i].prob, nmodel_vec[i].size, nmodel_vec[i].post_score, name.c_str());
        //}
        cump += nmodel_vec[i].prob;
        msize_mean += nmodel_vec[i].prob*nmodel_vec[i].size;
        msize_var  += nmodel_vec[i].prob*pow(nmodel_vec[i].size,2);

    }



    msize_var -= pow(msize_mean,2.0);
    if(msize_var < 0){
        msize_var = 0;
    }



    vector<NSNP> nsnp_vec_sort = nsnp_vec;
    std::sort(nsnp_vec_sort.begin(),nsnp_vec_sort.end(),sort_nsnp_dec_by_ip);

    map<string, int> snp2index;
    for(i=0;i<nsnp_vec_sort.size();i++){
        
        nsnp_vec_sort[i].cluster = -1;
        snp2index[nsnp_vec_sort[i].name] = i;
    }
    // estimate min_pip from BIC approximation
    int NN = 1000; //!/ needed?
    
    double min_pip = (1-null_prob)*prior_ratio/sqrt(NN);
    

    vector<double> cluster_pip;
    vector<double> cluster_r2;
    vector<int> cluster_count;
    vector<int> cluster_id;
    vector<vector<int> > grp_vec;
    map<int,int> grpr2_map;
    int cluster_index = 1;
    for(i=0;i<szm_vec.size();i++){

        double cluster_prob = 0;
        vector<int> member_vec;

        for(j=0;j<szm_vec[i].snp_cluster.size();j++){
            int snp = szm_vec[i].snp_cluster[j];
            string sname = geno_map[szm_vec[i].snp_cluster[j]];
            if(snp2index.find(sname) == snp2index.end())
                continue;
            int index = snp2index[sname];
            double prob = nsnp_vec_sort[index].incl_prob;  
            // we don't consider a SNP as a noteworthy cluster memeber if its pip < min_pip
            
            if(prob > min_pip){ //!/ deleted use_dap1
                member_vec.push_back(snp);
                cluster_prob += prob;
                nsnp_vec_sort[index].cluster = cluster_index;
            }
            
        }
        if(member_vec.size()>0){
            
            cluster_r2.push_back(compute_average_r2(member_vec));
            cluster_count.push_back(member_vec.size());
            cluster_pip.push_back(cluster_prob);
            cluster_id.push_back(cluster_index);

            if(cluster_prob >= cluster_pip_thresh){
                grp_vec.push_back(member_vec);
                grpr2_map[cluster_count.size()-1] = grp_vec.size()-1;
            }
            cluster_index++;
        }

    }       
	finemap_results fr;
	vector<string> variant_names;
	vector<double> incl_probs;
	vector<double> log_bfs;
	vector<double> z_scores;
	vector<int> cluster_ids;
	//!/
	
    for(i=0;i<nsnp_vec_sort.size();i++){
        if(nsnp_vec_sort[i].incl_prob < min_pip) //!/
            nsnp_vec_sort[i].incl_prob = min_pip;
       
        if(nsnp_vec_sort[i].cluster==-1 && !verbose)//!/
            continue; //!/ take out bhat and se
		//!/
		variant_names.push_back(nsnp_vec_sort[i].name.c_str());
		incl_probs.push_back(nsnp_vec_sort[i].incl_prob);
		log_bfs.push_back(single_log10_abfv[nsnp_vec_sort[i].name]);
		z_scores.push_back(single_z_score[nsnp_vec_sort[i].name]);
		cluster_ids.push_back(nsnp_vec_sort[i].cluster);
		
    }
	
	vector<int> cluster_counts;
	vector<double> cluster_pips;
	vector<double> cluster_r2s;
	if(grp_vec.size()>0){
      
        for(i=0;i<cluster_count.size();i++){
            if(cluster_pip[i]<cluster_pip_thresh)
                continue;
			cluster_counts.push_back(cluster_count[i]);
			cluster_pips.push_back(cluster_pip[i]);
			cluster_r2s.push_back(cluster_r2[i]);
		}
    }
	
	
	fr.variant_names=variant_names;
	fr.incl_probs=incl_probs;
	fr.log_bfs=log_bfs;
	fr.z_scores=z_scores;
	fr.cluster_ids=cluster_ids;
	
	fr.cluster_counts=cluster_counts;
	fr.cluster_pips=cluster_pips;
	fr.cluster_r2s=cluster_r2s;
	return fr;
    
}

///////////////////////////////////////////////////////////////////////////////////////

/*
void controller::parse_nmodel(Nmodel nmod){

    istringstream iss(nmod.id);
    string token;
    while (getline(iss, token, '&')){
        string snp_id = token;
        int index;
        if(nsnp_map.find(snp_id)!=nsnp_map.end()){
            index = nsnp_map[snp_id];
        }else{
            NSNP ns;
            ns.name = snp_id;
            ns.incl_prob = 0;
            nsnp_vec.push_back(ns);
            index = nsnp_vec.size()-1;
            nsnp_map[snp_id] = index;
        }
        nsnp_vec[index].incl_prob += nmod.prob;

    }
}*/
void controller::parse_nmodel(Nmodel nmod){

    istringstream iss(nmod.id);
    string token;
    while (getline(iss, token, '&')){
        string snp_id = token;
        int index;
        if(nsnp_map.find(snp_id)!=nsnp_map.end()){
            index = nsnp_map[snp_id];
        }else{
            NSNP ns;
            ns.name = snp_id;
            ns.incl_prob = 0;
            nsnp_vec.push_back(ns);
            index = nsnp_vec.size()-1;
            nsnp_map[snp_id] = index;
        }
        nsnp_vec[index].incl_prob += nmod.prob;

    }
}
//unchanged
// compared
///////////////////////////////////////////////////////////////////////////////////////

/*bool sort_nsnp_dec_by_ip(const NSNP &lhs, const NSNP &rhs){
    return lhs.incl_prob > rhs.incl_prob;
}*/
bool sort_nsnp_dec_by_ip(const NSNP &lhs, const NSNP &rhs){
    return lhs.incl_prob > rhs.incl_prob;
}
//unchanged
///////////////////////////////////////////////////////////////////////////////////////


/*bool sort_nmodel_dec(const Nmodel &lhs, const Nmodel &rhs){
    return lhs.prob > rhs.prob;
}*/
bool sort_nmodel_dec(const Nmodel &lhs, const Nmodel &rhs){
    return lhs.prob > rhs.prob;
}
//unchanged



