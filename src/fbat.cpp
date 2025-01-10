#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string> //!/
#include <map> //!/
#include "fbat.h"
#include "mystream.h"
#include "myerror.h"
#include "myutil.h"
#include "normbase.h"
#include "betabase.h"
#include "gammabas.h"
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"
#include "sufficient_stat.h"
#include "rand.h"
#include "dapg_fbat.h" //!/

const char *model_name[] = { "additive", "dominant", "recessive", "genotype", "cmm_estimate", "user defined" };
const char *mode_name[] = { "bi-allelic", "multi-allelic", "both" };
const double very_small_number = 1.0e-6;
double another_very_small_number = 1.0e-4;
double extreme_small_number = 1.0e-10; //!/ relative error

void FBAT::init() {
	
	char title[]="\n\
               *******************************************************\n\
               *                                                     *\n\
               *     *********  * * *          *       *********     *\n\
               *     *          *     *       * *          *         *\n\
               *     *******    *  * *       *   *         *         *\n\
               *     *          *     *     * *** *        *         *\n\
               *     *          *     *    *       *       *         *\n\
               *     *          * * *     *         *      *         *\n\
               *                                                     *\n\
               *          v2.0.9 	     			     *\n\
               *          Program for Population Genetics            *\n\
               *          Harvard School of Public Health            *\n\
               *                                                     *\n\
               *******************************************************\n\
               \n";
	
	static COMMAND cmds[] = {
 		{	"model", 
 			"model [a,d,r] -- select genetic models; (a)dditive, (d)ominant, (r)ecessive",
  			static_cast<funcp> (&FBAT::select_model)
 		},
 		{	"mode", 
 			"mode [a, b, m] -- select bi-allelic (b), multi-allelic (m) test, or both (a)",
  			static_cast<funcp> (&FBAT::select_mode)
 		},
 		{	"displayp", 
 			"displayp [p_value] -- report results with P<=p_value only",
  			static_cast<funcp> (&FBAT::select_alpha)
 		},
 		{	"minsize", 
 			"minsize [count] -- set minimum number of informative families to test",
  			static_cast<funcp> (&FBAT::select_minsize)
 		},
 		{	"minfreq", 
 			"minfreq [count] -- set minimum allele frequency to be included in test",
  			static_cast<funcp> (&FBAT::select_minfreq)
 		},
 		{	"offset",
 			"offset [value] -- set the trait offset",
  			static_cast<funcp> (&FBAT::select_offset)
 		},
 		{	"maxcmh", 
 			"maxcmh [count] -- set the maximum allowable # of compatible mating haplotypes",
  			static_cast<funcp> (&FBAT::select_maxcmh)
 		},
 		{	"rvminsize", 
 			"rvminsize [count] -- set minimum number of informative families to test for rare variant",
  			static_cast<funcp> (&FBAT::select_rvminsize)
 		},
  		{	"viewhap", 
 			"viewhap [-c] [-e] [-i] [-s] [marker...] -- view markers' haplotype configuration in each family\n     -c censored trait\n     -e empirical test\n     -i show informative families only\n     -s summary only", // removed -i for the offcial release
  			static_cast<funcp> (&FBAT::hapview_ss)
 		},
 		{	"hbat", 
 			"hbat [-c] [-e] [-o] [-p[#]] [marker...] -- haplotype version of FBAT test.\n     -c censored trait\n     -e empirical test\n     -o optimize offset\n     -p[#] permutation test, #=maximal number of runs",	// removed -w and -i 2/14/2004
  			static_cast<funcp> (&FBAT::hbat)
 		},
		{	"fbat", 
 			"fbat [-e] [-o] [-c] [-l] [-m] [-p#] [-r] [-s] [-t] [v#] [marker...] -- FBAT test.\n     -c censored trait\n     -e empirical test\n     -m FBAT_MM test\n     -l FBAT_LC test for multiple markers\n     -o optimize offset\n     -p# FBAT_minP test; # is permutation cycles (default = 10,000)\n     -s generate intermediate results\n     -t LC test for multiple traits\n     -v# FBAT test for rare variant; # 0 for unweighted test and 1 for weighted test",
  			static_cast<funcp> (&FBAT::fbat)
 		},
		{	"setafftrait", 
 			"setafftrait aff_t unaff_t unknown_t -- set trait value for affection status",
  			static_cast<funcp> (&FBAT::set_affection_trait)
		},
		{	"hapfreq", 
 			"hapfreq [-d] [-r] [markers...]  --  estimate haplotype frequency of specified markers, \n     -d: output LD(D')\n     -r: output LD(r^2)",
  			static_cast<funcp> (&FBAT::haplotypefrequency)
		}, //!/
		{	"nsimreg", //!/
 			"nsimreg [value] -- set the number of simulations for region-based analysis", //!/
  			static_cast<funcp> (&FBAT::select_nsim_region) //!/
 		}, //!/
		{	"maxsimreg", //!/
 			"maxsimreg [value] -- set the maximum number of simulations for region-based analysis", //!/
  			static_cast<funcp> (&FBAT::select_maxsim_region) //!/
 		}, //!/
		{	"gen_rv",  //!/
 			"gen_rv [marker...] -- perform rare variant testing for a region", //!/
  			static_cast<funcp> (&FBAT::gen_rvtest) //!/
		}, //!/
		{	"gen_cond_rv",  //!/
 			"gen_cond_rv", //!/
  			static_cast<funcp> (&FBAT::gen_cond_rvtest) //!/
		}, //!/
		{	"finemap",  //!/
 			"finemap [marker] -- perform DAP-G fine-mapping", //!/
  			static_cast<funcp> (&FBAT::finemap) //!/
		} //!/
/*		{	"pbat", 
 			"pbat [-o] [-c] [-i] [marker...] -- FBAT test.\n     -c censored trait\n     -i Use conditional mean medel\n     -o optimize offset",
  			static_cast<funcp> (&FBAT::pbat)
 		}
 		{	"viewstat", 
 			"viewstat [-e] [-o] [-c] [-s] marker -- view statistics.\n     -c censored trait\n     -e empirical test\n     -o optimize offset\n     -s summary only",
  			static_cast<funcp> (&FBAT::view_stat)
 		},
 		{	"sdt", 
 			"sdt [marker] -- sibship disequilibrium test",
  			static_cast<funcp> (&FBAT::sdt)
 		},
 		{	"transmit",
 			"transmit -- count of parental allele transmissions",
 			static_cast<funcp> (&FBAT::tdt)
		},
 		{	"fbatall", 
 			"fbatall  -- FBAT test for all traits, markers, genetic models and mode",
  			(funcp)&fbat_all
 		}
 		*/
	};

	addcmd(cmds, 18);
	
	MYOUTPUT(os, title)
	
}

void FBAT::select_model(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64];
		GetWord(s, 2, word);
		if (word[1]==0)
			switch(word[0]) {
				case 'a': gen_model = model_additive; break;
				case 'd': gen_model = model_dominant; break;
				case 'r': gen_model = model_recessive; break;
				case 'g': 
					MYOUTPUT(os, "The genotype model is not implemented\n");
					//gen_model = model_genotype;
					break;
				case 'e': gen_model = model_estimate; break;
				default:
					MYOUTPUT(os, "The user defined mode has not be implemented\n");
			}
	}
		
	MYOUTPUT(os, "current genetic model is " << model_name[gen_model] << "\n")
	
}

void FBAT::select_mode(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64];
		if (GetWord(s, 2, word))
			if (word[1]==0) {
				switch(word[0]) {
					case 'b': mode = by_allele; break;
					case 'm': mode = by_marker; break;
					case 'a': mode = by_both; break;
					default:
						throw Param_Error("unknown mode options");
				}
			} else
				throw Param_Error("unknown mode options");
	}
		
	MYOUTPUT(os, "current test mode is " << mode_name[mode] << "\n")
	
}

void FBAT::select_alpha(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		double dd = strtod(word, &c);
		if (c==word || dd<=0.0 || dd>1)
			throw ERROR("p_value must be number between (0..1]");
		alpha = dd;
	}
	
	MYOUTPUT(os, "current p_value is " << alpha << "\n")
}

void FBAT::select_offset(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		double dd = strtod(word, &c);
		if (c==word)
			throw ERROR("offset must be a number");
		tr_offset = dd;
	}
	
	MYOUTPUT(os, "current trait offset is " << setprecision(6) << tr_offset << "\n")
}

void FBAT::select_minsize(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		int cnt = strtol(word, &c, 10);
		if (c==word || cnt<0)
			throw ERROR("minimum informative family count must be a positive integer");
		min_size = cnt;
	}
	
	MYOUTPUT(os, "current minimum informative family count is " << min_size << "\n")
}

void FBAT::select_rvminsize(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		int cnt = strtol(word, &c, 10);
		if (c==word || cnt<0)
			throw ERROR("minimum informative family count for rare variant must be a positive integer");
		rv_min_size = cnt;
	}
	
	MYOUTPUT(os, "current minimum informative family count for rare varriant  is " << rv_min_size << "\n")
}

void FBAT::select_minfreq(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		double frq = strtod(word, &c);
		if (c==word || frq<0 || frq>=1)
			throw ERROR("minimum allele frequency must be within [0,1)");
		min_freq = frq;
	}
	
	MYOUTPUT(os, "minimum allele frequency to be included in test is " << min_freq << "\n")
}

void FBAT::select_maxcmh(char *s) {
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		int cnt = strtol(word, &c, 10);
		if (c==word || cnt<0)
			throw ERROR("minimum informative family count must be a positive integer");
		maxcmh = cnt;
	}
	
	MYOUTPUT(os, "current maximum allowable compatible mating haplotypes # is " << maxcmh << "\n")
}

void FBAT::print_settings() {
	
	if (n_sel_trt == 0) {
		n_sel_trt = 1;
		trait_id[0] = 0;
	}
	if (n_sel_trt>1) {
		MYOUTPUT(os, n_sel_trt << " traits selected: ")
		for (int i=0; i<n_sel_trt; i++)
			MYOUTPUT(os, trt[trait_id[i]].name << " ")
		MYOUTPUT(os, "\n")
	} else
		MYOUTPUT(os, "trait " << trt[trait_id[0]].name << "; offset " << setprecision(3) << tr_offset << "; ")
	MYOUTPUT(os, "model " << model_name[gen_model])
	MYOUTPUT(os, "; test " << mode_name[mode] << "; minsize " << min_size << "; min_freq " << min_freq << "; p " << alpha)
	MYOUTPUT(os, "; maxcmh " << maxcmh << "\n")
	
}
void FBAT::select_maxsim_region(char *s) { //!/
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		double dd = strtod(word, &c);
		if (c==word)
			throw ERROR("must be a number");
		unsigned long tmp = (unsigned long)dd;
		if(tmp>0 && tmp<=1000000000){ max_sim_region=tmp;}
		else{ max_sim_region=num_sim_region=10000; MYOUTPUT(os, "invalid input, setting max_sim_region and num_sim_region to 10000\n");}
		
	}
	if(num_sim_region>max_sim_region){
		MYOUTPUT(os, "WARNING: num_sim_region is greater than max_sim_region, setting max_sim_region = num_sim_region= " << num_sim_region<<"\n")
		max_sim_region=num_sim_region;
	}
	MYOUTPUT(os, "current maximum number of simulations for region-based tests is " << setprecision(10) << max_sim_region << ". num_sim_region is "<<num_sim_region<<"\n")
} //!/
void FBAT::select_nsim_region(char *s) { //!/
	
	int wcnt = WordCount(s);
	
	if (wcnt!=2 && wcnt!=1)
		throw Param_Error("wrong number of arguments");
	
	if (wcnt==2) {
		char word[64], *c;
		GetWord(s, 2, word);
		double dd = strtod(word, &c);
		if (c==word)
			throw ERROR("must be a number");
		unsigned long tmp = (unsigned long)dd;
		if(tmp>0 && tmp<=1000000000){ num_sim_region=tmp;}
		else{ num_sim_region=max_sim_region=10000; MYOUTPUT(os, "invalid input, setting num_sim_region and max_sim_region to 10000\n");}
		
	}
	
	MYOUTPUT(os, "current number of simulations for region-based tests is " << setprecision(10) << num_sim_region << ". max_sim_region is "<<max_sim_region<<"\n")
	if(num_sim_region>max_sim_region){
		MYOUTPUT(os, "WARNING: num_sim_region is greater than max_sim_region, setting max_sim_region = num_sim_region= " << num_sim_region<<"\n")
		max_sim_region=num_sim_region;
	}
} //!/
//!/
int FBAT::draw_index_pat(double* probs,int len) // draw index between 1 and len, assumes that sum(probs)=1
{
		double alpha=rand_uniform();
		int index_ctr=1;
		double cum_prob=0;
		while(cum_prob<alpha)
		{
			cum_prob+=probs[index_ctr];
			index_ctr++;
		}
		
		return index_ctr-1;
} //!/
int* FBAT::do_rand_sim_fam(Sufficient_Stat* ss) //!/
{
	Off_Genotype_Pattern *pat;
	int j,i;
	int count_pats=0;
	int ctr;
	int count_inds=0;
	for (pat=ss->ogpat; pat!=NULL; pat=pat->next){ count_pats++;} // count number of pats
	if(count_pats<1){ MYOUTPUT(os,"error count_pats\n"); exit(1);}
	double* probs=new double[count_pats+1];
	count_pats=0;
	for (pat=ss->ogpat; pat!=NULL; pat=pat->next) // extract corresponding probabilities
	{ 
		count_pats++;
		probs[count_pats]=pat->pg;
		
	}
	int pat_index=draw_index_pat(probs,count_pats); // choose random pat according to probabilities
	if(pat_index==0){MYOUTPUT(os,"error pat_index \n"); exit(1);}
	delete[] probs;
	count_pats=0;
	int* config=NULL;
	for (pat=ss->ogpat; pat!=NULL; pat=pat->next)
	{
		 count_pats++;
		 if(count_pats==pat_index)
		 { 
	        for (j=0; j<ss->n_common_ogt; j++) { count_inds+= pat->cnt[j];}
			config=new int[count_inds+1]; ctr=0;
			for (j=0; j<ss->n_common_ogt; j++) 
			{ 
		        for(i=1;i<=pat->cnt[j];i++)
				{
						ctr++;
						config[ctr]=j+1;
				}
			}
			break;
			
		 }
	}
	if(config==NULL){MYOUTPUT(os,"config\n"); exit(1);}
	int temp;
	for (i = count_inds; i >= 1; --i)
	{
		j =  (int)(rand_uniform()*i)+1;
		temp = config[i];
		config[i] = config[j];
		config[j] = temp;
	}
	
	return config;
}	//!/
bool FBAT::check_consistent_cond(Sufficient_Stat* ss, int index, int* random_index)
{
	int i;
	int tmp_counter=0;
	bool equal=true;
	for (i=0; i<ss->n_offspring; i++) 
	{
		if (ss->otr[i]!=NULL && ss->otr[i][trait_id[0]]<kVeryLargeNumber && (ss->which_comm_ogt(ss->ogt[i],ss->sex[i]))>=0 )
		{
			    tmp_counter++;
				if(ss->common_ogt[random_index[tmp_counter]-1][index].allele_cnt(1)!=ss->common_ogt[ss->which_comm_ogt(ss->ogt[i],ss->sex[i])][index].allele_cnt(1)){equal=false;}			
		}
	 }
	 return equal;
}

void FBAT::gen_cond_rvtest(char *s) {//!/
	
	char word[256], *c;
	int i, j, k, l;
	i = WordCount(s);
	double offset = tr_offset;
	int index; ///*/!!
	int nsel_n; ///*/!!
	//---------------------------------------------
	for (j=2; j<=i && GetWord(s, j, word) && *word=='-'; j++){ //!/ update 01/20/2021
		char msg[256]; //!/ update 01/20/2021
		sprintf(msg, "no options implemented for gen_rv test"); //!/ update 01/20/2021
		throw Param_Error(msg); //!/ update 01/20/2021
	} //!/ update 01/20/2021
	int nsel = j>=i? nLoci : i-j+1;
	int *idx = new int[nsel];
	if (j>i)
		for (nsel=0; nsel<nLoci; nsel++)
			idx[nsel] = nsel;
	else {
		for (nsel=0; j<=i; j++) {
			GetWord(s, j, word);
			if ((k=find_locus(word))>=0)
				idx[nsel++] = k;
			else {
				if (idx!=NULL)
					delete[] idx;
				char msg[256];
				sprintf(msg, "Marker %s not found", word);
				throw Param_Error(msg);
			}
		}
	}
	//--------------------------------------------------
	double gx[kMaxAlleles][3];
	if (gen_model==model_additive || gen_model==model_dominant || gen_model==model_recessive)
		for (i=0; i<kMaxAlleles; i++)
			gcode(&gx[i][0], gen_model);
	else
		throw Param_Error("please select additive, dominant or recessive genetic model.");
		
	Haplotype_List *emhl = hapfreq_ss(nsel, idx);
	if (emhl==NULL) 
		return;
    if (n_sel_trt>1) //!/ update 01/20/2021
		throw Param_Error("multivariate gen_rv test is not available"); //!/ update 01/20/2021
	// set trait to affection status if no other trait selected yet
	if (n_sel_trt == 0) {
		n_sel_trt = 1;
		trait_id[0] = 0;
 	}
    //------------------------------------------------------------	
	index= nsel - 1; ///*/!!
	nsel_n= nsel -1; ///*/!!
	
    MYOUTPUT(os, "region-based test for the "<< nsel_n << " markers between: ") ///*/!!
	MYOUTPUT(os, loc[idx[0]].name << " and " << loc[idx[nsel_n-1]].name<<" ") ///*/!!

	
	HaploidTable *htab = new HaploidTable(loc, nsel, idx, emhl);
	ColumnVector EX(nsel);
	ColumnVector EXtot(nsel);
	ColumnVector U(nsel);
	ColumnVector G(nsel);
	Matrix V(nsel,nsel);
	ColumnVector Utot(nsel);
	Matrix Vtot(nsel,nsel);
	
	Sufficient_Stat* ss;
	Sufficient_Stat** ss_array;

	Utot=0.0;
	Vtot=0.0;
	EXtot=0.0;
	
	
	int ctr=0;
	int inf_ctr=0;
	
	
	
	double t;
	double alpha;
	double beta;
	int* random_index=NULL; //!=!  ///*/!!
	int tmp_counter; //!=!
	
	ColumnVector weights(nsel);
	for(i=1;i<=nsel;i++) weights(i)=1.0;
	
	ColumnVector var_indicator(nsel);
	for(i=1;i<=nsel;i++) var_indicator(i)=0.0;
	
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node && flist->node->stat) { ctr++;}
		}
	}
	ss_array=new Sufficient_Stat*[ctr];
	
	ctr=0;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node && flist->node->stat) {
				ctr++;
				flist->node->stat->analyze(emhl);
				if(flist->node->stat->mhap!=NULL && flist->node->stat->mhap->len()>maxcmh){} 
			    else if(flist->node->stat->mhap==NULL){}
				else if(flist->node->stat->n_common_ogt<2){}
				else if(flist->node->stat->ogpat!=NULL) {
					inf_ctr++;
					ss=flist->node->stat;
					ss_array[inf_ctr]=ss;
					alpha=0.0;
					beta=0.0;
					ColumnVector P(ss->n_common_ogt);
					for (i=0; i<ss->n_common_ogt; i++)
						P(i+1) = ss->pg[i];
					DiagonalMatrix Q(ss->n_common_ogt);
					for (i=0; i<ss->n_common_ogt; i++)
						Q(i+1,i+1) = P(i+1);
					SymmetricMatrix PP(ss->n_common_ogt);
					for (i=0; i<ss->n_common_ogt; i++)
						for (j=i; j<ss->n_common_ogt; j++)
							PP(i+1,j+1) = ss->pgg[i][j];
							
					Matrix X(nsel,ss->n_common_ogt);
					for (i=0; i<nsel; i++){
						for (j=0; j<ss->n_common_ogt; j++){
							X(i+1,j+1) = gx[0][ss->common_ogt[j][i].allele_cnt(1)];//get_coding(ss->common_ogt[j][i].allele_cnt(2), gen_model);
						}
					}
					EX = X * P;
					U=0.0;
					V=0.0;
					for (i=0; i<ss->n_offspring; i++) {
						if (ss->otr[i]!=NULL && ss->otr[i][trait_id[0]]<kVeryLargeNumber && (j=ss->which_comm_ogt(ss->ogt[i],ss->sex[i]))>=0 ){
							t = ss->otr[i][trait_id[0]] - offset;
							beta += t;
							alpha += t*t;
							for(k=0;k<ss->nMarkers;k++){
							  G(k+1)= gx[0][ss->common_ogt[j][k].allele_cnt(1)];//get_coding(ss->common_ogt[j][k].allele_cnt(2), gen_model);
							}
							U+= t*(G-EX);
							EXtot+=t*EX;
						}
					}
					beta=beta*beta-alpha;
					V=alpha* X * Q * X.t()+beta* X * PP * X.t()- (alpha+beta)*X* P*P.t()*X.t(); // ok for missing phenotypes: nov29_2019: missing phenotype treated as t=0. 	
                    // 	X * PP * X.t() measures E[x_1 x2], the mean over choosing two genotypes randomly from the null, without ordering			
					for(i=1;i<=nsel;i++){ if(V(i,i)>very_small_number){ var_indicator(i)+=1.0;} }
					Utot+=U;
					Vtot+=V;
				
				}
			}
		}
	}
	
	double stat_max, stat_hc, stat_skat, stat_burden;
	
	int ctr_max, ctr_hc, ctr_skat, ctr_burden;
	ctr_max=ctr_hc=ctr_skat=ctr_burden=0;
	
	double stat, stat_tmp, stat_tmp1, stat_tmp2, stat_tmp3, stat_tmp4, tmp, ACAT_stat;
	double p_burden, p_skat, p_maxsv, p_hc, p_acat;
	
	// max-sv
	stat=0.0;
	for(j=1;j<=nsel_n;j++) ///*/!!
	{
		if(Vtot(j,j)>=very_small_number){
			tmp=fabs(Utot(j))/sqrt(Vtot(j,j));
			if(tmp>=stat){ stat=tmp;}
		}
	}
	stat_max=stat;
	// highcrit
	stat=0.0;
	ColumnVector q_vec(nsel_n); ///*/
	q_vec= 0.0;
	double x;
	double phi1;
	double phi2;
	int inf_var=0;
	for(i=1;i<=nsel_n;i++){ ///*/!!
		if(var_indicator(i)>=(double)min_size){//if(Vtot(i,i)>very_small_number){//!/ maybe implement option here ?
			inf_var++;
			x=fabs(Utot(i))/sqrt(Vtot(i,i));
			normbase(&x,&phi1);
			x=-fabs(Utot(i))/sqrt(Vtot(i,i));
			normbase(&x,&phi2);
			q_vec(inf_var)=fmin(1-phi1+phi2,1.0-very_small_number);
		}
	}
    
    for (i = 1; i <= inf_var-1; i++){       
		for (j = 1; j <= inf_var-i; j++)  {
				if (q_vec(j) > q_vec(j+1)) {
					tmp = q_vec(j+1); 
					q_vec(j+1) = q_vec(j); 
					q_vec(j) = tmp; 
		   }
        }			  
    }
	if(inf_var>0){
		stat=((double)(1/((double)inf_var))-q_vec(1))/sqrt(q_vec(1)*(1-q_vec(1)));
		for (i = 1; i <= (int)inf_var; i++){  //!/
		   tmp=((double)((double)i/((double)inf_var))-q_vec(i))/sqrt(q_vec(i)*(1-q_vec(i)));
		   if(tmp>=stat) stat=tmp;
		   
		}
		stat_hc=stat;
	}else{
		stat_hc=1.0;
	}
	//SKAT
	stat_skat=0.0; ///*/!!
	for(i=1; i<=nsel_n; i++) ///*/!!
	{
		  stat_skat+=Utot(i)*Utot(i); ///*/!!
	}
	
	// burden
	double U_sum=0.0;
	
	for (i=1; i<=nsel_n; i++) { ///*/!!
		U_sum+=weights(i)*Utot(i);
	}
	
    stat_burden=U_sum*U_sum;
    
	///// adaptive////
	unsigned long sim_counter=0;
	unsigned long max_sim=num_sim_region;
	unsigned long maximum_sim=max_sim_region;
	int running=1;
	bool valid_draw=false; ///*/!!
	while(sim_counter<maximum_sim && running==1)
	{
			U=0.0;
			for(j=1;j<=inf_ctr;j++){ // only U is updated, Vtot and EXtot keep unchanged
		  
					ss=ss_array[j];
					valid_draw=false; ///*/!!
					while(valid_draw==false) ///*/!!
					{
						if(random_index!=NULL){ delete[] random_index; random_index=NULL;}
						random_index=do_rand_sim_fam(ss);
						if(random_index==NULL){ MYOUTPUT(os,"error random_index\n"); exit(1);}
						valid_draw=check_consistent_cond(ss, index, random_index);
					} ///*/!!
					tmp_counter=0;
					for (i=0; i<ss->n_offspring; i++) {
							if (ss->otr[i]!=NULL && ss->otr[i][trait_id[0]]<kVeryLargeNumber && (ss->which_comm_ogt(ss->ogt[i],ss->sex[i]))>=0 ){
								t = ss->otr[i][trait_id[0]] - offset;
								tmp_counter++;
								for(l=0;l<ss->nMarkers;l++){
								    G(l+1)= gx[0][ss->common_ogt[random_index[tmp_counter]-1][l].allele_cnt(1)];//get_coding(ss->common_ogt[random_index[tmp_counter]-1][l].allele_cnt(2), gen_model);
								}
								U+= t*G;//-EX);
								
							}
					}
					if(random_index!=NULL){ delete[] random_index; random_index=NULL;}
			}
			U-=EXtot;
			/// MAX
			stat_tmp=0.0;
			for(i=1;i<=nsel_n;i++){ ///*/!!
				if(Vtot(i,i)>=very_small_number){
					tmp=fabs(U(i))/sqrt(Vtot(i,i));
					if(tmp>=stat_tmp){ stat_tmp=tmp;}
				}
			} ///*/!!
			
			if(stat_tmp > stat_max || fabs(stat_tmp-stat_max)<=extreme_small_number){ctr_max++;} //update 10/04/2020
			stat_tmp1=stat_tmp;
			
			/// HC
			inf_var=0;
			for(i=1;i<=nsel_n;i++){ ///*/!!
			   if(var_indicator(i)>=(double)min_size){
					inf_var++;
					x=fabs(U(i))/sqrt(Vtot(i,i));
					normbase(&x,&phi1);
					x=-fabs(U(i))/sqrt(Vtot(i,i));
					normbase(&x,&phi2);
					q_vec(inf_var)=fmin(1-phi1+phi2,1.0-very_small_number);
			   }
			}
			for (i = 1; i <= inf_var-1; i++){       
				for (l = 1; l <= inf_var-i; l++)  {
						if (q_vec(l) > q_vec(l+1)) {
							tmp = q_vec(l+1); 
							q_vec(l+1) = q_vec(l); 
							q_vec(l) = tmp; 
				   }
				}			  
			}
			if(inf_var>0){ //changed 02/18/2020 to catch case where no variants are informative for HC, resulting p-value should be 1.0
				stat_tmp=((double)(1/((double)inf_var))-q_vec(1))/sqrt(q_vec(1)*(1-q_vec(1)));
				for (i = 1; i <= (int)inf_var; i++){ 
				   tmp=((double)((double)i/((double)inf_var))-q_vec(i))/sqrt(q_vec(i)*(1-q_vec(i)));
				   if(tmp>=stat_tmp) stat_tmp=tmp;
				}
			}else{
				stat_tmp=1.0;
			}
			if(stat_tmp > stat_hc || fabs(stat_tmp-stat_hc)<=extreme_small_number){ctr_hc++;} //update 10/04/2020
			stat_tmp2=stat_tmp;
			
			
			/// SKAT
			stat_tmp=0.0; ///*/!!
			for(i=1; i<=nsel_n; i++) ///*/!!
			{
				stat_tmp+=U(i)*U(i); ///*/!!
			}
			
			if( stat_tmp > stat_skat || fabs(stat_tmp-stat_skat)<=extreme_small_number){ctr_skat++;}//update 10/04/2020
		    stat_tmp3=stat_tmp;
	
			/// Burden	
			U_sum=0.0;
			for(i=1;i<=nsel_n;i++){ ///*/!!
						 U_sum+=weights(i)*U(i);
			}
			stat_tmp=U_sum*U_sum;
			if( stat_tmp > stat_burden || fabs(stat_tmp-stat_burden)<=extreme_small_number){ctr_burden++;} //update 10/04/2020
			stat_tmp4=stat_tmp;
			
			///// adaptive////
			sim_counter++;
			if(sim_counter==max_sim){
				if(fmin(fmin(fmin(ctr_burden, ctr_skat),ctr_max),ctr_hc)<=5)
				{ 
					max_sim*=10;
				}
				else{ running=0;}
			}
			
			///// adaptive////
			
	
	}
	
	
	
	p_burden=(double)ctr_burden/(double)sim_counter;
	p_skat=(double)ctr_skat/(double)sim_counter;
	p_maxsv=(double)ctr_max/(double)sim_counter;
	p_hc=(double)ctr_hc/(double)sim_counter;
	
	//!/!/ update 10/09/2020: add ACAT statistic to FBAT software
	ACAT_stat=tan((0.5-p_burden)*M_PI)+tan((0.5-p_hc)*M_PI)+tan((0.5-p_maxsv)*M_PI)+tan((0.5-p_skat)*M_PI);
	p_acat=0.5-atan(ACAT_stat/4.0)/M_PI;
	
    MYOUTPUT(os,"--- N: "<<ctr<<" ")
    MYOUTPUT(os,"max_sv: "<< setprecision(10) << p_maxsv <<"	HC: "<< p_hc <<"	VC: " << p_skat <<"	Burden: " << p_burden <<"	ACAT: "<< p_acat<< "  max_sv_stat: " << stat_max)
    MYOUTPUT(os," nsim: "<<sim_counter<<"\n")
	if (idx)
		delete[] idx;
	if(ss_array)
		delete[] ss_array;
	
	
} //!/
void FBAT::gen_rvtest(char *s) {//!/
	
	char word[256], *c;
	int i, j, k, l;
	i = WordCount(s);
	double offset = tr_offset;
	
	//---------------------------------------------
	for (j=2; j<=i && GetWord(s, j, word) && *word=='-'; j++){ //!/ update 01/20/2021
		char msg[256]; //!/ update 01/20/2021
		sprintf(msg, "no options implemented for gen_rv test"); //!/ update 01/20/2021
		throw Param_Error(msg); //!/ update 01/20/2021
	} //!/ update 01/20/2021
	int nsel = j>=i? nLoci : i-j+1;
	int *idx = new int[nsel];
	if (j>i)
		for (nsel=0; nsel<nLoci; nsel++)
			idx[nsel] = nsel;
	else {
		for (nsel=0; j<=i; j++) {
			GetWord(s, j, word);
			if ((k=find_locus(word))>=0)
				idx[nsel++] = k;
			else {
				if (idx!=NULL)
					delete[] idx;
				char msg[256];
				sprintf(msg, "Marker %s not found", word);
				throw Param_Error(msg);
			}
		}
	}
	//--------------------------------------------------
	double gx[kMaxAlleles][3];
	if (gen_model==model_additive || gen_model==model_dominant || gen_model==model_recessive)
		for (i=0; i<kMaxAlleles; i++)
			gcode(&gx[i][0], gen_model);
	else
		throw Param_Error("please select additive, dominant or recessive genetic model.");
		
	Haplotype_List *emhl = hapfreq_ss(nsel, idx);
	if (emhl==NULL) 
		return;
    if (n_sel_trt>1) //!/ update 01/20/2021
		throw Param_Error("multivariate gen_rv test is not available"); //!/ update 01/20/2021
	// set trait to affection status if no other trait selected yet
	if (n_sel_trt == 0) {
		n_sel_trt = 1;
		trait_id[0] = 0;
	}
	//------------------------------------------------------------
	
	
    MYOUTPUT(os, "region-based test for the "<< nsel<< " markers between: ")
	MYOUTPUT(os, loc[idx[0]].name << " and " << loc[idx[nsel-1]].name<<" ")

	
	HaploidTable *htab = new HaploidTable(loc, nsel, idx, emhl);
	ColumnVector EX(nsel);
	ColumnVector EXtot(nsel);
	ColumnVector U(nsel);
	ColumnVector G(nsel);
	Matrix V(nsel,nsel);
	ColumnVector Utot(nsel);
	Matrix Vtot(nsel,nsel);

    Sufficient_Stat* ss;
	Sufficient_Stat** ss_array;
	
	Utot=0.0;
	Vtot=0.0;
	EXtot=0.0;
	
	int ninff=0; //!/!/
	int inf_ctr=0;
	int recomb_ctr=0;
	int geno1_ctr=0;
	int too_many_ctr=0;
	int ctr=0;
	int nonident_ctr=0;

	
	double t;
	double alpha;
	double beta;
	int* random_index=NULL; //!=!
	int tmp_counter; //!=!
	int info_available; //!/!
	ColumnVector weights(nsel);
	for(i=1;i<=nsel;i++) weights(i)=1.0;
	
	ColumnVector var_indicator(nsel);
	for(i=1;i<=nsel;i++) var_indicator(i)=0.0;
	
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node && flist->node->stat) { ctr++;}
		}
	}
	ss_array=new Sufficient_Stat*[ctr];
	
	ctr=0;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node && flist->node->stat) {
				ctr++;
				flist->node->stat->analyze(emhl);
				if(flist->node->stat->mhap!=NULL && flist->node->stat->mhap->len()>maxcmh){ too_many_ctr++;} 
			    else if(flist->node->stat->mhap==NULL){ recomb_ctr++; }
				else if(flist->node->stat->n_common_ogt<2){ geno1_ctr++; }
				else if(flist->node->stat->ogpat!=NULL) {
					inf_ctr++; info_available=0; //!/!
					ss=flist->node->stat;
					ss_array[inf_ctr]=ss;
					alpha=0.0;
					beta=0.0;
					ColumnVector P(ss->n_common_ogt);
					for (i=0; i<ss->n_common_ogt; i++)
						P(i+1) = ss->pg[i];
					DiagonalMatrix Q(ss->n_common_ogt);
					for (i=0; i<ss->n_common_ogt; i++)
						Q(i+1,i+1) = P(i+1);
					SymmetricMatrix PP(ss->n_common_ogt);
					for (i=0; i<ss->n_common_ogt; i++)
						for (j=i; j<ss->n_common_ogt; j++)
							PP(i+1,j+1) = ss->pgg[i][j];
							
					Matrix X(nsel,ss->n_common_ogt);
					for (i=0; i<nsel; i++){
						for (j=0; j<ss->n_common_ogt; j++){
							X(i+1,j+1) = gx[0][ss->common_ogt[j][i].allele_cnt(1)];
						}
					}
					EX = X * P;
					U=0.0;
					V=0.0;
					for (i=0; i<ss->n_offspring; i++) {
						if (ss->otr[i]!=NULL && ss->otr[i][trait_id[0]]<kVeryLargeNumber && (j=ss->which_comm_ogt(ss->ogt[i],ss->sex[i]))>=0 ){
							t = ss->otr[i][trait_id[0]] - offset;
							beta += t;
							alpha += t*t;
							for(k=0;k<ss->nMarkers;k++){
							  G(k+1)= gx[0][ss->common_ogt[j][k].allele_cnt(1)];
							}
							U+= t*(G-EX);
							EXtot+=t*EX;
						}
					}	
					beta=beta*beta-alpha;
					V=alpha* X * Q * X.t()+beta* X * PP * X.t()- (alpha+beta)*X* P*P.t()*X.t(); // ok for missing phenotypes: nov29_2019: missing phenotype treated as t=0. 	
                    // 	X * PP * X.t() measures E[x_1 x2], the mean over choosing two genotypes randomly from the null, without ordering			
					for(i=1;i<=nsel;i++){ if(V(i,i)>very_small_number){ var_indicator(i)+=1.0; info_available=1;} } //!/!
					Utot+=U;
					Vtot+=V;
					if(info_available==1){ninff++;} //!/!
				
				}
				else {
					nonident_ctr++;
				}
			}
		}
	}
	
	double stat_max, stat_hc, stat_skat, stat_burden;
	
	int ctr_max, ctr_hc, ctr_skat, ctr_burden;
	ctr_max=ctr_hc=ctr_skat=ctr_burden=0;
	
	double stat, stat_tmp, stat_tmp1, stat_tmp2, stat_tmp3, stat_tmp4, tmp, ACAT_stat;
	double p_burden, p_skat, p_maxsv, p_hc, p_acat;
	
	// max-sv
	stat=0.0;
	for(j=1;j<=nsel;j++)
	{
		if(Vtot(j,j)>=very_small_number){
			tmp=fabs(Utot(j))/sqrt(Vtot(j,j));
			if(tmp>=stat){ stat=tmp;}
		}
	}
	stat_max=stat;
	// highcrit
	stat=0.0;
	ColumnVector q_vec(nsel);
	q_vec= 0.0;
	double x;
	double phi1;
	double phi2;
	int inf_var=0;
	for(i=1;i<=nsel;i++){
		if(var_indicator(i)>=(double)min_size){//if(Vtot(i,i)>very_small_number){//!/ maybe implement option here ?
			inf_var++;
			x=fabs(Utot(i))/sqrt(Vtot(i,i));
			normbase(&x,&phi1);
			x=-fabs(Utot(i))/sqrt(Vtot(i,i));
			normbase(&x,&phi2);
			q_vec(inf_var)=fmin(1-phi1+phi2,1.0-very_small_number);
		}
	}
    
    for (i = 1; i <= inf_var-1; i++){       
		for (j = 1; j <= inf_var-i; j++)  {
				if (q_vec(j) > q_vec(j+1)) {
					tmp = q_vec(j+1); 
					q_vec(j+1) = q_vec(j); 
					q_vec(j) = tmp; 
		   }
        }			  
    }
	if(inf_var>0){
		stat=((double)(1/((double)inf_var))-q_vec(1))/sqrt(q_vec(1)*(1-q_vec(1)));
		for (i = 1; i <= (int)inf_var; i++){  //!/
		   tmp=((double)((double)i/((double)inf_var))-q_vec(i))/sqrt(q_vec(i)*(1-q_vec(i)));
		   if(tmp>=stat) stat=tmp;
		   
		}
		stat_hc=stat;
	}else{
		stat_hc=1.0;
	}
	//SKAT
	stat_skat= (Utot.t() * Utot).AsScalar(); 
	// burden
	double U_sum=0.0;
	
	for (i=1; i<=nsel; i++) {
		U_sum+=weights(i)*Utot(i);
	}
	
    stat_burden=U_sum*U_sum;
    
	///// adaptive////
	unsigned long sim_counter=0;
	unsigned long max_sim=num_sim_region;
	unsigned long maximum_sim=max_sim_region;
	int running=1;
	
	while(sim_counter<maximum_sim && running==1)
	{
			U=0.0;
			for(j=1;j<=inf_ctr;j++){ // only U is updated, Vtot and EXtot keep unchanged
		  
					ss=ss_array[j];
					
					random_index=do_rand_sim_fam(ss);
					if(random_index==NULL){ MYOUTPUT(os,"error random_index\n"); exit(1);}
					tmp_counter=0;
					for (i=0; i<ss->n_offspring; i++) {
							if (ss->otr[i]!=NULL && ss->otr[i][trait_id[0]]<kVeryLargeNumber && (ss->which_comm_ogt(ss->ogt[i],ss->sex[i]))>=0 ){
								t = ss->otr[i][trait_id[0]] - offset;
								tmp_counter++;
								for(l=0;l<ss->nMarkers;l++){
								    G(l+1)= gx[0][ss->common_ogt[random_index[tmp_counter]-1][l].allele_cnt(1)];
								}
								U+= t*G;//-EX);
							}
					}
					if(random_index!=NULL){ delete[] random_index; random_index=NULL;}
			}
			U-=EXtot;
			/// MAX
			stat_tmp=0.0;
			for(i=1;i<=nsel;i++){
				if(Vtot(i,i)>=very_small_number){
					tmp=fabs(U(i))/sqrt(Vtot(i,i));
					if(tmp>=stat_tmp){ stat_tmp=tmp;}
				}
			}
			
			if( stat_tmp > stat_max || fabs(stat_tmp-stat_max)<=extreme_small_number){ctr_max++;} //update 10/04/2020
			stat_tmp1=stat_tmp;
			
			/// HC
			inf_var=0;
			for(i=1;i<=nsel;i++){
			   if(var_indicator(i)>=(double)min_size){
					inf_var++;
					x=fabs(U(i))/sqrt(Vtot(i,i));
					normbase(&x,&phi1);
					x=-fabs(U(i))/sqrt(Vtot(i,i));
					normbase(&x,&phi2);
					q_vec(inf_var)=fmin(1-phi1+phi2,1.0-very_small_number);
			   }
			}
			for (i = 1; i <= inf_var-1; i++){       
				for (l = 1; l <= inf_var-i; l++)  {
						if (q_vec(l) > q_vec(l+1)) {
							tmp = q_vec(l+1); 
							q_vec(l+1) = q_vec(l); 
							q_vec(l) = tmp; 
				   }
				}			  
			}
			if(inf_var>0){ //changed 02/18/2020 to catch case where no variants are informative for HC, resulting p-value should be 1.0
				stat_tmp=((double)(1/((double)inf_var))-q_vec(1))/sqrt(q_vec(1)*(1-q_vec(1)));
				for (i = 1; i <= (int)inf_var; i++){ 
				   tmp=((double)((double)i/((double)inf_var))-q_vec(i))/sqrt(q_vec(i)*(1-q_vec(i)));
				   if(tmp>=stat_tmp) stat_tmp=tmp;
				}
			}else{
				stat_tmp=1.0;
			}
			if( stat_tmp > stat_hc || fabs(stat_tmp-stat_hc)<=extreme_small_number){ctr_hc++;} //update 10/04/2020
			stat_tmp2=stat_tmp;
			
			
			/// SKAT
			stat_tmp = (U.t() * U).AsScalar(); 
			if( stat_tmp > stat_skat || fabs(stat_tmp-stat_skat)<=extreme_small_number){ctr_skat++;}//update 10/04/2020
		    stat_tmp3=stat_tmp;
	
			/// Burden	
			U_sum=0.0;
			for(i=1;i<=nsel;i++){
						 U_sum+=weights(i)*U(i);
			}
			stat_tmp=U_sum*U_sum;
			if( stat_tmp > stat_burden || fabs(stat_tmp-stat_burden)<=extreme_small_number){ctr_burden++;} //update 10/04/2020
			stat_tmp4=stat_tmp;
			
			///// adaptive////
			sim_counter++;
			if(sim_counter==max_sim){
				if(fmin(fmin(fmin(ctr_burden, ctr_skat),ctr_max),ctr_hc)<=5)
				{ 
					max_sim*=10;
				}
				else{ running=0;}
			}
			///// adaptive////
	}
	
	p_burden=(double)ctr_burden/(double)sim_counter;
	p_skat=(double)ctr_skat/(double)sim_counter;
	p_maxsv=(double)ctr_max/(double)sim_counter;
	p_hc=(double)ctr_hc/(double)sim_counter;
	
	//!/!/ update 10/09/2020: add ACAT statistic to FBAT software
	ACAT_stat=tan((0.5-p_burden)*M_PI)+tan((0.5-p_hc)*M_PI)+tan((0.5-p_maxsv)*M_PI)+tan((0.5-p_skat)*M_PI);
	p_acat=0.5-atan(ACAT_stat/4.0)/M_PI;
	
    MYOUTPUT(os,"--- N: "<<ctr<<" NINFF: " <<ninff<<"	NINFF_g: "<<inf_ctr<<" other: "<< recomb_ctr<<"	" << geno1_ctr<<"	"<<too_many_ctr<<"	"<<nonident_ctr<<"	")
    MYOUTPUT(os,"max_sv: "<< setprecision(10) << p_maxsv <<"	HC: "<< p_hc <<"	VC: " << p_skat <<"	Burden: " << p_burden <<"	ACAT: "<< p_acat<< "  max_sv_stat: " << stat_max)
    MYOUTPUT(os," nsim: "<<sim_counter<<"\n")
	if (idx)
		delete[] idx;
	if(ss_array)
		delete[] ss_array;
} //!/

void FBAT::finemap(char *s) {
	
	char word[1024], *c;
	int i, k, j;
	/*########################*/
	double offset = tr_offset;
	if (n_sel_trt>1) 
		throw Param_Error("multivariate fine-mapping is not available"); 
	// set trait to affection status if no other trait selected yet
	if (n_sel_trt == 0) {
		n_sel_trt = 1;
		trait_id[0] = 0;
	}
	if(min_size<10){min_size=10;}
	/*########################*/
	int wcnt = WordCount(s);
	int flag=0;
	bool verbose=false;
	for (i=2; i<=wcnt && GetWord(s, i, word) && *word=='-'; i++) {
		if (word[2]==0)
			switch (word[1]) {
				case 'v': 
					verbose = true;
					break;
			}
	}

	
	if (gen_model == model_genotype)
		throw Param_Error("The requested genotype test has not been implemented.");
	
	MYOUTPUT(os, "Performing DAP-G fine-mapping...\n")
	
	print_settings();
	/*########################*/
	
	int nsel = i>=wcnt? nLoci : wcnt-i+1;
	int *idx = new int[nsel];
	if (i>wcnt)
		for (nsel=0; nsel<nLoci; nsel++)
			idx[nsel] = nsel;
	else {
		for (nsel=0; i<=wcnt; i++) {
				GetWord(s, i, word);
				if ((k=find_locus(word))>=0)
					idx[nsel++] = k;
				else {
					if (idx!=NULL)
						delete[] idx;
					char msg[256];
					sprintf(msg, "Marker %s not found", word);
					throw Param_Error(msg);
				}
		}
	}
	
	double gx[kMaxAlleles][3];
	if (gen_model==model_additive || gen_model==model_dominant || gen_model==model_recessive)
		for (i=0; i<kMaxAlleles; i++)
			gcode(&gx[i][0], gen_model);

	/*########################*/
	
	int na;
	double *frq; 

	for (i=0; i<nsel; i++) {
				na = loc[idx[i]].nAlleles;
				if(na>2){
					throw Param_Error("Fine-mapping implemented for bi-allelic variants only.");
				}
				if(loc[idx[i]].sexlinked()){
					throw Param_Error("Fine-mapping implemented for non-sexlinked variants only.");
				}
	}
		
	int loc_idx;
	Sufficient_Stat *ss;
	
	int geno_ctr=0;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
			for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				geno_ctr++;
			}
	}
	/*########################*/
	int dim_var=2*nsel;
	
	int *fcnt=new int[dim_var+1];
	if(fcnt==NULL){
		throw Param_Error("Memory allocation problem.");
	}
	for(i=0;i<=dim_var;i++){fcnt[i]=0;}
	
	ColumnVector Z_tmp(dim_var);
	Matrix R_work_tmp(dim_var, geno_ctr);	
	
	
	Z_tmp=0; R_work_tmp=-9; // -9 is missing value
	int ctr=0;
	std::vector<std::string> geno_map_all;
	std::string str_name;
	std::string a1="_1";
	std::string a2="_2";
	
	
	for(loc_idx=0;loc_idx<nsel;loc_idx++){
		
		str_name=loc[loc_idx].name;
		geno_map_all.push_back(str_name + a1);
		geno_map_all.push_back(str_name + a2);
		
		
		na=loc[loc_idx].nAlleles;
		
		if(na==2){
		
			ColumnVector S(na);	
			Matrix V(na,na);	//nuisance
			Matrix A(na,na);	
			Matrix B(na,na);
			ColumnVector X(na);
			ColumnVector FS(na);
			ColumnVector PS(na);	
			ColumnVector FES(na);
			Matrix FV(na,na);
			Matrix FA(na,na);
			Matrix FB(na,na);
			ColumnVector FX(na);
			
			X = 0;
			S = 0;
			A = 0;
			B = 0;
			V = 0;
			
			geno_ctr=0;
			
			for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
				
				for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
					geno_ctr++; //!/
					if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
						
						
						if (flist->node->stat!=NULL) {
							delete flist->node->stat;
							flist->node->stat = NULL;
						}
						ss = new Sufficient_Stat(flist->node, 1, &loc_idx, loc[loc_idx].sexlinked());
						flist->node->stat = ss;
						
						if(ss->ogt[0][0].defined()) R_work_tmp(ctr+1, geno_ctr) = gx[0][ss->ogt[0][0].allele_cnt(1)];
						if(ss->ogt[0][0].defined()) R_work_tmp(ctr+1, geno_ctr) = gx[1][ss->ogt[0][0].allele_cnt(2)];
						
						if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
							continue; 

						ss->get_common_ogt();
						ss->analyze();
						if (ss->informative() && ss->fbat_stat3(trait_id[0], offset, FES, FS, FV, FA, FB, FX, gx, flag)) {

							FS -= FES;
							PS += FS;
					
							S += FS;
							V += FV;
							
							if (FV(1,1)>very_small_number) {
										fcnt[ctr+1]++; //!/		
							}
							if (FV(2,2)>very_small_number) {
										fcnt[ctr+2]++; //!/		
							}
							
						}
					}
					
				}
	  
			}
			
			Z_tmp(ctr+1)=S(1)/sqrt(V(1,1));
			Z_tmp(ctr+2)=S(2)/sqrt(V(2,2));
		}
		ctr=ctr+2;
		
	}
	/*########################*/
	
	double m;
	double tmp;
	
	Matrix R_work(dim_var, dim_var);
	for(i=1;i<=dim_var;i++){
		m=0.0;
		ctr=0;
		for(k=1;k<=geno_ctr;k++){
				if(R_work_tmp(i,k)!=-9 ){
					m+=R_work_tmp(i,k);
					ctr++;
				}
		}
		m/=(double)ctr;
		for(k=1;k<=geno_ctr;k++){
				if(R_work_tmp(i,k)!=-9 ){
					R_work_tmp(i,k)=R_work_tmp(i,k)-m;
				}
				else{
					R_work_tmp(i,k)=0.0;
				}
		}
	}

	for(i=1;i<=dim_var;i++){
		for(j=1;j<=dim_var;j++){
			tmp=0.0;
			for(k=1;k<=geno_ctr;k++){
				tmp+=R_work_tmp(i,k)*R_work_tmp(j,k);
			}
			R_work(i,j)=tmp;
		}
	}
	
	// transform from 'covariance' to 'correlation'
	ColumnVector variances(dim_var); variances=0.0;
	for(i=1;i<=dim_var;i++){variances(i)=R_work(i,i);}
	for(i=1;i<=dim_var;i++){
		for(j=1;j<=dim_var;j++){
			R_work(i,j)=R_work(i,j)/sqrt(variances(i)*variances(j));
		}
	}
	
	/*########################*/
	int* include=new int[dim_var+1];
	if(include==NULL){
		throw Param_Error("Memory allocation problem.");
	}
	int incl_ctr;
	incl_ctr=0;
	for(i=1;i<=dim_var;i+=2){
		include[i]=0;
		include[i+1]=0;
		
		if(gen_model == model_additive && fcnt[i]>=min_size){
			include[i]=1;
			incl_ctr++;
			
		}
		if(gen_model == model_recessive && fcnt[i]>=min_size){
			include[i]=1;
			incl_ctr++;
		}
		if(gen_model == model_recessive && fcnt[i+1]>=min_size){
			include[i+1]=1;
			incl_ctr++;
		}
	}
	
	if(incl_ctr<5)
		throw Param_Error("After filtering, less than 5 variables left.");
    /* #######################################################*/
	
	nsel=incl_ctr;
	ColumnVector Z(nsel);
	Matrix R(nsel, nsel);
	R=0;
	
	int c1,c2; c1=c2=0;
	for(i=1;i<=dim_var;i++){
		if(include[i]==1){
			c2=0;
			c1++;
			Z(c1)=Z_tmp(i);
			for(j=1;j<=dim_var;j++){
				if(include[j]==1){
					c2++;
					R(c1,c2)=R_work(i,j);
				}
			}
		}
	}
	
	/*##############*/
	
	std::map<int, std::string> geno_map_p;
	incl_ctr=0;
    for(i=1;i<=dim_var;i++){
		if(include[i]==1){
			geno_map_p[incl_ctr] = geno_map_all[i-1];
			incl_ctr++;
		} 
    }

	/*########################################################*/
    controller con;
    con.initialize(Z, R, nsel, geno_map_p);
	
	
    con.fine_map();
    finemap_results fr=con.summarize_approx_posterior(verbose);
	
 
	MYOUTPUT(os, "\ncluster 	  #variants  	   PIP    	   r2\n")	   
    for(i=0;i<fr.cluster_counts.size();i++){
		MYOUTPUT(os, setw(3)<< ">>"<< setw(12)<< i+1<< setw(15) <<fr.cluster_counts[i]<< setw(18) << fr.cluster_pips[i]<< setw(19) <<fr.cluster_r2s[i]<<"\n")
    }
	
	MYOUTPUT(os, "\nvariant   	incl_probability 	    cluster    	   log_BF 		z-score\n")	   
    for(i=0;i<fr.variant_names.size();i++){
		MYOUTPUT(os, setw(4)<< ">>>" << setw(24) <<fr.variant_names[i]<< setw(25) <<fr.incl_probs[i]<< setw(18) << fr.cluster_ids[i]<< setw(19) <<fr.log_bfs[i]<< setw(15) <<fr.z_scores[i]<<"\n")
    }
	
	/*########################*/
	// free memory
	if (idx) delete[] idx;
	if(fcnt) delete fcnt;
	if(include) delete include;
}
void FBAT::hapview_ss(char *s) {
	
	int i = WordCount(s);
		
	char name[256];
	bool censored = false;
	
	bool informative_only = false;
	
	int j, k;
	int flag = 0;
	for (j=2; j<=i && GetWord(s,j,name) && *name=='-'; j++)
		switch (name[1]) {
			case 'i':
				informative_only = true;
				break;
			case 'c':
				censored = true;
				flag |= flag_censored;
				break;
			case 'e':
				throw Param_Error("statview with empirical variance option is not implemented");
				flag |= flag_empirical_var;
				break;
			case 'o':
				throw Param_Error("statview with offset optimization option is not implemented");
				break;
			default:
				throw Param_Error("unknown option");
		}
	
	if (n_sel_trt>1)
		throw Param_Error("statview for multivariate test is not implemented");\
		
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
	
	print_settings();
	
	Haplotype_List *emhl = hapfreq_ss(nsel, idx);
	
	if (emhl==NULL) 
		return;
	
	HaploidTable *htab = new HaploidTable(loc, nsel, idx, emhl);
	
	if (htab!=NULL)
		htab->print(os);
	
	MYOUTPUT(os, "\nFBAT Statistics\n------------------------------\n\n")
	
	char msg[255];
	int na = htab->size;
	
	if (na>kMaxAlleles)
		throw ERROR("too many alleles in haplotype test (max=80)");
		
	double gx[kMaxAlleles][3];
	if (gen_model==model_additive || gen_model==model_dominant || gen_model==model_recessive)
		for (i=0; i<kMaxAlleles; i++)
			gcode(&gx[i][0], gen_model);
	else
		throw Param_Error("Haplotype test under the current genetic model has not been implemented");

	ColumnVector FS(na);
	ColumnVector FES(na);
	ColumnVector X(na);
	Matrix FV(na,na);
	ColumnVector S(na);	// FS-FES
	Matrix V(na,na);
	Matrix A(na,na);
	Matrix B(na,na);
	S = 0;
	V = 0;
	int tot_cnt = 0;
	int *fcnt = new int[na+1];
	for (i=0; i<=na; i++)
		fcnt[i] = 0;
	bool good;
	
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next)
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
			if (flist->node && flist->node->stat) {
				tot_cnt++;
				flist->node->stat->analyze(emhl);
				sprintf(msg, "pedigree %s , nuclear family [%s x %s]", ped->node->name, (flist->node->fa()?flist->node->fa()->id:"-"), (flist->node->mo()?flist->node->mo()->id:"-"));
				if (!informative_only)
					flist->node->stat->output_info(os, htab, msg);
				if (flist->node->stat->informative() && flist->node->stat->fbat_stat3(trait_id[0], tr_offset, FES, FS, FV, A, B, X, gx, flag, htab)) {
					for (good=false,i=1; i<=na; i++)
						if (FV(i,i)>very_small_number) {
							fcnt[i]++;
							good = true;
						}
					if (informative_only && good)
						flist->node->stat->output_info(os, htab, msg);
					if (good) {	// informative for ar least one allele
						fcnt[0]++;
						S += FS-FES;
						V += FV;
						// print only informative alleles
						MYOUTPUT(os, "\nStatistic Summary:\n")
						MYOUTPUT(os, setw(10) << "allele" << setw(12) << "S" << setw(12) << "E(S)" << "Var(S)...\n")
						for (i=1; i<=na; i++)
							if (FV(i,i)>very_small_number) {
								MYOUTPUT(os, "a" << setw(9) << i << setw(12) << FS(i) << setw(12) << FES(i))
								for (j=1; j<=na; j++)
									if (FV(j,j)>very_small_number)
										MYOUTPUT(os, setw(10) << FV(i,j))
								MYOUTPUT(os, "\n")
							}
						MYOUTPUT(os, "\n------------------------------\n")						
					}
				}
			}
	
	MYOUTPUT(os, "total family count = " << tot_cnt << "; informative family count = " << fcnt[0] << "\n")
	MYOUTPUT(os, "Statistic summary:\n")
	MYOUTPUT(os, setw(10) << "allele" << setw(10) << "fam#" << setw(12) << "S-E(s)" << setw(12) << "Var(S)" << "\n")
	for (i=1; i<=na; i++) {
		MYOUTPUT(os, "a" << setw(9) << i << setw(10) << fcnt[i] << setw(12) << S(i))
		for (j=1; j<=na; j++)
			MYOUTPUT(os, setw(12) << V(i,j))
		MYOUTPUT(os, "\n")
	}
	MYOUTPUT(os, "\n")
	
	if (emhl)
		delete emhl;
	
	if (htab)
		delete htab;
		
	if (idx)
		delete[] idx;
	
	if (fcnt)
		delete[] fcnt;
}

// estimate the (diploid) haplotype frequency based on mhap of each nuclear family using EM, using  Sufficient_stat
Haplotype_List* FBAT::hapfreq_ss(int nmrk, int idx[]) {
	
	int i, j;
	double d, sum=0.0;
	
	bool xlinked = loc[idx[0]].sexlinked();
	for (i=1; i<nmrk; i++)
		if (loc[idx[i]].sexlinked() != xlinked)
			throw Param_Error("inconsistent sexlink code for the selected markers");

	Haplotype_List	*hl2, *hl = NULL; 
	Sufficient_Stat *ss = NULL;
	Mating_Haplotype *mh = NULL;
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	
	for (ped=peds; ped!=NULL; ped=ped->next)
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				if (flist->node->stat!=NULL) {
					delete flist->node->stat;
					flist->node->stat = NULL;
				}
				ss = new Sufficient_Stat(flist->node, nmrk, idx, xlinked);
				flist->node->stat = ss;
				
				if (ss->mhap==NULL || ss->mhap->len()>maxcmh)
					continue;
				ss->get_common_ogt();	// ss->mhap->pset is updated
				// first only pick the unique haplotypes in mh
				d = 1.0/ss->mhap->len();
				for (mh=ss->mhap; mh!=NULL; mh=mh->next) {
					for (i=0; i<2; i++) {
						if (!mh->unique(i))
							mh->hl[i] = NULL;
						else {
							if (hl==NULL || (hl2=hl->find(mh->h[2*i], mh->h[2*i+1]))==NULL)
								hl2 = hl = new Haplotype_List(nmrk, mh->h[2*i], mh->h[2*i+1], hl);
							hl2->p += d;
							mh->hl[i] = hl2;
							sum += d;
						}
					}
				}
				//sum += 2;
			}
	
	for (hl2=hl; hl2!=NULL; hl2=hl2->next)	// now hl->p is the freq from unique hl assuming equal likely mh.
		hl2->p /= sum;
	
	// now do EM based on hl->p
	double p[2];
	bool converged = false;
	while (!converged) {
		for (hl2=hl; hl2!=NULL; hl2=hl2->next)
			hl2->q = 0.0;
		for (ped=peds; ped!=NULL; ped=ped->next)
			for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				if (flist->node==NULL || flist->node->testflag(flag_uninformative)!=0 || (ss=flist->node->stat)==NULL || ss->mhap==NULL)
					continue;
				for (sum=0.0,mh=ss->mhap; mh!=NULL; mh=mh->next) {
					if (mh->hl[0]==NULL && mh->hl[1]==NULL)	// discard mh with both hl = non-unique
						continue;
					for (i=0; i<2; i++)
						p[i] = (mh->hl[i]!=NULL?mh->hl[i]->p:0.001); // fixed p=0.001 for all non-unique haplotypes
					sum += (mh->hl[0]==mh->hl[1]?1:2)*mh->pset*p[0]*p[1];
				}
				if (sum > 0.0) {
					for (mh=ss->mhap; mh!=NULL; mh=mh->next) {
						if (mh->hl[0]==NULL && mh->hl[1]==NULL)	// discard mh with both hl = non-unique
							continue;
						if (mh->hl[0]!=NULL && mh->hl[1]!=NULL) {
							d = (mh->hl[0]==mh->hl[1]?1:2)*mh->pset*mh->hl[0]->p*mh->hl[1]->p/sum;
							for (i=0; i<2; i++)
								mh->hl[i]->q += d;
						} else {
							j = (mh->hl[0]==NULL?0:1);	// mh->hj[j] is non-unique
							p[1-j] = mh->hl[1-j]->p;
							p[j] = 0.001;
							mh->hl[1-j]->q += 2*mh->pset*p[j]*p[1-j]/sum;
						}
					}
				}
			}
		
		
		for (sum=0.0,hl2=hl; hl2!=NULL; hl2=hl2->next)
			sum += hl2->q;
		
		if (sum>0.0)
			for (hl2=hl; hl2!=NULL; hl2=hl2->next)
				hl2->q /= sum;
		
		for (hl2=hl; hl2!=NULL && fabs(hl2->p-hl2->q)<kConvFreq; hl2=hl2->next)
			;
		
		converged = (hl2==NULL);	// converged
		
		for (hl2=hl; hl2!=NULL; hl2=hl2->next)
			hl2->p = hl2->q;
		 
	}

	return hl;
		
}

void FBAT::hbat(char *s) {
	
	bool empirical = false, censored=false;
	bool optimize_offset = false;
	bool permutation = false;
	char word[256], *c;
	int flag = 0;
	int i, j, k;
	int cycles = 0;
	
	i = WordCount(s);
	
	for (j=2; j<=i && GetWord(s, j, word) && *word=='-'; j++) {
		if (word[1]=='p') {
			permutation = true;
			if (word[2]!=0)
				cycles = strtol(word+2, &c, 10);
		}
		else if (word[2]==0)
			switch (word[1]) {
				case 'e': 
					empirical = true;
					flag |= flag_empirical_var;
					//throw Param_Error("Empirical test is oselete. Please use the pemutation test instead");
					break;
				case 'o': 
					optimize_offset = true;
					flag |= flag_optimize_offset;
					break;
				case 'c': 
					censored = true;
					flag |= flag_censored;
					if (trait_id[0]==0)
						throw Param_Error("affection status can't be the censored trait");
					break;
				default: throw Param_Error("unknown option");
			}
		else
			throw Param_Error("unknown option");
	}
	
	if (empirical && optimize_offset)
		throw Param_Error("Offset nuisance optimization is not available for empirical test");
	
	if (gen_model == model_genotype)
		throw Param_Error("The requested phased genotype test has not been implemented");
	
	if (n_sel_trt>1)
		throw Param_Error("The requested multivariate haplotype test has not been implemented");
	
	if (permutation && (empirical || optimize_offset || censored))
		throw Param_Error("the specified option is unavailable in the permutation test");
	
	print_settings();
	
	int nsel = j>=i? nLoci : i-j+1;
	int *idx = new int[nsel];
	if (j>i)
		for (nsel=0; nsel<nLoci; nsel++)
			idx[nsel] = nsel;
	else {
		for (nsel=0; j<=i; j++) {
			GetWord(s, j, word);
			if ((k=find_locus(word))>=0)
				idx[nsel++] = k;
			else {
				if (idx!=NULL)
					delete[] idx;
				char msg[256];
				sprintf(msg, "Marker %s not found", word);
				throw Param_Error(msg);
			}
		}
	}
	
	double gx[kMaxAlleles][3];
	if (gen_model==model_additive || gen_model==model_dominant || gen_model==model_recessive)
		for (i=0; i<kMaxAlleles; i++)
			gcode(&gx[i][0], gen_model);
	else
		throw Param_Error("Haplotype test under the current genetic model has not been implemented");
			
	if (permutation) {
		if (cycles>0)
			mc_hbat(nsel, idx, gx, flag, cycles);
		else
			mc_hbat(nsel, idx, gx, flag);
	} else 
		hdt(nsel, idx, gx, flag);
	
	MYOUTPUT(os, "\n");
	
	if (idx)
		delete[] idx;
	
}

void FBAT::hdt(int nmrk, int idx[], double gx[][3], int flag) {
	
	if (n_sel_trt>1)
		throw Param_Error("multivariate hbat is not available");
	
	double offset = ((flag&flag_censored)?uncensored_trait_mean(trait_id[0]) : tr_offset);
	if (flag&flag_censored)
		MYOUTPUT(os, "Uncensored sample trait mean = " << offset << "\n")
		
	Haplotype_List *emhl = hapfreq_ss(nmrk, idx);
	
	if (emhl==NULL) 
		return;
	
	MYOUTPUT(os, "\nhaplotype analysis for the following markers:\n")
	int i, j;
	for (i=0; i<nmrk; i++)
		MYOUTPUT(os, loc[idx[i]].name << " ")

	HaploidTable *htab = new HaploidTable(loc, nmrk, idx, emhl);
	
	MYOUTPUT(os, "\n\nhaplotypes and EM estimates of frequency:\n")
	i = 1;
	for (HaploidList *hpl2=htab->hpl; hpl2!=NULL; hpl2=hpl2->next,i++) {
		MYOUTPUT(os, "a" << setw(4) << i << setw(6) << ":")
		for (j=0; j<nmrk; j++)
			MYOUTPUT(os, setw(4) << loc[idx[j]].alleleName[hpl2->node->h[j]]);
		MYOUTPUT(os, setw(12) << hpl2->node->p << "\n")
	}

	int na = htab->size;

	if (na>kMaxAlleles)
	throw ERROR("too many alleles in haplotype test (max=80)");

	ColumnVector FS(na);
	ColumnVector PS(na);	// sum of S in whole pedigree, used in empirical var test
	ColumnVector FES(na);
	Matrix FV(na,na);

	Matrix A(na,na);	// these three are used for nuisance parameter etimation
	Matrix B(na,na);
	ColumnVector X(na);
	
	Matrix SA(na,na);
	Matrix SB(na,na);
	Matrix V(na,na);
	ColumnVector SX(na);
	ColumnVector S(na);	// FS-FES
	
	SA = 0;
	SB = 0;
	SX = 0;
	S = 0;
	V = 0;
	
	int tot_cnt = 0;
	int *fcnt = new int[na+1];
	for (i=0; i<=na; i++)
		fcnt[i] = 0;
	double *mu = new double[na+1];
	bool good;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		PS = 0;
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node && flist->node->stat) {
				tot_cnt++;
				flist->node->stat->analyze(emhl);
				if (flist->node->stat->informative() && flist->node->stat->fbat_stat3(trait_id[0], offset, FES, FS, FV, A, B, X, gx, flag, htab)) {
					FS -= FES;
					PS += FS;
					if (!(flag&flag_empirical_var)) {	// not empirical variance
						SA += A;
						SB += B;
						SX += X;
						S += FS;
						V += FV;
						if (!(flag&flag_optimize_offset)) {
							for (good=false,i=1; i<=na; i++)
								if (FV(i,i)>very_small_number) {
									fcnt[i]++;
									good = true;
								}
							if (good)	// informative for ar least one allele
								fcnt[0]++;
						}
					}
				}
			}
		}
		if (flag&flag_empirical_var) {
			FV = PS*PS.t();
			for (good=false,i=1; i<=na; i++)
				if (FV(i,i)>very_small_number) {
					fcnt[i]++;
					good = true;
				}
			if (good) {	// informative for ar least one allele
				fcnt[0]++;
				S += PS;
				V += FV;
			}
		}
	}
	
	if (flag&flag_optimize_offset) {	// update mu and fcnt
		double sa = 0;
		double sb = 0;
		for (i=1; i<=na; i++) {
			mu[i] = offset + SB(i,i)/SA(i,i);
			sa += SA(i,i);
			sb += SB(i,i);
		}
		mu[0] = offset + sb/sa;
		count_infomative_fam(mu, na, fcnt, gx, flag, htab);
		for (i=0; i<=na; i++)
			mu[i] -= offset;
	}

	if (mode==by_allele || mode==by_both) {
		MYOUTPUT(os, "\nAllele   afreq   fam#     S-E(S)      Var(S)       Z           P")
		if (flag&flag_optimize_offset)
			MYOUTPUT(os, "      Offset")
		MYOUTPUT(os, "\n----------------------------------------------------------------")
		if (flag&flag_optimize_offset)
			MYOUTPUT(os, "------------")
		MYOUTPUT(os, "\n")
		double z, p, st, var;
		char pstr[16];
		for (i=1; i<=na; i++) {
			double t = htab->hapfreq(i);
			if (fcnt[i]>=min_size && t>=min_freq) {
				if (flag&flag_optimize_offset) {
					st = S(i) - mu[i]*SX(i);
					var = V(i,i)-2*SB(i,i)*mu[i]+SA(i,i)*mu[i]*mu[i];
				} else {
					st = S(i);
					var = V(i,i);
				}
				z = st/sqrt(var);
				normbase(&z, &p);
				p = p>0.5?2*(1-p):2*p;
							
				if (alpha>0.99 ||  p<=alpha) {	//yes
					MYOUTPUT(os, "a" << setw(6) << i)
					os.setf(ios::right, ios::adjustfield);
					MYOUTPUT(os, setprecision(3) << setw(7) << t << setprecision(1) << setw(7) << fcnt[i])
					MYOUTPUT(os, setprecision(3) << setw(11) << st << setw(12) << var)
					pval2str(p, pstr);
					MYOUTPUT(os, setw(8) << z << setw(12) << pstr)
					if (flag&flag_optimize_offset) 
						MYOUTPUT(os, setw(12) << offset+mu[i])
					MYOUTPUT(os, "\n")
				}
				os.setf(ios::left, ios::adjustfield);
			}
		}
	}
	if (mode==by_marker || mode==by_both) {	// changed 5/13/03 to add by_both
		MYOUTPUT(os, "\nAllele#   Fam#     DF       CHISQ           P")
		if (flag&flag_optimize_offset)
			MYOUTPUT(os, "      Offset")
		MYOUTPUT(os, "\n---------------------------------------------")
		if (flag&flag_optimize_offset)
			MYOUTPUT(os, "------------")
		MYOUTPUT(os, "\n")		
		// calculate optimized_offset and update S and V
		if (flag&flag_optimize_offset) {
			S -= mu[0]*SX;
			V += mu[0]*mu[0]*SA - 2*mu[0]*SB;
		}
		int n = 0;
		int *index = new int[na+1];
		for (i=1; i<=na; i++)
			index[i] = (fcnt[i]>=min_size && htab->hapfreq(i)>=min_freq)? ++n : 0;		
		if (n<2)
			MYOUTPUT(os, "*** less than 2 major alleles ***\n")
		else {
			double df, chisq, p;
			bool err = false;
			try {
				p = calc_chisq(S, V, chisq, df, index);
			}
			catch (Exception error) {
				MYOUTPUT(os, "error in matrix inersion\n")
				err = true;
			}
			if (!err) {
				char pstr[16];
				os.setf(ios::right, ios::adjustfield);
				pval2str(p, pstr);
				MYOUTPUT(os, setw(7) << na << setw(7) << fcnt[0] << setw(7) << (int)(df+0.01) << setprecision(3) << setw(12) << chisq << setw(12) << pstr)
				if (flag&flag_optimize_offset)
					MYOUTPUT(os, setw(12) << offset+mu[0])
				MYOUTPUT(os, "\n")
				os.setf(ios::left, ios::adjustfield);
			}
		}
		if (index)
			delete[] index;
	} 
	
	if (emhl)
		delete emhl;
	
	if (mu)
		delete[] mu;
		
	if (fcnt)
		delete[] fcnt;

	if (htab)
		delete htab;
	
}

double FBAT::uncensored_trait_mean(int tr_idx) {
	
	double tsumaff = 0;
	int naff = 0;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next)
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0)
				for (PERSONLIST *sl=flist->node->sibs(); sl!=NULL; sl=sl->next)
					if (sl->node && sl->node->hastrait(tr_idx) && sl->node->affection==2) {
						tsumaff += sl->node->pt[tr_idx];
						naff++;
					}
		
	return tsumaff/naff;

}

	
// calculate UV(-)U in the submatrix indiced by idx (start from 1)
double FBAT::calc_chisq(ColumnVector &U, Matrix &V, double &chi, double &df, int idx[]) {	
		
	int na = V.Nrows();
	int i, j, n;
	for (i=1,n=0; i<=na; i++)
		if (idx[i]>0)
			n++;
	
	ColumnVector T(n);
	SymmetricMatrix S(n);			
	for (i=1; i<=na; i++)
		if (idx[i]>0) {
			T(idx[i]) = U(i);
			for (j=i; j<=na; j++)
				if (idx[j]>0)
					S(idx[i],idx[j]) = V(i,j);
		}
	
	Matrix R(n,n);
	DiagonalMatrix D(n);			
	Jacobi(S, D, R);
	df = 0;
	for (i=1; i<=n; i++)
		if (fabs(D(i,i))<another_very_small_number)
			D(i,i) = 0.0;
		else {
			D(i,i) = 1.0/D(i,i);
			df++;
		}

	Matrix Q = R * D * R.t();
	
	Matrix M = T.t() * Q * T;
	chi = M.AsScalar();
	chi *= 0.5;
	df *= 0.5;
	double p;
	gammabase(&chi, &df, &p);
	chi *= 2.0;
	df *= 2.0;
	
	return 1-p;

}

void FBAT::gcode(double *gx, GEN_MODEL model) {
	
	gx[0] = 0;
	gx[2] = 1;
	switch (model) {
		case model_additive:
			gx[1] = 1;
			gx[2] = 2;
			break;
		case model_dominant:
			gx[1] = 1;
			break;
		case model_recessive:
			gx[1] = 0;
			break;
		default:
			throw ERROR("specified genetic model is not recognized");
	}

}

void FBAT::count_infomative_fam(double mu[], int na, int *cnt, double gx[][3], int flag, HaploidTable *htab) {

	ColumnVector X(na);
	ColumnVector S(na);
	ColumnVector ES(na);
	Matrix V(na,na);	
	Matrix A(na,na);	// these six are used for nuisance parameter etimation
	Matrix B(na,na);
	
	int i;
	bool good;
	for (i=0; i<=na; i++)
		cnt[i] = 0;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				good = false;
				if ( flist->node->stat &&  flist->node->stat->informative() &&  flist->node->stat->fbat_stat3(trait_id[0], 0, ES, S, V, A, B, X, gx, flag, htab)) {
					for (i=1; i<=na; i++) {
						if (V(i,i)-2*mu[i]*B(i,i)+mu[i]*mu[i]*A(i,i) > very_small_number)
							cnt[i]++;
						if (!good && V(i,i)-2*mu[0]*B(i,i)+mu[0]*mu[0]*A(i,i) > very_small_number)
							good = true;
					}
					if (good)
						cnt[0]++;
				}
			}
		}
	}

}

void FBAT::fdt(int locidx, double gx[][3], int flag, int width, IRES* ires, bool rvtest, int rvoption, bool print_header) {

	if (locidx<0 || locidx>=nLoci)
		throw Param_Error("invalid locus id");
	
	int na = loc[locidx].nAlleles;
	double *frq = allelefrq(locidx);
	double offset = ((flag&flag_censored)?uncensored_trait_mean(trait_id[0]) : tr_offset);
		
	ColumnVector S(na);	// FS-FES
	Matrix V(na,na);	
	Matrix A(na,na);	// these six are used for nuisance parameter etimation
	Matrix B(na,na);
	ColumnVector X(na);

	int i;
	int *fcnt = new int[na+1];	// fam# that Var(S)!=0
	double *mu = new double[na+1];

	x_fbat_stat(locidx, trait_id[0], S, X, V, A, B, fcnt, gx, offset, flag, ires);
	
	if (flag&flag_optimize_offset) {	// update mu and fcnt
		double sa = 0;
		double sb = 0;
		for (i=1; i<=na; i++) {
			mu[i] = offset + B(i,i)/A(i,i);
			sa += A(i,i);
			sb += B(i,i);
		}
		mu[0] = offset + sb/sa;
		count_infomative_fam(mu, na, fcnt, gx, flag);
		for (i=0; i<=na; i++)
			mu[i] -= offset;
	}
	if (mode==by_allele || mode==by_both) {
		if (print_header) {
		  if (rvtest) {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "  afreq     fam#       weight     S-E(S)      Var(S)      Z        P   ")
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, "      Offset")
			if (gen_model==model_estimate)
				MYOUTPUT(os, "   Dominance")
			MYOUTPUT(os, "\n")
			for (i=0; i<72+width+((flag&flag_optimize_offset)?12:0)+((gen_model==model_estimate)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		  }
		  else {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "Allele   afreq   fam#     S-E(S)      Var(S)       Z           P")
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, "      Offset")
			if (gen_model==model_estimate)
				MYOUTPUT(os, "   Dominance")
			MYOUTPUT(os, "\n")
			for (i=0; i<64+width+((flag&flag_optimize_offset)?12:0)+((gen_model==model_estimate)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		  }
		}
		double z, p, st, var;
		char pstr[16];
		for (i=1; i<=na; i++) {
			if (fcnt[i]>=min_size && frq[i]>=min_freq) {
				if (flag&flag_optimize_offset) {
					st = S(i) - mu[i]*X(i);
					var = V(i,i)-2*B(i,i)*mu[i]+A(i,i)*mu[i]*mu[i];
				} else {
					st = S(i);
					var = V(i,i);
				}
				z = st/sqrt(var);
				normbase(&z, &p);
				p = p>0.5?2*(1-p):2*p;
							
				if ((alpha>0.99 ||  p<=alpha) && 
					((!rvtest) || (rvtest && (frq[i] < .5)))) {	//yes

					if (rvtest) 
					{
					  MYOUTPUT(os, setw(width) << loc[locidx].name )
					  os.setf(ios::right, ios::adjustfield);
					  MYOUTPUT(os, setprecision(3) << setw(7) << scientific << frq[i] << fixed << setprecision(1) << setw(7) << fcnt[i])
					}
					else
					{
					  MYOUTPUT(os, setw(width) << loc[locidx].name << setw(7) << loc[locidx].alleleName[i])
					  os.setf(ios::right, ios::adjustfield);
					  MYOUTPUT(os, setprecision(3) << setw(7) << frq[i] << setprecision(1) << setw(7) << fcnt[i])
					}
					if (rvtest) {
					  double weight = 1.0;
					  if (rvoption == $WT) 
					     weight = sqrt(ires->getfam_total()*frq[i]*(1-frq[i]));
					  ires->addwt(locidx, weight);
					  MYOUTPUT(os, setprecision(5) << setw(13) << 1/weight)
					}

					MYOUTPUT(os, setprecision(3) << setw(11) << st << setw(12) << var)
					pval2str(p, pstr);
					MYOUTPUT(os, setw(8) << z << setw(12) << pstr)
					if (flag&flag_optimize_offset) 
						MYOUTPUT(os, setw(12) << offset+mu[i])
					if (gen_model==model_estimate)
						MYOUTPUT(os, setw(12) << gx[i-1][1])
					MYOUTPUT(os, "\n")
				}
				os.setf(ios::left, ios::adjustfield);
			}
		}
	}
	if (mode==by_marker || mode==by_both) {	// changed 5/13/03 to add by_both
		if (print_header>0) {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "Allele#   Fam#     DF       CHISQ           P")
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, "      Offset")
			MYOUTPUT(os, "\n")
			for (i=0; i<45+width+((flag&flag_optimize_offset)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		}
		// calculate optimized_offset and update S and V
		if (flag&flag_optimize_offset) {
			S -= mu[0]*X;
			V += mu[0]*mu[0]*A - 2*mu[0]*B;
		}
		int n = 0;
		int *index = new int[na+1];
		for (i=1; i<=na; i++)
			index[i] = (fcnt[i]>=min_size && frq[i]>=min_freq)? ++n : 0;		
		if (n>1) {
			double df, chisq, p;
			bool err = false;
			try {
				p = calc_chisq(S, V, chisq, df, index);
			}
			catch (Exception error) {
				MYOUTPUT(os, "error in matrix inersion\n")
				err = true;
			}
			if (!err && (alpha>0.99 ||  p<=alpha)) {
				char pstr[16];
				MYOUTPUT(os, setw(width) << loc[locidx].name)
				os.setf(ios::right, ios::adjustfield);
				pval2str(p, pstr);
				MYOUTPUT(os, setw(7) << na << setw(7) << fcnt[0] << setw(7) << (int)(df+0.01) << setprecision(3) << setw(12) << chisq << setw(12) << pstr)
				if (flag&flag_optimize_offset)
					MYOUTPUT(os, setw(12) << offset+mu[0])
				MYOUTPUT(os, "\n")
				os.setf(ios::left, ios::adjustfield);
			}
		}
		if (index)
			delete[] index;
	}
	
	if (mu)
		delete[] mu;
		
	if (fcnt)
		delete[] fcnt;
	
}


void FBAT::fbat(char *s) {
	
	bool empirical = false, censored=false;
	bool optimize_offset = false;
	bool multimarker = false;
	bool lctest = false;
	bool mmtest = false;
	bool mptest = false;	// minP test
	bool lcmptest = false;	// lc test on multiple trait
	bool regression = false;
	bool rvtest = false;

	int rvoption = $UNW;	
	int genires = 0;	// 0- default: no extra output, 1- generate intermediate family statistics
	int old_minsize = 0;
	int ncycles = 100000;	// # of permutation for minP
	char word[1024], *c;
	int i, k;
	
	if (mode==by_both)
		throw Param_Error("by_both mode is only available for haplotype test");
		
	int wcnt = WordCount(s);
	
	int flag = 0;
	for (i=2; i<=wcnt && GetWord(s, i, word) && *word=='-'; i++) {
		if (word[2]==0)
			switch (word[1]) {
				case 'e': 
					empirical = true;
					flag |= flag_empirical_var;
					break;
				case 'o': 
					optimize_offset = true;
					flag |= flag_optimize_offset;
					break;
				case 'c': 
					flag |= flag_censored;
					censored = true; 
					if (trait_id[0]==0)
						throw Param_Error("affection can't be the censored trait");
					break;
				case 'm':
					multimarker = true;
					mmtest = true;
					break;
				case 'l':
					multimarker = true;
					lctest = true;
					break;
				case 'p':
					multimarker = true;
					mptest = true;
					break;
				case 'r':
					regression = true;
					break;
				case 's':
					genires = 1;
					break;
				case 't':
					multimarker = false;
					lcmptest = true;
					break;
				case 'v':
					/* rare variant collapsing method */
					rvtest = true;
					rvoption = $UNW;
					break;
				default: throw Param_Error("unknown option");
			}
		else if (word[1]=='p' && (k=strtol(word+2, &c, 10))>0) {
			multimarker = true;
			mptest = true;
			ncycles = k;
		}
		else if (word[1]=='v' && (k=strtol(word+2, &c, 10))>=0) {
			if (k == 0) rvoption = $UNW;
			else rvoption = $WT;
			rvtest = true;
		}
		else
			throw Param_Error("unknown option");
	}
	
	if (empirical && optimize_offset)
		throw Param_Error("Offset nuisance optimization is not available for empirical test");
	
	if (n_sel_trt>1 && (optimize_offset || censored))
		throw Param_Error("The requested multi-trait test has not been implemented");
	
	if (gen_model == model_genotype)
		throw Param_Error("The requested genotype test has not been implemented");
	
	if (gen_model == model_estimate && (multimarker || n_sel_trt>1))
		throw Param_Error("The model estimation is only available for single-marker univariate FBAT test");

	if (gen_model != model_additive && rvtest){
		MYOUTPUT(os, "current genetic model is " << model_name[gen_model] << "\n")
		throw Param_Error("Rare variant test can only be used for additive genetic model");
	}

	if ((lctest?1:0) + (mmtest?1:0) + (mptest?1:0) + (lcmptest?1:0) > 1)
		throw Param_Error("Different multi-marker/multi-trait tests cannot be combined");

	if (rvtest) {
	  old_minsize = this->min_size;
	  this->min_size = this->rv_min_size;
	}	
	print_settings();
	
	int nsel = i>=wcnt? nLoci : wcnt-i+1;
	int *idx = new int[nsel];
	if (i>wcnt)
		for (nsel=0; nsel<nLoci; nsel++)
			idx[nsel] = nsel;
	else {
		for (nsel=0; i<=wcnt; i++) {
				GetWord(s, i, word);
				if ((k=find_locus(word))>=0)
					idx[nsel++] = k;
				else {
					if (idx!=NULL)
						delete[] idx;
					char msg[256];
					sprintf(msg, "Marker %s not found", word);
					throw Param_Error(msg);
				}
		}
	}
	
	double gx[kMaxAlleles][3];
	if (gen_model==model_additive || gen_model==model_dominant || gen_model==model_recessive)
		for (i=0; i<kMaxAlleles; i++)
			gcode(&gx[i][0], gen_model);

	int mrkwidth = 0;
	for (i=0; i<nsel; i++)
		if ((k=strlen(loc[i].name)) > mrkwidth)
			mrkwidth = k;
	
	// If multiple trait are selected, do fbat-gee or lc for each of the selected marker
	if (n_sel_trt>1) {
		if (multimarker)
			throw Param_Error("The requested multi-trait multi-marker test has not been implemented");
		for (i=0; i<nsel; i++)
			if (!lcmptest)
				fbat_gee_mp(idx[i], gx, flag, i==0);
			else {
				for (k=1; k<=loc[idx[i]].nAlleles; k++)
					fbat_lc_test_mp(idx[i], k, gx, flag);
			}
	}
	else if (multimarker) {
		if (lctest)
			fbat_lc_test(nsel, idx, gx, flag);
		else if (mptest)
			fbat_minp_test(nsel, idx, gx, flag, ncycles);
		else
			fbat_mm_test(nsel, idx, gx, flag);
	}
	else if (nsel>0 && gen_model!=model_genotype) {
		double tval;
		IRES *ires=NULL;
		if (rvtest || genires !=0) {
		  int i=0;
		  for (PEDIGREELIST *ped = peds; ped != NULL; ped = ped->next) {
			i++;	
			//cout << i << ":" << ped->node->name << "\n";
		  }
		  ires = new IRES(this, empirical, peds, i);
		  if (rvtest) ires->save_minsize(old_minsize);
		}
		for (i=0; i<nsel; i++) {
			if (gen_model==model_estimate) {
				for (k=0; k<loc[idx[i]].nAlleles; k++) {
					gx[k][0] = 0;
					gx[k][2] = 1;
					gx[k][1] = estimate_dominance_cm(idx[i],k+1,trait_id[0],flag,tr_offset,tval);
				}
			}
			if (regression)
				fdt_reg(idx[i], gx, flag, mrkwidth+4, i==0);
			else 
				fdt(idx[i], gx, flag, mrkwidth+4, ires, rvtest, rvoption, i==0);
		}
		if (ires != NULL) {
			if (genires != 0) ires->print();
			if (rvtest) {
			  ires->rvtest(rvoption);
			  this->min_size = ires->restore_minsize();
			}
			delete ires;
		}
		MYOUTPUT(os, "\n")
	}
	
	if (idx)
		delete[] idx;
}

// estimate model parameter (dominance factor) using the conditional mean method
// EY = a + b1*EX + b2*P(X=1), where EX is the expected number of allele cnt.
// A b2 deviated from 0 suggests dominance component
double FBAT::estimate_dominance_cm(int loc_index, int a_idx, int trt_id, int flag, double offset, double &tval) {

	//if (loc[loc_index].nAlleles!=2)
	//	throw ("FBAT::estimate_dominance_cm works for bi-allelic markers only");
		
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	PERSONLIST *sl;
	
	bool censored = flag&flag_censored;
	
	tval = 0;
	
	int n = 0;
	// first look at the size of non-missing data
	for (ped=peds; ped!=NULL; ped=ped->next) {
		for (sl=ped->node->members; sl!=NULL; sl=sl->next)	// root members
			if (sl->node && sl->node->testroot() && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0) && sl->node->good_gt(loc_index, loc[loc_index].sexlinked()))
				n++;
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
			for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next)
				if (sl->node && sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0))
					n++;
	}
	
	ColumnVector Y(n);
	Matrix X(n,2); // X(,1) is the allele (allele a_idx) cnt, X(,2) = 1 if allele cnt = 1, 0 otherwise (or the probability of X=1)
	
	Y = 0;
	X = 0;
	
	int i = 0;
	int j;
	double x_prob[2];	// prob for allele count 1 and 2
	// assign value to Y and X
	for (ped=peds; ped!=NULL; ped=ped->next) {
		for (sl=ped->node->members; sl!=NULL; sl=sl->next)	// root members
			if (sl->node && sl->node->testroot() && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0) && sl->node->good_gt(loc_index, loc[loc_index].sexlinked())) {
				i++;
				Y(i) = sl->node->pt[trt_id] - ((censored && sl->node->affection!=2)? 0 : offset);
				j = sl->node->gt[loc_index].allele_cnt(a_idx, (loc[loc_index].sexlinked() && sl->node->sex==1));
				if (loc[loc_index].sexlinked() && sl->node->sex==1 && j==1)
					j = 2;
				X(i,1) = j;
				X(i,2) = (j==1?1:0);
			}
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			bool ss_infomative = flist->node && flist->node->stat && flist->node->stat->marker_prob(a_idx, x_prob);
			for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next)
				if (sl->node && sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
					i++;
					Y(i) = sl->node->pt[trt_id] - ((censored && sl->node->affection!=2)? 0 : offset);
					if (ss_infomative) {
						X(i,1) = x_prob[0] + 2*x_prob[1];
						X(i,2) = x_prob[0];
					} else {
						j = sl->node->gt[loc_index].allele_cnt(a_idx, (loc[loc_index].sexlinked() && sl->node->sex==1));
						if (loc[loc_index].sexlinked() && sl->node->sex==1 && j==1)
							j = 2;
						X(i,1) = j;
						X(i,2) = (j==1?1:0);
					}
				}
		}
	}
	
	X.Column(1) -= X.Column(1).Sum()/n;
	X.Column(2) -= X.Column(2).Sum()/n;
	
	Matrix T = X.Column(1).t() * Y;
	double b1 = T.AsScalar()/X.Column(1).SumSquare();
	
	Y -= Y.Sum()/n + b1*X.Column(1);	// Y is now residual from first regression
	T = X.Column(2).t() * Y;
	double xx2 = X.Column(2).SumSquare();
	double b2 = T.AsScalar()/xx2;
	Y -= b2*X.Column(1);
	double sse = Y.SumSquare();
	tval = b2/sqrt(sse/(xx2*(n-3)));
	
	if (b1==0)	// no effect, return additive moddel
		return 0.5;
		
	double d = (b1+b2)/(2*b1);
	// restrict d to [0,1]
	if (d<0)
		d = 0;
	else if (d>1)
		d = 1;
	
	return d;
}

// fitting Y_ij = u_i + e_ij, u_i=N(0,Vc), e_ij=N(0,Ve), return Vt=Vc+Ve, rho=Vc/Vt. Using all samples
double FBAT::tr_anova(int trt_id, double &rho) {
		
	double tsum, tssq, ssb, ssw, sst, t, y, yy;
	int n, k, m, mm;
	
	tsum = 0;
	tssq = 0;
	ssb = 0;
	n = 0;
	k = 0;
	mm = 0;
		
	PEDIGREELIST *ped;
	PERSONLIST *sl;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		m = 0;
		y = 0;
		yy = 0;
		for (sl=ped->node->members; sl!=NULL; sl=sl->next) {	// root members
			if (sl->node && sl->node->hastrait(trt_id)) {
				t = sl->node->pt[trt_id];
				y += t;
				yy += t*t;
				m++;
			}
		}
		if (m>0) {
			tsum += y;
			tssq += yy;
			ssb += y*y/m;
			n += m;
			mm += m*m;
			k++;
		}
	}
	sst = tssq - tsum*tsum/n;
	ssb -= tsum*tsum/n;
	ssw = sst - ssb;
	double q = (n-mm/n)/(k-1);
	double ve = ssw/(n-k);
	double vc = (ssb/(k-1)-ve)/q;
	rho = vc/(vc+ve);
	
	return vc+ve;
	
}

// variance weighted sample mean, where variance is the sum of variances of all markers 
// if a_idx is 0, use the trace of V as weight (sum of V for all alleles), otherwise use V for te particular allele specified
// return offset = vtsum/vsum

double FBAT::estimate_offset(int loc_index, int a_idx, double gscore[][3], double &vsum) {
		
	if (loc_index<0 || loc_index>=nLoci)
		throw Param_Error("invalid locus ID");
	if (a_idx<0 || a_idx>loc[loc_index].nAlleles)
		throw Param_Error("invalid allele ID");

	vsum =0;
	double vtsum = 0;
	
	int na = loc[loc_index].nAlleles;
			
	ColumnVector FS(na);
	ColumnVector FES(na);
	Matrix FV(na,na);
	Matrix V(na,na);
	
	Matrix A(na,na);	// these six are used for nuisance parameter etimation
	Matrix B(na,na);
	ColumnVector X(na);

	Sufficient_Stat *ss;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				if (flist->node->stat!=NULL) {
					delete flist->node->stat;
					flist->node->stat = NULL;
				}
				ss = new Sufficient_Stat(flist->node, 1, &loc_index, loc[loc_index].sexlinked());
				flist->node->stat = ss;
				if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
					continue;
				
				ss->get_common_ogt();
				ss->analyze();
				if (ss->informative() && ss->fbat_stat3(trait_id[0], 0, FES, FS, FV, A, B, X, gscore, 0)) {
					vsum += (a_idx>0?A(a_idx,a_idx):A.Trace());
					vtsum += (a_idx>0?B(a_idx,a_idx):B.Trace());						
				}
			}
		}
	}
					
	return vsum>very_small_number? vtsum/vsum : 0;

}

double FBAT::estimate_offset_mm(int nmrk, int loc_index[], int a_idx[], double gscore[][3]) {
	
	int i;
	double vsum = 0;
	double vtsum = 0;
	double offset, v;
	for (i=0; i<nmrk; i++) {
		offset = estimate_offset(loc_index[i], a_idx[i], gscore, v);
		if (v>very_small_number) {
			vsum += v;
			vtsum += v*offset;
		}
	}
	
	return vsum>very_small_number? vtsum/vsum : 0;

}

void FBAT::var2cor(Matrix &V) {	
	
	if (V.Nrows() != V.Ncols())
		throw ERROR("Variance matrix dimension error");
	
	int i, j;
	int nmrk = V.Nrows();
	for (i=1; i<=nmrk; i++)
		for (j=i+1; j<=nmrk; j++)
			if (V(i,i)>very_small_number && V(j,j)>very_small_number)
				V(i,j) = V(j,i) = V(i,j)/sqrt(V(i,i)*V(j,j));
	for (i=1; i<=nmrk; i++)
		V(i,i) = (V(i,i)>very_small_number? 1 : 0);

}

void FBAT::var2cor(SymmetricMatrix &V) {	
	
	if (V.Nrows() != V.Ncols())
		throw ERROR("Variance matrix dimension error");
	
	int i, j;
	int nmrk = V.Nrows();
	for (i=1; i<=nmrk; i++)
		for (j=i+1; j<=nmrk; j++)
			if (V(i,i)>very_small_number && V(j,j)>very_small_number)
				V(i,j) /= sqrt(V(i,i)*V(j,j));
	for (i=1; i<=nmrk; i++)
		V(i,i) = (V(i,i)>very_small_number? 1 : 0);

}

// multi-marker fbat test (T_mm);
// only additive model will be use and only allow for 2 alleles for each marker
// for single trait trait_id[0]
int FBAT::fbat_mm_test(int nmrk, int loc_idx[], double gscore[][3], int flag) {
	
	bool censored = flag&flag_censored;
	bool empirical = flag&flag_empirical_var;
	bool optimize_offset = flag&flag_optimize_offset;

	if (nmrk<1 || loc_idx==NULL)
		throw Param_Error("empty marker set");
	
	if (censored && optimize_offset)
		throw Param_Error("the request option has not been implemented");
	
	if (gen_model!=model_additive)
		throw Param_Error("T_mm works with additive genetic model only");
	
	int i;
	
	for (i=0; i<nmrk; i++) {
		if (loc_idx[i]<0 || loc_idx[i]>=nLoci)
			throw Param_Error("invalid locus index");
		if (loc[loc_idx[i]].nAlleles!=2)
			throw Param_Error("only bi-allelic markers are supported");
	}
	
	int *alle = new int[nmrk+1];		// allele to be tested, default to be the first allele
	int *info_cnt = new int[nmrk];	// informative fam#
	for (i=0; i<nmrk; i++) {
		info_cnt[i] = 0;
		alle[i] = 1;
	}
	
	double offset = censored? uncensored_trait_mean(trait_id[0]) : tr_offset;
	
	bool xlinked = loc[loc_idx[0]].sexlinked();
	for (i=1; i<nmrk; i++)
		if (loc[loc_idx[i]].sexlinked() != xlinked)
			throw Param_Error("inconsistent sexlink code for the selected markers");
			
	if (!censored && optimize_offset) {
		offset = estimate_offset_mm(nmrk, loc_idx, alle, gscore);
		MYOUTPUT(os, "offset estimated & used = " << offset << "\n")
	}

	Matrix VE(nmrk, nmrk); VE = 0.0;	// empirical variance matrix
	ColumnVector V(nmrk); V = 0.0;
	ColumnVector S(nmrk); S = 0.0;				// = S - ES
	ColumnVector FS(nmrk);
	
	double v, fs;
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	Sufficient_Stat *ss;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		FS = 0.0;
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				for (i=0;i<nmrk; i++) {
					if (flist->node->stat!=NULL) {
						delete flist->node->stat;
						flist->node->stat = NULL;
					}
					ss = new Sufficient_Stat(flist->node, 1, loc_idx+i, xlinked);
					flist->node->stat = ss;
					if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
						continue;
					ss->get_common_ogt();
					ss->analyze();
					if (ss->informative()) {
						v = ss->fbat_stat1(trait_id[0], loc[loc_idx[i]].nAlleles, alle[i], offset, gscore, flag, fs);
						if (v>very_small_number & !empirical)
							info_cnt[i]++;
						V(i+1) += v;
						FS(i+1) += fs;
					}
				}
			}
		S += FS;
		VE += FS*FS.t();
		
		if (empirical)
			for (i=0; i<nmrk; i++)
				if (fabs(FS(i+1))>very_small_number)
					info_cnt[i]++;
	}
	
	if (empirical)
		for (i=1; i<=nmrk; i++)
			V(i) = VE(i,i);
	int na = 0;
	double *frq;
	MYOUTPUT(os, "Marker  Allele  Freq  fam#       S-ES          Z\n")
	for (i=0; i<nmrk; i++) {
		MYOUTPUT(os, setw(10) << loc[loc_idx[i]].name)
		os.setf(ios::right, ios::adjustfield);
		MYOUTPUT(os, setw(4) << loc[loc_idx[i]].alleleName[1])
		frq = allelefrq(loc_idx[i]);
		MYOUTPUT(os, setprecision(3) << setw(6) << frq[1] << setw(6) << info_cnt[i] << setw(11) << S(i+1))
		if (V(i+1)>very_small_number && info_cnt[i]>=min_size && frq[1]>=min_freq) {
			S(i+1) /= sqrt(V(i+1));
			MYOUTPUT(os, setw(11) << S(i+1) << "\n")
			alle[i+1] = ++na;
		} else {
			MYOUTPUT(os, setw(11) << "---\n")
			alle[i+1] = -1;
		}
		os.setf(ios::left, ios::adjustfield);
	}

	if (na>0) {
		double df, chisq, p;
		bool err = false;
		var2cor(VE);
		try {
			p = calc_chisq(S, VE, chisq, df, alle);
		}
		catch (Exception error) {
			MYOUTPUT(os, "error in matrix inversion\n")
			err = true;
		}
		if (!err && (alpha>0.99 ||  p<=alpha)) {	
			char pstr[16];
			pval2str(p, pstr);
			MYOUTPUT(os, "FBAT_MM test: S_MM = " << setprecision(3) << chisq << "; df = " << df << "; p_value = " << pstr << "\n")
		}
	}

	if (info_cnt!=NULL)
		delete[] info_cnt;
	if (alle!=NULL)
		delete[] alle;
		
	return na;
}

// FBAT_LC test
// for single trait trait_id[0]
int FBAT::fbat_lc_test(int nmrk, int loc_idx[], double gscore[][3], int flag) {
	
	bool censored = flag&flag_censored;
	bool empirical = flag&flag_empirical_var;
	bool optimize_offset = flag&flag_optimize_offset;
	
	if (nmrk<1 || loc_idx==NULL)
		throw Param_Error("empty marker set");
	
	int i, j;
	
	for (i=0; i<nmrk; i++) {
		if (loc_idx[i]<0 || loc_idx[i]>=nLoci)
			throw Param_Error("invalid locus index");
		if (loc[loc_idx[i]].nAlleles!=2)
			throw Param_Error("only bi-allelic markers are supported");
	}
	
	bool xlinked = loc[loc_idx[0]].sexlinked();
	for (i=1; i<nmrk; i++)
		if (loc[loc_idx[i]].sexlinked() != xlinked)
			throw Param_Error("inconsistent sexlink code for the selected markers");

	double offset = censored? uncensored_trait_mean(trait_id[0]) : tr_offset;
	ColumnVector B(nmrk);		// betas for the allele in the condiitonal mean model, exclude random sample
	
	int *alle = new int[nmrk];
	double b2;
	for (i=0; i<nmrk; i++) {
		B(i+1) = pop_regression(loc_idx[i], 1, trait_id[0], flag, offset, gscore, true);
		alle[i] = 1;
		if (gen_model!=model_additive) {
			b2 = pop_regression(loc_idx[i], 2, trait_id[0], flag, offset, gscore, true);
			if (b2 > B(i+1)) {
				alle[i] = 2;
				B(i+1) = b2;
			}
		} 
		else if (B(i+1)<0) {
			B(i+1) = -B(i+1);
			alle[i] = 2;
		}
		
	}
		
	if (!censored && optimize_offset) {
		offset = estimate_offset_mm(nmrk, loc_idx, alle, gscore);
		MYOUTPUT(os, "offset estimated & used = " << offset << "\n")
	}

	int *info_cnt = new int[nmrk];	// informative fam#
	for (i=0; i<nmrk; i++)
		info_cnt[i] = 0;

	Matrix VE(nmrk, nmrk); VE = 0.0;	// empirical variance matrix
	ColumnVector V(nmrk); V = 0.0;
	ColumnVector S(nmrk); S = 0.0;				// = S - ES
	ColumnVector FS(nmrk);
		
	double v, fs;
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	Sufficient_Stat *ss;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		FS = 0.0;
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				for (i=0;i<nmrk; i++) {
					if (flist->node->stat!=NULL) {
						delete flist->node->stat;
						flist->node->stat = NULL;
					}
					ss = new Sufficient_Stat(flist->node, 1, loc_idx+i, xlinked);
					flist->node->stat = ss;
					if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
						continue;
					ss->get_common_ogt();
					ss->analyze();
					if (ss->informative()) {
						v = ss->fbat_stat1(trait_id[0], loc[loc_idx[i]].nAlleles, alle[i], offset, gscore, flag, fs);
						if (v>very_small_number & !empirical)
							info_cnt[i]++;
						V(i+1) += v;
						FS(i+1) += fs;
					}
				}
			}
		S += FS;
		VE += FS*FS.t();
		
		if (empirical)
			for (i=0; i<nmrk; i++)
				if (fabs(FS(i+1))>very_small_number)
					info_cnt[i]++;
	}
	
	if (empirical)
		for (i=1; i<=nmrk; i++)
			V(i) = VE(i,i);

	int na = 0;
	double *frq;
	MYOUTPUT(os, "Marker  Allele  Freq  fam#       S-ES          Z          W\n")
	for (i=0; i<nmrk; i++) {
		MYOUTPUT(os, setw(10) << loc[loc_idx[i]].name)
		os.setf(ios::right, ios::adjustfield);
		MYOUTPUT(os, setw(4) << loc[loc_idx[i]].alleleName[alle[i]])
		frq = allelefrq(loc_idx[i]);
		MYOUTPUT(os, setprecision(3) << setw(6) << frq[alle[i]] << setw(6) << info_cnt[i] << setw(11) << S(i+1))
		if (V(i+1)>very_small_number && info_cnt[i]>=min_size && frq[alle[i]]>=min_freq) {
			S(i+1) = S(i+1)/sqrt(V(i+1));
			MYOUTPUT(os, setw(11) << S(i+1) << setw(11) << B(i+1) << "\n")
			alle[i] = ++na;
		} else {
			MYOUTPUT(os, setw(11) << "---\n")
			alle[i] = -1;
		}
		os.setf(ios::left, ios::adjustfield);
	}
	
	// delete minor alleles
	ColumnVector SMM(na);
	ColumnVector W(na);
	SymmetricMatrix VMM(na);
	for (i=1; i<=nmrk; i++)
		if (alle[i-1]>0) {
			SMM(alle[i-1]) = S(i);
			W(alle[i-1]) = B(i);
			for (j=i; j<=nmrk; j++)
				if (alle[j-1]>0)
					VMM(alle[i-1],alle[j-1]) = VE(i,j);
		}
	
	var2cor(VMM);
	
	for (fs=0,i=1; i<=na; i++)
		fs += W(i)*SMM(i);
	Matrix VLC = W.t() * VMM * W;	// weighted variance
	double z = fs/sqrt(VLC.AsScalar());
	double plc;
	normbase(&z, &plc);
	plc = 1 - plc;
	char pstr[16];
	pval2str(plc, pstr);
	MYOUTPUT(os, "FBAT_LC test: Z_LC = " << setprecision(3) << z << "; p_value = " << pstr << "\n")
	
	if (info_cnt!=NULL)
		delete[] info_cnt;
	if (alle!=NULL)
		delete[] alle;
		
	return na;

}

double FBAT::pop_regression(int loc_index, int a_idx, int trt_id, int flag, double offset, double gscore[][3], bool use_cmm, double *b, double *v) {


	double y = 0.0;
	double x = 0.0;
	double xx = 0.0;
	double xy = 0.0;
	double nx = 0.0;
	double yy = 0.0;
	
	int n;
		
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	PERSONLIST *sl;
	
	Sufficient_Stat *ss;

	bool censored = flag&flag_censored;
	
	double t, xm, xms[2];
	for (ped=peds; ped!=NULL; ped=ped->next) {
		if (use_cmm) {
			for (sl=ped->node->members; sl!=NULL; sl=sl->next)
				sl->node->resetflag(flag_processed);
			for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
					if (flist->node->stat!=NULL) {
						delete flist->node->stat;
						flist->node->stat = NULL;
					}
					ss = new Sufficient_Stat(flist->node, 1, &loc_index, loc[loc_index].sexlinked());
					flist->node->stat = ss;
					if (ss==NULL || ss->mhap==NULL || ss->get_common_ogt()<2)
						continue;
					ss->analyze();
					if (!ss->informative())
						continue;
					// modified 11/16/2007 to allow sex-specific ex for X-markers
					ss->ex(xms, a_idx, gscore);
					for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next)
						if (sl->node && sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
							t = sl->node->trait(trt_id) - ((censored && sl->node->affection!=2)? 0 : offset);
							xm = xms[sl->node->sex==1?0:1];
							x += xm;
							y += t;
							xy += xm*t;
							xx += xm*xm;
							yy += t*t;
							nx += 1;
							sl->node->setflag(flag_processed);
						}
				}
			}
		}
		for (sl=ped->node->members; sl!=NULL; sl=sl->next)
			if (sl->node && (!use_cmm || !sl->node->testflag(flag_processed))) {
				if (sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
					t = sl->node->trait(trt_id) - ((censored && sl->node->affection!=2)? 0 : offset);
					n = sl->node->gt[loc_index].allele_cnt(a_idx,(loc[loc_index].sexlinked() && sl->node->sex==1));
					xm = (loc[loc_index].sexlinked() && sl->node->sex==1? (n>0?1:0) : gscore[a_idx-1][n]);
					x += xm;
					y += t;
					xy += xm*t;
					xx += xm*xm;
					yy += t*t;
					nx += 1;
				}
			}
	}
	
	if (nx*yy-y*y<very_small_number)
		throw ERROR("no variation of phenotype");
	
	double beta = (nx>1 && fabs(nx*xx-x*x)>very_small_number)? (xy-x*y/nx)/(xx-x*x/nx) : 0;
	double alfa = (y - beta*x)/nx;

	//double var = ((yy-y*y/nx)-beta*(xy-x*y/nx))/((nx-2)*(xx-x*x/nx));
	// calculate robust sandwitch variance
	double xmean = x/nx;
	double exex = 0;
	double fex;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		fex = 0;
		if (use_cmm) {
			for (sl=ped->node->members; sl!=NULL; sl=sl->next)
				sl->node->resetflag(flag_processed);
			for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0 && (ss=flist->node->stat)!=NULL && ss->informative()) {
					// modified 11/16/2007 to allow sex-specific ex for X-markers
					ss->ex(xms, a_idx, gscore);
					for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next)
						if (sl->node && sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
							t = sl->node->trait(trt_id) - ((censored && sl->node->affection!=2)? 0 : offset);
							xm = xms[sl->node->sex==1?0:1];
							fex += (t - alfa - beta*xm)*(xm-xmean);
							sl->node->setflag(flag_processed);
						}
				}
			}
		}
		for (sl=ped->node->members; sl!=NULL; sl=sl->next)
			if (sl->node && (!use_cmm || !sl->node->testflag(flag_processed))) {
				if (sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
					t = sl->node->trait(trt_id) - ((censored && sl->node->affection!=2)? 0 : offset);
					n = sl->node->gt[loc_index].allele_cnt(a_idx,(loc[loc_index].sexlinked() && sl->node->sex==1));
					xm = (loc[loc_index].sexlinked() && sl->node->sex==1? (n>0?1:0) : gscore[a_idx-1][n]);
					fex += (t - alfa - beta*xm)*(xm-xmean);
				}
			}
		
		exex += fex*fex;
	}		
			
	t = xx - nx*xmean*xmean;
	double var = exex/(t*t);
	double z = beta/sqrt(var);
	
	if (b)
		*b = beta;
	if (v)
		*v = var;
		
	return z;
}

// multi-marker minP fbat test (T_mp);
// only additive model will be use and only allow for 2 alleles for each marker
int FBAT::fbat_minp_test(int nmrk, int loc_idx[], double gscore[][3], int flag, int ncycles) {
	
	bool censored = flag&flag_censored;
	bool empirical = flag&flag_empirical_var;
	bool optimize_offset = flag&flag_optimize_offset;
	
	if (nmrk<1 || loc_idx==NULL)
		throw Param_Error("empty marker set");
	
	if (censored && optimize_offset)
		throw Param_Error("the request option has not been implemented");
	
	if (gen_model!=model_additive)
		throw Param_Error("T_mp works with additive genetic model only");
	
	int i, j;
	
	for (i=0; i<nmrk; i++) {
		if (loc_idx[i]<0 || loc_idx[i]>=nLoci)
			throw Param_Error("invalid locus index");
		if (loc[loc_idx[i]].nAlleles!=2)
			throw Param_Error("only bi-allelic markers are supported");
	}
	
	int *alle = new int[nmrk+1];		// allele to be tested, default to be the first allele
	int *info_cnt = new int[nmrk];	// informative fam#
	for (i=0; i<nmrk; i++) {
		info_cnt[i] = 0;
		alle[i] = 1;
	}

	double offset = censored? uncensored_trait_mean(trait_id[0]) : tr_offset;
		
	if (!censored && optimize_offset) {
		offset = estimate_offset_mm(nmrk, loc_idx, alle, gscore);
		MYOUTPUT(os, "offset estimated & used = " << offset << "\n")
	}

	Matrix VE(nmrk, nmrk); VE = 0.0;	// empirical variance matrix
	ColumnVector V(nmrk); V = 0.0;
	ColumnVector S(nmrk); S = 0.0;				// = S - ES
	ColumnVector FS(nmrk);
	
	double v, fs;
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	Sufficient_Stat *ss;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		FS = 0.0;
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next)
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				for (i=0;i<nmrk; i++) {
					if (flist->node->stat!=NULL) {
						delete flist->node->stat;
						flist->node->stat = NULL;
					}
					ss = new Sufficient_Stat(flist->node, 1, loc_idx+i, loc[loc_idx[i]].sexlinked());
					flist->node->stat = ss;
					if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
						continue;
					ss->get_common_ogt();
					ss->analyze();
					if (ss->informative()) {
						v = ss->fbat_stat1(trait_id[0], loc[loc_idx[i]].nAlleles, alle[i], offset, gscore, flag, fs);
						if (v>very_small_number & !empirical)
							info_cnt[i]++;
						V(i+1) += v;
						FS(i+1) += fs;
					}
				}
			}
		S += FS;
		VE += FS*FS.t();
		
		if (empirical)
			for (i=0; i<nmrk; i++)
				if (fabs(FS(i+1))>very_small_number)
					info_cnt[i]++;
	}
	
	if (empirical)
		for (i=1; i<=nmrk; i++)
			V(i) = VE(i,i);

	int na = 0;
	double *frq;
	MYOUTPUT(os, "Marker  Allele  Freq  fam#       S-ES          Z\n")
	for (i=0; i<nmrk; i++) {
		MYOUTPUT(os, setw(10) << loc[loc_idx[i]].name)
		os.setf(ios::right, ios::adjustfield);
		MYOUTPUT(os, setw(4) << loc[loc_idx[i]].alleleName[alle[i]])
		frq = allelefrq(loc_idx[i]);
		MYOUTPUT(os, setprecision(3) << setw(6) << frq[alle[i]] << setw(6) << info_cnt[i] << setw(11) << S(i+1))
		if (V(i+1)>very_small_number && info_cnt[i]>=min_size && frq[alle[i]]>=min_freq) {
			S(i+1) /= sqrt(V(i+1));
			MYOUTPUT(os, setw(11) << S(i+1) << "\n")
			alle[i] = ++na;
		} else {
			MYOUTPUT(os, setw(11) << "---\n")
			alle[i] = -1;
		}
		os.setf(ios::left, ios::adjustfield);
	}
	
	// delete minor alleles
	ColumnVector SMM(na);	// =Z
	SymmetricMatrix VMM(na);
	for (i=1; i<=nmrk; i++)
		if (alle[i-1]>0) {
			SMM(alle[i-1]) = S(i);
			for (j=i; j<=nmrk; j++)
				if (alle[j-1]>0)
					VMM(alle[i-1],alle[j-1]) = VE(i,j);
		}
	
	var2cor(VMM);
		
	double max_z = 0;
	for (i=1; i<=na; i++)
		if (fabs(SMM(i)) > max_z)
			max_z = fabs(SMM(i));
	
	double minp = empirical_min_p(max_z, VMM, ncycles);
	
	char pstr[16];
	pval2str(minp, pstr);
		
	MYOUTPUT(os, "FBAT_minP test: max_z = " << max_z << "; " << ncycles << " cycles; p_value = " << pstr << "\n")

	if (info_cnt!=NULL)
		delete[] info_cnt;
	if (alle!=NULL)
		delete[] alle;
		
	return na;
}

double FBAT::empirical_min_p(double maxz, SymmetricMatrix &V, int maxCycle = 1000000) {
	
	const int minCycle = 10000;
	const int minCnt = 100;
	
	int nz = V.Nrows();
	if (nz!=V.Ncols())
		throw("correlation matrix not othorgonal");
	
	ColumnVector Z(nz);
	
	int i;

	// derive the square root of V
	Matrix R(nz,nz);
	DiagonalMatrix D(nz);
	Jacobi(V, D, R);
	for (i=1; i<=nz; i++)
		if (D(i,i)>very_small_number)
			D(i,i) = sqrt(D(i,i));
	R = R*D*R.t();
	
	double mz;
	int ncycle;
	int pcnt = 0;
	for (ncycle=0; ncycle<minCycle || (pcnt<minCnt && ncycle<=maxCycle); ncycle++) {
		for (i=1; i<=nz; i++)
			Z(i) = rand_gaussian();
		Z = R*Z;
		for (mz=0,i=1; i<=nz; i++)
			if (fabs(Z(i))>mz)
				mz = fabs(Z(i));
		if (mz>=maxz)
			pcnt++;
	}
	
	return (double)(1+pcnt)/(double)(1+ncycle);
}

// FBAT_LC test
// for multiple trait defined in trait_id
int FBAT::fbat_lc_test_mp(int loc_idx, int a_idx, double gscore[][3], int flag) {
	
	bool censored = flag&flag_censored;
	bool empirical = flag&flag_empirical_var;
	bool optimize_offset = flag&flag_optimize_offset;
	
	if (loc_idx<0 || loc_idx>=nLoci)
		throw Param_Error("invalid loc_idx in FBAT::fbat_lc_test_mp");
	
	if (optimize_offset)
		throw Param_Error("T_LC_MVT currently does not support nuisance parameter estimation");
	
	if (censored)
		throw Param_Error("T_LC_MVT currently does not support censored traits");
		
//	if (print_header)
//		MYOUTPUT(os, "trait offset = 0 is used for all traits in T_LC_MVT\n")
	
	ColumnVector B(n_sel_trt);		// betas for the allele in the condiitonal mean model, exclude random sample
	int i, j;	
	for (i=0; i<n_sel_trt; i++)
		B(i+1) = pop_regression(loc_idx, a_idx, trait_id[i], flag, 0, gscore, true);
	
	int *info_cnt = new int[n_sel_trt];	// informative fam#
	for (i=0; i<n_sel_trt; i++)
		info_cnt[i] = 0;

	Matrix VE(n_sel_trt, n_sel_trt); VE = 0.0;	// empirical variance matrix
	Matrix FV(n_sel_trt, n_sel_trt);
	Matrix V(n_sel_trt, n_sel_trt); V = 0.0;
	ColumnVector S(n_sel_trt); S = 0.0;				// = S - ES
	ColumnVector FS(n_sel_trt);
	ColumnVector PS(n_sel_trt);
		
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	Sufficient_Stat *ss;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		PS = 0.0;
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				if (a_idx==1) {	// presume ss done for a_idx>1 alleles
					if (flist->node->stat!=NULL) {
						delete flist->node->stat;
						flist->node->stat = NULL;
					}
					ss = new Sufficient_Stat(flist->node, 1, &loc_idx, loc[loc_idx].sexlinked());
					flist->node->stat = ss;
				} else
					ss = flist->node->stat;
				if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
					continue;
				ss->get_common_ogt();
				ss->analyze();
				if (ss->informative() && ss->fbat_stat_mp(n_sel_trt, trait_id, a_idx, FS, FV, gscore, flag)) {
					PS += FS;
					V += FV;
					for (i=0; i<n_sel_trt; i++)
						if (FV(i+1,i+1)>very_small_number & !empirical)
							info_cnt[i]++;
				}
			}
		}
		S += PS;
		VE += PS*PS.t();
		
		if (empirical)
			for (i=0; i<n_sel_trt; i++)
				if (fabs(PS(i+1))>very_small_number)
					info_cnt[i]++;
	}
	
	if (empirical)
		V = VE;

	int nt = 0;
	double *frq;
	int *alle = new int[n_sel_trt];
	MYOUTPUT(os, "Marker  Allele  Freq       Trait  fam#       S-ES          Z          W\n")
	for (i=0; i<n_sel_trt; i++) {
		MYOUTPUT(os, setw(10) << loc[loc_idx].name)
		os.setf(ios::right, ios::adjustfield);
		MYOUTPUT(os, setw(4) << loc[loc_idx].alleleName[a_idx])
		frq = allelefrq(loc_idx);
		MYOUTPUT(os, setprecision(3) << setw(6) << frq[a_idx] << setw(12) << trt[trait_id[i]].name << setw(6) << info_cnt[i] << setw(11) << S(i+1))
		if (V(i+1,i+1)>very_small_number && info_cnt[i]>=min_size && frq[1]>=min_freq) {
			S(i+1) = S(i+1)/sqrt(V(i+1,i+1));
			MYOUTPUT(os, setw(11) << S(i+1) << setw(11) << B(i+1) << "\n")
			alle[i] = ++nt;
		} else {
			MYOUTPUT(os, setw(11) << "---\n")
			alle[i] = -1;
		}
		os.setf(ios::left, ios::adjustfield);
	}
	
	// delete minor alleles
	ColumnVector SMM(nt);
	ColumnVector W(nt);
	SymmetricMatrix VMM(nt);
	for (i=1; i<=n_sel_trt; i++)
		if (alle[i-1]>0) {
			SMM(alle[i-1]) = S(i);
			W(alle[i-1]) = B(i);
			for (j=i; j<=n_sel_trt; j++)
				if (alle[j-1]>0)
					VMM(alle[i-1],alle[j-1]) = V(i,j);
		}
	
	var2cor(VMM);
	
	double fs;
	for (fs=0,i=1; i<=nt; i++)
		fs += W(i)*SMM(i);
	Matrix VLC = W.t() * VMM * W;	// weighted variance
	double z = fs/sqrt(VLC.AsScalar());
	double plc;
	normbase(&z, &plc);
	plc = 1 - plc;
	char pstr[16];
	pval2str(plc, pstr);
	MYOUTPUT(os, "FBAT_LC_MVT test: Z_LC = " << setprecision(3) << z << "; p_value = " << pstr << "\n\n")
	
	if (info_cnt!=NULL)
		delete[] info_cnt;
	if (alle!=NULL)
		delete[] alle;
		
	return nt;

}

// monte-carlo test
int FBAT::mc_hbat(int nmrk, int idx[], double gx[][3], int flag, int max_runs) {

	const int flag_fam_info = 1<<17;
	const int kminobs = 100;

	if (nmrk<1 || idx==NULL)
		throw Param_Error("empty marker set");
	
	int i, j;

	double offset = ((flag&flag_censored)?uncensored_trait_mean(trait_id[0]) : tr_offset);
	if (flag&flag_censored)
		MYOUTPUT(os, "Uncensored sample trait mean = " << offset << "\n")

	Haplotype_List *emhl = hapfreq_ss(nmrk, idx);
	if (emhl==NULL) 
		return 0;
		
	HaploidTable *htab = new HaploidTable();
	htab->maketable(emhl);
	
	MYOUTPUT(os, "\nhaplotype permutation test for the following markers:\n")
	for (i=0; i<nmrk; i++)
		MYOUTPUT(os, loc[idx[i]].name << " ")
	MYOUTPUT(os, "\n")
	
	int na = htab->size;
	if (na>kMaxAlleles)
	throw ERROR("too many alleles in haplotype test (max=80)");

	ColumnVector FS(na);
	ColumnVector PS(na);	// sum of S in whole pedigree, used in empirical var test
	ColumnVector FES(na);
	ColumnVector X(na);
	Matrix FV(na,na);
	Matrix A(na,na);	// these three are used for nuisance parameter etimation
	Matrix B(na,na);
	Matrix V(na,na);
	ColumnVector S(na);	// FS-FES
	
	S = 0;
	V = 0;
	int *fcnt = new int[na+1];
	for (i=0; i<=na; i++)
		fcnt[i] = 0;
	bool good;
	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		PS = 0;
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			flist->node->resetflag(flag_fam_info);
			if (flist->node && flist->node->stat) {
				flist->node->stat->analyze(emhl);
				if (flist->node->stat->informative() && flist->node->stat->fbat_stat3(trait_id[0], offset, FES, FS, FV, A, B, X, gx, flag, htab)) {
					FS -= FES;
					PS += FS;
					if (!(flag&flag_empirical_var))	{	// not empirical variance
						V += FV;
						good = false;
						for (i=1; i<=na; i++)
							if (FV(i,i)>very_small_number) {
								fcnt[i]++;
								good = true;
							}
						if (good)
							fcnt[0]++;
					}
					for (i=1; i<=na; i++)
						if (FV(i,i)>very_small_number) {
							flist->node->setflag(flag_fam_info);
							flist->node->stat->init_permute_X(na, gx, htab);
							break;
						}
				}
			}
		}
		S += PS;
		if (flag&flag_empirical_var) {
			FV = PS*PS.t();
			V += FV;
			good = false;
			for (i=1; i<=na; i++)
				if (FV(i,i)>very_small_number) {
					fcnt[i]++;
					good = true;
				}
			if (good)
				fcnt[0]++;
		}
	}
	
	double *cnt = new double[na+1];
	for (i=0; i<=na; i++)
		cnt[i] = 0;
	
	double ts, tt;
	double maxchisq = 0;	// added 11/18/05 for min_p test (min_p is same as max Z^2)
	double minpcnt = 0;
	for (ts=0,i=1; i<=na; i++)
		if (V(i,i)>another_very_small_number) {
			tt = S(i)*S(i)/V(i,i);
			ts += tt;
			if (tt>maxchisq)
				maxchisq = tt;
		}
		
	double chisq, maxc;	// added 11/18/05 for minp
	ColumnVector T(na);
	int cycles=0;
	do {
		cycles++;
		T = 0;
		for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next)
			for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				if (flist->node==NULL || flist->node->testflag(flag_uninformative)!=0 || flist->node->testflag(flag_fam_info)==0)
					continue;
				flist->node->stat->permute(FS, trait_id[0], offset, flag);
				T += FS;
			}
		
		maxc = 0;
		for (tt=0,i=1; i<=na; i++) {
			if (V(i,i)>another_very_small_number) {
				chisq = T(i)*T(i)/V(i,i);
				tt += chisq;
				if (chisq>maxc)
					maxc = chisq;
			}
			if (S(i) == T(i))
				cnt[i] += 0.5;
			else if (S(i) > T(i))
				cnt[i] += 1;
//			if (cycles<100)
//				MYOUTPUT(os, T(1) << "\t" << T(2) << "\t" << T(1)+T(2) << "\n")
		}
		if (ts == tt)
			cnt[0] += 0.5;
		else if (ts>tt)
			cnt[0] += 1;
		
		// minp
		if (maxc==maxchisq)
			minpcnt += 0.5;
		else if (maxc<maxchisq)
			minpcnt += 1;
		
		for (i=0; i<=na && cnt[i]>=kminobs && cycles-cnt[i]>=kminobs; i++)
			;
		
	} while (i<=na && cycles<max_runs);
	
	double p;
	char pstr[16];
	MYOUTPUT(os, "Permutation cycles = " << cycles << "\n")
	MYOUTPUT(os, setw(11+4*nmrk) << "haplotype" << setw(12) << "afreq" << setw(10) << "fam#" << setw(16) << "rank(S_obs)" << "P_2side\n")
	i = 1;
	for (HaploidList *hpl2=htab->hpl; hpl2!=NULL; hpl2=hpl2->next,i++) {
		MYOUTPUT(os, "h" << setw(4) << i << setw(6) << ":")
		for (j=0; j<nmrk; j++)
			MYOUTPUT(os, setw(4) << loc[idx[j]].alleleName[hpl2->node->h[j]]);
		p = cnt[i]/cycles;
		p = p<0.5? 2*p : 2*(1-p);
		pval2str(p, pstr);
		if (fcnt[i]>0)
			MYOUTPUT(os,  setprecision(3) << setw(12) << hpl2->node->p << setw(10) << fcnt[i] << setw(16) << cycles-cnt[i] << pstr << "\n")
		else
			MYOUTPUT(os,  setprecision(3) << setw(12) << hpl2->node->p << setw(10) << "*** too few informative obs ***\n")
	}
	
	p = 1 - cnt[0]/cycles;
	
	for (i=1; i<=na && fcnt[i]<1; i++)
		;
	MYOUTPUT(os, "\n")
	if (i<=na) {
		pval2str(p, pstr);
		MYOUTPUT(os, "whole marker permutation test (chisq sum) P = " << pstr << "\n")
		pval2str(1.0-minpcnt/cycles, pstr);
		MYOUTPUT(os, "whole marker permutation test (minimal p) P = " << pstr << "\n")
	} else
		MYOUTPUT(os, setw(23+4*nmrk) << "whole marker" << setw(16) << "*** too few informative obs ***\n")
		
	delete[] cnt;
	
	delete[] fcnt;
	
	if (emhl)
		delete emhl;

	if (htab)
		delete htab;
	
	return 1;
}

// multi-allelic version for multiple traits, added under Nan's request, 02/15/2002
// modified 5/30/07, changed to fbat_gee
int FBAT::fbat_gee_mp(int loc_idx, double gscore[][3], int flag, bool print_header) {
	
	bool censored = flag&flag_censored;
	bool empirical = flag&flag_empirical_var;
	bool optimize_offset = flag&flag_optimize_offset;
	
	if (loc_idx<0 || loc_idx>=nLoci)
		throw Param_Error("invalid loc_idx in FBAT::fbat_lc_test_mp");
	
	if (optimize_offset)
		throw Param_Error("fbat_gee currently does not support nuisance parameter estimation");
	
	if (censored)
		throw Param_Error("fbat_gee currently does not support censored traits");
		
	if (print_header)
		MYOUTPUT(os, "trait offset = 0 is used for all traits in fbat_gee\n")
	
	int na = loc[loc_idx].nAlleles;	
	int ncols = n_sel_trt*na;
	int i, j, n;
	
	ColumnVector S(ncols), FS(ncols), PS(ncols), C(ncols), I(na);
	Matrix	V(ncols, ncols), FV(ncols, ncols);
	S = 0;
	V = 0;
	C = 0;
	I = 0;	//  fam# that are informative for >=1 traits
	
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	Sufficient_Stat *ss;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		PS = 0;
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				if (flist->node->stat!=NULL) {
					delete flist->node->stat;
					flist->node->stat = NULL;
				}
				ss = new Sufficient_Stat(flist->node, 1, &loc_idx, loc[loc_idx].sexlinked());
				flist->node->stat = ss;
				if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
					continue;
				ss->get_common_ogt();
				ss->analyze();
				if (ss->informative() && ss->fbat_gee_mp(na, n_sel_trt, trait_id, FS, FV, gscore, flag)) {
					PS += FS;
					if (!empirical) {
						V += FV;
						for (i=1; i<=ncols; i++)
							if (FV(i,i)>very_small_number)
								C(i) += 1;
						for (i=1; i<=na; i++) {
							for (j=0; j<n_sel_trt; j++)
								if (FV(j*na+i,j*na+i)>very_small_number) {
									I(i) += 1;
									break;
								}
						}
					}
				}
			}
		}
		S += PS;
		if (empirical) {
			V += PS*PS.t();
			for (i=1; i<=ncols; i++)
				if (fabs(PS(i))>very_small_number)
					C(i) += 1;
		}
	}
	/*
	cout << "C=\n" << C << "\n";
	cout << "I=\n" << I << "\n";
	cout << "V=\n" << V << "\n";
	*/
	
	int *idx = new int[ncols+1];
	double chi=0, df=0, p=0;
	char pstr[16];
	
	if (mode == by_marker || mode == by_both) {
		if (print_header)
			MYOUTPUT(os, "Marker   Allele#    DF       CHISQ           P\n")
		for (n=0,i=1; i<=ncols; i++)
			idx[i] = (C(i)>min_size-0.1? ++n : 0);
		if (n<2) {
			delete[] idx;
			MYOUTPUT(os, setw(10) << loc[loc_idx].name << "*** less than 2 major groups ***\n");
			return 0;
		}		
		calc_chisq(S, V, chi, df, idx);		
		chi *= 0.5;
		df *= 0.5;
		gammabase(&chi, &df, &p);
		p = 1.0 - p;
		chi *= 2.0;		
		if (p<=alpha) {
			MYOUTPUT(os, setw(10) << loc[loc_idx].name)
			os.setf(ios::right, ios::adjustfield);
			pval2str(p, pstr);
			MYOUTPUT(os, setw(6) << na << setw(6) << (int)(2*df+0.01) << setprecision(3) << setw(12) << chi << setw(12) << pstr << "\n")
			os.setf(ios::left, ios::adjustfield);
		}
	} else {
		if (print_header)
			MYOUTPUT(os, "Marker  Allele  afreq fam#       DF       CHISQ           P\n")
		double *frq = allelefrq(loc_idx);
		//double cnt = 0;
		for (j=1; j<=na; j++) {
			for (i=1; i<=ncols; i++)
				idx[i] = 0;
			for (n=0,i=j; i<=ncols; i+=na) {
				if (C(i)>min_size-0.1)
					idx[i] = ++n;
				//if (C(i)>cnt)
				//	cnt = C(i);
			}
			if (n>0) {
				calc_chisq(S, V, chi, df, idx);		
				chi *= 0.5;
				df *= 0.5;
				gammabase(&chi, &df, &p);
				p = 1.0 - p;
				chi *= 2.0;		
				if (p<=alpha) {
					pval2str(p, pstr);
					MYOUTPUT(os, setw(10) << loc[loc_idx].name)
					os.setf(ios::right, ios::adjustfield);
					MYOUTPUT(os, setw(4) << loc[loc_idx].alleleName[j])
					MYOUTPUT(os, setprecision(3) << setw(7) << frq[j] << setw(5) << (int)(I(j)+0.1))
					MYOUTPUT(os, setw(9) << (int)(2*df+0.01) << setw(12) << chi << setw(12) << pstr << "\n")
					os.setf(ios::left, ios::adjustfield);
				}
			}
		}
	}

		
	return 1;
	
}

// return of vector of X_mean for all the alleles of a locus
// use_cm: use conditional mean model
void FBAT::x_pop_mean(int loc_idx, int trt_id, ColumnVector &X, double gx[][3], int flag, bool use_cmm) {
	
	if (loc_idx<0 || loc_idx>=nLoci)
		throw Param_Error("Invalid locus_id in FBAT::x_pop_mean");
		
	int na = loc[loc_idx].nAlleles;
	
	if (X.Nrows()!=na)
		throw Param_Error("Inconsistent X vector dimension in FBAT::x_pop_mean");
	
	ColumnVector C(na), Y(na);
	ColumnVector EX(na);
	X = 0;
	C = 0;
	
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	PERSONLIST *sl;	
	Sufficient_Stat *ss;
	bool censored = flag&flag_censored;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		if (use_cmm) {
			for (sl=ped->node->members; sl!=NULL; sl=sl->next)
				sl->node->resetflag(flag_processed);
			for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
					if (flist->node->stat!=NULL) {
						delete flist->node->stat;
						flist->node->stat = NULL;
					}
					ss = new Sufficient_Stat(flist->node, 1, &loc_idx, loc[loc_idx].sexlinked());
					flist->node->stat = ss;
					if (ss==NULL || ss->mhap==NULL || ss->get_common_ogt()<2)
						continue;
					ss->analyze();
					if (ss->informative() && ss->fbat_ex(EX, gx,NULL)) {
						for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next) {
							sl->node->setflag(flag_processed);
							if (sl->node->good_gt(loc_idx, loc[loc_idx].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
								X += EX;
								C += 1;
							}
						}	
					}
				}
			}
		}
		for (sl=ped->node->members; sl!=NULL; sl=sl->next) {
			if (sl->node!=NULL && (!use_cmm || !sl->node->testflag(flag_processed)))
				if (sl->node->good_gt(loc_idx, loc[loc_idx].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
					gt2x(sl->node->gt[loc_idx], Y, sl->node->sex==2, loc[loc_idx].sexlinked(), gx);
					X += Y;
					C += 1;
				}
		}
	}
	
	for (int i=1; i<=na; i++)
		if (C(i)>0.1)
			X(i)/=C(i);

}

void FBAT::gt2x(Genotype &gt, ColumnVector &X, bool female, bool sexlink, double gx[][3]) {
	
	X = 0;
	
	if (gt.a[0]>0) {
		X(gt.a[0]) += ((female || !sexlink)?gx[gt.a[0]-1][gt.a[0]==gt.a[1]?2:1]:1);
			if (gt.a[1]>0 && gt.a[0]!=gt.a[1] && (female || !sexlink)) {
				X(gt.a[1]) += gx[gt.a[1]-1][1];
			}
	}

}
	

// return of vector of S and its covariance matrix V for all the alleles of a locus
// Statistic = S, Var = V(S) - u*Cov(S,X) + u^2Var(X) = V -2*u*B + u^2*A
// use_cm: use conditional mean model, fcnt: infomative count
void FBAT::x_pop_stat(int loc_idx, int trt_id, ColumnVector &S, Matrix &V, Matrix &A, Matrix &B, int *fcnt, double gx[][3], double offset, int flag, bool use_cmm) {
	
	if (loc_idx<0 || loc_idx>=nLoci)
		throw Param_Error("Invalid locus_id in FBAT::x_pop_mean");
		
	int na = loc[loc_idx].nAlleles;
	
	if (S.Nrows()!=na || V.Nrows()!=na || V.Ncols()!=na)
		throw Param_Error("Inconsistent matrix dimension in FBAT::x_pop_stat");
	
	int i;
	for (i=0; i<=na; i++)
		fcnt[i] = 0;
	S = 0;
	V = 0;
	A = 0;
	B = 0;
	
	ColumnVector PS(na);
	ColumnVector EX(na);
	ColumnVector Z(na), Y(na);	// E(X)
	x_pop_mean(loc_idx, trt_id, EX, gx, flag, use_cmm);
	
	bool censored = flag&flag_censored;
	bool good;
	
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	PERSONLIST *sl;	
	Sufficient_Stat *ss;
	double t;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		PS = 0;
		Z = 0;
		if (use_cmm) {
			for (sl=ped->node->members; sl!=NULL; sl=sl->next)
				sl->node->resetflag(flag_processed);
			for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
					ss = flist->node->stat;
					if (ss && ss->informative() && ss->fbat_ex(Y, gx,NULL)) {
						for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next) {
							sl->node->setflag(flag_processed);
							if (sl->node->good_gt(loc_idx, loc[loc_idx].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
								t = sl->node->pt[trt_id] - ( (!censored || sl->node->affection>0)? offset : 0 );	// for cencored data, t-u in affecteds, t in censored group
								PS += t*(Y-EX);
								Z += Y-EX;
							}
						}
					}
				}
			}
		}
		for (sl=ped->node->members; sl!=NULL; sl=sl->next) {
			if (sl->node && (!use_cmm || !sl->node->testflag(flag_processed)))
				if (sl->node->good_gt(loc_idx, loc[loc_idx].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
					t = sl->node->pt[trt_id] - ( (!censored || sl->node->affection>0)? offset : 0 );	// for cencored data, t-u in affecteds, t in censored group
					gt2x(sl->node->gt[loc_idx], Y, sl->node->sex==2, loc[loc_idx].sexlinked(), gx);
					PS += t*(Y-EX);
					Z += Y-EX;
				}
		}
		
		S += PS;
		V += PS*PS.t();
		B += Z*PS.t();
		A += Z*Z.t();
		for (good=false,i=1; i<=na; i++)
			if (fabs(PS(i))>very_small_number) {
				fcnt[i]++;
				good = true;
			}
		if (good)
			fcnt[0]++;
	}
	
	B = 0.5*(B+B.t());
}

void FBAT::pop_asso(int locidx, double gx[][3], int flag, bool use_cmm, int width, bool print_header) {
	
	if (locidx<0 || locidx>=nLoci)
		throw Param_Error("invalid locus id");
	
	int na = loc[locidx].nAlleles;
	double *frq = allelefrq(locidx);
	double *mu = new double[na+1];
	int *fcnt = new int[na+1]; 
	double offset = ((flag&flag_censored)?uncensored_trait_mean(trait_id[0]) : tr_offset);
	
	ColumnVector S(na);	
	Matrix V(na,na), A(na,na), B(na,na);
	
	x_pop_stat(locidx, trait_id[0], S, V, A, B, fcnt, gx, offset, flag, use_cmm);
	
	int i;
	if (flag&flag_optimize_offset) {	// update mu and fcnt
		mu[0] = B.Trace()/A.Trace();
		for (i=1; i<=na; i++)
			mu[i] = offset + B(i,i)/A(i,i);
		pop_count_infomative_fam(locidx, mu, fcnt, gx, flag, use_cmm);
		for (i=0; i<=na; i++)
			mu[i] -= offset;
	} else {
		for (i=1; i<=na; i++)
			mu[i] = 0;
		mu[0] = 0;
	}

	if (mode==by_allele || mode==by_both) {
		if (print_header) {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "Allele   afreq   fam#     S-E(S)      Var(S)       Z           P")
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, "      Offset")
			if (gen_model==model_estimate)
				MYOUTPUT(os, "   Dominance")
			MYOUTPUT(os, "\n")
			for (i=0; i<64+width+((flag&flag_optimize_offset)?12:0)+((gen_model==model_estimate)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		}
		double z, p, st, var;
		char pstr[16];
		for (i=1; i<=na; i++) {
			if (fcnt[i]>=min_size && frq[i]>=min_freq) {
				st = S(i);
				var = ((flag&flag_optimize_offset)? V(i,i)-2*B(i,i)*mu[i]+A(i,i)*mu[i]*mu[i] : V(i,i));
				z = st/sqrt(var);
				normbase(&z, &p);
				p = p>0.5?2*(1-p):2*p;
							
				if (alpha>0.99 ||  p<=alpha) {	//yes
					MYOUTPUT(os, setw(width) << loc[locidx].name << setw(7) << loc[locidx].alleleName[i])
					os.setf(ios::right, ios::adjustfield);
					MYOUTPUT(os, setprecision(3) << setw(7) << frq[i] << setprecision(1) << setw(7) << fcnt[i])
					MYOUTPUT(os, setprecision(3) << setw(11) << st << setw(12) << var)
					pval2str(p, pstr);
					MYOUTPUT(os, setw(8) << z << setw(12) << pstr)
					if (flag&flag_optimize_offset) 
						MYOUTPUT(os, setw(12) << offset+mu[i])
					if (gen_model==model_estimate)
						MYOUTPUT(os, setw(12) << gx[i-1][1])
					// debug: add regression Z
					z = pop_regression(locidx, i, trait_id[0], flag, offset, gx, use_cmm);
					MYOUTPUT(os, setprecision(3) << setw(8) << z)
					MYOUTPUT(os, "\n")
				}
				os.setf(ios::left, ios::adjustfield);
			}
		}
	}
	if (mode==by_marker || mode==by_both) {	// changed 5/13/03 to add by_both
		if (print_header>0) {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "Allele#   Fam#     DF       CHISQ           P")
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, "      Offset")
			MYOUTPUT(os, "\n")
			for (i=0; i<45+width+((flag&flag_optimize_offset)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		}
		// calculate optimized_offset and update S and V
		if (flag&flag_optimize_offset)
			V += mu[0]*mu[0]*A - 2*mu[0]*B;
		int n = 0;
		int *index = new int[na+1];
		for (i=1; i<=na; i++)
			index[i] = (fcnt[i]>=min_size && frq[i]>=min_freq)? ++n : 0;		
		if (n>1) {
			double df, chisq, p;
			bool err = false;
			try {
				p = calc_chisq(S, V, chisq, df, index);
			}
			catch (Exception error) {
				MYOUTPUT(os, "error in matrix inersion\n")
				err = true;
			}
			if (!err && (alpha>0.99 ||  p<=alpha)) {
				char pstr[16];
				MYOUTPUT(os, setw(width) << loc[locidx].name)
				os.setf(ios::right, ios::adjustfield);
				pval2str(p, pstr);
				MYOUTPUT(os, setw(7) << na << setw(7) << fcnt[0] << setw(7) << (int)(df+0.01) << setprecision(3) << setw(12) << chisq << setw(12) << pstr)
				if (flag&flag_optimize_offset)
					MYOUTPUT(os, setw(12) << offset+mu[0])
				MYOUTPUT(os, "\n")
				os.setf(ios::left, ios::adjustfield);
			}
		}
		if (index)
			delete[] index;
	}
	
	if (mu)
		delete[] mu;
		
	if (fcnt)
		delete[] fcnt;
	
}

void FBAT::pop_count_infomative_fam(int loc_idx, double mu[], int *cnt, double gx[][3], int flag, bool use_cmm) {

	if (loc_idx<0 || loc_idx>=nLoci)
		throw Param_Error("Invalid locus_id in FBAT::x_pop_mean");
		
	int na = loc[loc_idx].nAlleles;
	
	int i;
	for (i=0; i<=na; i++)
		cnt[i] = 0;
	
	ColumnVector PS(na);
	ColumnVector EX(na), Y(na), Z(na);
	x_pop_mean(loc_idx, trait_id[0], EX, gx, flag, use_cmm);

	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	PERSONLIST *sl;	
	Sufficient_Stat *ss;
	double t;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		PS = 0;
		Z = 0;
		if (use_cmm) {
			for (sl=ped->node->members; sl!=NULL; sl=sl->next)
				sl->node->resetflag(flag_processed);
			for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
				if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
					ss = flist->node->stat;
					if (ss && ss->mhap && ss->mhap->len()<maxcmh && ss->fbat_ex(Y, gx,NULL)) {
						for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next) {
							sl->node->setflag(flag_processed);
							if (sl->node->good_gt(loc_idx, loc[loc_idx].sexlinked()) && sl->node->hastrait(trait_id[0])) {
								t = sl->node->pt[trait_id[0]];
								PS += t*(Y-EX);
								Z += Y-EX;
							}
						}	
					}
				}
			}
		}
		for (sl=ped->node->members; sl!=NULL; sl=sl->next) {
			if (sl->node && (!use_cmm || !sl->node->testflag(flag_processed)))
				if (sl->node->good_gt(loc_idx, loc[loc_idx].sexlinked()) && sl->node->hastrait(trait_id[0])) {
					t = sl->node->pt[trait_id[0]];
					gt2x(sl->node->gt[loc_idx], Y, sl->node->sex==2, loc[loc_idx].sexlinked(), gx);
					PS += t*(Y-EX);
					Z += Y-EX;
				}
		}
		for (i=1; i<=na; i++)
			if (fabs(PS(i)-mu[i]*Z(i))>very_small_number)
				cnt[i]++;
		for (i=1; i<=na && fabs(PS(i)-mu[0]*Z(i))<very_small_number; i++)
			;
		if (i<=na)
			cnt[0]++;
	}
	
}

// population_base association test
void FBAT::pbat(char *s) {
	
	bool censored=false;
	bool empirical=false;
	bool optimize_offset = false;
	bool multimarker = false;
	bool use_cmm = false;	// use conditional mean model
	bool joint = false;
	
	char word[1024];
	int i, k;
	
	if (mode==by_both)
		throw Param_Error("by_both mode is only available for haplotype test");
		
	int wcnt = WordCount(s);
	
	int flag = 0;
	for (i=2; i<=wcnt && GetWord(s, i, word) && *word=='-'; i++) {
		if (word[2]==0)
			switch (word[1]) {
				case 'o': 
					optimize_offset = true;
					flag |= flag_optimize_offset;
					break;
				case 'e':
					empirical = true;
					flag |= flag_empirical_var;
					break;
				case 'c': 
					flag |= flag_censored;
					censored = true; 
					if (trait_id[0]==0)
						throw Param_Error("affection can't be the censored trait");
					break;
				case 'i':
					use_cmm = true;
					break;
				case 'j':
					joint = true;
					break;
				default: throw Param_Error("unknown option");
			}
		else
			throw Param_Error("unknown option");
	}
	
	if (n_sel_trt>1)
		throw Param_Error("Multivariate test is not implemented for pbat");
		
	if (optimize_offset && censored)
		throw Param_Error("Optimize offset option is not available for censored trait");
	
	if (empirical && optimize_offset)
		throw Param_Error("Offset nuisance optimization is not available for empirical test");
	
	if (gen_model == model_genotype)
		throw Param_Error("The requested genotype test has not been implemented");
	
	if (gen_model == model_estimate && (multimarker || n_sel_trt>1))
		throw Param_Error("The model estimation is only available for single-marker univariate FBAT test");
	
	print_settings();
	
	int nsel = i>=wcnt? nLoci : wcnt-i+1;
	int *idx = new int[nsel];
	if (i>wcnt)
		for (nsel=0; nsel<nLoci; nsel++)
			idx[nsel] = nsel;
	else {
		for (nsel=0; i<=wcnt; i++) {
			GetWord(s, i, word);
			if ((k=find_locus(word))>=0)
				idx[nsel++] = k;
			else {
				if (idx!=NULL)
					delete[] idx;
				char msg[256];
				sprintf(msg, "Marker %s not found", word);
				throw Param_Error(msg);
			}
		}
	}
	
	double gx[kMaxAlleles][3];
	if (gen_model==model_additive || gen_model==model_dominant || gen_model==model_recessive)
		for (i=0; i<kMaxAlleles; i++)
			gcode(&gx[i][0], gen_model);

	int mrkwidth = 0;
	for (i=0; i<nsel; i++)
		if ((k=strlen(loc[i].name)) > mrkwidth)
			mrkwidth = k;
	
	for (i=0; i<nsel; i++)
		if (joint)
			joint_asso(idx[i], gx, flag, mrkwidth+4, i==0);
		else
			pop_asso(idx[i], gx, flag, use_cmm, mrkwidth+4, i==0);
	
	MYOUTPUT(os, "\n")
	
	if (idx)
		delete[] idx;
}

void FBAT::x_fbat_stat(int loc_idx, int trt_id, ColumnVector &S, ColumnVector &X, Matrix &V, Matrix &A, Matrix &B, int *fcnt, double gx[][3], double offset, int flag, IRES* ires) {

	if (loc_idx<0 || loc_idx>=nLoci)
		throw Param_Error("invalid locus id");
	
	int na = loc[loc_idx].nAlleles;
		
	ColumnVector FS(na);
	ColumnVector PS(na);	// sum of S in whole pedigree, used in empirical var test
	ColumnVector FES(na);
	Matrix FV(na,na);
	Matrix FA(na,na);
	Matrix FB(na,na);
	ColumnVector FX(na);
	X = 0;
	S = 0;
	A = 0;
	B = 0;
	V = 0;
	int i;
	int tot_cnt = 0;
	for (i=0; i<=na; i++)
		fcnt[i] = 0;
	bool good;
	Sufficient_Stat *ss;
	bool noinfo;
	int  famcnt=0;

	for (PEDIGREELIST *ped=peds; ped!=NULL; ped=ped->next) {
		PS = 0;
		for (NUCFAMILYLIST *flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			noinfo = true;
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				//famcnt++;
				//MYOUTPUT(os, "\n[" << ped << "," << flist << "," << flist->node << "]\n")
				if (flist->node->stat!=NULL) {
					delete flist->node->stat;
					flist->node->stat = NULL;
				}
				ss = new Sufficient_Stat(flist->node, 1, &loc_idx, loc[loc_idx].sexlinked());
				flist->node->stat = ss;
				if (ss==NULL || ss->mhap==NULL || ss->mhap->len()>maxcmh)
					goto capture; //continue;

				tot_cnt++;
				ss->get_common_ogt();
				ss->analyze();
				if (ss->informative() && ss->fbat_stat3(trait_id[0], offset, FES, FS, FV, FA, FB, FX, gx, flag)) {
					noinfo = false;
					FS -= FES;
					PS += FS;
					if (!(flag&flag_empirical_var)) {	// not empirical variance
						A += FA;
						B += FB;
						X += FX;
						S += FS;
						V += FV;
						if (!(flag&flag_optimize_offset)) {	
							for (good=false,i=1; i<=na; i++)
								if (FV(i,i)>very_small_number) {
									fcnt[i]++;
									good = true;
								}
							if (good)	// informative for ar least one allele
								fcnt[0]++;
						}
					}
					
				}
			}
		capture:
		  if (ires != NULL){    
			famcnt++;
			if (noinfo)
				ires->additem(loc_idx, ped, flist, 0., 0., 0.);
			else {	
				double s = 0.;
				double es = 0.;
				double v = 0.;

				for (int i = 1; i <= na; i++){
					double *frq = allelefrq(loc_idx);
					if (frq[i] < 0.5) {					
					  s += FS(i) + FES(i);
					  es += FES(i);
					  v += FV(i,i);
					}
				}				
				ires->additem(loc_idx, ped, flist, s, es, v);
			}
		  }
		}
		if (flag&flag_empirical_var) {
			FV = PS*PS.t();
			for (good=false,i=1; i<=na; i++)
				if (FV(i,i)>very_small_number) {
					fcnt[i]++;
					good = true;
				}
			if (good) {	// informative for ar least one allele
				fcnt[0]++;
				S += PS;
				V += FV;
			}
		}
	}

	if (ires != NULL) {
		ires->setfam_total(famcnt);
		// MYOUTPUT(os, "[famcnt= " << famcnt << "\n")
	}

}

// joint analysis (pop + fbat)
void FBAT::joint_asso(int locidx, double gx[][3], int flag, int width, bool print_header) {
	
	if (locidx<0 || locidx>=nLoci)
		throw Param_Error("invalid locus id");
	
	int na = loc[locidx].nAlleles;
	double *frq = allelefrq(locidx);
	double *mu = new double[na+1];
	int *pop_fcnt = new int[na+1];
	int *fam_fcnt = new int[na+1];
	
	double offset = ((flag&flag_censored)?uncensored_trait_mean(trait_id[0]) : tr_offset);
	
	ColumnVector PS(na), FS(na), FX(na);	
	Matrix PV(na,na), PA(na,na), PB(na,na), FV(na,na), FA(na,na), FB(na,na);
	
	x_pop_stat(locidx, trait_id[0], PS, PV, PA, PB, pop_fcnt, gx, offset, flag, true);
	x_fbat_stat(locidx, trait_id[0], FS, FX, FV, FA, FB, fam_fcnt, gx, offset, flag);
	
	int i;
	if (flag&flag_optimize_offset) {	// update mu and fcnt
		mu[0] = (PB.Trace()+FB.Trace())/(PA.Trace()+FA.Trace());
		for (i=1; i<=na; i++)
			mu[i] = offset + (PB(i,i)+FB(i,i))/(PA(i,i)+FA(i,i));
		pop_count_infomative_fam(locidx, mu, pop_fcnt, gx, flag, true);
		count_infomative_fam(mu, na, fam_fcnt, gx, flag);
		for (i=0; i<=na; i++)
			mu[i] -= offset;
	} else {
		for (i=1; i<=na; i++)
			mu[i] = 0;
		mu[0] = 0;
	}

	if (mode==by_allele || mode==by_both) {
		if (print_header) {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "Allele   afreq  fam#p S-E(S)_pop  Var(S)_pop   Z_pop       P_pop  fam#f S-E(S)_fam  Var(S)_fam   Z_fam       P_fam Z_joint     P_joint")
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, "      Offset")
			if (gen_model==model_estimate)
				MYOUTPUT(os, "   Dominance")
			MYOUTPUT(os, "\n")
			for (i=0; i<64+70+width+((flag&flag_optimize_offset)?12:0)+((gen_model==model_estimate)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		}
		double z_pop, z_fam, z_joint, p_pop, p_fam, p_joint, st_pop, st_fam, st_joint, var_pop, var_fam, var_joint;
		char pstr[16];
		for (i=1; i<=na; i++) {
			if ((pop_fcnt[i]+fam_fcnt[i])>=min_size && frq[i]>=min_freq) {
				st_pop = PS(i);
				var_pop = ((flag&flag_optimize_offset)? PV(i,i)-2*PB(i,i)*mu[i]+PA(i,i)*mu[i]*mu[i] : PV(i,i));
				z_pop = st_pop/sqrt(var_pop);
				normbase(&z_pop, &p_pop);
				p_pop = p_pop>0.5?2*(1-p_pop):2*p_pop;
				
				st_fam = FS(i) - mu[i]*FX(i);
				var_fam = ((flag&flag_optimize_offset)? FV(i,i)-2*FB(i,i)*mu[i]+FA(i,i)*mu[i]*mu[i] : FV(i,i));
				z_fam = st_fam/sqrt(var_fam);
				normbase(&z_fam, &p_fam);
				p_fam = p_fam>0.5?2*(1-p_fam):2*p_fam;
				
				st_joint = st_pop + st_fam;
				var_joint = var_pop + var_fam;
				z_joint = st_joint/sqrt(var_joint);
				normbase(&z_joint, &p_joint);
				p_joint = p_joint>0.5?2*(1-p_joint):2*p_joint;
							
				if (alpha>0.99 ||  p_joint<=alpha) {	//yes
					MYOUTPUT(os, setw(width) << loc[locidx].name << setw(7) << loc[locidx].alleleName[i])
					os.setf(ios::right, ios::adjustfield);
					MYOUTPUT(os, setprecision(3) << setw(7) << frq[i])
					MYOUTPUT(os, setw(7) << pop_fcnt[i] << setprecision(3) << setw(11) << st_pop << setw(12) << var_pop)
					pval2str(p_pop, pstr);
					MYOUTPUT(os, setw(8) << z_pop << setw(12) << pstr)
					MYOUTPUT(os, setw(7) << pop_fcnt[i] << setprecision(3) << setw(11) << st_fam << setw(12) << var_fam)
					pval2str(p_fam, pstr);
					MYOUTPUT(os, setw(8) << z_fam << setw(12) << pstr)
					pval2str(p_joint, pstr);
					MYOUTPUT(os, setw(8) << z_joint << setw(12) << pstr)
					if (flag&flag_optimize_offset) 
						MYOUTPUT(os, setw(12) << offset+mu[i])
					if (gen_model==model_estimate)
						MYOUTPUT(os, setw(12) << gx[i-1][1])
					MYOUTPUT(os, "\n")
				}
				os.setf(ios::left, ios::adjustfield);
			}
		}
	}
	if (mode==by_marker || mode==by_both) {	// changed 5/13/03 to add by_both
		if (print_header>0) {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "Allele#  Fam#p DF_pop   CHISQ_pop       P_pop  Fam#f DF_fam   CHISQ_fam       P_fam DF_joi   CHISQ_joi     P_joint")
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, "      Offset")
			MYOUTPUT(os, "\n")
			for (i=0; i<45+69+width+((flag&flag_optimize_offset)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		}
		int *index = new int[na+1];
		int n;
		double df_pop=0, chisq_pop, p_pop, df_fam=0, chisq_fam, p_fam, df_joint=0, chisq_joint, p_joint;
		// calculate optimized_offset and update S and V
		if (flag&flag_optimize_offset) {
			FS -= mu[0]*FX;
			PV += mu[0]*mu[0]*PA - 2*mu[0]*PB;
			FV += mu[0]*mu[0]*FA - 2*mu[0]*FB;
		}
		for (n=0,i=1; i<=na; i++)
			index[i] = (pop_fcnt[i]>=min_size && frq[i]>=min_freq)? ++n : 0;		
		if (n>1) {
			try {
				p_pop = calc_chisq(PS, PV, chisq_pop, df_pop, index);
			}
			catch (Exception error) {
				throw("error in matrix inersion");
			}
		}
		for (n=0,i=1; i<=na; i++)
			index[i] = (fam_fcnt[i]>=min_size && frq[i]>=min_freq)? ++n : 0;		
		if (n>1) {
			try {
				p_fam = calc_chisq(FS, FV, chisq_fam, df_fam, index);
			}
			catch (Exception error) {
				throw("error in matrix inersion");
			}
		}
		FS += PS;
		FV += PV;
		for (n=0,i=1; i<=na; i++)
			index[i] = (fam_fcnt[i]>=min_size && frq[i]>=min_freq)? ++n : 0;		
		if (n>1) {
			try {
				p_joint = calc_chisq(FS, FV, chisq_joint, df_joint, index);
			}
			catch (Exception error) {
				throw("error in matrix inersion");
			}
		}
		if (df_joint>0 && (alpha>0.99 ||  p_joint<=alpha)) {
			char pstr[16];
			MYOUTPUT(os, setw(width) << loc[locidx].name)
			os.setf(ios::right, ios::adjustfield);
			MYOUTPUT(os, setw(7) << na << setw(7) << pop_fcnt[0])
			if (df_pop>0) {
				pval2str(p_pop, pstr);
				MYOUTPUT(os, setw(7) << (int)(df_pop+0.01) << setprecision(3) << setw(12) << chisq_pop << setw(12) << pstr)
			} else {
				MYOUTPUT(os, setw(34) << "")
			}
			MYOUTPUT(os, setw(7) << fam_fcnt[0])
			if (df_fam>0) {
				pval2str(p_fam, pstr);
				MYOUTPUT(os, setw(7) << (int)(df_fam+0.01) << setprecision(3) << setw(12) << chisq_fam << setw(12) << pstr)
			} else {
				MYOUTPUT(os, setw(34) << "")
			}
			
			//MYOUTPUT(os, setw(7) << (pop_fcnt[0]+fam_fcnt[0]))
			pval2str(p_joint, pstr);
			MYOUTPUT(os, setw(7) << (int)(df_joint+0.01) << setprecision(3) << setw(12) << chisq_joint << setw(12) << pstr)
			if (flag&flag_optimize_offset)
				MYOUTPUT(os, setw(12) << offset+mu[0])
			MYOUTPUT(os, "\n")
			os.setf(ios::left, ios::adjustfield);
		}
		if (index)
			delete[] index;
	}
	
	if (mu)
		delete[] mu;
		
	if (pop_fcnt)
		delete[] pop_fcnt;
	
	if (fam_fcnt)
		delete[] fam_fcnt;

}

double FBAT::fbat_regression(int loc_index, int a_idx, int trt_id, int flag, double offset, double gscore[][3], double *b, double *v) {


	double y = 0.0;
	double x = 0.0;
	double xx = 0.0;
	double xy = 0.0;
	double nx = 0.0;
	double yy = 0.0;
	
	int n;
		
	PEDIGREELIST *ped;
	NUCFAMILYLIST *flist;
	PERSONLIST *sl;
	
	Sufficient_Stat *ss;

	bool censored = flag&flag_censored;
	
	double t, xm, xs;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0) {
				if (flist->node->stat!=NULL) {
					delete flist->node->stat;
					flist->node->stat = NULL;
				}
				ss = new Sufficient_Stat(flist->node, 1, &loc_index, loc[loc_index].sexlinked());
				flist->node->stat = ss;
				if (ss==NULL || ss->mhap==NULL || ss->get_common_ogt()<2)
					continue;
				ss->analyze();
				if (!ss->informative())
					continue;					
				xm = ss->ex(a_idx, gscore);
				for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next)
					if (sl->node && sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
						t = sl->node->trait(trt_id) - ((censored && sl->node->affection!=2)? 0 : offset);
						n = sl->node->gt[loc_index].allele_cnt(a_idx,(loc[loc_index].sexlinked() && sl->node->sex==1));
						xs = (loc[loc_index].sexlinked() && sl->node->sex==1? (n>0?1:0) : gscore[a_idx-1][n]) - xm;
						x += xs;
						y += t;
						xy += xs*t;
						xx += xs*xs;
						yy += t*t;
						nx += 1;
					}
			}
		}
	}
	
	if (nx*yy-y*y<very_small_number)
		throw ERROR("no variation of phenotype");
	
	double beta = (nx>1 && fabs(nx*xx-x*x)>very_small_number)? (xy-x*y/nx)/(xx-x*x/nx) : 0;
	double alfa = (y - beta*x)/nx;

	//double var = ((yy-y*y/nx)-beta*(xy-x*y/nx))/((nx-2)*(xx-x*x/nx));
	// calculate robust sandwitch variance
	double xmean = x/nx;
	double exex = 0;
	double fex;
	for (ped=peds; ped!=NULL; ped=ped->next) {
		fex = 0;		
		for (flist=ped->node->nucfams; flist!=NULL; flist=flist->next) {
			if (flist->node!=NULL && flist->node->testflag(flag_uninformative)==0 && (ss=flist->node->stat)!=NULL && ss->informative()) {
				xm = ss->ex(a_idx, gscore);
				for (sl=flist->node->sibs(); sl!=NULL; sl=sl->next)
					if (sl->node && sl->node->good_gt(loc_index, loc[loc_index].sexlinked()) && sl->node->hastrait(trt_id) && (!censored || sl->node->affection>0)) {
						t = sl->node->trait(trt_id) - ((censored && sl->node->affection!=2)? 0 : offset);
						n = sl->node->gt[loc_index].allele_cnt(a_idx,(loc[loc_index].sexlinked() && sl->node->sex==1));
						xs = (loc[loc_index].sexlinked() && sl->node->sex==1? (n>0?1:0) : gscore[a_idx-1][n]) - xm;
						fex += (t - alfa - beta*xs)*(xs-xmean);
						sl->node->setflag(flag_processed);
					}
			}
		}
		
		exex += fex*fex;
	}		
			
	t = xx - nx*xmean*xmean;
	double var = exex/(t*t);
	double z = beta/sqrt(var);
	
	if (b)
		*b = beta;
	if (v)
		*v = var;
		
	return z;
}

// regression based fbat
void FBAT::fdt_reg(int locidx, double gx[][3], int flag, int width, bool print_header) {
	
	if (locidx<0 || locidx>=nLoci)
		throw Param_Error("invalid locus id");
	
	int i;
	int na = loc[locidx].nAlleles;
	double *frq = allelefrq(locidx);
	double offset = ((flag&flag_censored)?uncensored_trait_mean(trait_id[0]) : tr_offset);
		
	if (mode!=by_allele)
		throw ("the specified test mode is not available for regression based fbat");
	
	if (mode==by_allele || mode==by_both) {
		if (print_header) {
			MYOUTPUT(os, "\n" << setw(width) << "Marker" << "Allele   afreq        Beta   Var(Beta)       Z           P\n")
			for (i=0; i<64+width+((flag&flag_optimize_offset)?12:0)+((gen_model==model_estimate)?12:0); i++)
				MYOUTPUT(os, "-")
			MYOUTPUT(os, "\n")
		}
		double z, p, beta, var, z1, b1, v1, b2, v2, z2;
		char pstr[16];
		
		for (i=1; i<=na; i++) {
			z = fbat_regression(locidx, i, trait_id[0], flag, offset, gx, &beta, &var);
			normbase(&z, &p);
			p = p>0.5?2*(1-p):2*p;
						
			if (alpha>0.99 ||  p<=alpha) {	//yes
				os.setf(ios::left, ios::adjustfield);
				MYOUTPUT(os, setw(width) << loc[locidx].name << setw(7) << loc[locidx].alleleName[i])
				os.setf(ios::right, ios::adjustfield);
				MYOUTPUT(os, setprecision(3) << setw(7) << frq[i] << setprecision(3) << setw(12) << beta << setw(12) << var)
				pval2str(p, pstr);
				MYOUTPUT(os, setw(8) << z << setw(12) << pstr)
				
				z1 = pop_regression(locidx, i, trait_id[0], flag, offset, gx, true, &b1, &v1);
				MYOUTPUT(os, setprecision(3) << setw(12) << b1 << setw(12) << z1)
				b2 = (beta/var+b1/v1)/(1/var+1/v1);
				v2 = 1/(var*(1/var+1/v1)*(1/var+1/v1)) + 1/(v1*(1/var+1/v1)*(1/var+1/v1));
				z2 = b2/sqrt(v2);
				MYOUTPUT(os, setprecision(3) << setw(12) << b2 << setw(12) << z2)
				MYOUTPUT(os, "\n")
			}
		}
	}

}

void FBAT::set_affection_trait(char *s) {
	
	int i = WordCount(s);
	
	if (i!=4)
		throw Param_Error("wrong number of arguments");
	
	float xx[3];
	if (sscanf(GetWord(s,2,NULL), "%f%f%f", xx+2, xx+1, xx)!=3)
		throw Param_Error("number expected in arguments");
	
	MYOUTPUT(os, "conver affection status to trait " << xx[2] << "(aff), " << xx[1] << "(unaff), " << xx[0] << "(unknown)\n")
	
	for (i=0; i<3; i++)
		affection_trait[i] = xx[i];
	
	update_affection_trait();
	
}

// calculating pairwise LD and D' using the first allele of each marker
void FBAT::ld_matrix_print(int nmrk, int idx[], HaploidTable *htab, bool print_d) {
	
	if (htab==NULL)
		throw ERROR("no haplotype table available");
	
	if (nmrk<2 || idx==NULL)
		throw ERROR("haplotype markers not specified");
	
	int i, j;
	double **f = new double*[nmrk];	// f[i][i] is allele freq of marker[i], f[i][j] is haplotype 11 of marker i and j
	for (i=0; i<nmrk; i++) {
		f[i] = new double[nmrk];
		for (j=0; j<nmrk; j++)
			f[i][j] = 0.0;
	}
	
	for (HaploidList *hpl=htab->hpl; hpl!=NULL; hpl=hpl->next) {
		for (i=0; i<nmrk; i++)
			if (hpl->node->h[i]==1) {
				f[i][i] += hpl->node->p;
				for (j=i+1; j<nmrk; j++)
					if (hpl->node->h[j]==1)
						f[i][j] += hpl->node->p;
			}
	}
	
	if (print_d)
		MYOUTPUT(os, "\nPairwise LD(D') matrix using first allele of each marker\n")
	else
		MYOUTPUT(os, "\nPairwise LD(r^2) matrix using first allele of each marker\n")
	MYOUTPUT(os, setw(13) << " ")
	for (i=0; i<nmrk-1; i++)
		MYOUTPUT(os, setw(13) << loc[idx[i]].name)
	MYOUTPUT(os, "\n");
	double d, dp, dmax_neg, dmax_pos, rsq;
	for (i=1; i<nmrk; i++) {
		MYOUTPUT(os, setw(12) << loc[idx[i]].name)
		for (j=0; j<i; j++) {
			if (f[i][i]==0.0 || f[j][j]==0.0)
				MYOUTPUT(os, "     ----    ")
			else {
				d = f[j][i]-f[i][i]*f[j][j];
				/* corrected 10/16/03 for dmax
				dmax = fabs(f[i][i]-(f[i][i]<=1-f[j][j]?f[i][i]:1-f[j][j])-f[i][i]*f[j][j]);
				dp = (f[i][i]>f[j][j]?f[j][j]:f[i][i]) - f[i][i]*f[j][j];
				if (dp>dmax)
					dmax = dp;
				dp = fabs(d/dmax);
				*/
				// corrected 2/18/2004 to calculate D' by the sign of D
				dmax_neg = fabs(f[i][i]-(f[i][i]<=1-f[j][j]?f[i][i]:1-f[j][j])-f[i][i]*f[j][j]);
				dmax_pos = (f[i][i]>f[j][j]?f[j][j]:f[i][i]) - f[i][i]*f[j][j];
				dp = fabs(d/(d>0?dmax_pos:dmax_neg));
				rsq = d*d/(f[i][i]*(1-f[i][i])*f[j][j]*(1-f[j][j]));
				if (d>=0)
					MYOUTPUT(os, "  ")
				else
					MYOUTPUT(os, " ")
				MYOUTPUT(os, setprecision(3) << d << "(" << setprecision(2) << (print_d?dp:rsq) << ")")
			}
		}
		MYOUTPUT(os, "\n")
	}
	
	MYOUTPUT(os, "\n")
	
	for (i=0; i<nmrk; i++)
		if (f[i]!=NULL)
			delete[] f[i];
	
	if (f!=NULL)
		delete[] f;
		
}

void FBAT::haplotypefrequency(char *s) {

	int i, j, k;

	i = WordCount(s);
	
	char name[256];
	bool pair_ld = false;
	bool pair_rsq = false;
	
	for (j=2; j<=i && GetWord(s,j,name) && *name=='-'; j++)
		switch (name[1]) {
			case 'r':
				pair_rsq = true;
				break;
			case 'd':
				pair_ld = true;
				break;
			default:
				throw Param_Error("unknown option");
		}
	
	if (pair_ld && pair_rsq)
		throw Param_Error("can't output pairwise d' & r^2 together");
		
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
	
	Haplotype_List *emhl = hapfreq_ss(nsel, idx);
	
	if (emhl==NULL) 
		return;
	
	MYOUTPUT(os, "\nhaplotype analysis for the following markers:\n")
	for (i=0; i<nsel; i++)
		MYOUTPUT(os, loc[idx[i]].name << " ")

	HaploidTable *htab = new HaploidTable(loc, nsel, idx, emhl);
	
	MYOUTPUT(os, "\n\nhaplotypes and EM estimates of frequency:\n")
	i = 1;
	for (HaploidList *hpl2=htab->hpl; hpl2!=NULL; hpl2=hpl2->next,i++) {
		MYOUTPUT(os, "a" << setw(4) << i << setw(6) << ":")
		for (j=0; j<nsel; j++)
			MYOUTPUT(os, setw(4) << loc[idx[j]].alleleName[hpl2->node->h[j]]);
		MYOUTPUT(os, setw(12) << hpl2->node->p << "\n")
	}

	if (nsel>1 && (pair_ld || pair_rsq))
		ld_matrix_print(nsel, idx, htab, pair_ld);
		
	if (emhl)
		delete emhl;
	
	if (htab)
		delete htab;
		
	delete[] idx;

}


