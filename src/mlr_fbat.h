#pragma once
#include <vector>
#include <map>
#include "newmat.h"
#include "newmatap.h"
#include <math.h>

class MLR
{
    private:
        int n_var; // number of variants

        // summary data
		ColumnVector *Z;
		Matrix *R;
        
		vector<double> phi2_vec;

    public:
	
        MLR()
		{
            Z = NULL;
			R = NULL;
        }

        ~MLR(); //destructor

        // init
        void init(const ColumnVector *Z_, const Matrix *R_);  
        void copy(const MLR &mlr); // copy MLR
        void copy(const MLR & mlr, vector<int> & indicator); //copy MLR for indicator==1 variants

        void set_effect_vec(vector<double> &phi2); // set effect size vector

		// computes log10_ABF
		double compute_log10_ABF();
        double compute_log10_ABF(vector<int> &indicator);
		double get_z_score(int index);
       
};

