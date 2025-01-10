using namespace std;
#include "mlr_fbat.h"
#include <math.h>

// independent function, does not call other function, okay
double log10_w_sum(vector<double> &vec, vector<double> &wts){

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
double MLR::get_z_score(int index)
{
	if(index>=1 && index<=n_var)
	{
		return (*Z)(index);
	}
	else{
		return -1.0;
	}
}
void MLR::init(const ColumnVector *Z_, const Matrix *R_)
{
	int i,j;
    n_var = R_->Nrows();

	Z = new ColumnVector(n_var) ;
	R = new Matrix(n_var, n_var); 
	
	for(i=0;i<n_var;i++)
	{
		(*Z)(i+1)=(*Z_)(i+1); 
		for(j=0;j<n_var;j++)
		{
			(*R)(i+1,j+1)=(*R_)(i+1,j+1); 
		}
	}
	
}


void MLR::copy(const MLR & mlr)
{
	int i,j;
    
    n_var = mlr.n_var;
	if(Z!=NULL)
        delete Z;
    if(R!=NULL)
        delete R;
	
    if(mlr.R != NULL)
	{
		R = new Matrix(n_var, n_var);
		
		for(i=0;i<n_var;i++)
		{
			for(j=0;j<n_var;j++)
			{
				(*R)(i+1,j+1)=(*mlr.R)(i+1,j+1); 
			}
		}
    }
    if(mlr.Z != NULL)
	{
        Z = new ColumnVector(n_var);
	
		for(i=0;i<n_var;i++)
		{
			(*Z)(i+1)=(*mlr.Z)(i+1); 
		}
    }
 
}

void MLR::copy(const MLR & mlr, vector<int> & indicator)
{
    
    int nv = mlr.n_var;
	
    int ep = 0;
	int i,j;
    int count = 0;
    std::map<int,int> imap;


    for(i=0;i<indicator.size();i++)
	{
        if(indicator[i]==1)
		{
            ep++;
            imap[i] = count++;
        }
    }
    
    n_var = ep; 
	if(Z!=NULL)
        delete Z;
    if(R!=NULL)
        delete R;
	R = new Matrix(ep, ep); 
	Z = new ColumnVector(ep); 
   
    double val;

    for(i=0;i<nv;i++)
	{
        if(indicator[i] == 0)
            continue;
		(*Z)(imap[i]+1) = (*mlr.Z)(i+1); 
        for(j=0;j<nv;j++)
		{
            if(indicator[j] == 1)
			{
                val = (*mlr.R)(i+1,j+1); 
				(*R)(imap[i]+1, imap[j]+1)=val; 
            }
        }
    }
    phi2_vec = mlr.phi2_vec;
	
}


MLR::~MLR()
{  
    if(Z!=NULL)
        delete Z;
    if(R!=NULL)
        delete R;
}

void MLR::set_effect_vec(vector<double> &phi2_vec_)
{ 
    phi2_vec = phi2_vec_; //!/ why by reference?
}

double MLR::compute_log10_ABF()
{
   vector<int> indicator(n_var,1);
   return(compute_log10_ABF(indicator));
}

double MLR::compute_log10_ABF(vector<int> & indicator)
{

    // construct the sub-matrices

    vector<double> rstv;
	int i,j;
    int ep = 0;
    int count = 0;
    std::map<int,int> imap;


    for(i=0;i<indicator.size();i++)
	{
        if(indicator[i]==1)
		{
            ep++;
            imap[i] = count++;
        }
    }

	Matrix Rm(ep, ep);
	Matrix Zm(ep,1);
	Matrix V(ep, ep);
	DiagonalMatrix S(ep);
	Matrix work (ep,ep);
	Rm=0; Zm=0; V=0; S=0; work=0;
    double val;

    for(i=0;i<n_var;i++)
	{
        if(indicator[i] == 0)
            continue;
		Zm(imap[i]+1,1)=(*Z)(i+1);
        
        for(j=0;j<n_var;j++)
		{
            if(indicator[j] == 1)
			{
                val = (*R)(i+1,j+1);
				Rm(imap[i]+1, imap[j]+1)=val;
            }
        }
    }
    
    //SVD  gsl_linalg_SV_decomp (Rm, V, S,work); //!/
	//int gsl_linalg_SV_decomp(gsl_matrix *A, gsl_matrix *V, gsl_vector *S, gsl_vector *work)
	//void SVD(const Matrix& A, DiagonalMatrix& Q, Matrix& U, Matrix& V, bool withU, bool withV)
	SVD(Rm, S, work, V, true, true); //!/
	
	double v;
    for(i=0;i<phi2_vec.size();i++)
	{
        double kappa = phi2_vec[i];
		
        Matrix t1(ep, ep);
		Matrix t2(ep, ep);
		Matrix M_inv(ep, ep);
		Matrix t3(1,ep);
		Matrix t4(1,1);
		t1=0; t2=0; t3=0; t4=0; M_inv=0;
		
        double log_det = 0;
        for(j=0;j<ep;j++)
		{
            v = S(j+1);
			t1(j+1,j+1) = 1.0/(v+1/kappa);
            log_det += log(1+kappa*v);
        }

        
        //gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,t1,0,t2);
		t2 = V * t1; // ep x ep  xx ep x ep = ep x ep

        //gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,t2,V,0,M_inv);
		M_inv = t2 * V.t(); // ep x ep xx ep x ep = ep x ep
		
        //gsl_blas_dgemm(CblasTrans, CblasNoTrans,1, Zm,M_inv,0,t3);
		t3 = Zm.t() * M_inv; // 1 x ep xx ep x ep = 1 x ep
	
        //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,t3,Zm,0,t4);
		t4 = t3 * Zm; // 1 x ep xx ep x 1 = 1 x 1

        double log_BF = -0.5*log_det + 0.5*t4(1,1);
        
        rstv.push_back(log_BF/log(10));
        
    }
	
	
    vector<double> wv(phi2_vec.size(), 1.0/phi2_vec.size());
    double rst =  log10_w_sum(rstv,wv);
    return rst;
}




