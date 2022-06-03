// ires.cpp
//   implementation of IRES class
//   wai-ki yip (3/4/11)
//
#include <cmath>
#include "newmat.h"
#include "normbase.h"
#include "mystream.h"
#include "pedigree.h"
#include "ires.h"

IRES::IRES(FBAT *fbat, bool empirical, PEDIGREELIST *peds, int pedcount)
{
	this->fam_total = -1;
	this->fb = fbat;
	this->count = 0;
	this->wc = 0;
	this->empirical = empirical;
	this->peds = peds;
	this->ped_count = pedcount;
}


IRES::~IRES()
{
	this->ir.erase(this->ir.begin(), this->ir.end());
	this->weight.erase(this->weight.begin(), this->weight.end());
}

void IRES::print()
{
	if (this->fam_total <= 0){
		MYOUTPUT(fb->os, "Cannot generate intermediate results - Family Count not set\n");
	}
	else {
		int nsel;

		MYOUTPUT(fb->os, "\nContribution of FBAT statistics by family per pedigree" << '\n');
		nsel = this->getsize()/this->fam_total;

		/* generate V statistics */
		MYOUTPUT (fb->os, "V\n");
		MYOUTPUT (fb->os, "pid    fid    ");
		for (int i = 0; i < this->getsize(); i+=this->fam_total){
			MYOUTPUT(fb->os, fb->loc[ir[i].sel].name << " ");
		}
		MYOUTPUT (fb->os, '\n');

		/* enumerate through V statistics */
		for (int j = 0; j < this->fam_total; j++){
			for (int i = 0; i < nsel; i++){
				int k = i*this->fam_total + j;
				if (i == 0){
					MYOUTPUT(fb->os, ir[k].ped->node->name << " ");
	        			MYOUTPUT(fb->os, "["<<(ir[k].flist->node->fa() ? ir[k].flist->node->fa()->id:"-")<<"x");
	        			MYOUTPUT(fb->os, (ir[k].flist->node->mo() ? ir[k].flist->node->mo()->id:"-") << "] ");
				}
				MYOUTPUT(fb->os, ir[k].v << " ");
			}
			MYOUTPUT(fb->os, '\n');
		}	

		/* generate S statistics */
		MYOUTPUT (fb->os, "\nS\n");
		MYOUTPUT (fb->os, "pid    fid    ");
		for (int i = 0; i < this->getsize(); i+=this->fam_total){
			MYOUTPUT(fb->os, fb->loc[ir[i].sel].name << " ");
		}
		MYOUTPUT (fb->os, '\n');

		/* enumerate through S statistics */
		for (int j = 0; j < this->fam_total; j++){
			for (int i = 0; i < nsel; i++){
				int k = i*this->fam_total + j;
				if (i == 0){
					MYOUTPUT(fb->os, ir[k].ped->node->name << " ");
	        			MYOUTPUT(fb->os, "["<<(ir[k].flist->node->fa() ? ir[k].flist->node->fa()->id:"-")<<"x");
	        			MYOUTPUT(fb->os, (ir[k].flist->node->mo() ? ir[k].flist->node->mo()->id:"-") << "] ");
				}
				MYOUTPUT(fb->os, ir[k].s << " ");
			}
			MYOUTPUT(fb->os, '\n');
		}	

		/* generate ES statistics */
		MYOUTPUT (fb->os, "\nES\n");
		MYOUTPUT (fb->os, "pid    fid    ");
		for (int i = 0; i < this->getsize(); i+=this->fam_total){
			MYOUTPUT(fb->os, ir[i].sel << " ");
		}
		MYOUTPUT (fb->os, '\n');

		/* enumerate through ES statistics */
		for (int j = 0; j < this->fam_total; j++){
			for (int i = 0; i < nsel; i++){
				int k = i*this->fam_total + j;
				if (i == 0){
					MYOUTPUT(fb->os, ir[k].ped->node->name << " ");
	        			MYOUTPUT(fb->os, "["<<(ir[k].flist->node->fa() ? ir[k].flist->node->fa()->id:"-")<<"x");
	        			MYOUTPUT(fb->os, (ir[k].flist->node->mo() ? ir[k].flist->node->mo()->id:"-") << "] ");
				}
				MYOUTPUT(fb->os, ir[k].es << " ");
			}
			MYOUTPUT(fb->os, '\n');
		}	
	}
}

void IRES::additem(int idxsel, PEDIGREELIST* ped, NUCFAMILYLIST* flist, double s, double es, double v)
{
	struct IRESITEM t;
#if 0
	// for debugging only
	static int oldidxsel = 0;
	static PEDIGREELIST* oldped = NULL;
	static int fcount = 0;
	static int tcount = 0;
	static int ccount = 0;

	fcount++;
	tcount++;
	ccount++;
	if (idxsel != oldidxsel){
		if (oldidxsel != 0) MYOUTPUT(fb->os, "idxsel = " << oldidxsel << "\n")
		oldidxsel = idxsel;
		MYOUTPUT(fb->os, "total count = " << (tcount-1) << "," << (ccount -1) << "\n");
		tcount = 1;
	}

	if (ped != oldped){
		if (oldped != NULL) MYOUTPUT(fb->os, "[ped =" << oldped << " ][count = " << fcount << "]\n")
		oldped = ped;
		fcount = 1;
	}
#endif	
	t.sel 	 = idxsel;
	t.ped 	 = ped;
	t.flist  = flist;
	t.s 	 = s;
	t.es 	 = es;
	t.v 	 = v;

	this->ir.push_back(t);
	this->count++;
}

void IRES::addwt(int idxsel, double wt)
{
	struct WTITEM w;
	w.sel	= idxsel;
	w.wt	= wt;

	this->weight.push_back(w);
	this->wc++;
}


void IRES::rvtest(int option)
{
	double varW = 0.;
	double pv = 0.;

	int nsel;
	int nfam = this->fam_total;
	int nall = this->getsize()/this->fam_total;
	//PEDIGREELIST *parray[nfam];
	PEDIGREELIST **parray = new PEDIGREELIST *[nfam];

	if (this->wc == 0){
	  nsel = nall;
	}
	else {
	  nsel = this->wc;
	}

	// cout << "nfam=" << nfam << ",nsel=" << nsel << ",nall=" << nall << "\n"; 

	/* Perform the rare variant tests */
	Matrix S (nfam, nsel);
	Matrix ES(nfam, nsel);
	Matrix VS(nfam, nsel);
	Matrix VT(nfam, nsel);

	Matrix VE(nsel, nsel);
	Matrix VA(nsel, nsel);
	Matrix WT(nsel, 1);

	int l = 0;	// index in the weighted case
	for (int i=1; i<=nall; i++){
	  int k1 = (i-1)*nfam;
	  for (int j = 1; j <= nfam; j++) {
	    int k = k1 + j - 1;
	    parray[j-1] = this->ir[k].ped;

	    if (this->wc != 0) {
	      /* for the weighted case */
	      if (this->weight[l].sel == this->ir[k].sel) {
	        S(j,l+1)  = this->ir[k].s;
	        ES(j,l+1) = this->ir[k].es;
	        VS(j,l+1) = this->ir[k].v;
		if (j == nfam) {
	          WT(l+1, 1) = this->weight[l].wt;
	          l++;
	          if (l >= this->wc) goto cont;
		}
	      }     
	    }
	    else {
	      /* for the unweighted case */
	      S(j,i)  = this->ir[k].s;
	      ES(j,i) = this->ir[k].es;
	      VS(j,i) = this->ir[k].v;
	    }
	  }
	}
cont:
	VT = S - ES;
	if (this->wc !=0) {
	  DiagonalMatrix WD(nsel);
	  for (int i = 1; i <= nsel; i++){
	    WD(i) = 1.0/WT(i,1);
	  }
	  VT = VT * WD;
	}

	// compute W 
	double W = Sum(VT);

#if 1 
	if (this->empirical)
	// Compute varW
	{
	  /* count no. of families within each pedigree */
	  int ip=0;
	  //int pedarray[this->ped_count];
	  int *pedarray = new int[this->ped_count]; 
	  for (int i=0; i < this->ped_count; i++) pedarray[i] = 0;
	  PEDIGREELIST *curped = parray[0];

	  for (int j=0; j < nfam; j++)
	  {
	    // count the number of families within each pedigree
	    if (curped == parray[j]) pedarray[ip]++;
	    else
	    {
	      // move onto the next pedigree
	      ip++;
	      pedarray[ip]++;
	      curped = parray[j];
	    }
	  }

	  /* compute the variance */
	  int famindex = 1;
	  varW = 0;
	  double ttotal = 0.;
	  for (int i = 0; i < this->ped_count; i++)
	  {
	    int nf = pedarray[i];
	    double temp = 0.;
	    ColumnVector P(nsel); P = 0.0;
	    double pedvar;
	    for (int k = 1; k <= nsel; k++)
	    {
	      for (int j = famindex; j < famindex + nf; j++)
	      {
		temp += VT(j, k);
		P(k) += VT(j, k);
	      }
	    } 
	    Matrix PT(nsel, nsel);
	    PT = P*P.t();
	    pedvar = Sum(PT);
	    ttotal += temp;
	    // MYOUTPUT(fb->os, "Contribution from " << i << "th ped = " << temp << "(" << nf << ")\n");
	    // MYOUTPUT(fb->os, "Var Contribution from " << i << "th ped = " << pedvar << "\n");
	    famindex = famindex + nf;
	    varW += pedvar;
	  }

	  delete[] pedarray;

	  // MYOUTPUT(fb->os, "\nW (computed independently) = " << ttotal << "\n");
	  // MYOUTPUT(fb->os, "\nvarW (computed independently) = " << varW << "\n");
	}
	else
#endif
	{
	  VE = VT.t() * VT;

	  // Compute VA 
	  DiagonalMatrix DE(nsel);
	  for (int i = 1; i <= nsel; i++){
	    if (VE(i,i) == 0)
	      DE(i,i) = 0.;
	    else
	      DE(i,i) = 1./sqrt(VE(i,i));
	  }  

	  DiagonalMatrix D(nsel);
	  for (int i = 1; i <= nsel; i++) {
	    D(i,i) = 0.;
	    for (int j = 1; j <= nfam; j++)
	      D(i,i) += VS(j,i);
	  }  

	  for (int i = 1; i <= nsel; i++) {
	    if (this->wc == 0)
	      D(i,i) = sqrt(D(i,i));
	    else
	      D(i,i) = sqrt(D(i,i))/WT(i,1);
	  }
	
	  VA = D * DE * VE * DE * D;
	  varW = Sum(VA);
	}

	// Compute Z
	double Z = W/sqrt(varW);
	double Z1 = fabs(Z);

	/* look up the p-value */
	normbase(&Z1, &pv);
	pv = 2*(1-pv);

	/* Output the results */
	MYOUTPUT(fb->os, "\n");
	if (option == $UNW){ 
	  MYOUTPUT(fb->os, "Unweighted ");
	}
	else {
	  MYOUTPUT(fb->os, "Weighted ");
	}
	MYOUTPUT(fb->os, "FBAT rare variant statistics for the SNPs:\n");
	MYOUTPUT(fb->os, "\n");
	MYOUTPUT(fb->os, 
	"W           Var(W)      Z           p-value(2-sided)\n"); 
	MYOUTPUT(fb->os,
	"----------------------------------------------------\n");
	MYOUTPUT(fb->os, setw(11) << W << " ");
	MYOUTPUT(fb->os, setw(11) << varW << " ");
	MYOUTPUT(fb->os, setw(11) << Z << " ");
	MYOUTPUT(fb->os, setw(15) << setprecision(8) << scientific << pv << fixed << "\n");
	MYOUTPUT(fb->os,
	"----------------------------------------------------\n");
	MYOUTPUT(fb->os, "\n");

	// cleanup 
	delete [] parray;

}
