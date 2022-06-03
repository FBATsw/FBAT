/*
*(C)* The file is part of the source distribution of MacAnova
*(C)* version 4.00 or later
*(C)*
*(C)* Copyright (c) 1996 by Gary Oehlert and Christopher Bingham
*(C)* unless indicated otherwise
*(C)*
*(C)* You may give out copies of this software;  for conditions see the
*(C)* file COPYING included with this distribution
*(C)*
*(C)* This file is distributed WITHOUT ANY WARANTEE; without even
*(C)* the implied warantee of MERCHANTABILITY or FITNESS FOR
*(C)* A PARTICULAR PURPOSE
*(C)*
*/



#include <math.h>
#include <stdio.h>
#ifndef HUGEDBL
#define HUGEDBL 1e38
#endif

#include "normbase.h"

/* Normal CDF */
/*
   Function normbase() implements approximations from
   Cody, W. J. Rational Chebyshev approximations for the error function
     Math. Comp. 23 (1969) 631-637
   Also in W. J. Kennedy and J. E. Gentle, Statistical Computing,
      Marcel Dekker 1980, pp. 90-92
*/

/* coefficients from Cody Table II, n=3*/
#define P1_0 242.66795523053175
#define P1_1 21.979261618294152
#define P1_2 6.9963834886191355
#define P1_3 -.035609843701815385
#define Q1_0 215.05887586986120
#define Q1_1 91.164905404514901
#define Q1_2 15.082797630407787
#define Q1_3 1.0

/* coefficients from Cody Table III, n=7*/
#define P2_0 300.4592610201616005
#define P2_1 451.9189537118729422
#define P2_2 339.3208167343436870
#define P2_3 152.9892850469404039
#define P2_4 43.16222722205673530
#define P2_5 7.211758250883093659
#define P2_6 .5641955174789739711
#define P2_7 -.0000001368648573827167067
#define Q2_0 300.4592609569832933
#define Q2_1 790.9509253278980272
#define Q2_2 931.3540948506096211
#define Q2_3 638.9802644656311665
#define Q2_4 277.5854447439876434
#define Q2_5 77.00015293522947295
#define Q2_6 12.78272731962942351
#define Q2_7 1.0

/* coefficients from Cody Table IV, n=4*/
#define P3_0 -.00299610707703542174
#define P3_1 -.0494730910623250734
#define P3_2 -.226956593539686930
#define P3_3 -.278661308609647788
#define P3_4 -.0223192459734184686
#define Q3_0 .0106209230528467918
#define Q3_1 .191308926107829841
#define Q3_2 1.05167510706793207
#define Q3_3 1.98733201817135256
#define Q3_4 1.0

#define SQRT2 1.414213562373095049
#define SQRTPI 1.772453850905516027

void normbase(double * x, double * phi)
{
	long            sn;
	double          r1, r2, y, y2, erfVal, erfcVal, z;

	y = *x / SQRT2;
	if (y < 0)
	{
		y = -y;
		sn = -1;
	}
	else
	{
		sn = 1;
	}

	y2 = y * y;
	if (y < 0.46875)
	{
		r1 = ((P1_3 * y2 + P1_2) * y2 + P1_1) * y2 + P1_0;
		r2 = ((Q1_3 * y2 + Q1_2) * y2 + Q1_1) * y2 + Q1_0;
		erfVal = y * r1 / r2;
		if (sn == 1)
		{
			*phi = 0.5 + 0.5 * erfVal;
		}
		else
		{
			*phi = 0.5 - 0.5 * erfVal;
		}
	} /*if (y < 0.46875)*/
	else
	{		
		if (y < 4.0)
		{
			r1 = ((((((P2_7 * y + P2_6) * y + P2_5) * y + P2_4) * y + P2_3) * 
				   y + P2_2) * y + P2_1) * y + P2_0;
			r2 = ((((((Q2_7 * y + Q2_6) * y + Q2_5) * y + Q2_4) * y + Q2_3) *
				   y + Q2_2) * y + Q2_1) * y + Q2_0;
			erfcVal = exp(-y2) * r1 / r2;
		} /*if (y < 4.0)*/
		else
		{
			z = y2*y2;
			r1 = (((P3_4 * z + P3_3) * z + P3_2) * z + P3_1) * z + P3_0;
			r2 = (((Q3_4 * z + Q3_3) * z + Q3_2) * z + Q3_1) * z + Q3_0;
			erfcVal = (exp(-y2) / y) * (1.0 / SQRTPI + r1 / (r2 * y2));
		} /*if (y < 4.0){}else{}*/
		if (sn == 1)
		{
			*phi = 1.0 - 0.5 * erfcVal;
		}
		else
		{
			*phi = 0.5 * erfcVal;
		}
	} /*if (y < 0.46875){}else{}*/	
} /*normbase()*/

/* Normal inverse */

#define SPLIT 0.42
#define A0 2.50662823884
#define A1 -18.61500062529
#define A2 41.39119773534
#define A3 -25.44106049637

#define B1 -8.47351093090
#define B2 23.08336743743
#define B3 -21.06224101826
#define B4 3.13082909833

#define C0 -2.78718931138
#define C1 -2.29796479134
#define C2 4.85014127135
#define C3 2.32121276858

#define D1 3.54388924762
#define D2 1.63706781897
/*
  
  Algorithm as 111 Applied Statistics (1977), Vol 26 No 1 Page 121
  Produces normal deviate corresponding to lower tail area of p.
  The hash sums are the sums of the moduli of the coefficients;
  they have no inherent meanings but are incuded for use in
  checking transcriptions.  Functions fabs, log and sqrt are used.  

  Derived from AS111 Fortran version
*/

double          ppnd(double p, long  *ifault)
{
	double          q, r, ppn;

	*ifault = 0;
	q = p - 0.5;
	if (fabs(q) <= SPLIT)
	{
		r = q * q;
		ppn = q * (((A3 * r + A2) * r + A1) * r + A0) /
				((((B4 * r + B3) * r + B2) * r + B1) * r + 1.0);
	}
	else
	{
		r = p;
		if (q > 0.0)
		{
			r = 1.0 - p;
		}
		if (r <= 0.0)
		{
			*ifault = 1;
			ppn = 0.0;
		}
		else
		{			
			r = sqrt(-log(r));
			ppn = (((C3 * r + C2) * r + C1) * r + C0) /
				((D2 * r + D1) * r + 1.0);
			if (q < 0.0)
			{
				ppn = -ppn;
			}
		}
	}
	
	return (ppn);
} /*ppnd()*/
