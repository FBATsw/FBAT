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
//#include "dbug.h"
#include "gammabas.h"

#define lgamma mylgamma /*960220*/

#ifndef HUGEDBL
#define HUGEDBL 1e38
#endif

#undef TRUE
#undef FALSE
#define TRUE 1
#define FALSE 0

/* external function */
/* old
extern double macheps(), lgamma();
*/

#include "gammabas.h"
#include "betabase.h"

double macheps(void);
void betabase(double */*x*/, double */*a*/, double */*b*/, long */*gia*/, long */*gib*/, double */*cdf*/);

#if  FALSE
double          macheps(void)
{
	return (1e-14);
}
#else
#define macheps()  (1e-14)
#endif

/* forward declarations */
/*
  Static routines
*/

static double   logbeta(double p, double q)
{
	return (lgamma(p) + lgamma(q) - lgamma(p + q));
}

/*
  xinbta.f -- translated by f2c and modified

  algorithm as 109 appl. statist. (1977), vol.26, no.1
  (replacing algorithm as 64  appl. statist. (1973), vol.22, no.3)

  Remark AS R83 has been incorporated in this version.

  Computes inverse of the incomplete beta function
  ratio for given positive values of the arguments
  p and q, alpha between zero and one.
  log of complete beta function, beta, is assumed to be known.

  Auxiliary function required: betai

  SAE below is the most negative decimal exponent which does not
  cause an underflow; a value of -308 or thereabouts will often be
*/

static double   xinbta(double p, double q, double beta, double alpha,
					   long * ifault)
{
	/* Initialized data */
	double   sae = -30.0;	/* this should be sufficient */

	/* System generated locals */
	double          ret_val, d_1, d_2;

	/* Local variables */
	long     indx;
	double   prev, a, g, h, r, s, t, w, y, yprev, pp, qq;
	double   sq, tx, adj, acu;
	long     iex;
	double   fpu, xin;

	/* Define accuracy and initialise. */
	fpu = sae * 10.;
	ret_val = alpha;

	/* test for admissibility of parameters */
	if (p <= 0.0 || q <= 0.0)
	{
		*ifault = 1;
		return ret_val;
	}
	if (alpha < 0.0 || alpha > 1.0)
	{
		*ifault = 2;
		return ret_val;
	}
	*ifault = 0;
	if (alpha == 0.0 || alpha == 1.0)
	{
		return ret_val;
	}
	
	/* change tail if necessary */
	if (alpha <= .5)
	{
		a = alpha;
		pp = p;
		qq = q;
		indx = FALSE;
	}
	else
	{
		a = 1.0 - alpha;
		pp = q;
		qq = p;
		indx = TRUE;
	}

	/* calculate the initial approximation */
	r = sqrt(-log(a * a));
	y = r - (r * .27061 + 2.30753) / (1.0 + (r * .04481 + .99229) * r);
	if (pp > 1.0 && qq > 1.0)
	{
		r = (y * y - 3.0) / 6.0;
		s = 1.0 / (pp + pp - 1.0);
		t = 1.0 / (qq + qq - 1.0);
		h = 2.0 / (s + t);
		d_1 = y * sqrt(h + r) / h;
		d_2 = (t - s) * (r + 5.0 / 6.0 - 2.0 / (3.0 * h));
		w = d_1 - d_2;
		ret_val = pp / (pp + qq * exp(w + w));
	} /*if (pp > 1.0 && qq > 1.0)*/
	else
	{
		r = qq + qq;
		t = 1.0 / (qq * 9.);
		/* Computing 3rd power */
		d_1 = 1.0 - t + y * sqrt(t);
		t = r * (d_1 * d_1 * d_1);
		if (t <= 0.0)
		{
			ret_val = 1.0 - exp((log((1.0 - a) * qq) + beta) / qq);
		}
		else
		{
			t = (4.0 * pp + r - 2.0) / t;
			if (t <= 1.0)
			{
				ret_val = exp((log(a * pp) + beta) / pp);
			}
			else
			{
				ret_val = 1.0 - 2.0 / (t + 1.0);
			}
		}
	} /*if (pp > 1.0 && qq > 1.0){}else{}*/


	/*
      solve for x by a modified newton-raphson method, using the function betai
    */
	r = 1.0 - pp;
	t = 1.0 - qq;
	yprev = 0.0;
	sq = 1.0;
	prev = 1.0;
	if (ret_val < 1e-4)
	{
		ret_val = 1e-4;
	}
	else if (ret_val > .9999)
	{
		ret_val = .9999;
	}
	/* Computing MAX, two 2nd powers */
	d_1 = -5.0 / (pp * pp) - 1.0 / (a * a) - 13.0;
	iex = (sae > d_1) ? sae : d_1;
	acu = pow(10.0, (double) iex);
	do
	{
		y = betai(ret_val, pp, qq);
		xin = ret_val;
		y = (y - a) * exp(beta + r * log(xin) + t * log(1.0 - xin));
		if (y * yprev <= 0.0)
		{
			prev = (sq > fpu) ? sq : fpu;
		}
		g = 1.0;
		do
		{
			adj = g * y;
			sq = adj * adj;
			if (sq < prev)
			{
				tx = ret_val - adj;
				if (tx >= 0.0 && tx <= 1.0)
				{
					if (prev <= acu || y * y <= acu)
					{
						if (indx)
						{
							ret_val = 1.0 - ret_val;
						}
						return ret_val;
					}
					if (tx != 0.0 && tx != 1.0)
					{
						break;
					}
				}
			}
			g /= 3.0;
		} while (TRUE);
		if (tx == ret_val)
		{
			if (indx)
			{
				ret_val = 1.0 - ret_val;
			}
			break;
		}
		ret_val = tx;
		yprev = y;
	} while (TRUE);
	return (ret_val);
} /*xinbta()*/

#if (0) /*betabase calls changed to call betai() directly*/
void betabase(x, a, b, gia, gib, cdf)
long           *gia, *gib;
double         *x, *a, *b, *cdf;
{
	*cdf = betai(*x, *a, *b);
}
#endif

double          ppbeta(double p, double a, double b, long * ifault)
{
	return (xinbta(a, b, logbeta(a, b), p, ifault));
}

#define Min(x,y) (((x) < (y)) ? (x) : (y))
#define Max(x,y) (((x) > (y)) ? (x) : (y))

double   betai(double x, double pin, double qin)
{
	/* Translated from FORTRAN
           july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
           based on bosten and battiste, remark on algorithm 179, comm. acm,
           v 17, p 153, (1974).
      
           input arguments --
           x      upper limit of integration.  x must be in (0,1) inclusive.
           p      first beta distribution parameter.  p must be gt 0.0.
           q      second beta distribution parameter.  q must be gt 0.0.
           betai  the incomplete beta function ratio is the probability that a
                  random variable from a beta distribution having parameters
                  p and q will be less than or equal to x.
        */
	double          c, finsum, p, ps, q, term, xb, xi, y, dbetai, p1;
	static double   eps = 0.0, alneps = 0.0, sml = 0.0, alnsml = 0.0;
	long            i, n, ib;

	/* I'm not sure these tolerances are optimal */
	if (eps == 0.0)
	{
		eps = macheps();
		alneps = log(eps);
		sml = macheps();
		alnsml = log(sml);
	}

	y = x;
	p = pin;
	q = qin;
	if (q > p || x >= 0.8)
	{
		if (x >= 0.2)
		{
			y = 1.0 - y;
			p = qin;
			q = pin;
		}
	} /*if (q > p || x >= 0.8)*/
	
	if ((p + q) * y / (p + 1.0) < eps)
	{
		dbetai = 0.0;
		xb = p * log(Max(y, sml)) - log(p) - logbeta(p, q);
		if (xb > alnsml && y != 0.0)
		{
			dbetai = exp(xb);
		}
		if (y != x || p != pin)
		{
			dbetai = 1.0 - dbetai;
		}
	} /*if ((p + q) * y / (p + 1.0) < eps)*/
	else
	{
		/*
	         *  evaluate the infinite sum first.  term will equal
	         *  y**p/beta(ps,p) * (1.-ps)-sub-i * y**i / fac(i) .
	         */
		ps = q - floor(q);
		if (ps == 0.0)
		{
			ps = 1.0;
		}
		xb = p * log(y) - logbeta(ps, p) - log(p);
		dbetai = 0.0;
		if (xb >= alnsml)
		{
			dbetai = exp(xb);
			term = dbetai * p;
			if (ps != 1.0)
			{
				n = (long) (alneps / log(y));
				if(n < 4)
				{
					n = 4;
				}
				for (i = 1; i <= n; i++)
				{
					xi = i;
					term *= (xi - ps) * y / xi;
					dbetai += term / (p + xi);
				} /*for (i = 1; i <= n; i++)*/
			} /*if (ps != 1.0)*/
		} /*if (xb >= alnsml)*/
		/*
	         * now evaluate the finite sum, maybe.
         */
		if (q > 1.0)
		{

			xb = p * log(y) + q * log(1.0 - y) - logbeta(p, q) - log(q);
			ib = (long) (xb / alnsml);
			if(ib < 0)
			{
				ib = 0;
			}
			term = exp(xb - ((double) ib) * alnsml);
			c = 1.0 / (1.0 - y);
			p1 = q * c / (p + q - 1.0);

			finsum = 0.0;
			n = q;
			if (q == (double) n)
			{
				n--;
			}
			for (i = 1; i <= n; i++)
			{
				if (p1 <= 1.0 && term / eps <= finsum)
				{
					break;
				}
				xi = i;
				term = (q - xi + 1.0) * c * term / (p + q - xi);

				if (term > 1.0)
				{
					ib--;
					term *= sml;
				}
				
				if (ib == 0)
				{
					finsum += term;
				}
			} /*for (i = 1; i <= n; i++)*/

			dbetai += finsum;
		} /*if (q > 1.0)*/
		if (y != x || p != pin)
		{
			dbetai = 1.0 - dbetai;
		}
		dbetai = Max(Min(dbetai, 1.0), 0.0);
	} /*if ((p + q) * y / (p + 1.0) < eps){}else{}*/
	return (dbetai);
} /*betai()*/


/*
	Algorithm AS226 Appl. Statist. (1987) Vol. 36, No. 2
	Incorporates modification as R84 from as Vol. 39, pp311-2, 1990
	Incorporates bug fix by C. Bingham.

	Returns the cumulative probability of x for the non-central beta
	distribution with parameters a, b and non-centrality lambda

	Auxiliary routines required: alogam - log-gamma function (acm
	291 or as 245), and betain - incomplete-beta function (AS 63)
*/

/*      change ERRMAX and ITRMAX if desired ... */

#define ERRMAX   1e-6
#define ITRMAX   100

double betanc(double x, double a, double b, double lambda, long * ifault)
{
	double       ax, c, errbd, gx, q, sumq, temp, xj;
	double       a0, x0;
	double       top, bottom;
	double       ualpha = 5.0;
	double       retVal = 0.0;
	long         iter;

	*ifault = 0;
	if(lambda < 0.0 || a <= 0.0 || b <= 0.0)
	{
		*ifault = 2;
	}
	else if(x < 0.0 || x > 1.0)
	{
		*ifault = 3;
	}
	else if (x == 0.0 || x == 1.0)
	{
		retVal = x;
	}
	else if(lambda == 0.0)
	{
		retVal = betai(x, a, b);
	}
	else
	{/* lambda != 0.0 */
		c = lambda*0.5;

		/*      initialize the series ... */

		x0 = floor(c - ualpha*sqrt(c));
		if(x0 < 0.0)
		{
			x0 = 0.0;
		}
			
		a0 = a + x0;
		temp = betai(x, a0, b);
		gx = exp(a0*log(x) + b*log(1.0 - x) -
				 logbeta(a0, b) - log(a0));
		if (a0 <= a)
		{
			q = exp(-c);
		}
		else
		{
			q = exp(-c + x0*log(c) - lgamma(x0 + 1.0));
		}
		xj = x0 + 1.0;			/* corrected from xj = 1.0; by kb */

		ax = q*temp;
		sumq = 1.0 - q;
		retVal = ax;
		top = a0 + b;
		bottom = a0 + 1.0;
		for (iter = 0;iter < ITRMAX;iter++)
		{
			temp -= gx;
			gx = x*(top++)*gx/(bottom++);
			q *= c/xj++;
			sumq -= q;
			ax = temp*q;
			retVal += ax;

			/*      check for convergence and act accordingly... */

			if((errbd = (temp - gx)*sumq) <= ERRMAX)
			{
				break;
			}
		} /*for (iter = 0;iter < ITRMAX;iter++)*/
		if (errbd > ERRMAX)
		{
			*ifault = iter;
		}
	}

	return (retVal);
} /*betanc()*/

