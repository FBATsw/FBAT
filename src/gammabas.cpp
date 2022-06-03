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
#include "normbase.h"
#include "gammabas.h"

#define lgamma mylgamma

#ifndef HUGEDBL
#define HUGEDBL 1e38
#endif


#undef TRUE 
#undef FALSE 
#define TRUE 1
#define FALSE 0

//static double lgamma(double /*x*/);
static double gammp(double /*x*/, double /*a*/);
static double gnorm(double /*x*/, double /*a*/);
static double gser(double /*x*/, double /*a*/, double /*gln*/);
static double gcf(double /*x*/, double /*a*/, double /*gln*/);
static double ppchi2(double /*p*/, double /*v*/,double /*g*/,long * /*ifault*/);
static double gammad(double /*x*/,double /*a*/);

void gammabase(double * x, double * a, double * p)
{
	*p = gammp(*a, *x);
} /*gammabase()*/

double          ppgamma(double p, double a, long * ifault)
{
	return (.5 * ppchi2(p, 2.0*a, lgamma(a), ifault));
} /*ppgamma()*/

/*
  Static Routines
*/

/*
  From Numerical Recipes, with normal approximation from Appl. Stat. 239
*/

#define EPSILON 1.0e-14
#define LARGE_A 10000.0
#define ITMAX 1000


static double   ln2pio2 = 0.9189385332046727417803297;

static double   clocal[] =
{
	-0.577191652,
	0.988205891,
	-0.897056937,
	0.918206857,
	-0.756704078,
	0.482199394,
	-0.193527818,
	0.035868343,
};

static double   casymp[] =
{
	12.,
	-360.,
	1260.,
	-1680.,
};


double mylgamma(double x)
{
	long            i;
	double          out, z, z2, value;

	if (x > 8.)
	{			/* asymptotic, 6.1.41 of Abramowitz and Stegun */
		z = 1 / x;
		z2 = z * z;
		out = 0.0;
		for (i = 0; i < 4; i++)
		{
			out += z / casymp[i];
			z *= z2;
		}
		return ((x - .5) * log(x) - x + ln2pio2 + out);
	}

	/* small x, 6.1.36 of Abramowitz and Stegun */
	value = 0.0;
	while (x < 1.0)
	{
		value -= log(x++);
	}
	while(x > 2.0)
	{
		value += log(--x);
	}
	
	z = x - 1.;
	z2 = z;
	out = 1.;
	for (i = 0; i < 8; i++)
	{
		out += clocal[i] * z2;
		z2 *= z;
	}
	return (log(out) + value);
} /*mylgamma()*/


static double   gammp(double a, double x)
{
	double          p;

	if (x <= 0.0 || a <= 0.0)
	{
		p = 0.0;
	}
	else if (a > LARGE_A)
	{
		p = gnorm(a, x);
	}
	else
	{
		p = (x < (a + 1.0)) ?
			gser(a, x, lgamma(a)) : 1.0 - gcf(a, x, lgamma(a));
	}
	return (p);
} /*gammp()*/

/* compute gamma cdf by a normal approximation */
static double   gnorm(double a, double x)
{
	double          p, sx;

	if (x > 0.0 && a > 0.0)
	{
		sx = sqrt(a) * 3.0 * (pow(x/a, 1.0/3.0) + 1.0/(a * 9.0) - 1.0);
		normbase(&sx, &p);
	}
	else
	{
		p = 0.0;
	}

	return (p);
} /*gnorm()*/

/* compute gamma cdf by its series representation */
static double   gser(double a, double x, double gln)
{
		double          p, sum, del, ap;
	long            n;

	if (x > 0.0 && a > 0.0)
	{
		ap = a;
		del = sum = 1.0/a;
		for (n = 1; n < ITMAX; n++)
		{
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < EPSILON)
			{
				break;
			}
		} /*for (n = 1; n < ITMAX; n++)*/
		p = sum * exp(-x + a * log(x) - gln);
	} /*if (x > 0.0 && a > 0.0)*/
	else
	{
		p = 0.0;
	}
	
	return (p);
} /*gser()*/

/* compute complementary gamma cdf by its continued fraction expansion */
static double   gcf(double a, double x, double gln)
{
	double          gold = 0.0, g, fac = 1.0, b1 = 1.0;
	double          b0 = 0.0, anf, ana, an, a1, a0 = 1.0;
	double          p;

	a1 = x;
	p = 0.0;
	for (an = 1.0; an <= ITMAX; an += 1.0)
	{
		ana = an - a;
		a0 = (a1 + a0 * ana) * fac;
		b0 = (b1 + b0 * ana) * fac;
		anf = an * fac;
		a1 = x * a0 + anf * a1;
		b1 = x * b0 + anf * b1;
		if (a1 != 0.0)
		{
			fac = 1.0/a1;
			g = b1 * fac;
			if (fabs((g - gold)/g) < EPSILON)
			{
				p = exp(-x + a * log(x) - gln) * g;
				break;
			}
			gold = g;
		} /*if (a1 != 0.0)*/
	} /*for (an = 1.0; an <= ITMAX; an += 1.0)*/
	return (p);
} /*gcf()*/

static double   gammad(double x,double a)
{
	double          cdf;

	gammabase(&x, &a, &cdf);
	return (cdf);
} /*gammad()*/

/*
  ppchi2.f -- translated by f2c and modified

  Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35
  To evaluate the percentage points of the chi-squared
  probability distribution function.

  p must lie in the range 0.000002 to 0.999998,
    (but I am using it for 0 < p < 1 - seems to work)
  v must be positive,
  g must be supplied and should be equal to ln(gamma(v/2.0))

  Auxiliary routines required: ppnd = AS 111 (or AS 241) and gammad.
*/
#define      LOG2        0.6931471805599453
#define      TWONINETHS  0.2222222222222222

static double   ppchi2(double p, double v, double g, long * ifault)
{
	/* Initialized data */
	/* coefficients are now built in to code */
	/*
	  static double pmin = 2e-6;
	  static double pmax = .999998;
	*/
	double   pmin = 0.0;
	double   pmax = 1.0;
	double   epsilon = 5e-7;
	/* System generated locals */
	double          d1, d2;

	/* Local variables */
	double   a, b, c, q, t, x, p1, p2, s1, s2, s3, s4, s5, s6, ch;
	double   xx;
	long     if1, iter = 0, maxit = 20;


	/* test arguments and initialise */
	ch = -1.0;
	if (p <= pmin || p >= pmax)
	{
		*ifault = 1;
		goto errorExit;
	}
	if (v <= 0.0)
	{
		*ifault = 2;
		goto errorExit;
	}
	*ifault = 0;
	xx = 0.5 * v;
	c = xx - 1.0;

	if (v < -1.24 * log(p))
	{
		/* starting approximation for small chi-squared */
		ch = pow(p * xx * exp(g + xx * LOG2), 1.0/xx);
		if (ch < epsilon)
		{
			goto normalExit;
		}
	}
	else if (v > .32)
	{
		/* call to algorithm AS 111 - note that p has been tested above. */
		/* AS 241 could be used as an alternative. */
		x = ppnd(p, &if1);

		/* starting approximation using Wilson and Hilferty estimate */
		p1 = TWONINETHS/v;
		/* Computing 3rd power */
		d1 = x * sqrt(p1) + 1.0 - p1;
		ch = v * (d1 * d1 * d1);

		/* starting approximation for p tending to 1 */
		if (ch > 2.2 * v + 6.0)
		{
			ch = -2.0 * (log(1.0 - p) - c * log(0.5 * ch) + g);
		}
	}
	else
	{
		/* starting approximation for v less than or equal to 0.32 */
		ch = .4;
		a = log(1.0 - p);
		do
		{
			q = ch;
			p1 = 1.0 + ch * (4.67 + ch);
			p2 = ch * (6.73 + ch * (6.66 + ch));
			d1 = -0.5 + (4.67 + 2.0 * ch)/p1;
			d2 = (6.73 + ch * (13.32 + 3.0 * ch))/p2;
			t = d1 - d2;
			ch -= (1.0 - exp(a + g + 0.5 * ch + c * LOG2) * p2/p1)/t;
		} while (fabs(q/ch - 1.0) > .01);
	}

	do
	{
		/* call to gammad and calculation of seven term Taylor series */
		q = ch;
		p1 = 0.5 * ch;
		p2 = p - gammad(p1, xx);
		t = p2 * exp(xx * LOG2 + g + p1 - c * log(ch));
		b = t/ch;
		a = 0.5*t - b*c;
		s1 = (210. + a*(140. + a*(105. + a*(84. + a*(70. + 60.*a)))))/420.;
		s2 = (420. + a*(735. + a*(966. + a*(1141. + 1278.*a))))/2520.;
		s3 = (210. + a*(462. + a*(707. + 932.*a)))/2520.;
		s4 = (252. + a*(672. + 1182.*a) + c*(294. + a*(889. + 1740.*a)))/5040.;
		s5 = (84. + 2264.*a + c*(1175. + 606.*a))/2520.;
		s6 = (120. + c*(346. + 127.*c))/5040.;
		d1 = (s3 - b*(s4 - b*(s5 - b*s6)));
		d1 = (s1 - b*(s2 - b*d1));
		ch += t*(1.0 + 0.5*t*s1 - b*c*d1);
	} while (++iter < maxit && fabs(q/ch - 1.0) > epsilon);
	if (iter < maxit)
	{
		goto normalExit;
	}
	*ifault = 4;

  errorExit:
  normalExit:
	return ch;
} /*ppchi2()*/
