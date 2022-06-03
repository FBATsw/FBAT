/* Long period random number generator (>2.0E18)
   generate a random number between 0 and 1 uniformly
   adapted from Numberic Receipt in C ran2(long *)
*/

#include "rand.h"
#include <math.h>

const double EPS = 1.2e-7;
const double RNMX = 1.0 - EPS;

static long idum = 0;

void rand_seed(long seed) {
	
	idum = seed;
	for (seed=0; seed<32; seed++)
		rand_gaussian();
}

double rand_uniform() {
	
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[32];
	double temp;
	
	if (idum<=0) {
		idum2 = idum = idum==0?1:-idum;
		for (j=32+7; j>=0; j--) {
			k = idum/53668;
			idum = 40014*(idum-k*53668) - k*12211;
			if (idum<0) idum += 2147483563;
			if (j<32) iv[j] = idum;
		}
		iy = iv[0];
	}
	k = idum/53668;
	idum = 40014*(idum-k*53668) - k*12211;
	if (idum<0) idum += 2147483563;
	k = idum2/52774;
	idum2 = 40692*(idum2-k*52774)-k*3791;
	if (idum2<0) idum2 += 2147483399;
	j = iy/67108862;
	iy = iv[j] - idum2;
	iv[j] = idum;
	if (iy<1) iy += 2147483562;
	temp = (double)iy/2147483563.0;
	return temp>RNMX? RNMX : temp;

}

double rand_gaussian() {
	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;
	
	if (iset==0) {
		do {
			v1 = 2.0*rand_uniform() - 1.0;
			v2 = 2.0*rand_uniform() - 1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1.0 || rsq == 0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	} else {
		iset = 0;
		return gset;
	}
}