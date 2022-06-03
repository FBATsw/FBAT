#include <math.h>
#include "xfunc.h"

void gtindex2allele(int gtidx, char a[]) {
	
	int i;
	
	gtidx++;
	for (i=1; (i*(i+1))/2<gtidx; i++)
		;
	a[1] = i;
	a[0] = gtidx-(i*(i-1))/2;
		
}

int gtindex(int a0, int a1) {

	return a0<=a1? (a1*(a1-1))/2+a0-1 : (a0*(a0-1))/2+a1-1; 

}

double xfunc(int g, int a, GEN_MODEL model) {
	
	if (g<0 || a<0)
		return 0;
	char b[2];
	gtindex2allele(g, b);
	
	return xfunc(b[0]-1,b[1]-1,a, model);

}

double xfunc(int a1, int a2, int a, enum GEN_MODEL model) {
	
	if (a1==a && a2==a)
		return model==model_additive?2.0:1.0;
	else if ((a1==a || a2==a) && model!=model_recessive)
		return 1.0;
	
	return 0.0;
}

double pdf_norm(double x, double m, double sd) {
	
	double pi = 3.1415926536;
	
	return exp(-(x-m)*(x-m)/(2*sd*sd))/(sd*sqrt(2*pi));

}
