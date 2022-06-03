// xfunc.h
#pragma once

enum GEN_MODEL { model_additive, model_dominant, model_recessive, model_genotype, model_estimate, model_user_defined };

double xfunc(int a1, int a2, int a, GEN_MODEL model);

double xfunc(int g, int a, GEN_MODEL model);

int gtindex(int a0, int a1);

void gtindex2allele(int gtidx, char a[]);

double pdf_norm(double x, double m, double sd);	// return the pdf of x for N(m, sd^2)