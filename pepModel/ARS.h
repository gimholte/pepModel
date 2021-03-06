/*
 * ARS.h
 *
 *  Created on: Jan 24, 2012
 *      Author: hoblitz
 */

#ifndef ARS_H_
#define ARS_H_

#define NMAX 100

struct ARS_WORKSPACE{
	double hwv[NMAX];
	double hpwv[NMAX];
	double scum[NMAX];
	double scum_norm[NMAX];
	double s[NMAX];
	double z[NMAX];
};

typedef struct ARS_WORKSPACE ARS_workspace;

#endif /* ARS_H_ */

#include <R.h>
#include <Rmath.h>
#include "RngStream.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Rversion.h>
#include <float.h>

#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif

double sample_conditional(double* restrict x,
		int* restrict num_x,
		double* restrict argvec,
		int* restrict arglen,
		ARS_workspace *ws,
		RngStream rng,
		double (*h)(double, double *, int *),
		double (*h_prime)(double , double *, int *));

int update_hull(double* restrict x,
		ARS_workspace *ws,
		double* restrict argvec,
		int* restrict arglen,
		int* restrict num_x,
		double xnew,
		double hnew,
		int l_section,
		double* restrict huzmax,
		double (*h)(double, double *, int*), double (*h_prime)(double , double *, int*));

double sample_hull(double* restrict x,
		ARS_workspace *ws,
		int* restrict num_x,
		int* restrict section,
		double p,
		double huzmax);

void initialize_hull(double* restrict x,
		ARS_workspace *ws,
		int num_x,
		double huzmax);

double log_apb(double loga, double logb);
