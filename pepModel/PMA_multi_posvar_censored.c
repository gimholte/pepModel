/*
 * PMA_multi_posvar_censored.c
 *
 *  Created on: Feb 14, 2012
 *      Author: Gregory Imholte
 */

#include <omp.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rversion.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "norm_gamma_generation.h"
#include "RngStream.h"
#include "ARS.h"

#if (R_VERSION >= R_Version(2,3,0))
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>
#endif

#define ZZERO DBL_MIN

//
double lc_AB(double x, double *argvec);

double lcp_AB(double x, double *argvec);

double lc_alpha(double x, double *argvec);

double lcp_alpha(double x, double *argvec);

inline double log1m_from_logit(double x);

inline double log_from_logit(double x);

inline double expit(double x);

void update_u_pars(double *P, double *a_0, double *b_0,
		double *xA0, double *xB0, int n_position, double psi_a,
		double psi_b, RngStream rng, int nmax, double eps);

void update_var_pars(double *Sig2, double *alpha, double *beta,
		double alpha_0, double beta_0,  double *xAlpha,
		int n_max, int ind_pep, int n_peptide,
		int n_position, RngStream rng, double eps);

void update_mu_pep(double *Mu, double **Exprs, double **Alpha,
		int **Gamma, double Sig2, double kappa, int n_indiv, int pep,
		RngStream rng);

void update_indiv_j1(double *Exprs, double *Alpha, double *mu_j, int* Omega_Ind,
		double *Omega_Logit, int *Gamma,
		double *Sig2, double m, double c_var,
		double tau_0, int *pnum, int *pstart,
		int *n_position, int *n_peptide, int ind_pep, RngStream rng);

int update_position_p1(double **Exprs, int* Omega_Ind, double *Omega_Logit,
		double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P,
		double a_0, double b_0, double lambda_a, double lambda_b,
		int p, int *n_indiv, int *pnum, int *pstart, RngStream rng, int nmax,
		double *xA, double *xB, ARS_workspace *workspace, double eps,
		int ind_pep, int MRF, double alpha, double beta);

void InitializeChain1(int *Omega_Ind, double *Omega_Logit, double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P, double *Theta,
		double *c_var, double *m, double *alpha, double *beta,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		int nmax, double **xA, double **xB,
		int MRF, int ind_pep);

void update_global_parms(double **Alpha, int **Gamma, double *m, double *c_var,
		double *A, double *B, double *Sig2, double m_0, double v_0, double alpha_0,
		double beta_0, double *alpha, double *beta, int lambda_prior,
		double *lambda_a, double *lambda_b, double r_a, double r_b,
		int n_position, int n_peptide, int n_indiv, int *accept_m,
		int *accept_c, RngStream rng, double adptm, double adptv, int ind_pep,
		int var_prior, double *xAlpha, int n_max, double eps);

inline double expit(double x);

void tnorm_test(double* x, int *n, double *m, double *sigmasqr);

void update_lambdas(double *lambda_a, double *lambda_b, double r_a,
		double r_b, double *A, double *B, int n_position,
		RngStream rng);

void update_MRF(double *P, double *Theta, double *Omega,
		double *kappa, int *pnum, int *pstart, int *n_position,
		int *accept, RngStream rng);

void update_prob_include(int *n_peptide, int *n_indiv, int **Gamma, int **ProbSum, double *mean_fitted,
		double **Alpha, double *Mu, double n);

void update_kappa(double *Mu, double *kappa, int n_peptide, RngStream rng);

double truncNorm_parallel(double mean, double sigmasqr, RngStream rng);

double truncNorm(double mean, double sigmasqr);

void update_data(double *D, int cen_num, int* cen_ind, int* cen_pep, double *Y, double **Exprs,
		double **Alpha, double *Mu, int n_peptide, double *Sig2, int *cen_pos);

void store_mcmc_output1(double *Mu, double *A, double *B, double *P, double *Sig2, double *D,
		double *Theta, double *Omega_Logit, int* Omega_Ind, double kappa, double alpha, double beta,
		double m, double c_var, int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		double lambda_a, double lambda_b, double a_0, double b_0, int MRF, int ind_pep,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *THETAFILE, FILE *OFILE);

void finalize_prob_include(int *n_iter, int *n_peptide, int *n_indiv,
		double *OutProbs, int **ProbSum);

/*
 * Y: expressions column major vector, rows are peptides, columns are individuals
 * hyper_parms: model hyperparameters
 * pstart: starting indices of positions
 * pnum: number of clades per position
 * cen_ind: indicies of individual with censoring
 * cen_pep: indicies denoting which peptide censored
 * cen_num: total number of censored peptides
 *
 * nP: number of threads for parallel processing
 * mean_fitted: mean fitted value of all peptides
 */

void PMA_mcmc_MS(double *Y, double *hyper_parms, int *pstart,
		int *pnum, int *n_position, int *n_peptide, int *n_indiv, int *nP,
		int *cen_ind, int *cen_pep, int *cen_num, int *cen_pos,
		int *n_iter, int *n_sweep, int *n_burn,
		int *MRF, int *lambda_prior, int *ind_pep, int *var_prior, int *update_u,
		double *OutProbs, double *mean_fitted, int *write, double *eps, int *nmax,
		int *accept, int *silent)
{
	R_CStackLimit=(uintptr_t)-1;
	int p, i, n = 0, j, k;
	*nmax = (*nmax < NMAX) ? *nmax : NMAX;

	double **Exprs;
	double **Alpha;
	double **xA, **xB;
	double *Omega_Logit, *Theta;
	int *Omega_Ind;
	int **ProbSum;
	int **Gamma;
	double *Mu;
	double *A, *B, *P, *Sig2, *D;
	double *xAlpha, *xA0, *xB0;

	int accept_m = 0, accept_c = 0;

	double a_0, b_0, lambda_a, lambda_b, tau, alpha_0, beta_0, m_0, v_0;
	double r_a, r_b, kappa, psi_a, psi_b;

	a_0 = hyper_parms[0];
	b_0 = hyper_parms[1];
	psi_a = hyper_parms[0];
	psi_b = hyper_parms[1];

	r_a = hyper_parms[2];
	r_b = hyper_parms[3];
	lambda_a = hyper_parms[2];
	lambda_b = hyper_parms[3];

	tau = hyper_parms[4];
	alpha_0 = hyper_parms[5];
	beta_0 = hyper_parms[6];
	m_0 = hyper_parms[7];
	v_0 = hyper_parms[8];

	kappa = hyper_parms[9];

	Exprs = (double**) malloc(*n_indiv*sizeof(double*));
	Alpha = (double**) malloc(*n_indiv*sizeof(double*));

	// xA and xB hold starting values for ARS of a_p and b_p
	xA = (double**) malloc(*n_position*sizeof(double*));
	xB = (double**) malloc(*n_position*sizeof(double*));
	xAlpha = (double*) malloc(*nmax*sizeof(double));
	xAlpha[0] = 1.0;
	xAlpha[1] = 2.0;

	xA0 = (double*) malloc(*nmax*sizeof(double));
	xB0 = (double*) malloc(*nmax*sizeof(double));
	xA0[0] = .1;
	xB0[0] = .1;
	xB0[1] = 2.0;
	xA0[1] = 2.0;

	Omega_Ind = (int*) malloc(*n_peptide*sizeof(int));
	Omega_Logit = (double*) malloc(*n_peptide*sizeof(double));

	Gamma = (int**) malloc(*n_indiv*sizeof(int*));
	ProbSum = (int**) malloc(*n_indiv*sizeof(int*));

	for(p = 0; p < *n_position; p++)
	{
		xA[p] = (double*) malloc(*nmax*sizeof(double));
		xB[p] = (double*) malloc(*nmax*sizeof(double));
	}

	for(i = 0; i < *n_indiv; i++)
	{
		Exprs[i] = (double*) malloc(*n_peptide*sizeof(double));
		Alpha[i] = (double*) malloc(*n_peptide*sizeof(double));
		Gamma[i] = (int*) malloc(*n_peptide*sizeof(int));
		ProbSum[i] = (int*) malloc(*n_peptide*sizeof(int));
	}

	if(*ind_pep == 1)
	{
		Sig2 = (double*) malloc(*n_peptide*sizeof(double));
	}
	else
	{
		Sig2 = (double*) malloc(*n_position*sizeof(double));
	}

	if(*MRF == 1)
	{
		Theta = (double*) malloc(*n_position*sizeof(double));
	}

	A = (double*) malloc(*n_position*sizeof(double));
	B = (double*) malloc(*n_position*sizeof(double));
	P = (double*) malloc(*n_position*sizeof(double));
	Mu = (double*) malloc(*n_peptide*sizeof(double));

	// check whether our data is censored at all,
	// if so prepare for augmentation.
	if(*cen_num > 0)
	{
		D = (double*) malloc(*cen_num*sizeof(double));
		for(i = 0; i < *cen_num; i++)
		{
			D[i] = 0.0;
		}
	}

	double m = 0.0;
	double c_var = 1.0;
	double alpha = alpha_0;
	double beta = beta_0;
	double adptm = 5.0;
	double adptv = 1.0;

	FILE *AFILE, *BFILE, *PFILE, *VARFILE, *Sig2FILE, *MUFILE, *DFILE, *THETAFILE, *OFILE;

	if(*write == 1)
	{
		AFILE = fopen("afile.txt", "w");
		BFILE = fopen("bfile.txt", "w");
		PFILE = fopen("pfile.txt", "w");
		OFILE = fopen("ofile.txt", "w");
		VARFILE = fopen("varfile.txt", "w");
		Sig2FILE = fopen("sig2file.txt", "w");
		MUFILE = fopen("mufile.txt", "w");
		if(*MRF == 1)
		{
			THETAFILE = fopen("thetafile.txt", "w");
		}
		if(*cen_num > 0)
		{
			DFILE = fopen("dfile.txt", "w");
		}
	}

	//statically allocated struct for Adaptive Rejection Sampling
	ARS_workspace workspace;
	*nmax = (*nmax < NMAX) ? *nmax : NMAX;

	// initialize the probsum values
	for(i = 0; i < *n_indiv; i++)
	{
		for(p = 0; p < *n_peptide; p++)
		{
			ProbSum[i][p] = 0;
			Exprs[i][p] = Y[*n_peptide*i + p];
		}
	}

	int th_id, pep;
	RngStream rng[*nP];
	double sigma;
	for(i = 0; i < *nP; i++)
	{
		rng[i] = RngStream_CreateStream("");
		//RngStream_IncreasedPrecis(rngs[i], 1);
	}
	GetRNGstate();
	InitializeChain1(Omega_Ind, Omega_Logit, Alpha,
			Gamma, Sig2, Mu, A, B, P, Theta,
			&c_var, &m, &alpha, &beta,
			n_position, pstart, pnum, n_peptide, n_indiv,
			*nmax, xA, xB,
			*MRF, *ind_pep);
	PutRNGstate();

	Rprintf("parameters initialized \n");
	for(i = 0; i <= (*n_iter)*(*n_sweep) + *n_burn; i++)
	{
		if((i%1000 == 0) & (i > 0) & (*silent == 0))
		{
			Rprintf("***** iteration %d ***** \n", i);
			if(i > *n_burn)
			{
				Rprintf("m acceptance prob = %f\n", (double)(accept_m)/(double)(i - *n_burn));
				Rprintf("c acceptance prob = %f\n", (double)(accept_c)/(double)(i - *n_burn));
			}
		}

		R_CheckUserInterrupt();

		update_global_parms(Alpha, Gamma, &m, &c_var, A, B, Sig2,
				 m_0, v_0, alpha_0, beta_0, &alpha, &beta,
				*lambda_prior, &lambda_a, &lambda_b, r_a, r_b,
				*n_position, *n_peptide, *n_indiv, &accept_m,
				&accept_c, rng[0], adptm, adptv, *ind_pep, *var_prior,
				xAlpha, *nmax, *eps);

		update_kappa(Mu, &kappa, *n_peptide, rng[0]);

#pragma omp parallel private(th_id, workspace, sigma) num_threads(*nP)
		{
			th_id = omp_get_thread_num();
#pragma omp for nowait
			for(p = 0; p < *n_position; p++)
			{
				update_position_p1(Exprs, Omega_Ind, Omega_Logit, Alpha, Gamma,
						Sig2, Mu, A, B, P,
						a_0, b_0, lambda_a, lambda_b,
						p, n_indiv, pnum, pstart, rng[th_id], *nmax,
						xA[p], xB[p], &workspace, *eps, *ind_pep, *MRF,
						alpha, beta);
			}

#pragma omp for nowait
			for(j = 0; j < *n_indiv; j++)
			{
				update_indiv_j1(Exprs[j], Alpha[j], Mu, Omega_Ind,
						Omega_Logit, Gamma[j],
						Sig2, m, c_var, tau, pnum, pstart,
						n_position, n_peptide, *ind_pep, rng[th_id]);
			}
		}
#pragma omp for nowait
			for(pep = 0; pep < *n_peptide; pep++)
			{
				update_mu_pep(Mu, Exprs, Alpha,
						Gamma, Sig2[pep], kappa, *n_indiv, pep,
						rng[th_id]);
			}

		if(*update_u == 1)
		{
			update_u_pars(P, &a_0, &b_0, xA0, xB0, *n_position,
					psi_a, psi_b, rng[0], *nmax, *eps);
		}

		// check whether we need to update complete data
		if(*cen_num > 0)
		{
			update_data(D, *cen_num, cen_ind, cen_pep, Y, Exprs,
					Alpha, Mu, *n_peptide, Sig2, cen_pos);
		}

		if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ))
		{
			n = n + 1.0;
			update_prob_include(n_peptide, n_indiv, Gamma, ProbSum, mean_fitted,
					Alpha, Mu, n);
			if(*write == 1)
			{
				store_mcmc_output1(Mu, A, B, P, Sig2, D, Theta, Omega_Logit,
						Omega_Ind, kappa, alpha, beta,
						m, c_var, n_peptide, n_indiv, n_position, cen_num,
						lambda_a, lambda_b, a_0, b_0, *MRF, *ind_pep,
						AFILE, BFILE, PFILE, VARFILE, Sig2FILE, MUFILE, DFILE,
						THETAFILE, OFILE);
			}
		}

		// adaptation increment
		if(i < *n_burn & (i % 50 == 0))
		{
			if(accept_c > 15)
			{
				adptv = adptv + 1.0;
			}
			else if(accept_c < 5)
			{
				adptv = adptv/2.0;
			}
			accept_c = 0;
		}

		if(i < *n_burn & (i % 50 == 0))
		{
			if(accept_m > 15)
			{
				adptm = adptm + .1;
			}
			else if(accept_m < 5)
			{
				adptm = adptm - fmin(.1, adptm/2.0);
			}
			accept_m = 0;
		}
	}

	finalize_prob_include(n_iter, n_peptide, n_indiv, OutProbs, ProbSum);
	Rprintf("closing files\n");
	if(*write == 1)
	{
		fclose(AFILE);
		fclose(BFILE);
		fclose(PFILE);
		fclose(OFILE);
		fclose(VARFILE);
		fclose(Sig2FILE);
		fclose(MUFILE);
		if(*cen_num > 0)
		{
			fclose(DFILE);
		}
		if(*MRF == 1)
		{
			fclose(THETAFILE);
		}
	}

	free(A);
	free(B);
	free(P);
	free(Sig2);
	free(Mu);
	free(xAlpha);
	free(xA0);
	free(xB0);
	free(Omega_Ind);
	free(Omega_Logit);

	for(i = 0; i < *nP; i++)
	{
		RngStream_DeleteStream(rng[i]);
		//RngStream_IncreasedPrecis(rngs[i], 1);
	}


	if(*MRF == 1)
	{
		free(Theta);
	}

	if(*cen_num > 0)
	{
		free(D);
	}

	for(i = 0; i < *n_indiv; i++)
	{
		free(Exprs[i]);
		free(Alpha[i]);
		free(ProbSum[i]);
		free(Gamma[i]);
	}

	for(p = 0; p < *n_position; p++)
	{
		free(xA[p]);
		free(xB[p]);
	}

	free(xA);
	free(xB);
	free(Exprs);
	free(Alpha);
	free(ProbSum);
	free(Gamma);


	return;
}

void finalize_prob_include(int *n_iter, int *n_peptide, int *n_indiv, double *OutProbs, int **ProbSum)
{
	int p, i;
	for(i = 0; i < *n_indiv; i++)
	{
		for(p = 0; p < *n_peptide; p++)
		{
			OutProbs[*n_peptide*i + p] = ((double)(ProbSum[i][p]))/((double)(*n_iter));
		}
	}
	return;
}

void store_mcmc_output1(double *Mu, double *A, double *B, double *P, double *Sig2, double *D,
		double *Theta, double *Omega_Logit, int* Omega_Ind, double kappa, double alpha, double beta,
		double m, double c_var, int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		double lambda_a, double lambda_b, double a_0, double b_0, int MRF, int ind_pep,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *THETAFILE, FILE *OFILE)
{
	int p, i, k;

	for(p = 0; p < *n_position; p++)
	{
		fprintf(AFILE, "%.6lf \t", A[p]);
		fprintf(BFILE, "%.6lf \t", B[p]);
		fprintf(PFILE, "%.6lf \t", P[p]);
		if(ind_pep == 0)
		{
			fprintf(Sig2FILE, "%.6lf \t", Sig2[p]);
		}
		if(MRF == 1)
		{
			fprintf(THETAFILE, "%.6lf \t", Theta[p]);
		}
	}

	if(ind_pep == 1)
	{
		for(p = 0; p < *n_peptide; p++)
		{
			fprintf(Sig2FILE, "%.6lf \t", Sig2[p]);
		}
	}

	for(p = 0; p < *n_peptide; p++)
	{
		fprintf(OFILE, "%.6lf \t", (Omega_Ind[p] == 1) ? expit(Omega_Logit[p]):0.0);
		fprintf(MUFILE, "%.6lf \t", Mu[p]);
	}

	if(*cen_num > 0)
	{
		for(k = 0; k < *cen_num; k++)
		{
			fprintf(DFILE, "%.6lf \t", D[k]);

		}
		fprintf(DFILE, "\n");
	}

	if(MRF == 1)
	{
		fprintf(THETAFILE,"\n");
	}
	fprintf(AFILE, "\n");
	fprintf(BFILE, "\n");
	fprintf(PFILE, "\n");
	fprintf(OFILE, "\n");
	fprintf(Sig2FILE, "\n");
	fprintf(MUFILE, "\n");
	fprintf(VARFILE, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \n", m, c_var,
			lambda_a, lambda_b, kappa, alpha, beta, a_0, b_0);
	return;
}

void update_mu_pep(double *Mu, double **Exprs, double **Alpha,
		int **Gamma, double Sig2, double kappa, int n_indiv, int pep,
		RngStream rng)
{
	double S_dot = 0.0, V_dot = 0.0;
	double zi;
	int i;

	for(i = 0; i < n_indiv; i++)
	{
		zi = (Gamma[i][pep] == 1) ? (Exprs[i][pep] - Alpha[i][pep]) : Exprs[i][pep];
		S_dot += zi;
	}
	S_dot = S_dot/((double)(n_indiv) + kappa*Sig2); // con. mean
	V_dot = Sig2/((double)(n_indiv) + kappa*Sig2); // con. var

	Mu[pep] = RngStream_N01(rng)*sqrt(V_dot) + S_dot;
	return;
}

void update_kappa(double *Mu, double *kappa, int n_peptide, RngStream rng)
{
	int pep;
	double SS_mu = 0.0;

	for(pep = 0; pep < n_peptide; pep++)
	{
		SS_mu += gsl_pow_2(Mu[pep]);
	}

	SS_mu += 1.0;
	SS_mu = SS_mu/2.0;
	*kappa = RngStream_GA1(((double)n_peptide + 1)/2.0,rng)/SS_mu;

	return;
}

void update_data(double *D, int cen_num, int* cen_ind, int* cen_pep, double *Y, double **Exprs,
		double **Alpha, double *Mu, int n_peptide, double *Sig2, int *cen_pos)
{
	int i, p, k, pos;
	double tmp;
	for(k = 0; k < cen_num; k++)
	{
		// retrieve indices
		i = cen_ind[k];
		p = cen_pep[k];
		pos = cen_pos[k];

		tmp = truncNorm(Mu[p] + Alpha[i][p] - Y[i*n_peptide + p], Sig2[pos]);
		Exprs[i][p] = Y[i*n_peptide + p] + tmp;
		D[k] = tmp;
		//Rprintf("censored data %d updated \n", k);
	}
	return;
}

void update_prob_include(int *n_peptide, int *n_indiv, int **Gamma, int **ProbSum, double *mean_fitted,
		double **Alpha, double *Mu, double n)
{
	int p, i;

	for(i = 0; i < *n_indiv; i++)
	{
		for(p = 0; p < *n_peptide; p++)
		{
			ProbSum[i][p] += Gamma[i][p];
			mean_fitted[*n_peptide*i + p] = mean_fitted[*n_peptide*i + p]*((n-1.0)/n) + (Alpha[i][p] + Mu[p])/n;
		}
	}
	return;
}


void InitializeChain1(int *Omega_Ind, double *Omega_Logit,
		double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P, double *Theta,
		double *c_var, double *m, double *alpha, double *beta,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		int nmax, double **xA, double **xB,
		int MRF, int ind_pep)
{
	int c, p, i, g, pep;
	double u;
	for(p = 0; p < *n_position; p++)
	{
		A[p] = 1.0;
		B[p] = 2.0;
		if(MRF == 1)
		{
			Theta[p] = 0;
		}
		P[p] = 0.0;

		if(ind_pep == 0)
		{
			Sig2[p] = 1.0;
		}

		for(g = 0; g < nmax; g++)
			{
				xA[p][g] = 0.0;
				xB[p][g] = 0.0;
			}
			xA[p][0] = .1;
			xB[p][0] = .1;
			xA[p][1] = 	2.0;
			xB[p][1] =  2.0;
	}

	if(MRF == 1)
	{
		Theta[0] = -4.0;
		Theta[*n_position - 1] = -4.0;

		P[0] = expit(-4.0);
		P[*n_position - 1] = expit(-4.0);
	}

	*c_var = 1.0;
	*m = rnorm(0.0, 1.0);

	for(p = 0; p < *n_position; p++)
	{
		for(c = 0; c < pnum[p]; c++)
		{
			pep = pstart[p] + c;
			Mu[pep] = 0.0;
			if(ind_pep == 1)
			{
				Sig2[pep] = 1.0;
			}

			if(runif(0.0,1.0) <= expit(P[p]))
			{
				Omega_Ind[pep] = 0;
				Omega_Logit[pep] = 0.0;
			}

			else
			{
				Omega_Ind[pep] = 1;
				Omega_Logit[pep] = log(rgamma(A[p], 1.0)/rgamma(B[p], 1.0));
			}

			for(i = 0; i < *n_indiv; i++)
			{
				u = runif(0.0, 1.0);
				if(Omega_Ind[pep] == 1)
				{
					Gamma[i][pep] = (int)(u <= expit(Omega_Logit[pep]));
				}

				else
				{
					Gamma[i][pep] = 0;
				}

				if(Gamma[i][pep] == 0)
				{
					Alpha[i][pep] = 0.0;
				}
				else
				{
					Alpha[i][pep] = truncNorm(*m, *c_var);
				}
			}
		}
	}

	return;
}

void update_MRF(double *P, double *Theta, double *Omega,
		double *kappa, int *pnum, int *pstart, int *n_position, int *accept, RngStream rng)
{
	//Rprintf("entering MRF update \n");
	int S_p = 0;
	int p,c;
	double theta_prop, log_ratio, M, S = 0, shape;

	for(p = 0; p < (*n_position); p++)
	{
		for(c = 0; c < pnum[p]; c++)
		{
			if(Omega[pstart[p] + c] > 0.0)
			{
				S_p ++;
			}
		}

		theta_prop = Theta[p] + 2.0*(RngStream_RandU01(rng)*2.0 - 1.0)/sqrt(*kappa);

		if(p == 0)
		{
			M = Theta[1];
		}
		else if(p == (*n_position - 1))
		{
			M = Theta[*n_position - 2];
		}
		else
		{
			M = (Theta[p-1] + Theta[p+1])/2.0;
		}

		log_ratio = (*kappa/2.0)*(gsl_pow_2(Theta[p]) - gsl_pow_2(theta_prop)) +
				*kappa*M*(theta_prop - Theta[p]) +
				(double)(pnum[p] - S_p)*(theta_prop - Theta[p]) +
				(double)(pnum[p])*(log1p(exp(Theta[p])) - log1p(exp(theta_prop)));
		if(log(RngStream_RandU01(rng)) <= log_ratio)
		{
			Theta[p] = theta_prop;
			P[p] = expit(Theta[p]);
			accept[p] = accept[p] + 1;
		}
		S_p = 0;
	}

	for(p = 0; p < *n_position - 1; p++)
	{
		S += gsl_pow_2(Theta[p]) - Theta[p]*Theta[p+1];
	}
	S += gsl_pow_2(Theta[*n_position]);
	S += .5;
	shape = ((double)(*n_position) + 1.0)/2.0;

	*kappa = RngStream_GA1(shape, rng)/S;

	return;
}
void update_global_parms(double **Alpha, int **Gamma, double *m, double *c_var,
		double *A, double *B, double *Sig2, double m_0, double v_0, double alpha_0,
		double beta_0, double *alpha, double *beta, int lambda_prior,
		double *lambda_a, double *lambda_b, double r_a, double r_b,
		int n_position, int n_peptide, int n_indiv, int *accept_m,
		int *accept_c, RngStream rng, double adptm, double adptv, int ind_pep,
		int var_prior, double *xAlpha, int n_max, double eps)
{
	double n_alpha = 0.0;
	double sum_alpha = 0.0;
	double v_prime, d_prime;
	double SS_alpha_resid = 0.0;
	double m_prop, c_prop;
	double log_DF_ratio;

	if((var_prior != 0))
	{
		update_var_pars(Sig2, alpha, beta, alpha_0, beta_0,
				xAlpha, n_max, ind_pep, n_peptide, n_position,
				rng, eps);
	}

	int i, p;

	// find total number of active peptides
	for(i = 0; i < n_indiv; i++)
	{
		for(p = 0; p < n_peptide; p++)
		{
			n_alpha += (double)(Gamma[i][p]);
			sum_alpha += Alpha[i][p];
		}
	}

	//update m
	v_prime = (*c_var*v_0)/(*c_var + v_0*n_alpha);
	m_prop = *m + (RngStream_RandU01(rng)*2.0 - 1.0)*sqrt(v_prime)*adptm; //sqrt(v_prime)*rnorm(0.0,1.0);
	log_DF_ratio = pnorm(0.0, (*m), sqrt(*c_var), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, m_prop, sqrt(*c_var), 0, 1);
	log_DF_ratio = (double)(n_alpha)*log_DF_ratio;
	log_DF_ratio += -(.5/v_prime)*(gsl_pow_2(m_prop) - gsl_pow_2(*m));
	log_DF_ratio += ((*c_var*m_0 + v_0*sum_alpha)/(*c_var*v_0))*(m_prop - *m);

	if(log(RngStream_RandU01(rng)) <= log_DF_ratio)
	{
		*m = m_prop;
		*accept_m = *accept_m + 1;
	}

	// update c
	for(i = 0; i < n_indiv; i++)
	{
		for(p = 0; p < n_peptide; p++)
		{
			if(Gamma[i][p] == 1)
			{
				SS_alpha_resid += gsl_pow_2((*m) - Alpha[i][p]);
			}
		}
	}

	SS_alpha_resid += 1.0;

	// log-scale random walk proposal
	d_prime = n_alpha/2.0 + 0.5;
	c_prop = *c_var*exp(adptv*(RngStream_RandU01(rng)*1.0 - .5)/
			sqrt((double)(n_alpha)));

	log_DF_ratio = pnorm(0.0, (*m), sqrt(*c_var), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, (*m), sqrt(c_prop), 0, 1);
	log_DF_ratio = (double)(n_alpha)*log_DF_ratio;
	log_DF_ratio += .5*(1.0/(*c_var) - 1.0/c_prop)*SS_alpha_resid;
	log_DF_ratio += d_prime*(log(*c_var) - log(c_prop));

	// accept or reject step
	if(log(RngStream_RandU01(rng)) <= log_DF_ratio)
	{
		*c_var = c_prop;
		*accept_c = *accept_c + 1;
	}

	// update lambda_a, lambda_b
	if(lambda_prior == 1)
	{
		update_lambdas(lambda_a, lambda_b, r_a, r_b, A, B,
				n_position, rng);
	}

	return;
}


void update_var_pars(double *Sig2, double *alpha, double *beta,
		double alpha_0, double beta_0,  double *xAlpha,
		int n_max, int ind_pep, int n_peptide,
		int n_position, RngStream rng, double eps)
{
	double s_inv_sig = beta_0;
	double s_log_sig = 0.0;
	double N_p;
	double argvec[2];

	ARS_workspace workspace;

	int p;
	int num_x = 2;

	if(ind_pep != 0) //peptide specific variances
	{
		for(p = 0; p < n_peptide; p++)
		{
			s_inv_sig += 1.0/Sig2[p];
			s_log_sig += log(Sig2[p]);
		}
		N_p = (double)(n_peptide);
	}

	else // position specific variances
	{
		for(p = 0; p < n_position; p++)
			{
				s_inv_sig += 1.0/Sig2[p];
				s_log_sig += log(Sig2[p]);
			}
			N_p = (double)(n_position);
	}

	*beta = RngStream_GA1(*alpha*N_p, rng)/(s_inv_sig);
	argvec[0] = alpha_0 + s_log_sig - log(*beta)*N_p;
	argvec[1] = N_p;
	//Rprintf("updating alpha \n");
	*alpha = sample_conditional(xAlpha, &num_x, n_max, argvec,
			&workspace, rng, eps, lc_alpha, lcp_alpha);
	//Rprintf("beta = %lf, alpha = %lf \n", *beta, *alpha);
	return;
}

void update_lambdas(double *lambda_a, double *lambda_b, double r_a,
		double r_b, double *A, double *B, int n_position,
		RngStream rng)
{

	int p;
	double S_a = 0.0, S_b = 0.0;
	double r, s, rate, shape;

	for(p = 0; p < n_position; p++)
	{
		S_a += A[p];
		S_b += B[p];
	}
/*
	*lambda_a = RngStream_GA1((double)(n_position) + 1.0, rng)/(S_a + r_a);
	*lambda_b = RngStream_GA1((double)(n_position) + 1.0, rng)/(S_b + r_b);
*/
	s = *lambda_b;

	rate = S_a*s + s*r_a;
	shape = (double)(n_position + 1);
	r = RngStream_GA1(shape, rng)/rate;

	rate = S_b + S_a*r + r_a*r + r_b;
	shape = 2.0*(double)(n_position + 1);
	s = RngStream_GA1(shape, rng)/rate;

	*lambda_b = s;
	*lambda_a = r*s;

	return;
}

void update_indiv_j1(double *Exprs, double *Alpha, double *Mu, int* Omega_Ind,
		double *Omega_Logit, int *Gamma,
		double *Sig2, double m, double c_var,
		double tau_0, int *pnum, int *pstart,
		int *n_position, int *n_peptide, int ind_pep, RngStream rng)
{
	int p, c, pep;

	// update Alpha_icp and Gamma_icp, the individual fitted effects.
	int cur;
	double m_prime, v_prime, R;
	double var, omega, u;
	for(p = 0; p < *n_position; p++)
	{
		if(ind_pep == 0)
		{
			var = Sig2[p];
		}

		for(c = 0; c < pnum[p]; c++)
		{
			pep = pstart[p] + c;
			if(Omega_Ind[pep] == 1)
			{
				if(ind_pep == 1)
				{
					var = Sig2[pep];
				}
				v_prime = (c_var*var)/(var + c_var);
				m_prime = ((Exprs[pep] - Mu[pep])*c_var + m*(var))/(c_var + var);
				R = .5*log(v_prime/c_var);
				R = R + pnorm(0.0, m_prime, sqrt(v_prime), 0, 1);
				R = R - pnorm(0.0, m, sqrt(c_var), 0, 1);

				R = exp(gsl_pow_2(m_prime)/(2.0*v_prime)
						- gsl_pow_2(m)/(2.0*(c_var)) + R);

				omega = expit(Omega_Logit[pep]);
				u = RngStream_RandU01(rng);

				cur = (u <= ((R*omega)/(1.0 - omega + omega*R))) ? 1 : 0;

				if(cur == 1)
				{
					Gamma[pep] = 1;
					Alpha[pep] = truncNorm_parallel(m_prime, v_prime, rng);
				}
				else if(cur == 0)
				{
					Gamma[pep] = 0;
					Alpha[pep] = 0.0;
				}
				else
				{
					Rprintf("Something weird! Cur = %d, R = %lf \n", cur, R);
				}
			}
			else
			{
				Gamma[pep] = 0;
				Alpha[pep] = 0.0;
			}
		}
	}
	return;
}

int update_position_p1(double **Exprs, int* Omega_Ind, double *Omega_Logit,
		double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P,
		double a_0, double b_0, double lambda_a, double lambda_b,
		int p, int *n_indiv, int *pnum, int *pstart, RngStream rng, int nmax,
		double *xA, double *xB, ARS_workspace *workspace, double eps,
		int ind_pep, int MRF, double alpha, double beta)
{
	double s_log_w = 0.0, s_log_1minus_w = 0.0;
	double frac, R;

	int c, S_p = 0, i;
	int num_x = 2;
	int p_begin = pstart[p];
	double N_cp = 0.0;
	double argvec[3];

	// update A's and B's
	for(c = 0; c < pnum[p]; c++)
	{
		if(Omega_Ind[c] == 1)
		{
			s_log_w += log_from_logit(Omega_Logit[p_begin + c]);
			s_log_1minus_w += log1m_from_logit(Omega_Logit[p_begin + c]);
			S_p += 1;
		}
	}

	if(S_p == 0)
	{
		A[p] = RngStream_GA1(1.0, rng)/(lambda_a);
		B[p] = RngStream_GA1(1.0, rng)/(lambda_b);
	}

	else
	{
		R = (lambda_a - s_log_w);
		argvec[2] = R;
		argvec[0] = (double)S_p;
		argvec[1] = B[p];

		A[p] = sample_conditional(xA, &num_x, nmax, argvec,
				workspace, rng, eps, lc_AB, lcp_AB);

		if(A[p] == -1.0)
		{
			return(0);
		}

		num_x = 2;
		R = (lambda_b - s_log_1minus_w);
		argvec[2] = R;
		argvec[1] = A[p];

		B[p] = sample_conditional(xB, &num_x, nmax, argvec,
				workspace, rng, eps, lc_AB, lcp_AB);
		if(B[p] == -1.0)
		{
			return(0);
		}
	}
	//update Ps
	if(MRF == 0)
	{
		P[p] = RngStream_LogitBeta(a_0 + (double)(pnum[p] - S_p), b_0 + (double)(S_p), rng);
	}

	frac = lbeta(A[p], (double)(*n_indiv) + B[p])- lbeta(A[p], B[p]);
	frac = exp(frac);

	int pep;
	double u;
	u = expit(P[p]);
	// update Omegas
	for(c = 0; c < pnum[p]; c++)
	{
		pep = p_begin + c;
		for(i = 0; i < *n_indiv; i++)
		{
			N_cp += (double)Gamma[i][pep];
		}

		R = (double)(N_cp == 0.0)*u/(u + (1.0 - u)*frac);
		if(RngStream_RandU01(rng) <= R)
		{
			Omega_Ind[pep] = 0;
		}
		else
		{
			Omega_Ind[pep] = 1;
			Omega_Logit[pep] = RngStream_LogitBeta(A[p] + N_cp,
					B[p] + (double)(*n_indiv) - N_cp, rng);
		}
		N_cp = 0.0;
	}

	double SS_i = 0.0;
	//update variances
	if(ind_pep == 0)
	{
		for(i = 0; i < *n_indiv; i++)
		{
			for(c = 0; c < pnum[p]; c++)
			{
				pep = p_begin + c;
				SS_i += gsl_pow_2(Exprs[i][pep] - Mu[pep] - Alpha[i][pep]);
			}
		}
		SS_i = SS_i/2.0;
		SS_i += beta;
		Sig2[p] = SS_i/RngStream_GA1(((double)(*n_indiv*pnum[p]))/2.0 + alpha, rng);
	}

	else
	{
		for(c = 0; c < pnum[p]; c++)
		{
			pep = p_begin + c;
			for(i = 0; i < *n_indiv; i++)
			{
				SS_i += gsl_pow_2(Exprs[i][pep] - Mu[pep] - Alpha[i][pep]);
			}
			SS_i = SS_i/2.0;
			SS_i += beta;
			Sig2[pep] = SS_i/RngStream_GA1(((double)(*n_indiv))/2.0 + alpha, rng);
			SS_i = 0.0;
		}

	}
	return(1);
}

void update_u_pars(double *P, double *a_0, double *b_0,
		double *xA0, double *xB0, int n_position, double psi_a,
		double psi_b, RngStream rng, int nmax, double eps)
{
	int p;
	int num_x = 2;

	double s_log_u = psi_a;
	double s_log_1mu = psi_b;
	double argvec[3];

	ARS_workspace workspace;
	argvec[0] = (double)(n_position);
	argvec[1] = *b_0;

	for(p = 0; p < n_position; p++)
	{
		s_log_u -= log_from_logit(P[p]);
		s_log_1mu -= log1m_from_logit(P[p]);
	}

	argvec[2] = s_log_u;

	*a_0 = sample_conditional(xA0, &num_x, nmax, argvec,
			&workspace, rng, eps, lc_AB, lcp_AB);

	num_x = 2;

	argvec[1] = *a_0;
	argvec[2] = s_log_1mu;
	*b_0 = sample_conditional(xB0, &num_x, nmax, argvec,
			&workspace, rng, eps, lc_AB, lcp_AB);
	return;
}


double truncNorm_parallel(double mean, double sigmasqr, RngStream rng)
{
	double sigma;
	double left;
	double x;

	sigma = sqrt(sigmasqr);
	left = -1.0*mean/sigma;

	if(left < 0) // don't bother with exponential proposal
	{
		while(1)
		{
			x = RngStream_N01(rng);
			if(x > left)
			{
				x = x*sigma + mean;
				break;
			}
		}
	}

	else // we use the exponential proposal.
	{
		double a_star, z, pz;
		a_star = (left + sqrt(left*left + 4.0))/2.0;
		while(1)
		{
			z = RngStream_GA1(1.0, rng)*a_star + left;
			pz = exp(-.5*gsl_pow_2(z - a_star));
			if(RngStream_RandU01(rng) < pz) // accept z
			{
				x = sigma*z + mean;
				break;
			}
		}
	}
	return(x);
}

double truncNorm(double mean, double sigmasqr)
{
	double sigma;
	double left;
	double x;

	sigma = sqrt(sigmasqr);
	left = -1.0*mean/sigma;

	if(left < 0) // don't bother with exponential proposal
	{
		while(1)
		{
			x = rnorm(0.0, 1.0);
			if(x > left)
			{
				x = x*sigma + mean;
				break;
			}
		}
	}

	else // we use the exponential proposal.
	{
		double a_star, z, pz;
		a_star = (left + sqrt(left*left + 4.0))/2.0;
		while(1)
		{
			z = rexp(1.0)*a_star + left;
			pz = exp(-.5*gsl_pow_2(z - a_star));
			if(runif(0.0, 1.0) < pz) // accept z
			{
				x = sigma*z + mean;
				break;
			}
		}
	}
	return(x);
}

inline double expit(double x)
{
	double p;
	p = 1.0/(1.0 + exp(-x));
	return(p);
}

inline double log_from_logit(double x)
{
	if(x > 0.0)
	{
		return(-log1p(exp(-x)));
	}
	else
	{
		return(x - log1p(exp(x)));
	}
}

inline double log1m_from_logit(double x)
{
	if(x > 0.0)
	{
		return(-x - log1p(exp(-x)));
	}
	else
	{
		return(-log1p(exp(x)));
	}
}

void tnorm_test(double* x, int *n, double *m, double *sigmasqr)
{
	int i;
	for(i = 0; i < *n; i++)
	{
		GetRNGstate();
		x[i] = truncNorm(*m, *sigmasqr);
		PutRNGstate();
	}
	return;
}

double lc_AB(double x, double *argvec)
{
	double out;
	out = -1.0*x*argvec[2] + argvec[0]*lgamma(x + argvec[1]) - argvec[0]*lgamma(x);
	return(out);
}

double lcp_AB(double x, double *argvec)
{
	double out;
	out = -1.0*argvec[2] + argvec[0]*gsl_sf_psi(x + argvec[1]) - argvec[0]*gsl_sf_psi(x);
	return(out);
}

double lc_alpha(double x, double *argvec)
{
	double out;
	out = -x*argvec[0] - argvec[1]*lgamma(x);
	return(out);
}

double lcp_alpha(double x, double *argvec)
{
	double out;
	out = -argvec[0] - argvec[1]*gsl_sf_psi(x);
	return(out);
}

/********************************************************8 */
