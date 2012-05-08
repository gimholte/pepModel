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
double lc_AB(double x, double *argvec, int *arglen);

double lcp_AB(double x, double *argvec, int *arglen);

double lc_alpha(double x, double *argvec, int *arglen);

double lcp_alpha(double x, double *argvec, int *arglen);

inline double log1m_from_logit(double x);

inline double log_from_logit(double x);

inline double expit(double x);

inline double m1expit(double x);

void update_u_pars(double *P, double *a_0, double *b_0,
		double *xA0, double *xB0, int n_position, double psi_a,
		double psi_b, RngStream rng, int nmax, double eps);

void update_var_pars(double *Sig2, double *alpha, double *beta,
		double alpha_0, double beta_0,  double *xAlpha,
		int n_max, int ind_pep, int n_peptide,
		int n_position, RngStream rng, double eps);

void update_indiv_j1(double *Exprs, double *Alpha, double *mu_j,
		int *Gamma, double *Sig2,
		double tau_0, int *pnum, int *pstart,
		int *n_position, int *n_peptide, int ind_pep, RngStream rng);

int update_position_p1(double **Exprs, int* Omega_Ind, double *Omega_Logit,
		double *Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P,
		double a_0, double b_0, double lambda_a, double lambda_b,
		int p, int *n_indiv, int *pnum, int *pstart, RngStream rng, int nmax,
		double *xA, double *xB, ARS_workspace *workspace, double eps,
		int ind_pep, int MRF, double alpha, double beta);

void update_effect_means(double *M, double *Alpha, double m_0, double v_0,
		double v, double adptm,
		int p, int pnum, int pbegin, int n_position, int *accept_m,
		RngStream rng);

void InitializeChain1(int *Omega_Ind, double *Omega_Logit, double *Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P, double *Theta,
		double *M,
		double *c_var, double m_0, double v_0, double *alpha, double *beta,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		int nmax, double **xA, double **xB, double **xEffect,
		int MRF, int ind_pep);

void update_global_parms(double *Alpha, int **Gamma, double *M, double *c_var,
		double *A, double *B, double *Sig2, double m_0, double v_0, double alpha_0,
		double beta_0, double *alpha, double *beta, int lambda_prior,
		double *lambda_a, double *lambda_b, double r_a, double r_b,
		int n_position, int n_peptide, int *pnum, int *pstart,
		int *accept_c, RngStream rng, double adptv, int ind_pep,
		int var_prior, double *xAlpha, int n_max, double eps);

void update_peptide_p(double **Exprs, double *Alpha, int **Gamma, int *Omega_Ind,
		double *Omega_Logit, double *Mu, double Sig2, double m, double v,
		int pep, int n_indiv, int pos, double *xEffect,
		double *argvec, int nmax, ARS_workspace *workspace, double eps, RngStream rng);

void tnorm_test(double* x, int *n, double *m, double *sigmasqr);

void update_lambdas(double *lambda_a, double *lambda_b, double r_a,
		double r_b, double *A, double *B, int n_position,
		RngStream rng);

void update_MRF(double *P, double *Theta, double *Omega,
		double *kappa, int *pnum, int *pstart, int *n_position,
		int *accept, RngStream rng);

void update_prob_include(int *n_peptide, int *n_indiv, int **Gamma, int **ProbSum, double *mean_fitted,
		double *Alpha, double *Mu, double n);

double truncNorm_parallel(double mean, double sigmasqr, RngStream rng);

double truncNorm(double mean, double sigmasqr);

void update_data(double *D, int cen_num, int* cen_ind, int* cen_pep, double *Y, double **Exprs,
		int **Gamma, double *Alpha, double *Mu, int n_peptide, double *Sig2, int *cen_pos);

void store_mcmc_output1(double *M, double *Alpha, double *Mu, double *A, double *B, double *P, double *Sig2, double *D,
		double *Theta, double *Omega_Logit, int* Omega_Ind, double kappa, double alpha, double beta,
		double c_var, int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		double lambda_a, double lambda_b, double a_0, double b_0, int MRF, int ind_pep,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *THETAFILE, FILE *OFILE, FILE *ALPHAFILE, FILE *MFILE);

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
	double *Alpha;
	double **xA, **xB;
	double *Omega_Logit, *Theta;
	int *Omega_Ind, *pos_ind_by_pep;
	int *accept_m;
	double *adptm;
	int **ProbSum;
	int **Gamma;
	double *Mu;
	double *A, *B, *P, *Sig2, *D, *M;
	double *xAlpha, *xA0, *xB0;

	int accept_c = 0;

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
	double **xEffect;
	double **Alpha_argvec;


	Exprs = (double**) malloc(*n_indiv*sizeof(double*));
	Alpha = (double*) malloc(*n_peptide*sizeof(double));
	xEffect = (double**) malloc(*n_peptide*sizeof(double*));
	Alpha_argvec = (double**) malloc(*nP * sizeof(double*));
	accept_m = (int*) malloc(*n_position*sizeof(int));
	adptm = (double*) malloc(*n_position*sizeof(double));

	for(i = 0; i < *nP; i++)
	{
		Alpha_argvec[i] = (double*) malloc((*n_indiv + 5)*sizeof(double));
	}

	for(p = 0; p < *n_peptide; p++)
	{
		xEffect[p] = (double*) malloc(*nmax*sizeof(double));
	}

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
	pos_ind_by_pep = (int*) malloc(*n_peptide*sizeof(int));
	Omega_Logit = (double*) malloc(*n_peptide*sizeof(double));

	Gamma = (int**) malloc(*n_indiv*sizeof(int*));
	ProbSum = (int**) malloc(*n_indiv*sizeof(int*));

	int c;
	k = 0;
	for(p = 0; p < *n_position; p++)
	{
		xA[p] = (double*) malloc(*nmax*sizeof(double));
		xB[p] = (double*) malloc(*nmax*sizeof(double));
		accept_m[p] = 0;
		adptm[p] = 5.0;

		for(c = 0; c < pnum[p]; c++)
		{
			pos_ind_by_pep[pstart[p] + c] = k;
		}
		k++;
	}

	for(i = 0; i < *n_indiv; i++)
	{
		Exprs[i] = (double*) malloc(*n_peptide*sizeof(double));
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
	Mu = (double*) malloc(*n_indiv*sizeof(double));
	M = (double*) malloc(*n_position*sizeof(double));

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

	double c_var = 1.0;
	double alpha = alpha_0;
	double beta = beta_0;
	double adptv = 1.0;
	double sigma;

	FILE *AFILE, *BFILE, *PFILE, *VARFILE, *Sig2FILE, *MUFILE, *DFILE, *THETAFILE, *OFILE;
	FILE *ALPHAFILE, *MFILE;
	if(*write == 1)
	{
		MFILE = fopen("mfile.txt", "w");
		AFILE = fopen("afile.txt", "w");
		BFILE = fopen("bfile.txt", "w");
		PFILE = fopen("pfile.txt", "w");
		OFILE = fopen("ofile.txt", "w");
		ALPHAFILE = fopen("alphafile.txt", "w");
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

	int th_id, pos;
	RngStream rng[*nP];
	for(i = 0; i < *nP; i++)
	{
		rng[i] = RngStream_CreateStream("");
		//RngStream_IncreasedPrecis(rngs[i], 1);
	}
	GetRNGstate();
	InitializeChain1(Omega_Ind, Omega_Logit, Alpha,
			Gamma, Sig2, Mu, A, B, P, Theta, M,
			&c_var, m_0, v_0, &alpha, &beta,
			n_position, pstart, pnum, n_peptide, n_indiv,
			*nmax, xA, xB, xEffect,
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
				Rprintf("c acceptance prob = %f\n", (double)(accept_c)/(double)(i - *n_burn));
			}
		}

		R_CheckUserInterrupt();

		update_global_parms(Alpha, Gamma, M, &c_var, A, B, Sig2,
				 m_0, v_0, alpha_0, beta_0, &alpha, &beta,
				*lambda_prior, &lambda_a, &lambda_b, r_a, r_b,
				*n_position, *n_peptide, pnum, pstart,
				&accept_c, rng[0], adptv, *ind_pep, *var_prior,
				xAlpha, *nmax, *eps);


#pragma omp parallel private(th_id, workspace, pos, sigma) num_threads(*nP)
		{
			th_id = omp_get_thread_num();
#pragma omp for nowait
			for(p = 0; p < *n_peptide; p++)
			{
				pos = pos_ind_by_pep[p];
				sigma = (*ind_pep == 1) ? Sig2[p] : Sig2[pos];

				update_peptide_p(Exprs, Alpha, Gamma, Omega_Ind, Omega_Logit,
						Mu, sigma, M[pos], c_var, p, *n_indiv, pos, xEffect[p],
						Alpha_argvec[th_id], *nmax, &workspace, *eps, rng[th_id]);
			}

#pragma omp for nowait
			for(p = 0; p < *n_position; p++)
			{
				update_position_p1(Exprs, Omega_Ind, Omega_Logit, Alpha, Gamma,
						Sig2, Mu, A, B, P,
						a_0, b_0, lambda_a, lambda_b,
						p, n_indiv, pnum, pstart, rng[th_id], *nmax,
						xA[p], xB[p], &workspace, *eps, *ind_pep, *MRF,
						alpha, beta);

				update_effect_means(M, Alpha, m_0, v_0,
						c_var, adptm[p],
						p, pnum[p], pstart[p], *n_position, &(accept_m[p]),
						rng[th_id]);
			}

#pragma omp for nowait
			for(j = 0; j < *n_indiv; j++)
			{
				update_indiv_j1(Exprs[j], Alpha, &Mu[j], Gamma[j],
						Sig2, tau, pnum, pstart,
						n_position, n_peptide, *ind_pep, rng[th_id]);
			}

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
					Gamma, Alpha, Mu, *n_peptide, Sig2, cen_pos);
		}

		if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ))
		{
			n = n + 1.0;
			update_prob_include(n_peptide, n_indiv, Gamma, ProbSum, mean_fitted,
					Alpha, Mu, n);
			if(*write == 1)
			{
				store_mcmc_output1(M, Alpha, Mu, A, B, P, Sig2, D, Theta, Omega_Logit,
						Omega_Ind, kappa, alpha, beta,
						c_var, n_peptide, n_indiv, n_position, cen_num,
						lambda_a, lambda_b, a_0, b_0, *MRF, *ind_pep,
						AFILE, BFILE, PFILE, VARFILE, Sig2FILE, MUFILE, DFILE,
						THETAFILE, OFILE, ALPHAFILE, MFILE);
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
			for(p = 0; p < *n_position; p++)
			{
				if(accept_m[p] > 15)
				{
					adptm[p] = adptm[p] + .1;
				}
				else if(accept_m[p] < 5)
				{
					adptm[p] = adptm[p] - fmin(.1, adptm[p]/2.0);
				}
				accept_m[p] = 0;
			}
		}
	}

	finalize_prob_include(n_iter, n_peptide, n_indiv, OutProbs, ProbSum);
	Rprintf("closing files\n");
	if(*write == 1)
	{
		fclose(MFILE);
		fclose(AFILE);
		fclose(BFILE);
		fclose(PFILE);
		fclose(OFILE);
		fclose(VARFILE);
		fclose(ALPHAFILE);
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
	free(M);
	free(B);
	free(P);
	free(Sig2);
	free(Mu);
	free(xAlpha);
	free(xA0);
	free(xB0);
	free(Omega_Ind);
	free(pos_ind_by_pep);
	free(Omega_Logit);
	free(accept_m);
	free(adptm);

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
		free(ProbSum[i]);
		free(Gamma[i]);
	}

	for(p = 0; p < *n_position; p++)
	{
		free(xA[p]);
		free(xB[p]);
	}
	for(p = 0; p < *n_peptide; p++)
	{
		free(xEffect[p]);
	}
	for(i = 0; i < *nP; i++)
	{
		free(Alpha_argvec[i]);
	}
	free(xEffect);
	free(Alpha_argvec);
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

void store_mcmc_output1(double *M, double *Alpha, double *Mu, double *A, double *B, double *P, double *Sig2, double *D,
		double *Theta, double *Omega_Logit, int* Omega_Ind, double kappa, double alpha, double beta,
		double c_var, int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		double lambda_a, double lambda_b, double a_0, double b_0, int MRF, int ind_pep,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *THETAFILE, FILE *OFILE, FILE *ALPHAFILE, FILE *MFILE)
{
	int p, i, k;

	for(p = 0; p < *n_position; p++)
	{
		fprintf(AFILE, "%.6lf \t", A[p]);
		fprintf(MFILE, "%.6lf \t", M[p]);
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
		fprintf(ALPHAFILE, "%.6lf \t", Alpha[p]);
	}

	for(i = 0; i < *n_indiv; i++)
	{
		fprintf(MUFILE, "%.6lf \t", Mu[i]);
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
	fprintf(ALPHAFILE, "\n");
	fprintf(MFILE, "\n");
	fprintf(AFILE, "\n");
	fprintf(BFILE, "\n");
	fprintf(PFILE, "\n");
	fprintf(OFILE, "\n");
	fprintf(Sig2FILE, "\n");
	fprintf(MUFILE, "\n");
	fprintf(VARFILE, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \n", c_var,
			lambda_a, lambda_b, kappa, alpha, beta, a_0, b_0);
	return;
}

void update_data(double *D, int cen_num, int* cen_ind, int* cen_pep, double *Y, double **Exprs,
		int **Gamma, double *Alpha, double *Mu, int n_peptide, double *Sig2, int *cen_pos)
{
	int i, p, k, pos;
	double tmp;
	for(k = 0; k < cen_num; k++)
	{
		// retrieve indices
		i = cen_ind[k];
		p = cen_pep[k];
		pos = cen_pos[k];

		tmp = truncNorm(Mu[i] + Alpha[p]*(double)Gamma[i][p] - Y[i*n_peptide + p], Sig2[pos]);
		Exprs[i][p] = Y[i*n_peptide + p] + tmp;
		D[k] = tmp;
		//Rprintf("censored data %d updated \n", k);
	}
	return;
}

void update_prob_include(int *n_peptide, int *n_indiv, int **Gamma, int **ProbSum, double *mean_fitted,
		double *Alpha, double *Mu, double n)
{
	int p, i;

	for(i = 0; i < *n_indiv; i++)
	{
		for(p = 0; p < *n_peptide; p++)
		{
			ProbSum[i][p] += Gamma[i][p];
			mean_fitted[*n_peptide*i + p] = mean_fitted[*n_peptide*i + p]*((n-1.0)/n) + (Alpha[p] + Mu[i])/n;
		}
	}
	return;
}


void InitializeChain1(int *Omega_Ind, double *Omega_Logit, double *Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P, double *Theta,
		double *M,
		double *c_var, double m_0, double v_0, double *alpha, double *beta,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		int nmax, double **xA, double **xB, double **xEffect,
		int MRF, int ind_pep)
{
	int c, p, i, g, pep;
	double u;
	for(p = 0; p < *n_position; p++)
	{
		A[p] = 1.0;
		B[p] = 2.0;
		M[p] = rnorm(m_0, sqrt(v_0));
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

	for(i = 0; i < *n_indiv; i++)
	{
		Mu[i] = 0.0;
	}

	for(p = 0; p < *n_position; p++)
	{
		for(c = 0; c < pnum[p]; c++)
		{
			pep = pstart[p] + c;
			xEffect[pep][0] = .5;
			xEffect[pep][1] = 2.5;
			Alpha[pep] = truncNorm(M[p], *c_var);
			if(ind_pep == 1)
			{
				Sig2[pep] = .6324555;
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

void update_effect_means(double *M, double *Alpha, double m_0, double v_0,
		double v, double adptm,
		int p, int pnum, int pbegin, int n_position, int *accept_m,
		RngStream rng)
{
	int pep = 0, c = 0;
	double sum_alpha = 0.0, c_p = (double)pnum;
	double m_prop, v_prime;
	double log_DF_ratio = 0.0;
	double m_cur = M[p];

	for(c = 0; c < pnum; c++)
	{
		pep = pbegin + c;
		sum_alpha += Alpha[pep];
	}

	v_prime = v/(c_p + v/v_0);
	m_prop = RngStream_N01(rng)*sqrt(v_prime)*adptm + m_cur;

	log_DF_ratio = pnorm(0.0, m_cur, sqrt(v), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, m_prop, sqrt(v), 0, 1);
	log_DF_ratio = c_p*log_DF_ratio;
	log_DF_ratio += .5*(m_cur*m_cur - m_prop*m_prop)/v_prime;
	log_DF_ratio += (sum_alpha/v + m_0/v_0)*(m_prop - m_cur);

	if(log(RngStream_RandU01(rng)) <= log_DF_ratio)
	{
		M[p] = m_prop;
		*accept_m = *accept_m + 1;
	}

	return;
}

void update_global_parms(double *Alpha, int **Gamma, double *M, double *c_var,
		double *A, double *B, double *Sig2, double m_0, double v_0, double alpha_0,
		double beta_0, double *alpha, double *beta, int lambda_prior,
		double *lambda_a, double *lambda_b, double r_a, double r_b,
		int n_position, int n_peptide, int *pnum, int *pstart,
		int *accept_c, RngStream rng, double adptv, int ind_pep,
		int var_prior, double *xAlpha, int n_max, double eps)
{
	double n_alpha = 0.0;
	double sum_alpha = 0.0;
	double d_prime, c_p, n_p = 0.0;
	double SS_alpha_resid = 0.0;
	double c_prop, c_cur = *c_var;
	double log_DF_cur = 0.0, log_DF_prop = 0.0, log_ratio;

	if((var_prior != 0))
	{
		update_var_pars(Sig2, alpha, beta, alpha_0, beta_0,
				xAlpha, n_max, ind_pep, n_peptide, n_position,
				rng, eps);
	}

	int i, p, pep, c;

	c_prop = c_cur*exp(adptv*(RngStream_RandU01(rng)*1.0 - .5));

	// update c
	for(p = 0; p < n_position; p++)
	{
		for(c = 0; c < pnum[p]; c++)
		{
			pep = pstart[p] + c;
			SS_alpha_resid += gsl_pow_2(Alpha[pep] - M[p]);
		}
		c_p = (double)(pnum[p]);
		log_DF_cur += c_p*pnorm(0.0, M[p], sqrt(c_cur), 0, 1);
		log_DF_prop += c_p*pnorm(0.0, M[p], sqrt(c_prop), 0, 1);

	}
	n_p = (double)n_peptide;
	SS_alpha_resid += 1.0;

	// log-scale random walk proposal
	d_prime = n_p/2.0 + 0.5;

	log_ratio = log_DF_cur - log_DF_prop;
	log_ratio += .5*(1.0/(c_cur) - 1.0/c_prop)*SS_alpha_resid;
	log_ratio += d_prime*(log(c_cur) - log(c_prop));

	// accept or reject step
	if(log(RngStream_RandU01(rng)) <= log_ratio)
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
	int arglen = 2;
	//Rprintf("updating alpha \n");
	*alpha = sample_conditional(xAlpha, &num_x, n_max, argvec, &arglen,
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

void update_indiv_j1(double *Exprs, double *Alpha, double *mu_j,
		int *Gamma, double *Sig2,
		double tau_0, int *pnum, int *pstart,
		int *n_position, int *n_peptide, int ind_pep, RngStream rng)
{
	// update Mu_j
	double S_dot = 0.0, V_dot = 1.0/tau_0, zi;
	int p, c, pep;

	// peptide variances common across position
	if(ind_pep == 0)
	{
		for(p = 0; p < *n_position; p++)
		{
			for(c = 0; c < pnum[p]; c++)
			{
				pep = pstart[p] + c;
				zi = (Gamma[pep] == 1) ? (Exprs[pep] - Alpha[pep]) : Exprs[pep];
				zi = zi/Sig2[p];
				S_dot += zi;
			}
			V_dot += pnum[p]/Sig2[p];
		}
		*mu_j = RngStream_N01(rng)/sqrt(V_dot) + S_dot/V_dot;
	}
	// peptide specific variances
	else
	{
		for(p = 0; p < *n_peptide; p++)
		{
			zi = (Gamma[p] == 1) ? (Exprs[p] - Alpha[p]) : Exprs[p];
			zi = zi/Sig2[p];

			S_dot += zi;
			V_dot += 1.0/Sig2[p];
		}
		*mu_j = RngStream_N01(rng)/sqrt(V_dot) + S_dot/V_dot;
	}
	return;
}

int update_position_p1(double **Exprs, int* Omega_Ind, double *Omega_Logit,
		double *Alpha, int **Gamma,
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

		int arglen = 3;

		A[p] = sample_conditional(xA, &num_x, nmax, argvec, &arglen,
				workspace, rng, eps, lc_AB, lcp_AB);

		if(A[p] == -1.0)
		{
			return(0);
		}

		num_x = 2;
		R = (lambda_b - s_log_1minus_w);
		argvec[2] = R;
		argvec[1] = A[p];

		B[p] = sample_conditional(xB, &num_x, nmax, argvec, &arglen,
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

	double SS_i = 0.0, zi;
	//update variances

	if(ind_pep == 0)
	{
		for(i = 0; i < *n_indiv; i++)
		{
			for(c = 0; c < pnum[p]; c++)
			{
				pep = p_begin + c;
				zi = Exprs[i][pep] - Mu[i];
				zi = (Gamma[i][pep] == 1) ? zi - Alpha[pep] : zi;
				SS_i += gsl_pow_2(zi);
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
				zi = Exprs[i][pep] - Mu[i];
				zi = (Gamma[i][pep] == 1) ? zi - Alpha[pep] : zi;
				SS_i += gsl_pow_2(zi);
			}
			SS_i = SS_i/2.0;
			SS_i += beta;
			Sig2[pep] = SS_i/RngStream_GA1(((double)(*n_indiv))/2.0 + alpha, rng);
			SS_i = 0.0;
		}
	}
	return(1);
}

void update_peptide_p(double **Exprs, double *Alpha, int **Gamma, int *Omega_Ind,
		double *Omega_Logit, double *Mu, double Sig2, double m, double v,
		int pep, int n_indiv, int pos, double *xEffect,
		double *argvec, int nmax, ARS_workspace *workspace, double eps, RngStream rng)
{
	int i, cur;
	double zi, N_cp = 0.0, M_cp = 0.0;
	double alpha = Alpha[pep];
	double alpha_sqd_sig2 = .5*alpha*alpha/Sig2;
	double m_prime, v_prime;
	double one_prob, denom, num;
	const double log_w = log_from_logit(Omega_Logit[pep]);
	const double w = expit(Omega_Logit[pep]);
	const double m1w = m1expit(Omega_Logit[pep]);

	if(Omega_Ind[pep] == 1)
	{
		for(i = 0; i < n_indiv; i++)
		{
			zi = (Exprs[i][pep] - Mu[i]);
			denom = log1p(w*expm1(zi*alpha/Sig2 - alpha_sqd_sig2));
			num = zi*alpha/Sig2 - alpha_sqd_sig2 + log_w;
			one_prob = num - denom;

			cur = (log(RngStream_RandU01(rng)) <= one_prob) ? 1 : 0;
			Gamma[i][pep] = cur;
			if(cur == 1)
			{
				N_cp += 1.0;
				M_cp += zi;
			}
		}
	}

	else
	{
		for(i = 0; i < n_indiv; i++)
		{
			Gamma[i][pep] = 0;
		}
	}

	m_prime = (M_cp + Sig2*m/v)/(N_cp + Sig2/v);
	v_prime = (Sig2) / (N_cp + Sig2/v);

	alpha = truncNorm_parallel(m_prime, v_prime, rng);
	Alpha[pep] = alpha;
	return;
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
	int arglen = 3;

	*a_0 = sample_conditional(xA0, &num_x, nmax, argvec, &arglen,
			&workspace, rng, eps, lc_AB, lcp_AB);

	num_x = 2;

	argvec[1] = *a_0;
	argvec[2] = s_log_1mu;
	*b_0 = sample_conditional(xB0, &num_x, nmax, argvec, &arglen,
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

inline double m1expit(double x)
{
	double p;
	p = 1.0/(1.0 + exp(x));
	return(p);
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

double lc_AB(double x, double *argvec, int *arglen)
{
	double out;
	out = -1.0*x*argvec[2] + argvec[0]*lgamma(x + argvec[1]) - argvec[0]*lgamma(x);
	return(out);
}

double lcp_AB(double x, double *argvec, int *arglen)
{
	double out;
	out = -1.0*argvec[2] + argvec[0]*gsl_sf_psi(x + argvec[1]) - argvec[0]*gsl_sf_psi(x);
	return(out);
}

double lc_alpha(double x, double *argvec, int *arglen)
{
	double out;
	out = -x*argvec[0] - argvec[1]*lgamma(x);
	return(out);
}

double lcp_alpha(double x, double *argvec, int *arglen)
{
	double out;
	out = -argvec[0] - argvec[1]*gsl_sf_psi(x);
	return(out);
}

/*
 * 			N_minus_i -= Gamma[i][pep];
			//compute transition probability
			if(N_minus_i == 0)
			{
				frac = gsl_sf_lnbeta((double)n_indiv + B, A) - gsl_sf_lnbeta(A, B);
				frac = exp(frac);
				theta_i_one = (1 - U)*frac / (U + (1 - U)*frac);
			}
			else
			{
				theta_i_one = ((double)N_minus_i + A)/((double)n_indiv + A + B - 1.0);
			}

			//compute update
			if(Gamma[i][pep] == 1)
			{
				flip = (RngStream_RandU01(rng) <= (1.0 - theta_i_one)) ? 1 : 0;
			}
			else
			{
				flip = (RngStream_RandU01(rng) <= theta_i_one) ? 1 : 0;
			}

			//compute acceptance probability if we flip

			if(flip == 1)
			{
				S_0 = y_minus_mu_gamma - (Exprs[i][pep] - Mu[i])*((double)Gamma[i][pep]);
				S_1 = y_minus_mu_gamma + (Exprs[i][pep] - Mu[i])*(1.0 - (double)Gamma[i][pep]);

				//Rprintf("s0 = %.4lf, s1 = %.4lf\n", S_0, S_1);
				v_prime_1 = (v*Sig2)/(((double)N_minus_i + 1.0)*v + Sig2);
				v_prime_0 = (v*Sig2)/(((double)N_minus_i)*v + Sig2);

				m_prime_1 = (v*S_1 + m*Sig2)/(((double)N_minus_i + 1.0)*v + Sig2);
				m_prime_0 = (v*S_0 + m*Sig2)/(((double)N_minus_i)*v + Sig2);

				L_0 = -.5*gsl_pow_2(m_prime_0)/v_prime_0 + log(v_prime_0) +
						pnorm(0.0, m_prime_0, sqrt(v_prime_0), 0, 1);
				L_1 = -.5*gsl_pow_2(m_prime_1)/v_prime_1 + log(v_prime_1) +
						pnorm(0.0, m_prime_1, sqrt(v_prime_1), 0, 1);

				// log acceptance ratio;
				log_alpha = (Gamma[i][pep] == 1) ? (L_0 - L_1):(L_1 - L_0);

				// check acceptance of the flip move:
				if(log_alpha <= log(RngStream_RandU01(rng)))
				{
					//flip!
					Gamma[i][pep] = (Gamma[i][pep] == 1) ? 0 : 1;
					//must update y_minus_mu_gamma if we flip
					// if we flipped from on to off, we take away the expression diff
					y_minus_mu_gamma += (Gamma[i][pep] == 1) ?
							(Exprs[i][pep] - Mu[i]) : -1.0*(Exprs[i][pep] - Mu[i]);
				}
			}
			// add back our Gamma
			N_minus_i += Gamma[i][pep];
		}

		//then we update alpha!
		if(N_minus_i > 0)
		{
			m_prime = (y_minus_mu_gamma*v + Sig2*m)/(v*(double)N_minus_i + Sig2);
			v_prime = (Sig2*v)/(v*(double)N_minus_i + Sig2);

			alpha_new = truncNorm_parallel(m_prime, v_prime, rng);
			Alpha[pep] = alpha_new;
		}

		else Alpha[pep] = 0.0;
	}
	else Alpha[pep] = 0.0;
 * *******************************************************8 */
