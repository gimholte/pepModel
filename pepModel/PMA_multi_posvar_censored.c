/*
 * PMA_multi_posvar_censored.c
 *
 *  Created on: Feb 14, 2012
 *      Author: Gregory Imholte
 */

#include <omp.h>
#include <math.h>
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

#define ZZERO 2.0e-308

//
void update_indiv_j1(double *Exprs, double *Alpha, double *mu_j, double *Omega, int *Gamma,
		double *Sig2, double m, double c_var, double tau_0, int *pnum, int *pstart,
		int *n_position, int *n_peptide, int ind_pep, RngStream rng);

int update_position_p1(double **Exprs, double *Omega, double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P,
		double a_0, double b_0, double lambda_a, double lambda_b,
		int p, int *n_indiv, int *pnum, int *pstart, RngStream rng, int nmax,
		double *xA, double *xB, ARS_workspace *workspace, double eps,
		int ind_pep, int MRF);

void InitializeChain1(double *Omega, double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P, double *Theta,
		double *c_var, double *m, double *alpha, double *beta,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		int nmax, double **xA, double **xB,
		int MRF, int ind_pep);

void update_global_parms(double **Alpha, int **Gamma, double *m, double *c_var,
		double *A, double *B, double m_0, double v_0, int lambda_prior,
		double *lambda_a, double *lambda_b, double r_a, double r_b,
		int n_position, int n_peptide, int n_indiv, int *accept_m,
		int *accept_c, RngStream rng, double adptm, double adptv);

double expit(double x);

void tnorm_test(double* x, int *n, double *m, double *sigmasqr);

void update_lambdas(double *lambda_a, double *lambda_b, double r_a,
		double r_b, double *A, double *B, int n_position,
		RngStream rng);

void update_MRF(double *P, double *Theta, double *Omega,
		double *kappa, int *pnum, int *pstart, int *n_position, int *accept, RngStream rng);

void update_prob_include(int *n_peptide, int *n_indiv, int **Gamma, int **ProbSum, double *mean_fitted,
		double **Alpha, double *Mu, double n);

double truncNorm_parallel(double mean, double sigmasqr, RngStream rng);

double truncNorm(double mean, double sigmasqr);

void update_data(double *D, int cen_num, int* cen_ind, int* cen_pep, double *Y, double **Exprs,
		double **Alpha, double *Mu, int n_peptide, double *Sig2, int *cen_pos);

void store_mcmc_output1(double *Mu, double *A, double *B, double *P, double *Sig2, double *D,
		double *Theta, double kappa,
		double m, double c_var, int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		double lambda_a, double lambda_b, int MRF, int ind_pep,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *THETAFILE);

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
		int *MRF, int *lambda_prior, int *ind_pep,
		double *OutProbs, double *mean_fitted, int *write, double *eps, int *nmax,
		int *accept, int *silent)
{
	R_CStackLimit=(uintptr_t)-1;
	int p, i, n = 0, j;

	double **Exprs;
	double **Alpha;
	double **xA, **xB;
	double *Omega, *Theta;
	int **ProbSum;
	int **Gamma;
	double *Mu;
	double *A, *B, *P, *Sig2, *D;

	int accept_m = 0, accept_c = 0;

	double a_0, b_0, lambda_a, lambda_b, tau, alpha_0, beta_0, m_0, v_0;
	double r_a, r_b, kappa;

	a_0 = hyper_parms[0];
	b_0 = hyper_parms[1];

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
	Omega = (double*) malloc(*n_peptide*sizeof(double));
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
	Mu = (double*) malloc(*n_indiv*sizeof(double));

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
	double alpha = 1.0;
	double beta = 1.0;
	double adptm = 5.0;
	double adptv = 1.0;

	FILE *AFILE, *BFILE, *PFILE, *VARFILE, *Sig2FILE, *MUFILE, *DFILE, *THETAFILE;

	if(*write == 1)
	{
		AFILE = fopen("afile.txt", "w");
		BFILE = fopen("bfile.txt", "w");
		PFILE = fopen("pfile.txt", "w");
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
	*nmax = (*nmax < 10) ? *nmax : 10;

	// initialize the probsum values
	for(i = 0; i < *n_indiv; i++)
	{
		for(p = 0; p < *n_peptide; p++)
		{
			ProbSum[i][p] = 0;
			Exprs[i][p] = Y[*n_peptide*i + p];
		}
	}

	int th_id;
	RngStream rng[*nP];
	for(i = 0; i < *nP; i++)
	{
		rng[i] = RngStream_CreateStream("");
		//RngStream_IncreasedPrecis(rngs[i], 1);
	}
	GetRNGstate();
	InitializeChain1(Omega, Alpha, Gamma, Sig2, Mu, A, B, P, Theta,
			&c_var, &m, &alpha, &beta,
			n_position, pstart, pnum, n_peptide, n_indiv,
			*nmax, xA, xB,
			*MRF, *ind_pep);

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

		update_global_parms(Alpha, Gamma, &m, &c_var, A, B,
				m_0, v_0,
				*lambda_prior, &lambda_a, &lambda_b, r_a, r_b,
				*n_position, *n_peptide, *n_indiv, &accept_m,
				&accept_c, rng[0], adptm, adptv);


#pragma omp parallel private(th_id, workspace) num_threads(*nP)
		{
			th_id = omp_get_thread_num();
#pragma omp for nowait
			for(p = 0; p < *n_position; p++)
			{
				update_position_p1(Exprs, Omega, Alpha, Gamma,
						Sig2, Mu, A, B, P,
						a_0, b_0, lambda_a, lambda_b,
						p, n_indiv, pnum, pstart, rng[th_id], *nmax,
						xA[p], xB[p], &workspace, *eps, *ind_pep, *MRF);
			}

#pragma omp for nowait
			for(j = 0; j < *n_indiv; j++)
			{
				update_indiv_j1(Exprs[j], Alpha[j], &Mu[j], Omega, Gamma[j],
						Sig2, m, c_var, tau, pnum, pstart,
						n_position, n_peptide, *ind_pep, rng[th_id]);
			}
		}

		// check whether we need to update complete data
		if(*cen_num > 0)
		{
			update_data(D, *cen_num, cen_ind, cen_pep, Y, Exprs,
					Alpha, Mu, *n_peptide, Sig2, cen_pos);
		}
		// update MRF if applicable
		if(*MRF == 1)
		{
			update_MRF(P, Theta, Omega,
					&kappa, pnum, pstart, n_position, accept, rng[0]);
		}

		if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ))
		{
			n = n + 1.0;
			update_prob_include(n_peptide, n_indiv, Gamma, ProbSum, mean_fitted,
					Alpha, Mu, n);
			if(*write == 1)
			{
				store_mcmc_output1(Mu, A, B, P, Sig2, D, Theta, kappa,
						m, c_var, n_peptide, n_indiv, n_position, cen_num,
						lambda_a, lambda_b, *MRF, *ind_pep,
						AFILE, BFILE, PFILE, VARFILE, Sig2FILE, MUFILE, DFILE,
						THETAFILE);
			}
		}

		// adaptation increment
		if(i < *n_burn & (i % 50 == 0))
		{
			if(accept_c > 24)
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
			if(accept_m > 24)
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

	PutRNGstate();
	free(A);
	free(B);
	free(P);
	free(Sig2);
	free(Mu);

	if(*MRF == 1)
	{
		fclose(THETAFILE);
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
	free(Omega);
	free(Exprs);
	free(Alpha);
	free(ProbSum);
	free(Gamma);


	fclose(AFILE);
	fclose(BFILE);
	fclose(PFILE);
	fclose(VARFILE);
	fclose(Sig2FILE);
	fclose(MUFILE);
	if(*cen_num > 0)
	{
		fclose(DFILE);
	}

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
		double *Theta, double kappa,
		double m, double c_var, int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		double lambda_a, double lambda_b, int MRF, int ind_pep,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *THETAFILE)
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
	fprintf(AFILE, "\n");
	fprintf(BFILE, "\n");
	fprintf(PFILE, "\n");
	fprintf(Sig2FILE, "\n");
	fprintf(MUFILE, "\n");
	fprintf(VARFILE, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \n", m, c_var,
			lambda_a, lambda_b, kappa);
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

		tmp = truncNorm(Mu[i] + Alpha[i][p] - Y[i*n_peptide + p], Sig2[pos]);
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
			mean_fitted[*n_peptide*i + p] = mean_fitted[*n_peptide*i + p]*((n-1.0)/n) + (Alpha[i][p] + Mu[i])/n;
		}
	}
	return;
}


void InitializeChain1(double *Omega, double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P, double *Theta,
		double *c_var, double *m, double *alpha, double *beta,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		int nmax, double **xA, double **xB,
		int MRF, int ind_pep)
{
	int c, p, i, g;
	for(p = 0; p < *n_position; p++)
	{
		A[p] = 1.0;
		B[p] = 2.0;
		if(MRF == 1)
		{
			Theta[p] = 0;
		}
		P[p] = 0.5;

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

	*alpha = 1.0;
	*beta = 1.0;
	*c_var = 1.0;
	*m = rnorm(0.0, 1.0);

	for(i = 0; i < *n_indiv; i++)
	{
		Mu[i] = 0.0;
		for(p = 0; p < *n_position; p++)
		{
			for(c = 0; c < pnum[p]; c++)
			{
				if(ind_pep == 1)
				{
					Sig2[pstart[p] + c] = 1.0;
				}

				if(runif(0.0,1.0) <= P[p])
				{
					Omega[pstart[p] + c] = 0.0;
				}
				else
				{
					Omega[pstart[p] + c] = rbeta(A[p], B[p]);
				}

				Gamma[i][pstart[p] + c] =
						(int)(runif(0.0, 1.0) <= Omega[pstart[p] + c]);

				if(Gamma[i][pstart[p] + c] == 0)
				{
					Alpha[i][pstart[p] + c] = 0.0;
				}
				else
				{
					Alpha[i][pstart[p] + c] = truncNorm(*m, *c_var);
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
		double *A, double *B, double m_0, double v_0, int lambda_prior,
		double *lambda_a, double *lambda_b, double r_a, double r_b,
		int n_position, int n_peptide, int n_indiv, int *accept_m,
		int *accept_c, RngStream rng, double adptm, double adptv)
{
	double n_alpha = 0.0;
	double sum_alpha = 0.0;
	double v_prime, d_prime;
	double SS_alpha_resid = 0.0;
	double m_prop, c_prop;
	double log_DF_ratio;

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

void update_lambdas(double *lambda_a, double *lambda_b, double r_a,
		double r_b, double *A, double *B, int n_position,
		RngStream rng)
{
	int p;
	double S_a = 0, S_b = 0;

	for(p = 0; p < n_position; p++)
	{
		S_a += A[p];
		S_b += B[p];
	}

	*lambda_a = RngStream_GA1((double)(n_position) + 1.0, rng)/(S_a + r_a);
	*lambda_b = RngStream_GA1((double)(n_position) + 1.0, rng)/(S_b + r_b);

	return;
}

void update_indiv_j1(double *Exprs, double *Alpha, double *mu_j, double *Omega, int *Gamma,
		double *Sig2, double m, double c_var, double tau_0, int *pnum, int *pstart,
		int *n_position, int *n_peptide, int ind_pep,  RngStream rng)
{
	// update Mu_j
	double S_dot = 0.0, V_dot = 1.0/tau_0;
	int p, c;

	// peptide variances common across position
	if(ind_pep == 0)
	{
		for(p = 0; p < *n_position; p++)
		{
			for(c = 0; c < pnum[p]; c++)
			{
				S_dot += (Exprs[pstart[p] + c] - Alpha[pstart[p] + c])/Sig2[p];
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
			S_dot += (Exprs[p] - Alpha[p])/Sig2[p];
			V_dot += 1.0/Sig2[p];
		}
		*mu_j = RngStream_N01(rng)/sqrt(V_dot) + S_dot/V_dot;
	}

	// update Alpha_icp and Gamma_icp, the individual fitted effects.
	int cur;
	double m_prime, v_prime, R;
	double var;
	for(p = 0; p < *n_position; p++)
	{
		if(ind_pep == 0)
		{
			var = Sig2[p];
		}

		for(c = 0; c < pnum[p]; c++)
		{
			if(Omega[pstart[p] + c] > 0.0)
			{
				if(ind_pep == 1)
				{
					var = Sig2[pstart[p] + c];
				}
				v_prime = (c_var*var)/(var + c_var);
				m_prime = ((Exprs[pstart[p] + c] - *mu_j)*c_var + m*(var))/(c_var + var);
				R = .5*log(v_prime/c_var);
				R = R + pnorm(0.0, m_prime, sqrt(v_prime), 0, 1);
				R = R - pnorm(0.0, m, sqrt(c_var), 0, 1);

				R = exp(gsl_pow_2(m_prime)/(2.0*v_prime)
						- gsl_pow_2(m)/(2.0*(c_var)) + R);

				cur = (int)(RngStream_RandU01(rng) <= ((R*Omega[pstart[p] + c])/(1.0 - Omega[pstart[p] + c] + Omega[pstart[p] + c]*R)));

				if(cur == 1)
				{
					Gamma[pstart[p] + c] = 1;
					Alpha[pstart[p] + c] = truncNorm_parallel(m_prime, v_prime, rng);
				}
				else if(cur == 0)
				{
					Gamma[pstart[p] + c] = 0;
					Alpha[pstart[p] + c] = 0.0;
				}
				else
				{
					Rprintf("Something weird! Cur = %d, R = %lf \n", cur, R);
				}
			}
			else
			{
				Gamma[pstart[p] + c] = 0;
				Alpha[pstart[p] + c] = 0.0;
			}
		}
	}
	return;
}

int update_position_p1(double **Exprs, double *Omega, double **Alpha, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *P,
		double a_0, double b_0, double lambda_a, double lambda_b,
		int p, int *n_indiv, int *pnum, int *pstart, RngStream rng, int nmax,
		double *xA, double *xB, ARS_workspace *workspace, double eps,
		int ind_pep, int MRF)
{
	double s_log_w = 0.0, s_log_1minus_w = 0.0;
	double frac, R;

	int c, S_p = 0, i;
	int num_x = 2;
	int p_begin = pstart[p];
	double N_cp = 0.0;

	// update A's and B's
	for(c = 0; c < pnum[p]; c++)
	{
		if(Omega[p_begin + c] > 1.0 - ZZERO)
		{
			Omega[p_begin + c] = Omega[p_begin + c] - .0000001;
			//Rprintf("Omega essentially 1. Slightly reducing for numerical considerations\n");
		}
		if((Omega[p_begin + c] > ZZERO) & ((1.0 - Omega[p_begin + c]) > ZZERO))
		{
			s_log_w += log(Omega[p_begin + c]);
			s_log_1minus_w += log1p(-1.0*Omega[p_begin + c]);
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
		//Rprintf("%lf", R);
		A[p] = sample_conditional(xA, &num_x, nmax, (double)S_p, B[p], R,
				workspace, rng, eps);

		if(A[p] == -1.0)
		{
			return(0);
		}

		num_x = 2;
		R = (lambda_b - s_log_1minus_w);
		//Rprintf("%lf\n", R);
		B[p] = sample_conditional(xB, &num_x, nmax, (double)S_p, A[p], R,
				workspace, rng, eps);
		if(B[p] == -1.0)
		{
			return(0);
		}
	}
	//update Ps
	if(MRF == 0)
	{
		P[p] = RngStream_Beta(a_0 + (double)(pnum[p] - S_p), b_0 + (double)(S_p), rng);
	}
	frac = beta(A[p], (double)(*n_indiv) + B[p])/beta(A[p], B[p]);

	// update Omegas
	for(c = 0; c < pnum[p]; c++)
	{
		for(i = 0; i < *n_indiv; i++)
		{
			N_cp += (double)Gamma[i][p_begin + c];
		}

		R = (double)(N_cp == 0.0)*P[p]/(P[p] + (1.0 - P[p])*frac);
		if(RngStream_RandU01(rng) <= R)
		{
			Omega[p_begin + c] = 0.0;
		}
		else
		{
			Omega[p_begin + c] = RngStream_Beta(A[p] + N_cp,
					B[p] + (double)(*n_indiv) - N_cp, rng);
		}
		N_cp = 0.0;
	}

	double SS_i = 1.0;
	//update variances
	if(ind_pep == 0)
	{
		for(i = 0; i < *n_indiv; i++)
		{
			for(c = 0; c < pnum[p]; c++)
			{
				SS_i += gsl_pow_2(Exprs[i][p_begin + c] - Mu[i] - Alpha[i][p_begin + c]);
			}
		}
		SS_i = SS_i/2.0;
		Sig2[p] = SS_i/RngStream_GA1(((double)(*n_indiv*pnum[p]) + 1.0)/2.0, rng);
	}

	else
	{
		for(c = 0; c < pnum[p]; c++)
		{
			for(i = 0; i < *n_indiv; i++)
			{
				SS_i += gsl_pow_2(Exprs[i][p_begin + c] - Mu[i] - Alpha[i][p_begin + c]);
			}
			SS_i = SS_i/2.0;
			Sig2[p_begin + c] = SS_i/RngStream_GA1(((double)(*n_indiv) + 1.0)/2.0, rng);
			SS_i = 1.0;
		}

	}
	return(1);
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

double expit(double x)
{
	double p;
	p = 1.0/(1.0 + exp(-x));
	return(p);
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

