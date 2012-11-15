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

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "RngFunctions.h"
#include "RngStream.h"
#include "ARS.h"

#define	TEST_INT_LENGTH 200
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>

#include "pepBayes_v1.h"

/*
 * individual specific mean, with prior N(0, kappa)
 * peptide specific scale with IG(A,B) prior
 * peptide specific alpha effects with global mean/var
 * t-distributed errors with estimated DoF.
 */

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

void pepbayes_v1(double *Y, double *hyper_param, int *pstart,
		int *pnum, int *n_position, int *n_peptide, int *n_indiv, int *nP,
		int *cen_ind, int *cen_pep, int *cen_num, int *cen_pos,
		int *n_iter, int *n_sweep, int *n_burn,
		double *OutProbs, int *write)
{
	R_CStackLimit = (uintptr_t)-1;
	// miscellaneous useful quantities
	int p, i, j, k = 0, c, d;
	int pep, pos, th_id;
	int *pos_ind_by_pep;
	double *ProbSum;
	int total_iterations = *n_burn+(*n_sweep)*(*n_iter);
	int percent_complete = 10;

	adpt zeta_adpt;
	adpt m_adpt;

	zeta_adpt.total_count = 0;
	zeta_adpt.count = 0;
	zeta_adpt.tune = 2.5;
	m_adpt.total_count = 0;
	m_adpt.count = 0;
	m_adpt.tune = 5.0;

	RngStream rng[*nP];

	//statically allocated struct for Adaptive Rejection Sampling
	ARS_workspace workspace;

	// missing data variables
	double *Exprs;			// censor-completed expressions
	int *Gamma;			// cluster membership indicator
	double *W;				// t-distribution weights
	double dof = 4.0;				// t-distribution degrees of freedom
	double *D;				// imputed censored tails

	// component membership probability parameters
	double *Omega_logit;	// membership probabilities, on logit scale
	int *Omega_ind;			// indicator for whether omega_{cp} > 0
	double *A, *B, *U; 		// prior parameters for omega probabilities
	double lambda_a, lambda_b; //hyperparameters for A, B
	double a_0, b_0;		// hyperparameters for U
	double **xA, **xB;		// vectors required for ARS

	// location-related parameters
	double *Mu;				// individual specific means
	double kappa = 10.0;	// prior mean precision
	double *Alpha_pep;		// peptide effects
	double m = 2.0, zeta = 1.0;			// prior mean and variance of peptide effects
	double m_0, v_0;		// hyperparameters for m

	// variance parameters
	double *Sig2;			// peptide response variance
	double alpha = 10.0;	// prior parameter
	double beta = 10.0; 	// prior parameter
	double alpha_0, beta_0; // rates on alpha, beta (hyperparameters)
	double *xAlpha;			// vector needed for ARS

	// file pointers
	FILE *AFILE, *BFILE, *PFILE, *VARFILE, *Sig2FILE, *MUFILE, *DFILE, *OFILE;
	FILE *ALPHAFILE;

	// retreive and initialize hyperparameter values
	a_0 = hyper_param[0];
	b_0 = hyper_param[1];
	lambda_a = hyper_param[2];
	lambda_b = hyper_param[3];
	alpha_0 = hyper_param[4];
	beta_0 = hyper_param[5];
	m_0 = hyper_param[6];
	v_0 = hyper_param[7];

	// begin memory allocation
	Gamma = (int*) R_alloc(*n_peptide*(*n_indiv), sizeof(int));
	Exprs = (double*) R_alloc(*n_peptide*(*n_indiv), sizeof(double));
	W = (double *) R_alloc(*n_peptide*(*n_indiv), sizeof(double));
	double *RB = (double *) R_alloc(*n_peptide*(*n_indiv), sizeof(double));

	// xA and xB hold starting values for ARS of a_p, b_p, and alpha.
	xA = (double**) R_alloc(*n_position, sizeof(double*));
	xB = (double**) R_alloc(*n_position, sizeof(double*));
	xAlpha = (double*) R_alloc(NMAX, sizeof(double));
	// initial values for hull quantiles.
	xAlpha[0] = 1.0;
	xAlpha[1] = 2.0;

	Omega_ind = (int*) R_alloc(*n_peptide, sizeof(int));
	Omega_logit = (double*) R_alloc(*n_peptide, sizeof(double));
	pos_ind_by_pep = (int*) R_alloc(*n_peptide, sizeof(int));
	ProbSum = (double*) R_alloc(*n_peptide*(*n_indiv), sizeof(double));

	for(p = 0; p < *n_position; p++)
	{
		xA[p] = (double*) R_alloc(NMAX, sizeof(double));
		xB[p] = (double*) R_alloc(NMAX, sizeof(double));
		for(c = 0; c < pnum[p]; c++)
		{
			pos_ind_by_pep[pstart[p] + c] = k;
		}
		k++;
	}

	Alpha_pep = (double*) R_alloc(*n_peptide, sizeof(double));
	Sig2 = (double*) R_alloc(*n_peptide, sizeof(double));
	A = (double*) R_alloc(*n_position, sizeof(double));
	B = (double*) R_alloc(*n_position, sizeof(double));
	U = (double*) R_alloc(*n_position, sizeof(double));
	Mu = (double*) R_alloc(*n_indiv, sizeof(double));
	double* likworkspace = (double *) R_alloc((*n_peptide)*(*n_indiv), sizeof(double));

	// check whether our data is censored at all,
	// if so prepare for augmentation.

	if(*cen_num > 0)
	{
		D = (double*) R_alloc(*cen_num, sizeof(double));
		for(i = 0; i < *cen_num; i++)
		{
			D[i] = 0.0;
		}
	}

	if(*write == 1)
	{
		AFILE = fopen("afile.txt", "w");
		BFILE = fopen("bfile.txt", "w");
		PFILE = fopen("pfile.txt", "w");
		OFILE = fopen("ofile.txt", "w");
		ALPHAFILE = fopen("alphafile.txt", "w");
		VARFILE = fopen("varfile.txt", "w");
		Sig2FILE = fopen("sig2file.txt", "w");
		MUFILE = fopen("mufile.txt", "w");
		if(*cen_num > 0)
		{
			DFILE = fopen("dfile.txt", "w");
		}
	}
	for(i = 0; i < *nP; i++)
	{
		rng[i] = RngStream_CreateStream("");
	}

	GetRNGstate();
	initialize_chain(ProbSum, Exprs, Y, W, Omega_ind, Omega_logit, Alpha_pep,
			Gamma, Sig2, Mu, A, B, U,
			n_position, pstart, pnum, n_peptide, n_indiv,
			xA, xB, RB);
	PutRNGstate();


	Rprintf("parameters initialized \n");
	for(i = 1; i <= total_iterations; i++)
	{
		R_CheckUserInterrupt();
		update_dof_integrated(&dof, Exprs, W, Alpha_pep, Gamma,
				Sig2, Mu, likworkspace, *n_indiv, *n_peptide, rng[0]);

#pragma omp parallel private(th_id, workspace, pos) num_threads(*nP)
		{
			th_id = omp_get_thread_num();
#pragma omp for
			for(pep = 0; pep < *n_peptide; pep++)
			{
				pos = pos_ind_by_pep[pep];
				update_peptide(Exprs, Mu, Alpha_pep + pep,
						W, Omega_ind + pep, Omega_logit + pep, Gamma,
						Sig2 + pep, alpha, beta,
						U[pos], A[pos], B[pos],
						m, zeta, dof, *n_indiv, *n_peptide, pep, rng[th_id], RB);
			}
			// update individual means
#pragma omp for
			for(j = 0; j < *n_indiv; j++)
			{
				update_indiv_mu(Exprs, W, Alpha_pep, Mu + j, Gamma,
						Sig2, kappa, *n_peptide, j, rng[th_id]);
			}
#pragma omp for
			for(p = 0; p < *n_position; p++)
			{
				update_position_p(Omega_ind, Omega_logit,
						A + p, B + p, U + p,
						a_0, b_0, lambda_a, lambda_b,
						p, *n_indiv, pnum, pstart, rng[th_id],
						xA[p], xB[p], &workspace);
			}
			// update gammas, alphas, omegas

			if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ))
			{
#pragma omp for
				for(d = 0; d < (*n_indiv)*(*n_peptide); d++)
				{
					ProbSum[d] += RB[d];
				}
			}
		}
		update_global_params(&m, &zeta, &alpha, &beta, &kappa,
				Mu, Sig2, Alpha_pep,
				m_0, v_0, alpha_0, beta_0,
				*n_indiv, *n_peptide, rng[0], &m_adpt, &zeta_adpt,
				&workspace, xAlpha, *nP);

		if(i % TEST_INT_LENGTH == 0)
		{
			update_tuning(&m_adpt, i, *n_burn);
			update_tuning(&zeta_adpt, i, *n_burn);
		}

		// check whether we need to update complete data
		if((*cen_num > 0) & (i > *n_burn/2))
		{
			update_censoring(W, D, *cen_num, cen_ind, cen_pep, Y, Exprs,
					Gamma, Alpha_pep, Mu, *n_peptide, Sig2, rng[0]);
		}

		if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ) && *write == 1)
		{
			store_mcmc_output(Alpha_pep, Mu, A, B, U, Sig2, D, Omega_logit,
					Omega_ind, kappa, alpha, beta, m, zeta,
					dof, *n_peptide, *n_indiv, *n_position, cen_num,
					AFILE, BFILE, PFILE, VARFILE, Sig2FILE, MUFILE, DFILE,
					OFILE, ALPHAFILE);
		}

		if( 100*((double) i)/total_iterations >= percent_complete)
		{
			Rprintf("MCMC %d percent complete\n", percent_complete);
			percent_complete += 10;
		}
	}

	// finalize marginal PPA
	for(d = 0; d < (*n_indiv)*(*n_peptide); d++)
	{
		OutProbs[d] = ProbSum[d]/(*n_iter);
	}

	Rprintf("M acceptance rate: %.3lf, Zeta acceptance rate: %.3lf\n", (double)(m_adpt.total_count)/(*n_iter*(*n_sweep)),
			(double)(zeta_adpt.total_count)/(*n_iter*(*n_sweep)));

	Rprintf("closing files\n");
	if(*write == 1)
	{
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
	}
	return;
}

void update_global_params(double *m, double *zeta, double *alpha, double *beta, double *kappa,
		double *Mu, double *Sig2, double *Alpha_pep,
		double m0, double v0, double alpha0, double beta0,
		int n_indiv, int n_peptide, RngStream rng, adpt *m_tune, adpt *z_tune,
		ARS_workspace *ws, double *xAlpha, int nP)
{
	double v_prime;
	double m_prop, zeta_prop;
	double log_DF_ratio;
	double sum_alpha = 0.0;
	int pep;
	//update m
#pragma omp parallel for reduction(+:sum_alpha)
	for(pep = 0; pep < n_peptide; pep++)
	{
		sum_alpha = Alpha_pep[pep] + sum_alpha;
	}

	v_prime = (*zeta*v0)/(*zeta + v0*n_peptide);
	m_prop = *m + (RngStream_RandU01(rng)*2.0 - 1.0)*sqrt(v_prime)*m_tune->tune; //sqrt(v_prime)*rnorm(0.0,1.0);
	log_DF_ratio = pnorm(0.0, *m, sqrt(*zeta), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, m_prop, sqrt(*zeta), 0, 1);
	log_DF_ratio = n_peptide*log_DF_ratio;
	log_DF_ratio += -(.5/v_prime)*(gsl_pow_2(m_prop) - gsl_pow_2(*m));
	log_DF_ratio += ((*zeta*m0 + v0*sum_alpha)/(*zeta*v0))*(m_prop - *m);

	if(log(RngStream_RandU01(rng)) <= log_DF_ratio)
	{
		*m = m_prop;
		m_tune->count = m_tune->count + 1;
	}

	// update zeta
	double SS_alpha = 0.0;
	double d_prime;
#pragma omp parallel for reduction(+:SS_alpha)
	for(pep = 0; pep < n_peptide; pep++)
	{
		SS_alpha = SS_alpha + gsl_pow_2(*m - Alpha_pep[pep]);
	}
	SS_alpha += 10.0;

	// log-scale random walk proposal
	d_prime = n_peptide/2.0 + 5.0;
	zeta_prop = *zeta*exp(z_tune->tune*(RngStream_RandU01(rng) - .5));

	log_DF_ratio = pnorm(0.0, (*m), sqrt(*zeta), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, (*m), sqrt(zeta_prop), 0, 1);
	log_DF_ratio = n_peptide*log_DF_ratio;
	log_DF_ratio += .5*(1.0/(*zeta) - 1.0/zeta_prop)*SS_alpha;
	log_DF_ratio += d_prime*(log(*zeta) - log(zeta_prop));

	// accept or reject step
	if(log(RngStream_RandU01(rng)) <= log_DF_ratio)
	{
		*zeta = zeta_prop;
		z_tune->count++;
	}

	//update alpha and kappa
	int i;
	double sum_u_sqr = 0.0;
	double sum_sig2_inv = 0.0;
	double sum_log_sig2 = 0.0;
	double argvec[4];
#pragma omp parallel num_threads(nP)
	{
#pragma omp for reduction(+:sum_u_sqr)
		for(i = 0; i < n_indiv; i++)
		{
			sum_u_sqr = sum_u_sqr + gsl_pow_2(Mu[i]);
		}
#pragma omp for reduction(+:sum_sig2_inv)
		for(pep = 0; pep < n_peptide; pep++)
		{
			sum_sig2_inv = sum_sig2_inv + 1.0/Sig2[pep];
		}
#pragma omp for reduction(+:sum_log_sig2)
		for(pep = 0; pep < n_peptide; pep++)
		{
			sum_log_sig2 = sum_log_sig2 + log(Sig2[pep]);
		}
	}
	sum_u_sqr += .5;
	sum_sig2_inv += beta0;

	*kappa = RngStream_GA1(n_indiv/2.0 + .5, rng)/sum_u_sqr;

	argvec[0] = alpha0;
	argvec[1] = (double)(n_peptide);
	argvec[2] = log(sum_sig2_inv);
	argvec[3] = sum_log_sig2;
	int arglen = 4;
	int num_x = 2;
	*alpha = sample_conditional(xAlpha, &num_x, argvec, &arglen, ws, rng,
			lc_alpha_int, lcp_alpha_int);

	//update beta
	*beta = RngStream_GA1(*alpha*n_peptide + 1.0, rng)/sum_sig2_inv;
	return;
}

void update_peptide(double *Exprs, double *Mu, double *Alpha_pep, double *W,
		int *Omega_ind, double *Omega_logit, int *Gamma,
		double *Sig2_pep, double alpha, double beta,
		double u_logit, double a, double b,
		double m, double zeta, double dof, int n_indiv, int n_peptide, int pep,
		RngStream rng, double *RB)
{
	int i, S = 0;
	double s_init = 0.0, w_init = 0.0;
	double w, C;
	int cur;
	double s2_pep = *Sig2_pep;
	double a_pep = *Alpha_pep;
	int d_ind;
	double SS = 0.0;
	double delta0, delta1;
	double w1 = .5 + dof/2.0;
	double w2 = dof/2.0;
	double tmp;

	d_ind = pep; // start indexing through large arrays
	if(*Omega_ind == 1)
	{
		for(i = 0; i < n_indiv; i++)
		{
			delta0 = Exprs[d_ind] - Mu[i];
			delta1 = delta0 - a_pep;

			C = (dof/2.0 + .5)*(log1p(gsl_pow_2(delta0)/(s2_pep*dof)) -
					log1p(gsl_pow_2(delta1)/(s2_pep*dof)));
			RB[d_ind] = expit(C + *Omega_logit);

			cur = (RngStream_RandU01(rng) < RB[d_ind]) ? 1 : 0;
			Gamma[d_ind] = cur;

			tmp = gsl_pow_2(delta0 - a_pep*cur)/s2_pep;
			w = RngStream_GA1(w1, rng)/(w2 + tmp/2.0);
			W[d_ind] = w;
			//quantities for alpha, omega, sigma updates
			S += cur;
			s_init += cur*delta0*w;
			w_init += cur*w;
			SS += w*tmp*s2_pep;

			d_ind += n_peptide;
		}
	}
	// if here, all gamma's for peptide p are zero
	else
	{
		for(i = 0; i < n_indiv; i++)
		{
			delta0 = Exprs[d_ind] - Mu[i];
			tmp = gsl_pow_2(delta0);
			w = RngStream_GA1(w1, rng)/(w2 + tmp/(2.0*s2_pep));
			W[d_ind] = w;
			RB[d_ind] = 0.0;

			//quantities for sigma updates
			SS += w*tmp;
			d_ind += n_peptide;
		}
	}

	// update scale sigma^2_{cp}
	SS = SS*.5 + beta;
	s2_pep = SS/RngStream_GA1(alpha + n_indiv/2.0, rng);
	*Sig2_pep = s2_pep;

	// update omega, first check whether omega_{cp} > 0:
	if(S == 0)
	{
		double log_frac = gsl_sf_lnbeta(a, b + n_indiv) - gsl_sf_lnbeta(a, b);
		*Omega_ind = (logit(RngStream_RandU01(rng)) < (u_logit - log_frac)) ? 0 : 1;
	}

	// if omega_{cp} > 0, update this.
	if(*Omega_ind == 1)
	{
		*Omega_logit = RngStream_LogitBeta(a + S, b + n_indiv - S, rng);
	}

	// update alpha
	double s_alpha, v_alpha;
	s_alpha = s_init/s2_pep + m/zeta;
	v_alpha = w_init/s2_pep	+ 1.0/zeta;
	a_pep = truncNorm(s_alpha/v_alpha, 1.0/v_alpha, rng);
	*Alpha_pep = a_pep;

	return;
}

void store_mcmc_output(double *Alpha, double *Mu, double *A, double *B, double *U,
		double *Sig2, double *D, double *Omega_Logit, int* Omega_Ind, double kappa,
		double alpha, double beta, double m, double zeta, double dof,
		int n_peptide, int n_indiv, int n_position, int *cen_num,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *OFILE, FILE *ALPHAFILE)
{
	int p, i, k;

	for(p = 0; p < n_position; p++)
	{
		fprintf(AFILE, "%.6lf \t", A[p]);
		fprintf(BFILE, "%.6lf \t", B[p]);
		fprintf(PFILE, "%.6lf \t", U[p]);
	}

	for(i = 0; i < n_indiv; i++)
	{
		fprintf(MUFILE, "%.6lf \t", Mu[i]);
	}

	for(p = 0; p < n_peptide; p++)
	{
		fprintf(OFILE, "%.6lf \t", (Omega_Ind[p] == 1) ? expit(Omega_Logit[p]):0.0);
		fprintf(ALPHAFILE, "%.6lf \t", Alpha[p]);
		fprintf(Sig2FILE, "%.6lf \t", Sig2[p]);
	}

	if(*cen_num > 0)
	{
		for(k = 0; k < *cen_num; k++)
		{
			fprintf(DFILE, "%.6lf \t", D[k]);
		}
		fprintf(DFILE, "\n");
	}

	fprintf(ALPHAFILE, "\n");
	fprintf(AFILE, "\n");
	fprintf(BFILE, "\n");
	fprintf(PFILE, "\n");
	fprintf(OFILE, "\n");
	fprintf(Sig2FILE, "\n");
	fprintf(MUFILE, "\n");
	fprintf(VARFILE, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t\n",
			m, zeta, kappa, alpha, beta, dof);
	return;
}

void update_censoring(double *W, double *D, int cen_num, int* cen_ind,
		int* cen_pep, double *Y, double *Exprs,
		int *Gamma, double *Alpha, double *Mu, int n_peptide, double *Sig2,
		RngStream rng)
{
	int i, p, k, ind;
	double tmp, weight;
	for(k = 0; k < cen_num; k++)
	{
		// retrieve indices
		i = cen_ind[k];
		p = cen_pep[k];
		ind = i*n_peptide + p;
		weight = W[ind];

		tmp = truncNorm(Mu[i] + Alpha[p]*Gamma[ind] - Y[ind], Sig2[p]/weight,
				rng);
		Exprs[ind] = Y[ind] + tmp;
		D[k] = tmp;
	}
	return;
}


void initialize_chain(double *ProbSum, double *Exprs, double *Y,
		double *W, int *Omega_Ind, double *Omega_Logit,
		double *Alpha_pep, int *Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *U,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		double **xA, double **xB, double *RB)
{
	int c, p, i, g, pep;
	double prob;
	for(i = 0; i < *n_indiv; i++)
	{
		Mu[i] = 0.0;
		for(pep = 0; pep < *n_peptide; pep++)
		{
			Exprs[(*n_peptide)*i + pep] = Y[(*n_peptide)*i + pep];
			ProbSum[(*n_peptide)*i + pep] = 0.0;
			RB[(*n_peptide)*i + pep] = 0.0;
		}
	}

	for(p = 0; p < *n_position; p++)
	{
		A[p] = 1.0;
		B[p] = 3.0;
		U[p] = .25;

		for(g = 0; g < NMAX; g++)
		{
			xA[p][g] = 0.0;
			xB[p][g] = 0.0;
		}
		xA[p][0] = .1;
		xB[p][0] = .1;
		xA[p][1] = 	2.0;
		xB[p][1] =  2.0;

		for(c = 0; c < pnum[p]; c++)
		{
			pep = pstart[p] + c;

			Alpha_pep[pep] = 2.0;
			Sig2[pep] = 1.0;
			Omega_Ind[pep] = 1;
			Omega_Logit[pep] = log(rgamma(A[p], 1.0)/rgamma(B[p], 1.0));

			for(i = 0; i < *n_indiv; i++)
			{
				prob = gsl_pow_2(Y[(*n_peptide)*i + pep] - Mu[i])/
						(gsl_pow_2(Y[(*n_peptide)*i + pep] - Mu[i]) +
						gsl_pow_2(Y[(*n_peptide)*i + pep] - Mu[i] - Alpha_pep[pep]));
				W[(*n_peptide)*i + pep] = 1.0;
				Gamma[(*n_peptide)*i + pep] = (int)(unif_rand() <= prob);
			}
		}
	}
	return;
}

void update_dof_integrated(double *dof, double *Exprs, double *W,
		double *Alpha, int *Gamma,
		double *Sig2, double *Mu, double *workspace,
		int n_indiv, int n_peptide,
		RngStream rng)
{
	double V[] = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 32.0, 64.0};

	int i, pep, k, nv = 9;
	const int n_obs = n_indiv*n_peptide;
	double log_cprob[nv];
	double log_lik[nv];
	double lik_sum[nv];
	double v;
	double ds;
	double lik_temp = 0.0;

	for(k = 0; k < nv; k++)
	{
		log_lik[k] = 0.0;
	}

	int d_ind;
#pragma omp parallel for private(ds, d_ind, pep)
	// compute the vicious log-sum
	for(i = 0; i < n_indiv; i++)
	{
		d_ind = i*n_peptide;
		const double mu = Mu[i];
		for(pep = 0; pep < n_peptide; pep++)
		{
			const double y = Exprs[d_ind];
			const int gam = Gamma[d_ind];
			const double sig2 = Sig2[pep];
			const double alpha = gam*Alpha[pep];

			ds = gsl_pow_2(y - mu - alpha);
			ds = ds / sig2;

			workspace[d_ind] = ds;
			d_ind++;
		}
	}

#pragma omp parallel for private(lik_temp, i, v)
	for(k = 0; k < nv; k++)
	{
		v = V[k];
		lik_temp = 0.0;
		for(i = 0; i < (n_indiv*n_peptide); i++)
		{
			lik_temp = lik_temp + log1p(workspace[i]/v);
		}
		log_lik[k] = lik_temp;
	}

	double lik_max;
	// compute normalized probabilities on log scale
	for(i = 0; i < nv; i++)
	{
		v = V[i];
		log_lik[i] = -n_obs*log(v)/2.0 + n_obs*lgamma(.5 + v/2.0) -
				n_obs*lgamma(v/2.0) - (v/2.0 + .5)*log_lik[i];
		lik_max = (i == 0) ? log_lik[i] : lik_max;
		lik_max = (lik_max > log_lik[i]) ? lik_max : log_lik[i];

	}

	for(i = 0; i < nv; i++)
	{
		v = V[i];
		log_lik[i] = log_lik[i] - lik_max;
		lik_sum[i] = (i == 0) ? log_lik[i] : log_apb(log_lik[i], lik_sum[i-1]);
	}

	for(i = 0; i < nv; i++)
	{
		log_cprob[i] = lik_sum[i] - lik_sum[nv - 1];
	}
	// sample from V
	int j = 0;
	double u = log(RngStream_RandU01(rng));
	while(1)
	{
		if(log_cprob[j] >= u)
		{
			break;
		}
		j++;
	}

	*dof = V[j];
	return;
}

void update_indiv_mu(double *Exprs, double *W,
		double *Alpha_pep, double *mu_j, int *Gamma,
		double *Sig2, double kappa,
		int n_peptide, int j, RngStream rng)
{
	// update Mu_j
	double S = 0.0, V = 0.0, w, ssq_pep, r;
	int pep, d_ind = n_peptide*j;
	for(pep = 0; pep < n_peptide; pep++)
	{
		ssq_pep = Sig2[pep];
		w = W[d_ind];
		r = w/ssq_pep;
		S += (Exprs[d_ind] - Gamma[d_ind]*Alpha_pep[pep])*r;
		V += r;
		d_ind++;
	}
	V += kappa;
	*mu_j = RngStream_N01(rng)/sqrt(V) + S/V;
	return;
}

void update_position_p(int* Omega_Ind, double *Omega_Logit,
		double *A_p, double *B_p, double *U_p,
		double a_0, double b_0, double lambda_a, double lambda_b,
		int p, int n_indiv, int *pnum, int *pstart, RngStream rng,
		double *xA, double *xB, ARS_workspace *workspace)
{
	double s_log_w = 0.0, s_log_1minus_w = 0.0;
	double R;

	int c, S_p = 0;
	int num_x = 2;
	int p_begin = pstart[p];
	double argvec[3];

	// update A's and B's
	for(c = 0; c < pnum[p]; c++)
	{
		if(Omega_Ind[p_begin + c])
		{
			s_log_w += log_from_logit(Omega_Logit[p_begin + c]);
			s_log_1minus_w += log1m_from_logit(Omega_Logit[p_begin + c]);
			S_p += 1;
		}
	}

	if(S_p == 0)
	{
		*A_p = RngStream_GA1(1.0, rng)/(lambda_a);
		*B_p = RngStream_GA1(1.0, rng)/(lambda_b);
	}

	else
	{
		R = (lambda_a - s_log_w);
		argvec[2] = R;
		argvec[0] = (double)S_p;
		argvec[1] = *B_p;

		int arglen = 3;

		*A_p = sample_conditional(xA, &num_x, argvec, &arglen,
				workspace, rng, lc_AB, lcp_AB);

		if(*A_p == -1.0)
		{
			error("invalid ARS sample");
		}

		num_x = 2;
		R = (lambda_b - s_log_1minus_w);
		argvec[2] = R;
		argvec[1] = *A_p;

		*B_p = sample_conditional(xB, &num_x, argvec, &arglen,
				workspace, rng, lc_AB, lcp_AB);
		if(*B_p == -1.0)
		{
			error("invalid ARS sample");
		}
	}
	//update Ps
	*U_p = RngStream_LogitBeta(a_0 + (double)(pnum[p] - S_p), b_0 + (double)(S_p), rng);

	return;
}

double truncNorm(double mean, double sigmasqr, RngStream rng)
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

inline double m1expit(double x)
{
	return(1.0/(1.0 + exp(x)));
}


inline double expit(double x)
{
	return(1.0/(1.0 + exp(-1.0*x)));
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

double lc_AB(double x, double *argvec, int *arglen)
{
	return(-1.0*x*argvec[2] + argvec[0]*lgamma(x + argvec[1]) - argvec[0]*lgamma(x));
}

double lcp_AB(double x, double *argvec, int *arglen)
{
	return(-1.0*argvec[2] + argvec[0]*gsl_sf_psi(x + argvec[1])
			- argvec[0]*gsl_sf_psi(x));
}

double lc_alpha_int(double x, double *argvec, int *arglen)
{
	const double Np = argvec[1];
	return(-argvec[0]*x + lgamma(Np*x + 1.0) - Np*lgamma(x) - Np*x*argvec[2] -
			x*argvec[3]);
}

double lcp_alpha_int(double x, double *argvec, int *arglen)
{

	const double Np = argvec[1];
	return(-argvec[0] + Np*gsl_sf_psi(Np*x + 1.0) - Np*gsl_sf_psi(x) - Np*argvec[2] -
			argvec[3]);
}

double lc_alpha(double x, double *argvec, int *arglen)
{
	return(-x*argvec[0] - argvec[1]*lgamma(x));
}

double lcp_alpha(double x, double *argvec, int *arglen)
{
;
	return (-argvec[0] - argvec[1]*gsl_sf_psi(x));
}

double logit(double x)
{
	if(x == 0.0)
	{
		return(-DBL_MAX);
	}
	else if(x == 1.0)
	{
		return(DBL_MAX);
	}
	else
	{
		return(log(x) - log1p(-x));
	}
}

void update_tuning(adpt *tune_set, int i, int n_burn)
{
	double accept_prop = ((double) tune_set->count)/TEST_INT_LENGTH;
	if(i <= n_burn)
	{
		if(accept_prop > .28)
		{
			tune_set->tune = tune_set->tune + .1*(n_burn + 1.0 - i)/(n_burn);
		}
		else if(accept_prop < .18)
		{
			tune_set->tune = (tune_set->tune)*
					((i)/(n_burn + 1.0) + .8*(n_burn + 1.0 - i)/(1.0 + n_burn));
		}
		tune_set->count = 0;
	}
	else
	{
		tune_set->total_count += tune_set->count;
		tune_set->count = 0;
	}
	return;
}


// initialize sums for updates
/*if(marginalize == 1)
	{
		double m0, m1, v0, v1;
		double w0 = 0.0, w1 = 0.0, s0 = 0.0, s1 = 0.0;
		d_ind = pep; // start indexing through large arrays
		for(i = 0; i < n_indiv; i++)
		{
			s_init += (Exprs[d_ind] - Mu[i])*W[d_ind]*(Gamma[d_ind]);
			w_init += W[d_ind]*Gamma[d_ind];
			d_ind += n_peptide;
		}
		// loop through subjects, updating cluster membership
		d_ind = pep; // start indexing through large arrays
		for(i = 0; i < n_indiv; i++)
		{
			delta = Exprs[d_ind] - Mu[i];
			w = W[d_ind];
			w0 = w_init - w*(Gamma[d_ind]);
			s0 = s_init - delta*w*(Gamma[d_ind]);
			w1 = w0 + w;
			s1 = s0 + delta*w;

			m0 = (m*s2_pep + s0*zeta)/(w0*zeta + s2_pep);
			m1 = (m*s2_pep + s1*zeta)/(w1*zeta + s2_pep);
			v0 = (zeta*s2_pep)/(w0*zeta + s2_pep);
			v1 = (zeta*s2_pep)/(w1*zeta + s2_pep);

			C = .5*log(v1/v0) + .5*m1*m1/v1 - .5*m0*m0/v0 +
					pnorm(0.0, m1, sqrt(v1), 0, 1) - pnorm(0.0, m0, sqrt(v0), 0, 1);
			// update gamma and get ready for next iteration
			RB[d_ind] = expit(C + *Omega_logit);
			cur = (int)(RngStream_RandU01(rng) < RB[d_ind]);
			//cur = (int)(logit(RngStream_RandU01(rng)) < (C + *Omega_logit));
			Gamma[d_ind] = cur;

			S += cur;
			w_init = (cur) ? w1 : w0;
			s_init = (cur) ? s1 : s0;
			d_ind += n_peptide;
		}
	}*/


/*
	else
	{
		d_ind = pep; // start indexing through large arrays
		for(i = 0; i < n_indiv; i++)
		{
			delta = Exprs[d_ind] - Mu[i];
			w = W[d_ind];

			C = -.5*w*a_pep*a_pep/s2_pep + w*a_pep*delta/s2_pep;
			RB[d_ind] = expit(C + *Omega_logit);

			cur = (int)(RngStream_RandU01(rng) < RB[d_ind]);
			Gamma[d_ind] = cur;

			S += cur;
			s_init += cur*delta*w;
			w_init += cur*w;
			d_ind += n_peptide;
		}
	}
*/

