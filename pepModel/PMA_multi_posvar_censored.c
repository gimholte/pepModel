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
#include <Rinterface.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "PMA_multi_posvar_censored.h"
#include "norm_gamma_generation.h"
#include "RngStream.h"
#include "ARS.h"


/*
 * individual specific mean, with prior N(0, kappa)
 * peptide specific scale with IG(A,B) prior
 * peptide specific alpha effects with global mean/var
 * t-distributed errors with estimated DoF.
 */

double lc_AB(double x, double *argvec, int *arglen);
double lcp_AB(double x, double *argvec, int *arglen);
double lc_alpha(double x, double *argvec, int *arglen);
double lcp_alpha(double x, double *argvec, int *arglen);
double logit(double x);

inline double log1m_from_logit(double x);
inline double log_from_logit(double x);
inline double expit(double x);
inline double m1expit(double x);

void tnorm_test(double* x, int *n, double *m, double *sigmasqr);

void update_prob_include(int *n_peptide, int *n_indiv, int **Gamma, int **ProbSum, double *mean_fitted,
		double *Alpha, double *Mu, double n);

double truncNorm(double mean, double sigmasqr, RngStream rng);

void update_data(double **W, double *D, int cen_num, int* cen_ind, int* cen_pep, double *Y, double **Exprs,
		int **Gamma, double *Alpha, double *Mu, int n_peptide, double *Sig2, int *cen_pos,
		RngStream rng);

void store_mcmc_output(double *Alpha, double *Mu, double *A, double *B, double *U, double *Sig2, double *D,
		double *Theta, double *Omega_Logit, int* Omega_Ind, double kappa, double alpha, double beta,
		double m, double c_var, double dof, int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		double lambda_a, double lambda_b, double a_0, double b_0, int MRF, int ind_pep,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *THETAFILE, FILE *OFILE, FILE *ALPHAFILE);

void finalize_prob_include(int *n_iter, int *n_peptide, int *n_indiv,
		double *OutProbs, double **ProbSum);

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

void PMA_mcmc_MS(double *Y, double *hyper_param, int *pstart,
		int *pnum, int *n_position, int *n_peptide, int *n_indiv, int *nP,
		int *cen_ind, int *cen_pep, int *cen_num, int *cen_pos,
		int *n_iter, int *n_sweep, int *n_burn,
		double *OutProbs, double *mean_fitted, int *write,
		int *silent, int *cladePos, int *n_clade, int *cladeCounts)
{
	R_CStackLimit=(uintptr_t)-1;

	// miscellaneous useful quantities
	int p, i, n = 0, j, k = 0, c;
	int *pos_ind_by_pep;
	int **ProbSum;
	adpt zeta_adpt;
	adpt m_adpt;

	zeta_adpt.total_count = 0;
	zeta_adpt.count = 0;
	zeta_adpt.tune = 1.0;
	m_adpt.total_count = 0;
	m_adpt.count = 0;
	m_adpt.tune = 1.0;

	// missing data variables
	double **Exprs;			// censor-completed expressions
	int **Gamma;			// cluster membership indicator
	double **W;				// t-distribution weights
	double dof;				// t-distribution degrees of freedom
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
	double m, zeta;			// prior mean and variance of peptide effects
	double m_0, v_0;		// hyperparameters for m

	// variance parameters
	double *Sig2;			// peptide response variance
	double alpha = 10.0;	// prior parameter
	double beta = 10.0; 	// prior parameter
	double alpha_0, beta_0; // rates on alpha, beta (hyperparameters)
	double *xAlpha;			// vector needed for ARS

	// retreive and initialize hyperparameter values
	a_0 = hyper_param[0];
	b_0 = hyper_param[1];
	lambda_a = hyper_param[2];
	lambda_b = hyper_param[3];
	alpha_0 = hyper_param[4];
	beta_0 = hyper_param[5];
	m_0 = hyper_param[6];
	v_0 = hyper_param[7];
	dof = hyper_param[8];

	// begin memory allocation
	Gamma = (int**) R_alloc(*n_peptide*sizeof(int*));
	Exprs = (double**) R_alloc(*n_peptide*sizeof(double*));
	W = (double **) R_alloc(*n_peptide*sizeof(double*));
	for(p = 0; p < *n_peptide; p++)
	{
		Gamma[p] = (int*) R_alloc(*n_indiv*sizeof(int));
		Exprs[p] = (double *) R_alloc(*n_indiv*sizeof(double));
		W[p] = (double *) R_alloc(*n_indiv*sizeof(double));
	}

	// xA and xB hold starting values for ARS of a_p, b_p, a_0, b_0, and alpha.
	xA = (double**) R_alloc(*n_position*sizeof(double*));
	xB = (double**) R_alloc(*n_position*sizeof(double*));
	xAlpha = (double*) R_alloc(NMAX*sizeof(double));
	// initial values for hull quantiles.
	xAlpha[0] = 1.0;
	xAlpha[1] = 2.0;

	Omega_ind = (int*) R_alloc(*n_peptide*sizeof(int));
	Omega_logit = (double*) R_alloc(*n_peptide*sizeof(double));
	pos_ind_by_pep = (int*) R_alloc(*n_peptide*sizeof(int));

	ProbSum = (double**) R_alloc(*n_indiv*sizeof(double*));

	for(p = 0; p < *n_position; p++)
	{
		xA[p] = (double*) R_alloc(NMAX*sizeof(double));
		xB[p] = (double*) R_alloc(NMAX*sizeof(double));
		for(c = 0; c < pnum[p]; c++)
		{
			pos_ind_by_pep[pstart[p] + c] = k;
		}
		k++;
	}

	for(i = 0; i < *n_indiv; i++)
	{
		Gamma[i] = (int*) R_alloc(*n_peptide*sizeof(int));
		ProbSum[i] = (double*) R_alloc(*n_peptide*sizeof(double));
	}

	Alpha_pep = (double*) R_alloc(*n_peptide*sizeof(double));
	Sig2 = (double*) R_alloc(*n_peptide*sizeof(double));
	A = (double*) R_alloc(*n_position*sizeof(double));
	B = (double*) R_alloc(*n_position*sizeof(double));
	U = (double*) R_alloc(*n_position*sizeof(double));
	Mu = (double*) R_alloc(*n_indiv*sizeof(double));
	double* likworkspace = R_alloc((*n_peptide)*(*n_indiv)*sizeof(double));

	// check whether our data is censored at all,
	// if so prepare for augmentation.

	if(*cen_num > 0)
	{
		D = (double*) R_alloc(*cen_num*sizeof(double));
		for(i = 0; i < *cen_num; i++)
		{
			D[i] = 0.0;
		}
	}

	FILE *AFILE, *BFILE, *PFILE, *VARFILE, *Sig2FILE, *MUFILE, *DFILE, *OFILE;
	FILE *ALPHAFILE;
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

	//statically allocated struct for Adaptive Rejection Sampling
	ARS_workspace workspace;

	int th_id, pos;
	RngStream rng[*nP];
	for(i = 0; i < *nP; i++)
	{
		rng[i] = RngStream_CreateStream("");
	}

	GetRNGstate();
	initialize_chain(ProbSum, Exprs, Y, W, Omega_ind, Omega_logit, Alpha_pep,
			Gamma, Sig2, Mu, A, B, U,
			n_position, pstart, pnum, n_peptide, n_indiv,
			xA, xB);
	PutRNGstate();
	int pep;
	Rprintf("parameters initialized \n");
	for(i = 0; i <= (*n_iter)*(*n_sweep) + *n_burn; i++)
	{
		R_CheckUserInterrupt();

#pragma omp parallel private(th_id, workspace, pos, sigma) num_threads(*nP)
		{
			th_id = omp_get_thread_num();
			// update gammas, alphas, omegas
#pragma omp for
			for(pep = 0; pep < *n_peptide; pep++)
			{
				p = pos_ind_by_pep[pep];
				update_peptide(Exprs[pep], Mu, Alpha_pep + pep,
						W[pep], Omega_ind + pep, Omega_logit + pep, Gamma[pep],
						Sig2 + pep, alpha, beta,
						U[p], A[p], B[p],
						m, zeta, *n_indiv, rng[th_id]);
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
						p, n_indiv, pnum, pstart, rng[th_id],
						xA[p], xB[p], &workspace);
			}

			update_global_parms(Alpha_pep, Gamma, &m, &zeta, &kappa, Mu, A, B, Sig2,
					m_0, v_0, alpha_0, beta_0, &alpha, &beta,
					*n_position, *n_peptide, *n_indiv,
					rng[0], xAlpha);

			update_dof_integrated(&dof, Exprs, Alpha_pep, Gamma,
					Sig2, Mu, likworkspace, *n_indiv, *n_peptide, rng[0]);

		}

		// check whether we need to update complete data
		if(*cen_num > 0 & i > *n_burn/2)
		{
			update_data(W, D, *cen_num, cen_ind, cen_pep, Y, Exprs,
					Gamma, Alpha_pep, Mu, *n_peptide, Sig2, cen_pos, rng[0]);
		}

		if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ))
		{
			n = n + 1.0;
			update_prob_include(n_peptide, n_indiv, Gamma, ProbSum, mean_fitted,
					Alpha_pep, Mu, n);
			if(*write == 1)
			{
				store_mcmc_output(Alpha_pep, Mu, A, B, U, Sig2, D, Omega_logit,
						Omega_ind, kappa, alpha, beta, m, zeta,
						dof, n_peptide, n_indiv, n_position, cen_num,
						AFILE, BFILE, PFILE, VARFILE, Sig2FILE, MUFILE, DFILE,
						OFILE, ALPHAFILE);
			}
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
		double Mu, double Sig2, double Alpha_pep,
		double m0, double v0, double alpha0, double beta0,
		int n_indiv, int n_peptide, RngStream rng, adpt m_tune, adpt c_tune,
		ARS_workspace ws)
{
	//update m

}

void update_peptide(double *Y, double *Mu, double *Alpha_pep, double *W,
		int *Omega_ind, double *Omega_logit, int *Gamma,
		double *Sig2_pep, double alpha, double beta,
		double u, double a, double b,
		double m, double zeta, int n_indiv, RngStream rng)
{
	int i, S = 0;
	double w0 = 0.0, w1 = 0.0, s0 = 0.0, s1 = 0.0;
	double s_init = 0.0, w_init = 0.0;
	double m0, m1, v0, v1;
	double delta, w, C, cur;
	double s2_pep = *Sig2_pep;

	if(*Omega_ind == 1)
	{
		// initialize sums for updates
		for(i = 0; i < n_indiv; i++)
		{
			s_init += (Gamma[i] == 1) ? (Y[i] - Mu[i])*W[i] : 0.0;
			w_init += (Gamma[i] == 1) ? W[i] : 0.0;
		}

		// loop through subjects, updating cluster membership
		for(i = 0; i < n_indiv; i++)
		{
			delta = Y[i] - Mu[i];
			w = W[i];
			w0 = (Gamma[i] == 1) ? w_init - w : w_init;
			w1 = w0 + w;
			s0 = (Gamma[i] == 1) ? s_init - delta*w : w_init;
			s1 = s0 + delta*w;

			m0 = (m*s2_pep + s0*zeta)/(w0*zeta + s2_pep);
			m1 = (m*s2_pep + s1*zeta)/(w1*zeta + s2_pep);
			v0 = (zeta*s2_pep)/(w0*zeta + s2_pep);
			v1 = (zeta*s2_pep)/(w1*zeta + s2_pep);

			C = .5*log(v1/v0) + .5*m1*m1/v1 - .5*m0*m0/v0 +
					pnorm(0.0, m1, sqrt(v1), 0, 1) - pnorm(0.0, m0, sqrt(v0), 0, 1);

			// update gamma and get ready for next iteration
			cur = logit(RngStream_U01(rng)) < C + *Omega_logit;
			Gamma[i] = cur;
			S += cur;
			w_init = (cur == 1) ? w1 : w0;
			s_init = (cur == 1) ? s1 : s0;
		}
	}

	// update omega, first check whether omega_{cp} > 0:
	if(S == 0)
	{
		double frac = exp(gsl_sf_lnbeta(a, b + n_indiv) - gsl_sf_lnbeta(a, b));
		*Omega_ind = (RngStream_U01(rng) < u/(u + (1 - u)*frac)) ? 0 : 1;
	}
	// if omega_{cp} > 0, update this.
	if(*Omega_ind == 1)
	{
		*Omega_logit = RngStream_LogitBeta(a + S, b + n_indiv - S, rng);
	}

	// update alpha
	double s_alpha, v_alpha;
	s_alpha = s_init/s2_pep + m/zeta;
	v_alpha = w_init/s2_pep	 + 1.0/zeta;
	*Alpha_pep = truncNorm(s_alpha/v_alpha, 1/v_alpha, rng);

	// update scale sigma^2_{cp}
	double SS = 0.0;
	for(i = 0; i < n_indiv; i++)
	{
		SS += W[i]*gsl_pow_2(Y[i] - Mu[i] - Gamma[i]*(*Alpha_pep));
	}
	SS = SS*.5 + beta;
	*Sig2_pep = beta/RngStream_GA1(alpha + n_indiv/2.0, rng);

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

void store_mcmc_output(double *Alpha, double *Mu, double *A, double *B, double *U,
		double *Sig2, double *D, double *Omega_Logit, int* Omega_Ind, double kappa,
		double alpha, double beta, double m, double zeta, double dof,
		int *n_peptide, int *n_indiv, int *n_position, int *cen_num,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *OFILE, FILE *ALPHAFILE)
{
	int p, i, k;

	for(p = 0; p < *n_position; p++)
	{
		fprintf(AFILE, "%.6lf \t", A[p]);
		fprintf(BFILE, "%.6lf \t", B[p]);
		fprintf(PFILE, "%.6lf \t", U[p]);
		fprintf(MUFILE, "%.6lf \t", Mu[i]);
	}

	for(p = 0; p < *n_peptide; p++)
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
	fprintf(VARFILE, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \n",
			m, zeta, kappa, alpha, beta, dof);
	return;
}

void update_data(double **W, double *D, int cen_num, int* cen_ind,
		int* cen_pep, double *Y, double **Exprs,
		int **Gamma, double *Alpha, double *Mu, int n_peptide, double *Sig2, int *cen_pos,
		RngStream rng)
{
	int i, p, k, pos;
	double tmp, weight;
	for(k = 0; k < cen_num; k++)
	{
		// retrieve indices
		i = cen_ind[k];
		p = cen_pep[k];
		weight = W[p][i];

		tmp = truncNorm(Mu[p] + Alpha[p]*Gamma[i][p] - Y[i*n_peptide + p], Sig2[p]/weight,
				rng);
		Exprs[i][p] = Y[i*n_peptide + p] + tmp;
		D[k] = tmp;
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
			mean_fitted[*n_peptide*i + p] = mean_fitted[*n_peptide*i + p]*((n-1.0)/n) + (Alpha[p]*Gamma[i][p] + Mu[p])/n;
		}
	}
	return;
}


void initialize_chain(double **ProbSum, double **Exprs, double *Y,
		double **W, int *Omega_Ind, double *Omega_Logit,
		double *Alpha_pep, int **Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *U,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		double **xA, double **xB)
{
	int c, p, i, g, pep;
	double u;
	for(i = 0; i < *n_indiv; i++)
	{
		Mu[i] = 0.0;
		for(pep = 0; pep < *n_peptide; pep++)
		{
			Exprs[i][pep] = Y[i*(*n_peptide) + pep];
			ProbSum[i][pep] = 0.0;
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
			Sig2[pep] = .6324555;
			Omega_Ind[pep] = 1;
			Omega_Logit[pep] = log(rgamma(A[p], 1.0)/rgamma(B[p], 1.0));

			for(i = 0; i < *n_indiv; i++)
			{
				W[pep][i] = 1.0;
				u = runif(0.0, 1.0);
				Gamma[i][pep] = 0;
			}
		}
	}
	return;
}

void update_dof_integrated(double *dof, double **Exprs, double *Alpha, int **Gamma,
		double *Sig2, double *Mu, double *workspace,
		int n_indiv, int n_peptide,
		RngStream rng)
{
	double V[] = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 32.0, 64.0};

	int i,p,k, nv = 9;
	const double n_obs = (double)(n_indiv*n_peptide);
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

#pragma omp parallel for private(p, ds)
	// compute the vicious log-sum

	for(i = 0; i < n_indiv; i++)
	{
		int *gam_i = Gamma[i];
		double *y_i = Exprs[i];
		for(p = 0; p < n_peptide; p++)
		{
			const double y = y_i[p];
			const double mu = Mu[p];
			const int gam = gam_i[p];
			const double alpha = (gam == 1) ? Alpha[p]:0.0;
			const double sig2 = Sig2[p];

			ds = gsl_pow_2(y - mu - alpha);
			ds = ds / (2.0 * sig2);

			workspace[i*n_peptide + p] = ds;
		}
	}

#pragma omp parallel for private(lik_temp, i, v, ds)
	for(k = 0; k < nv; k++)
	{
		v = V[k];
		lik_temp = 0.0;
		for(i = 0; i < (n_indiv*n_peptide); i++)
		{
			lik_temp = lik_temp + log(workspace[i] + v/2.0);
		}
		log_lik[k] = lik_temp;
	}

	double lik_max;
	// compute normalized probabilities on log scale
	for(i = 0; i < nv; i++)
	{
		v = V[i];
	//	Rprintf("%.0lf ", log_lik[i]);
		log_lik[i] = n_obs*(v/2.0)*log(v/2.0) + n_obs*lgamma(.5 + v/2.0) -
				n_obs*lgamma(v/2.0) - (v/2.0 + .5)*log_lik[i];
		lik_max = (i == 0) ? log_lik[i] : lik_max;
		lik_max = (lik_max > log_lik[i]) ? lik_max : log_lik[i];

	}
	//Rprintf("\n");

	for(i = 0; i < nv; i++)
	{
		v = V[i];
		log_lik[i] = log_lik[i] - lik_max;
		lik_sum[i] = (i == 0) ? log_lik[i] : log_apb(log_lik[i], lik_sum[i-1]);
		//Rprintf("%.5lf ", log_lik[i]);
	}
	//Rprintf("\n");

	for(i = 0; i < nv; i++)
	{
		log_cprob[i] = lik_sum[i] - lik_sum[nv - 1];
		//Rprintf("%.2lf ", log_cprob[i]);
	}
	//Rprintf("\n \n");
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


void update_var_pars(double *Sig2, double *alpha, double *beta,
		double alpha_0, double beta_0,  double *xAlpha,
		int n_max, int ind_pep, int n_peptide,
		int n_position, RngStream rng)
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
			&workspace, rng, lc_alpha, lcp_alpha);
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

void update_indiv_mu(double **Exprs, double **W,
		double *Alpha_pep, double *mu_j, int **Gamma,
		double *Sig2, double kappa,
		int n_peptide, int j, RngStream rng)
{
	// update Mu_j
	double S = 0.0, V = 0.0, w, ssq_pep, r;
	int pep;
	for(pep = 0; pep < n_peptide; pep++)
	{
		ssq_pep = Sig2[pep];
		w = W[pep][j];
		r = w/ssq_pep;
		S += (Exprs[pep][j] - Gamma[pep][j]*Alpha_pep[pep])*r;
		V += r;
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
	double frac, R;

	int c, S_p = 0;
	int num_x = 2;
	int p_begin = pstart[p];
	double argvec[3];

	// update A's and B's
	for(c = 0; c < pnum[p]; c++)
	{
		if(Omega_Ind[p_begin + c] == 1)
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
	U_p = RngStream_LogitBeta(a_0 + (double)(pnum[p] - S_p), b_0 + (double)(S_p), rng);
	frac = gsl_sf_lnbeta(*A_p, n_indiv + *B_p)- gsl_sf_lnbeta(*A_p, *B_p);
	frac = exp(frac);

	return;
}

void update_peptide_p(double *W, double **Exprs, double *Alpha, int **Gamma, int *Omega_Ind,
		double *Omega_Logit, double *Mu, double Sig2,
		double lambda, double kappa, double dof,
		int pep, int n_indiv, int pos, RngStream rng)
{
	int i, cur;
	double zi, M_cp = 0.0;
	double alpha = Alpha[pep];
	double alpha_sqd = alpha*alpha;
	double wi, weightsum = 0.0, weightsum_alpha = 0.0;
	double wi_sig2;
	double m_prime, v_prime;
	double one_prob, denom, num;
	const double log_w = log_from_logit(Omega_Logit[pep]);
	const double w = expit(Omega_Logit[pep]);
	const double m1w = m1expit(Omega_Logit[pep]);
	double mu = Mu[pep];

	if(Omega_Ind[pep] == 1)
	{
		for(i = 0; i < n_indiv; i++)
		{
			zi = (Exprs[i][pep] - mu);
			wi = W[i];
			wi_sig2 = wi/Sig2;
			weightsum += wi;

			denom = log1p(w*expm1(zi*alpha*wi_sig2 - .5*alpha_sqd*wi_sig2));
			num = wi_sig2*zi*alpha -.5*alpha_sqd*wi_sig2 + log_w;
			one_prob = num - denom;

			cur = (log(RngStream_RandU01(rng)) <= one_prob) ? 1 : 0;
			Gamma[i][pep] = cur;
			if(cur == 1)
			{
				M_cp += wi*zi;
				weightsum_alpha += wi;
			}
		}
	}

	else
	{
		for(i = 0; i < n_indiv; i++)
		{
			wi = W[i];
			Gamma[i][pep] = 0;
			weightsum += wi;
		}
	}

	double S = 0.0;
	double V = Sig2/(weightsum + Sig2*kappa);

	for(i = 0; i < n_indiv; i++)
	{
		S += (Exprs[i][pep] - Gamma[i][pep]*alpha)*W[i];
	}
	S = S/(weightsum + Sig2*kappa);

	mu = RngStream_N01(rng)*sqrt(V) + S;
	Mu[pep] = mu;

	for(i = 0; i < n_indiv; i++)
	{
		zi = (Exprs[i][pep] - alpha*Gamma[i][pep] - mu);
		zi = gsl_pow_2(zi)/Sig2;
		W[i] = RngStream_GA1((dof + 1.0)/2.0, rng)/(dof/2.0 + zi/2.0);
	}
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

double dof_likelihood(double **W, int n_indiv, int n_peptide)
{
	double S = 0.0;
	double wip;
	int i, p;

#pragma omp parallel for private(wip, i) reduction(+:S) schedule(static, 1)
		for(p = 0; p < n_peptide; p++)
		{
			double *W_pep = W[p];
			for(i = 0; i < n_indiv; i++)
			{
				wip = W_pep[i];
				S = S + log(wip) - wip;
			}
		}

	return(.5*S);
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
		double out = log(x) - log1p(-x);
		return(out);
	}
}
