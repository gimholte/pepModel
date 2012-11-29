/*
 * pepBayes_v5.c
 *
 *  Created on: Nov 29, 2012
 *      Author: hoblitz
 */

#include <omp.h>
#include <float.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "RngFunctions.h"
#include "RngStream.h"
#include "ARS.h"

#define	TEST_INT_LENGTH 200
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#include <Rinterface.h>

#include "pepBayes_v5.h"

/*
 * individual specific mean, with prior N(0, kappa)
 * also includes peptide specific means
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

void pepbayes_v5(double *Y, double *hyper_param, int *pstart,
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
	double *Mu_pep;				// peptide specific means
	double *Mu_ind;				// individual specific means
	double kappa = 10.0;	// prior mean precision
	double psi = 10.0;		// prior ind-mean precision
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
	FILE *ALPHAFILE, *MUINDFILE;

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
	Mu_pep = (double*) R_alloc(*n_peptide, sizeof(double));
	Mu_ind = (double*) R_alloc(*n_indiv, sizeof(double));
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
		MUINDFILE = fopen("muindfile.txt", "w");
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
			Gamma, Sig2, Mu_pep, Mu_ind, A, B, U,
			n_position, pstart, pnum, n_peptide, n_indiv,
			xA, xB, RB);
	PutRNGstate();


	Rprintf("parameters initialized \n");
	for(i = 1; i <= total_iterations; i++)
	{
		R_CheckUserInterrupt();
		update_dof_integrated(&dof, Exprs, W, Alpha_pep, Gamma,
				Sig2, Mu_pep, Mu_ind, likworkspace, *n_indiv, *n_peptide, rng[0]);

#pragma omp parallel private(th_id, workspace, pos) num_threads(*nP)
		{
			th_id = omp_get_thread_num();
#pragma omp for
			for(pep = 0; pep < *n_peptide; pep++)
			{
				pos = pos_ind_by_pep[pep];
				update_peptide(Exprs, Mu_pep + pep, Mu_ind, Alpha_pep + pep,
						W, Omega_ind + pep, Omega_logit + pep, Gamma,
						Sig2 + pep, alpha, beta,
						U[pos], A[pos], B[pos],
						m, zeta, dof, kappa, *n_indiv, *n_peptide, pep, rng[th_id], RB);
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
#pragma omp for
			for(j = 0; j < *n_indiv; j++)
			{
					update_indiv_mu(Exprs,
					W, Mu_pep, Alpha_pep,
					Mu_ind + j, Gamma, Sig2, psi,
					*n_peptide, *n_indiv, j, rng[th_id]);
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

		update_global_params(&m, &zeta, &alpha, &beta, &kappa, &psi,
				Mu_pep, Mu_ind, Sig2, Alpha_pep, Omega_ind,
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
					Gamma, Alpha_pep, Mu_pep, *n_indiv, Sig2, rng[0]);
		}

		if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ) && *write == 1)
		{
			store_mcmc_output(Alpha_pep, Mu_pep, Mu_ind, A, B, U, Sig2, D, Omega_logit,
					Omega_ind, kappa, psi, alpha, beta, m, zeta,
					dof, *n_peptide, *n_indiv, *n_position, cen_num,
					AFILE, BFILE, PFILE, VARFILE, Sig2FILE, MUFILE,
					MUINDFILE, DFILE,
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
		fclose(MUINDFILE);
		if(*cen_num > 0)
		{
			fclose(DFILE);
		}
	}
	return;
}

void update_indiv_mu(double* restrict Exprs,
		double* restrict W,
		double* restrict Mu_pep,
		double* restrict Alpha_pep,
		double* restrict mu_j,
		int* restrict Gamma,
		double* restrict Sig2,
		double psi,
		int n_peptide, int n_indiv, int j, RngStream rng)
{
	// update Mu_j
	double S = 0.0, V = 0.0, w, ssq_pep, r;
	int pep, d_ind = j;
	for(pep = 0; pep < n_peptide; pep++)
	{
		ssq_pep = Sig2[pep];
		w = W[d_ind];
		r = w/ssq_pep;
		S += (Exprs[d_ind] - Mu_pep[pep] - Gamma[d_ind]*Alpha_pep[pep])*r;
		V += r;
		d_ind += n_indiv;
	}
	V += psi;
	*mu_j = RngStream_N01(rng)/sqrt(V) + S/V;
	return;
}

void update_global_params(double* restrict m,
		double* restrict zeta,
		double* restrict alpha,
		double* restrict beta,
		double* restrict kappa,
		double* restrict psi,
		double* restrict Mu,
		double* restrict Mu_ind,
		double* restrict Sig2,
		double* restrict Alpha_pep,
		int* restrict Omega_ind,
		double m0, double v0, double alpha0, double beta0,
		int n_indiv, int n_peptide, RngStream rng, adpt *m_tune, adpt *z_tune,
		ARS_workspace *ws, double *xAlpha, int nP)
{
	double v_prime;
	double m_prop, zeta_prop;
	double log_DF_ratio;

	int pep, i;
	double sum_alpha = 0.0;
	int n_alpha = 0;
	double sum_u_sqr = 0.0;
	double sum_sig2_inv = 0.0;
	double sum_log_sig2 = 0.0;
	double SS_alpha = 0.0;
	double sum_muind_sqr = 0.0;
	int o_ind;
	double a_pep, s2_pep, mu_pep;

	// compute numerous quanities
#pragma omp parallel for num_threads(nP) \
		default(shared) private(pep, a_pep, s2_pep, mu_pep, o_ind) \
		reduction(+:sum_u_sqr, sum_sig2_inv, sum_log_sig2, sum_alpha, n_alpha, SS_alpha)
	for(pep = 0; pep < n_peptide; pep++)
	{
		mu_pep = Mu[pep];
		a_pep = Alpha_pep[pep];
		s2_pep = Sig2[pep];
		o_ind = Omega_ind[pep];

		sum_u_sqr = sum_u_sqr + mu_pep*mu_pep;
		sum_sig2_inv = sum_sig2_inv + 1.0/s2_pep;
		sum_log_sig2 = sum_log_sig2 + log(s2_pep);
		n_alpha += o_ind;
		sum_alpha = o_ind*a_pep + sum_alpha;
		SS_alpha = SS_alpha + o_ind*R_pow_di(*m - a_pep, 2);
	}

	for(i = 0; i < n_indiv; i++)
	{
		sum_muind_sqr += Mu_ind[i]*Mu_ind[i];
	}

	sum_muind_sqr += .5;
	sum_u_sqr += .5;
	sum_sig2_inv += beta0;
	SS_alpha += 10.0;

	// update zeta
	double d_prime;
	// log-scale random walk proposal
	d_prime = n_alpha/2.0 + 5.0;
	zeta_prop = *zeta*exp(z_tune->tune*(RngStream_RandU01(rng) - .5));
	log_DF_ratio = pnorm(0.0, (*m), sqrt(*zeta), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, (*m), sqrt(zeta_prop), 0, 1);
	log_DF_ratio = n_alpha*log_DF_ratio;
	log_DF_ratio += .5*(1.0/(*zeta) - 1.0/zeta_prop)*SS_alpha;
	log_DF_ratio += d_prime*(log(*zeta) - log(zeta_prop));
	// accept or reject step
	if(log(RngStream_RandU01(rng)) <= log_DF_ratio)
	{
		*zeta = zeta_prop;
		z_tune->count++;
	}

	// update m
	v_prime = (*zeta*v0)/(*zeta + v0*n_alpha);
	m_prop = *m + (RngStream_RandU01(rng)*2.0 - 1.0)*sqrt(v_prime)*m_tune->tune; //sqrt(v_prime)*rnorm(0.0,1.0);
	log_DF_ratio = pnorm(0.0, *m, sqrt(*zeta), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, m_prop, sqrt(*zeta), 0, 1);
	log_DF_ratio = n_alpha*log_DF_ratio;
	log_DF_ratio += -(.5/v_prime)*(m_prop*m_prop - (*m)*(*m));
	log_DF_ratio += ((*zeta*m0 + v0*sum_alpha)/(*zeta*v0))*(m_prop - *m);

	if(log(RngStream_RandU01(rng)) <= log_DF_ratio)
	{
		*m = m_prop;
		m_tune->count = m_tune->count + 1;
	}

	//update alpha and kappa and psi
	*kappa = RngStream_GA1(n_peptide/2.0 + .5, rng)/sum_u_sqr;
	*psi = RngStream_GA1(n_indiv/2.0 + .5, rng)/sum_muind_sqr;
	double argvec[4];
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

void update_peptide(double* restrict Exprs,
		double* restrict Mu_pep,
		double* restrict Mu_ind,
		double* restrict Alpha_pep,
		double* restrict W,
		int* restrict Omega_ind,
		double* restrict Omega_logit,
		int* restrict Gamma,
		double* restrict Sig2_pep,
		double alpha, double beta,
		double u_logit, double a, double b,
		double m, double zeta, double dof, double kappa,
		int n_indiv, int n_peptide, int pep,
		RngStream rng, double* restrict RB)
{
	int i, S = 0;
	double s_init = 0.0, w_init = 0.0;
	double w, C;
	int cur;
	double s2_pep = *Sig2_pep;
	double a_pep = *Alpha_pep;
	double mu_pep = *Mu_pep;
	int d_ind;
	double SS = 0.0;
	double delta0, delta1;
	const double w1 = .5 + dof/2.0;
	const double w2 = dof/2.0;
	double tmp;

	d_ind = pep*n_indiv; // start indexing through large arrays

	if(*Omega_ind == 1)
	{
		for(i = 0; i < n_indiv; i++)
		{
			delta0 = Exprs[d_ind] - mu_pep - Mu_ind[i];
			delta1 = delta0 - a_pep;

			C = (dof/2.0 + .5)*(log(1.0 + delta0*delta0/(s2_pep*dof)) -
					log(1.0 + delta1*delta1/(s2_pep*dof)));
			RB[d_ind] = expit(C + *Omega_logit);

			cur = (RngStream_RandU01(rng) < RB[d_ind]) ? 1 : 0;
			Gamma[d_ind] = cur;

			tmp = (delta0 - a_pep*cur)*(delta0 - a_pep*cur)/s2_pep;
			w = RngStream_GA1(w1, rng)/(w2 + tmp/2.0);
			W[d_ind] = w;

			//quantities for alpha, omega, sigma updates
			S += cur;
			s_init += cur*delta0*w;
			w_init += cur*w;
			SS += w*tmp*s2_pep;

			d_ind++;
		}
	}
	// if here, all gamma's for peptide p are zero
	else
	{
		for(i = 0; i < n_indiv; i++)
		{
			delta0 = Exprs[d_ind] - mu_pep - Mu_ind[i];
			tmp = R_pow_di(delta0, 2);
			w = RngStream_GA1(w1, rng)/(w2 + tmp/(2.0*s2_pep));
			W[d_ind] = w;
			RB[d_ind] = 0.0;

			//quantities for sigma updates
			SS += w*tmp;
			d_ind++;
		}
	}
	// update scale sigma^2_{cp}
	SS = SS*.5 + beta;
	s2_pep = SS/RngStream_GA1(alpha + n_indiv/2.0, rng);
	*Sig2_pep = s2_pep;

	// update omega, first check whether omega_{cp} > 0:
	if(S == 0)
	{
		double log_frac = lbeta(a, b + n_indiv) - lbeta(a, b);
		*Omega_ind = (logit(RngStream_RandU01(rng)) < (u_logit - log_frac)) ? 0 : 1;
	}

	// if omega_{cp} > 0, update this.
	if(*Omega_ind == 1)
	{
		*Omega_logit = RngStream_LogitBeta(a + S, b + n_indiv - S, rng);
		// update alpha
		double s_alpha, v_alpha;
		s_alpha = s_init/s2_pep + m/zeta;
		v_alpha = w_init/s2_pep	+ 1.0/zeta;
		a_pep = truncNorm(s_alpha/v_alpha, 1.0/v_alpha, rng);
		*Alpha_pep = a_pep;
	}

	double M = 0.0, V = 0.0;
	d_ind = pep*n_indiv;
	for(i = 0; i < n_indiv; i++)
	{
		w = W[d_ind];
		M += (Exprs[d_ind] - Mu_ind[i] - Gamma[d_ind]*a_pep)*w;
		V += w;
		d_ind ++;
	}
	V = (V/s2_pep) + kappa;
	M = M/s2_pep;
	*Mu_pep = RngStream_N01(rng)/sqrt(V) + M/V;

	return;
}

void store_mcmc_output(double *Alpha, double *Mu, double *Mu_ind,
		double *A, double *B, double *U,
		double *Sig2, double *D, double *Omega_Logit, int* Omega_Ind, double kappa,
		double psi, double alpha, double beta, double m, double zeta, double dof,
		int n_peptide, int n_indiv, int n_position, int *cen_num,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *MUINDFILE, FILE *DFILE, FILE *OFILE, FILE *ALPHAFILE)
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
		fprintf(MUINDFILE, "%.6lf\t", Mu_ind[i]);
	}
	for(p = 0; p < n_peptide; p++)
	{
		fprintf(MUFILE, "%.6lf \t", Mu[p]);
		fprintf(OFILE, "%.6lf \t", (Omega_Ind[p] == 1) ? expit(Omega_Logit[p]):0.0);
		if(Omega_Ind[p] == 1)
		{
			fprintf(ALPHAFILE, "%.6lf \t", Alpha[p]);
		}
		else
		{
			fprintf(ALPHAFILE, "NA \t");
		}
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
	fprintf(MUINDFILE, "\n");
	fprintf(OFILE, "\n");
	fprintf(Sig2FILE, "\n");
	fprintf(MUFILE, "\n");
	fprintf(VARFILE, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t\n",
			m, zeta, kappa, psi, alpha, beta, dof);
	return;
}

void update_censoring(double *W, double *D, int cen_num, int* cen_ind,
		int* cen_pep, double *Y, double *Exprs,
		int *Gamma, double *Alpha, double *Mu, int n_indiv, double *Sig2,
		RngStream rng)
{
	int i, p, k, ind;
	double tmp, weight;
	for(k = 0; k < cen_num; k++)
	{
		// retrieve indices
		i = cen_ind[k];
		p = cen_pep[k];
		ind = i + p*n_indiv;
		weight = W[ind];

		tmp = truncNorm(Mu[p] + Alpha[p]*Gamma[ind] - Y[ind], Sig2[p]/weight,
				rng);
		Exprs[ind] = Y[ind] + tmp;
		D[k] = tmp;
	}
	return;
}


void initialize_chain(double *ProbSum, double *Exprs, double *Y,
		double *W, int *Omega_Ind, double *Omega_Logit,
		double *Alpha_pep, int *Gamma,
		double *Sig2, double *Mu_pep, double *Mu_ind, double *A, double *B, double *U,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		double **xA, double **xB, double *RB)
{
	int c, p, i, g, pep;
	double prob;
	for(pep = 0; pep < *n_peptide; pep++)
	{
		for(i = 0; i < *n_indiv; i++)
		{
			Exprs[(*n_indiv)*pep + i] = Y[(*n_indiv)*pep + i];
			ProbSum[(*n_indiv)*pep + i] = 0.0;
			RB[(*n_indiv)*pep + i] = 0.0;
		}
	}

	for(i = 0; i < *n_indiv; i++)
	{
		Mu_ind[i] = 0.0;
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
			Mu_pep[pep] = 0.0;
			Omega_Ind[pep] = 1;
			Omega_Logit[pep] = log(rgamma(A[p], 1.0)/rgamma(B[p], 1.0));

			for(i = 0; i < *n_indiv; i++)
			{
				prob = R_pow_di(Y[(*n_indiv)*pep + i] - Mu_pep[pep],2)/
						(R_pow_di(Y[(*n_indiv)*pep + i] - Mu_pep[pep],2) +
						R_pow_di(Y[(*n_indiv)*pep + i] - Mu_pep[pep] - Alpha_pep[pep],2));
				W[(*n_indiv)*pep + i] = 1.0;
				Gamma[(*n_indiv)*pep + i] = (int)(unif_rand() <= prob);
			}
		}
	}
	return;
}

void update_dof_integrated(double* restrict dof,
		double* restrict Exprs,
		double* restrict W,
		double* restrict Alpha,
		int* restrict Gamma,
		double* restrict Sig2,
		double* restrict Mu_pep,
		double* restrict Mu_ind,
		double* restrict workspace,
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

	for(k = 0; k < nv; k++)
	{
		log_lik[k] = 0.0;
	}

	int d_ind;
	double tmp2 = 0.0, tmp4 = 0.0, tmp6 = 0.0, tmp8 = 0.0, tmp10 = 0.0,
			tmp12 = 0.0, tmp16 = 0.0, tmp32 = 0.0, tmp64 = 0.0;

#pragma omp parallel for private(ds, d_ind, i) \
	reduction(+:tmp2,tmp4,tmp6,tmp8,tmp10,tmp12,tmp16,tmp32,tmp64)
	// compute the vicious log-sum
	for(pep = 0; pep < n_peptide; pep++)
	{
		d_ind = pep*n_indiv;
		const double mu = Mu_pep[pep];
		const double sig2 = Sig2[pep];
		const double alpha = Alpha[pep];
		for(i = 0; i < n_indiv; i++)
		{
			const double y = Exprs[d_ind];
			const double mu_ind = Mu_ind[i];
			const int gam = Gamma[d_ind];

			ds = (y - mu - mu_ind - alpha*gam)*(y - mu - mu_ind - alpha*gam)/sig2;

			tmp2 += log(1.0 + ds/2);
			tmp4 += log(1.0 + ds/4);
			tmp6 += log(1.0 + ds/6);
			tmp8 += log(1.0 + ds/8);
			tmp10 += log(1.0 + ds/10);
			tmp12 += log(1.0 + ds/12);
			tmp16 += log(1.0 + ds/16);
			tmp32 += log(1.0 + ds/32);
			tmp64 += log(1.0 + ds/64);

			//workspace[d_ind] = ds;
			d_ind ++;
		}
	}

	log_lik[0] = tmp2;
	log_lik[1] = tmp4;
	log_lik[2] = tmp6;
	log_lik[3] = tmp8;
	log_lik[4] = tmp10;
	log_lik[5] = tmp12;
	log_lik[6] = tmp16;
	log_lik[7] = tmp32;
	log_lik[8] = tmp64;

	double lik_max;
	// compute normalized probabilities on log scale
	for(i = 0; i < nv; i++)
	{
		v = V[i];
		log_lik[i] = -n_obs*log(v)/2.0 + n_obs*lgammafn(.5 + v/2.0) -
				n_obs*lgammafn(v/2.0) - (v/2.0 + .5)*log_lik[i];
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
			pz = exp(-.5*R_pow_di(z - a_star,2));
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

double lc_AB(double x, double* restrict argvec, int *arglen)
{
	return(-1.0*x*argvec[2] + argvec[0]*lgammafn(x + argvec[1]) - argvec[0]*lgammafn(x));
}

double lcp_AB(double x, double* restrict argvec, int *arglen)
{
	return(-1.0*argvec[2] + argvec[0]*digamma(x + argvec[1])
			- argvec[0]*digamma(x));
}

double lc_alpha_int(double x, double* restrict argvec, int *arglen)
{
	const double Np = argvec[1];
	return(-argvec[0]*x + lgammafn(Np*x + 1.0) - Np*lgammafn(x) - Np*x*argvec[2] -
			x*argvec[3]);
}

double lcp_alpha_int(double x, double* restrict argvec, int *arglen)
{

	const double Np = argvec[1];
	return(-argvec[0] + Np*digamma(Np*x + 1.0) - Np*digamma(x) - Np*argvec[2] -
			argvec[3]);
}

double lc_alpha(double x, double *argvec, int *arglen)
{
	return(-x*argvec[0] - argvec[1]*lgammafn(x));
}

double lcp_alpha(double x, double *argvec, int *arglen)
{
;
	return (-argvec[0] - argvec[1]*digamma(x));
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
		if(accept_prop > .25)
		{
			tune_set->tune = tune_set->tune + .1*(n_burn + 1.0 - i)/(n_burn);
		}
		else if(accept_prop < .18)
		{
			tune_set->tune = (tune_set->tune)*
					((i)/(n_burn + 1.0) + .7*(n_burn + 1.0 - i)/(1.0 + n_burn));
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
