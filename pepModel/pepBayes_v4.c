/*
 * pepBayes_v4.c
 *
 *  Created on: Nov 19, 2012
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
#include <R_ext/Lapack.h>

#define	TEST_INT_LENGTH 200
#define R_INTERFACE_PTRS 1
#define CSTACK_DEFNS 1
#define CHUNK 32
#include <Rinterface.h>

#include "pepBayes_v4.h"

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

void pepbayes_v4(double *Y, double *hyper_param, int *pstart,
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
	double *pnum_double;
	double *Q_work, *L_cur, *b_cur, *mu_cur, *L_prop, *b_prop, *mu_prop;
	double *ProbSum;
	int total_iterations = *n_burn+(*n_sweep)*(*n_iter);
	int percent_complete = 10;

	adpt zeta_adpt;
	adpt m_adpt;
	adpt phi_adpt;
	adpt *omega_adpt = (adpt*) R_alloc(*n_peptide, sizeof(adpt));

	zeta_adpt.total_count = 0;
	zeta_adpt.count = 0;
	zeta_adpt.tune = 1.0;
	m_adpt.total_count = 0;
	m_adpt.count = 0;
	m_adpt.tune = 5.0;
	phi_adpt.total_count = 0;
	phi_adpt.count = 0;
	phi_adpt.tune = 1.0;

	for(p = 0; p < *n_peptide; p++)
	{
		omega_adpt[p].total_count = 0;
		omega_adpt[p].count = 0;
		omega_adpt[p].tune = 4.0;
	}
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
	double *Nu;				// mean of position membership probabilities
	double *Omega_logit;	// membership probabilities, on logit scale (theta)
	double nu_0;			// prior mean for position 1.
	double phi = 1.0;		// between-position variance of thetas
	double tau = 1.0; 		// within-position variance of thetas

	// location-related parameters
	double *Mu;				// peptide specific means
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
	FILE *NUFILE, *VARFILE, *Sig2FILE, *MUFILE, *DFILE, *OFILE;
	FILE *ALPHAFILE;

	// retreive and initialize hyperparameter values
	alpha_0 = hyper_param[4];
	beta_0 = hyper_param[5];
	m_0 = hyper_param[6];
	v_0 = hyper_param[7];
	nu_0 = hyper_param[8];

	// begin memory allocation
	Gamma = (int*) R_alloc(*n_peptide*(*n_indiv), sizeof(int));
	Exprs = (double*) R_alloc(*n_peptide*(*n_indiv), sizeof(double));
	W = (double *) R_alloc(*n_peptide*(*n_indiv), sizeof(double));
	double *RB = (double *) R_alloc(*n_peptide*(*n_indiv), sizeof(double));

	// starting values for ARS of alpha.
	xAlpha = (double*) R_alloc(NMAX, sizeof(double));
	// initial values for hull quantiles.
	xAlpha[0] = 1.0;
	xAlpha[1] = 2.0;

	Omega_logit = (double*) R_alloc(*n_peptide, sizeof(double));
	pos_ind_by_pep = (int*) R_alloc(*n_peptide, sizeof(int));
	pnum_double = (double*) R_alloc(*n_position, sizeof(double));
	ProbSum = (double*) R_alloc(*n_peptide*(*n_indiv), sizeof(double));

	Q_work = (double*) R_alloc(*n_position*2, sizeof(double));
	L_cur = (double*) R_alloc(*n_position*2, sizeof(double));
	b_cur = (double*) R_alloc(*n_position, sizeof(double));
	mu_cur = (double*) R_alloc(*n_position, sizeof(double));
	L_prop = (double*) R_alloc(*n_position*2, sizeof(double));
	b_prop = (double*) R_alloc(*n_position, sizeof(double));
	mu_prop = (double*) R_alloc(*n_position, sizeof(double));

	for(p = 0; p < *n_position; p++)
	{
		b_cur[p] = 0.0;
		b_prop[p] = 0.0;
		L_cur[2*p] = 0.0;
		L_cur[2*p + 1] = 0.0;
		L_prop[2*p] = 0.0;
		L_prop[2*p + 1] = 0.0;
		mu_cur[p] = 0.0;
		mu_prop[p] = 0.0;
		Q_work[2*p] = 2.0;
		Q_work[2*p + 1] = -1.0;
		pnum_double[p] = (double)(pnum[p]);
		for(c = 0; c < pnum[p]; c++)
		{
			pos_ind_by_pep[pstart[p] + c] = k;
		}
		k++;
	}
	Q_work[2*(*n_position - 1)] = 1.0;

	Alpha_pep = (double*) R_alloc(*n_peptide, sizeof(double));
	Sig2 = (double*) R_alloc(*n_peptide, sizeof(double));
	Mu = (double*) R_alloc(*n_peptide, sizeof(double));
	Nu = (double*) R_alloc(*n_position, sizeof(double));

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
		NUFILE = fopen("nufile.txt", "w");
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
	initialize_chain(ProbSum, Exprs, Y, W, Omega_logit, Alpha_pep,
			Gamma, Sig2, Mu, Nu,
			n_position, pstart, pnum, n_peptide, n_indiv, RB);
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
				update_peptide(Exprs, Mu + pep, Alpha_pep + pep,
						W, Omega_logit + pep, Gamma,
						Sig2 + pep, alpha, beta,
						Nu[pos], tau,
						m, zeta, dof, kappa, *n_indiv, *n_peptide, pep, rng[th_id], RB,
						omega_adpt + pep);
			}
			if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ))
			{
#pragma omp for schedule(static, CHUNK)
				for(d = 0; d < (*n_indiv)*(*n_peptide); d++)
				{
					ProbSum[d] += RB[d];
				}
			}
		}

		update_Nu(Omega_logit,
				Nu, nu_0, tau, &phi,
				*n_position, pnum_double, pnum,
				Q_work, L_cur, b_cur, mu_cur,
				L_prop, b_prop, mu_prop, rng[th_id],
				&phi_adpt);

		update_global_params(&m, &zeta, &alpha, &beta, &kappa, &tau, &phi,
				Mu, Nu, Omega_logit, Sig2, Alpha_pep,
				m_0, v_0, alpha_0, beta_0, nu_0,
				*n_indiv, *n_peptide, *n_position, pos_ind_by_pep,
				rng[0], &m_adpt, &zeta_adpt,
				&workspace, xAlpha, *nP);

		if(i % TEST_INT_LENGTH == 0)
		{
			update_tuning(&m_adpt, i, *n_burn);
			update_tuning(&zeta_adpt, i, *n_burn);
			update_tuning(&phi_adpt, i, *n_burn);
			for(p = 0; p < *n_peptide; p++)
			{
				update_tuning(omega_adpt + p, i, *n_burn);
			}
		}

		// check whether we need to update complete data
		if((*cen_num > 0) & (i > *n_burn/2))
		{
			update_censoring(W, D, *cen_num, cen_ind, cen_pep, Y, Exprs,
					Gamma, Alpha_pep, Mu, *n_indiv, Sig2, rng[0]);
		}

		if((i > *n_burn) && ((i - *n_burn) % (*n_sweep) == 0 ) && *write == 1)
		{
			store_mcmc_output(Alpha_pep, Mu, Nu, Sig2, D, Omega_logit,
					kappa, alpha, beta, m, zeta, phi, tau,
					dof, *n_peptide, *n_indiv, *n_position, cen_num,
					NUFILE, VARFILE, Sig2FILE, MUFILE, DFILE,
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
	Rprintf("Phi acceptance rate: %.3lf\n",
			(double)(phi_adpt.total_count)/(*n_iter*(*n_sweep)), phi_adpt.tune);

	if(*write == 1)
	{
		Rprintf("closing files\n");
		fclose(NUFILE);
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

void update_global_params(double* restrict m,
		double* restrict zeta,
		double* restrict alpha,
		double* restrict beta,
		double* restrict kappa,
		double* restrict tau,
		double* restrict phi,
		double* restrict Mu,
		double* restrict Nu,
		double* restrict Omega_logit,
		double* restrict Sig2,
		double* restrict Alpha_pep,
		double m0, double v0, double alpha0, double beta0, double nu0,
		int n_indiv, int n_peptide, int n_position, int* pos_ind_by_pep,
		RngStream rng, adpt *m_tune, adpt *z_tune,
		ARS_workspace *ws, double *xAlpha, int nP)
{
	double v_prime;
	double m_prop, zeta_prop;
	double log_DF_ratio;

	int pep;
	double sum_alpha = 0.0;
	int n_alpha = n_peptide;
	double sum_u_sqr = 0.0;
	double sum_sig2_inv = 0.0;
	double sum_log_sig2 = 0.0;
	double SS_alpha = 0.0;
	double tau_SS = 0.0;

	int pos;
	double a_pep, s2_pep, mu_pep, omega_pep;

	// compute numerous quanities
#pragma omp parallel for num_threads(nP) \
		default(shared) private(pep, a_pep, s2_pep, mu_pep, omega_pep, pos) \
		reduction(+:sum_u_sqr, sum_sig2_inv, sum_log_sig2, sum_alpha, SS_alpha, tau_SS)
	for(pep = 0; pep < n_peptide; pep++)
	{
		mu_pep = Mu[pep];
		a_pep = Alpha_pep[pep];
		s2_pep = Sig2[pep];
		omega_pep = Omega_logit[pep];
		pos = pos_ind_by_pep[pep];

		sum_u_sqr = sum_u_sqr + mu_pep*mu_pep;
		sum_sig2_inv = sum_sig2_inv + 1.0/s2_pep;
		sum_log_sig2 = sum_log_sig2 + log(s2_pep);
		sum_alpha = a_pep + sum_alpha;
		SS_alpha = SS_alpha + R_pow_di(*m - a_pep, 2);
		tau_SS = tau_SS + (omega_pep - Nu[pos])*(omega_pep - Nu[pos]);
	}
	sum_u_sqr += .5;
	sum_sig2_inv += beta0;
	SS_alpha += 10.0;

	// update tau
	*tau = RngStream_GA1(n_peptide/2.0 + 2.0, rng)/(2.0 + .5*tau_SS);

	// update zeta
	double d_prime;
	double log_zeta = log(*zeta);

	// log-scale random walk proposal
	d_prime = n_alpha/2.0 + 5.0;
	double delta = z_tune->tune*(RngStream_RandU01(rng) - .5);
	zeta_prop = exp(log_zeta + delta);
	log_DF_ratio = pnorm(0.0, (*m), sqrt(*zeta), 0, 1);
	log_DF_ratio = log_DF_ratio - pnorm(0.0, (*m), sqrt(zeta_prop), 0, 1);
	log_DF_ratio = n_alpha*log_DF_ratio;
	log_DF_ratio += .5*(1.0/(*zeta) - 1.0/zeta_prop)*SS_alpha;
	log_DF_ratio += d_prime*(-1.0*delta);
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

	//update alpha and kappa
	*kappa = RngStream_GA1(n_peptide/2.0 + .5, rng)/sum_u_sqr;
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
		double* restrict Alpha_pep,
		double* restrict W,
		double* restrict Omega_logit,
		int* restrict Gamma,
		double* restrict Sig2_pep,
		double alpha, double beta,
		double nu_pos, double tau,
		double m, double zeta, double dof, double kappa,
		int n_indiv, int n_peptide, int pep,
		RngStream rng, double* restrict RB, adpt *omega_adpt)
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

	for(i = 0; i < n_indiv; i++)
	{
		delta0 = Exprs[d_ind] - mu_pep;
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

	// update scale sigma^2_{cp}
	SS = SS*.5 + beta;
	s2_pep = SS/RngStream_GA1(alpha + n_indiv/2.0, rng);
	*Sig2_pep = s2_pep;

	// update omega
	double omega_prop, omega_cur;
	double logratio, scale = omega_adpt->tune;

	// precision of proposal from theta -> theta*
	omega_cur = *Omega_logit;
	omega_prop = RngStream_N01(rng)*scale + omega_cur;

	logratio = .5*tau*((omega_cur - nu_pos)*(omega_cur - nu_pos) -
			(omega_prop - nu_pos)*(omega_prop - nu_pos));
	logratio += S*(omega_prop - omega_cur);
	logratio += n_indiv*(log1pexp(omega_cur) - log1pexp(omega_prop));

	if(log(RngStream_RandU01(rng)) <= logratio)
	{
		*Omega_logit = omega_prop;
		omega_adpt->count = omega_adpt->count + 1;
	}

	// update alpha
	double s_alpha, v_alpha;
	s_alpha = s_init/s2_pep + m/zeta;
	v_alpha = w_init/s2_pep	+ 1.0/zeta;
	a_pep = truncNorm(s_alpha/v_alpha, 1.0/v_alpha, rng);
	*Alpha_pep = a_pep;

	double M = 0.0, V = 0.0;
	d_ind = pep*n_indiv;
	for(i = 0; i < n_indiv; i++)
	{
		w = W[d_ind];
		M += (Exprs[d_ind] - Gamma[d_ind]*a_pep)*w;
		V += w;
		d_ind ++;
	}
	V = (V/s2_pep) + kappa;
	M = M/s2_pep;
	*Mu_pep = RngStream_N01(rng)/sqrt(V) + M/V;

	return;
}

void store_mcmc_output(double *Alpha, double *Mu, double *Nu,
		double *Sig2, double *D, double *Omega_Logit, double kappa,
		double alpha, double beta, double m, double zeta, double phi,
		double tau, double dof,
		int n_peptide, int n_indiv, int n_position, int *cen_num,
		FILE *NUFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *OFILE, FILE *ALPHAFILE)
{
	int p, i, k;

	for(p = 0; p < n_position; p++)
	{
		fprintf(NUFILE, "%.6lf \t", Nu[p]);
	}

	for(p = 0; p < n_peptide; p++)
	{
		fprintf(MUFILE, "%.6lf \t", Mu[p]);
		fprintf(OFILE, "%.6lf \t", Omega_Logit[p]);
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
	fprintf(NUFILE, "\n");
	fprintf(OFILE, "\n");
	fprintf(Sig2FILE, "\n");
	fprintf(MUFILE, "\n");
	fprintf(VARFILE, "%.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf \n",
			m, zeta, kappa, alpha, beta, phi, tau, dof);
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
		double *W, double *Omega_Logit,
		double *Alpha_pep, int *Gamma,
		double *Sig2, double *Mu, double *Nu,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		double *RB)
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

	for(p = 0; p < *n_position; p++)
	{
		Nu[p] = 0.0;
		for(c = 0; c < pnum[p]; c++)
		{
			pep = pstart[p] + c;

			Alpha_pep[pep] = 5.0 + unif_rand();
			Sig2[pep] = 1.0;
			Mu[pep] = 0.0;
			Omega_Logit[pep] = -.5 + unif_rand();
			double C;
			for(i = 0; i < *n_indiv; i++)
			{
				C = exp(-R_pow_di(Y[(*n_indiv)*pep + i] - Mu[pep],2)
						+ R_pow_di(Y[(*n_indiv)*pep + i] - Mu[pep] - Alpha_pep[pep],2));
				prob = C/(C + 1.0);
				W[(*n_indiv)*pep + i] = 1.0;
				Gamma[(*n_indiv)*pep + i] = (int)(unif_rand() >= prob);
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
		double* restrict Mu,
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
		const double mu = Mu[pep];
		const double sig2 = Sig2[pep];
		const double alpha = Alpha[pep];
		for(i = 0; i < n_indiv; i++)
		{
			const double y = Exprs[d_ind];
			const int gam = Gamma[d_ind];
			ds = (y - mu - alpha*gam)*(y - mu - alpha*gam)/sig2;

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

void update_Nu(double *Omega_Logit,
		double *Nu,	double nu0, double tau, double *phi,
		int n_position, double *pnum_double, int *pnum,
		double *Q_work, double *L_cur, double *b_cur, double *mu_cur,
		double *L_prop, double *b_prop, double *mu_prop,
		RngStream rng, adpt *phi_tune)
{
	/*
	 * QB initialized as an (lda, n_position) array of ones/twos
	 * L must be an (lda, n_position) array of zeros
	 * b is a n_position length array
	 */

	char low = 'L', trans = 'T', diag = 'N';
	char *UPLO = &low, *TRANS = &trans, *DIAG = &diag;
	int p = 0, c = 0, k = 0;
	const int inc1 = 1, inc2 = 2;
	const int kd = 1;
	const int ld = 2;
	const int length_L = ld*n_position;
	int status0, status1;
	const double ZERO = 0.0, ONE = 1.0;
	double scale = phi_tune->tune;

	double lik_prop, phi_prop, lik_cur;
	double half_ldetQ_prop = 0.0, half_ldetQ_cur = 0.0;

	// log scale random walk proposal
	phi_prop = (*phi)*exp(scale*(RngStream_RandU01(rng) - .5));

	for(p = 0; p < n_position; p++)
	{
		b_cur[p] = 0.0;
		for(c = 0; c < pnum[p]; c++)
		{
			b_cur[p] += Omega_Logit[k];
			k++;
		}
		b_cur[p] = b_cur[p]*tau;
		b_prop[p] = b_cur[p];
		Nu[p] = RngStream_N01(rng);
	}
	// get b vectors for proposal and current value of phi
	b_prop[0] = b_prop[0] + phi_prop*nu0;
	b_cur[0] = b_cur[0] + (*phi)*nu0;

	F77_NAME(dcopy)(&n_position, b_cur, &inc1, mu_cur, &inc1);
	F77_NAME(dcopy)(&n_position, b_prop, &inc1, mu_prop, &inc1);

	// prepare band matrices L_cur and L_prop for current values
	F77_NAME(daxpy)(&n_position, &tau, pnum_double, &inc1, L_cur, &inc2);
	F77_NAME(dcopy)(&length_L, L_cur, &inc1, L_prop, &inc1);

	F77_NAME(daxpy)(&length_L, phi, Q_work, &inc1, L_cur, &inc1);
	F77_NAME(daxpy)(&length_L, &phi_prop, Q_work, &inc1, L_prop, &inc1);

	/* L matricies now hold the packed band representation of Q for the current
	 * and proposal parameter values. Next compute the cholesky decomposition
	 * of the banded matrices Q and store in L matrices */
	F77_NAME(dpbtrf)(UPLO, &n_position, &kd, L_cur, &ld, &status0);
	F77_NAME(dpbtrf)(UPLO, &n_position, &kd, L_prop, &ld, &status1);

	if(status0 != 0 || status1 != 0)
	{
		error("an argument in dpbtrf had illegal value\n");
	}

	// solves Q*mu = b. On input mu is b, on output mu is, in fact, mu.
	F77_NAME(dpbtrs)(UPLO, &n_position, &kd, &inc1, L_prop,
			&ld, mu_prop, &n_position, &status0);
	F77_NAME(dpbtrs)(UPLO, &n_position, &kd, &inc1, L_cur,
				&ld, mu_cur, &n_position, &status1);
	if(status0 != 0 || status1 != 0)
	{
		error("an argument in dpdtrs had illegal value\n");
	}

	// get .5*log(det(Q)), exploiting the lower triangle factorization
	for(p = 0; p < n_position; p++)
	{
		half_ldetQ_cur += log(L_cur[p*2]);
		half_ldetQ_prop += log(L_prop[p*2]);
	}

	// get quantities mu**T*Q*mu through b**T*mu
	double btmu_cur, btmu_prop;
	btmu_cur = F77_NAME(ddot)(&n_position, b_cur, &inc1, mu_cur, &inc1);
	btmu_prop = F77_NAME(ddot)(&n_position, b_prop, &inc1, mu_prop, &inc1);

	lik_cur = (2.0 + .5*n_position)*log((*phi)) - half_ldetQ_cur - (*phi)*2.0 + .5*btmu_cur;
	lik_prop = (2.0 + .5*n_position)*log(phi_prop) - half_ldetQ_prop - phi_prop*2.0 + .5*btmu_prop;

	double logU = log(RngStream_RandU01(rng));
	if(logU < (lik_prop - lik_cur))
	{
		*phi = phi_prop;
		phi_tune->count++;

		/* efficiently solves L**T*y = z for y, giving a normal sample with
		 * mean zero and precision Q */
		F77_NAME(dtbsv)(UPLO, TRANS, DIAG, &n_position, &kd,
				L_prop, &ld, Nu, &inc1);

		/* On entry, Nu holds a normal vector with
		 * mean zero and variance Q**-1. On exit,
		 * Nu holds a conditional distribution
		 * sample for the current iteration. */
		F77_NAME(daxpy)(&n_position, &ONE, mu_prop, &inc1, Nu, &inc1);
	}
	else
	{
		F77_NAME(dtbsv)(UPLO, TRANS, DIAG, &n_position, &kd,
				L_cur, &ld, Nu, &inc1);
		F77_NAME(daxpy)(&n_position, &ONE, mu_cur, &inc1, Nu, &inc1);
	}

	// zero out L's for the next iteration
	F77_NAME(dscal)(&length_L, &ZERO, L_cur, &inc1);
	F77_NAME(dscal)(&length_L, &ZERO, L_prop, &inc1);
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

double p1mp(double x)
{
	return(1.0/(2.0 + exp(-x) + exp(x)));
}

inline double log_from_logit(double x)
{
	if(x > 0.0)
	{
		return(-log1pexp(-x));
	}
	else
	{
		return(x - log1pexp(x));
	}
}

inline double log1m_from_logit(double x)
{
	if(x > 0.0)
	{
		return(-x - log1pexp(-x));
	}
	else
	{
		return(-log1pexp(x));
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
			tune_set->tune = (tune_set->tune)*(1.0 + (n_burn + 1.0 - i)/(n_burn));
		}
		else if(accept_prop < .21)
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





