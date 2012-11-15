/*
 * PMA_multi_posvar_censored.h
 *
 *  Created on: Nov 2, 2012
 *      Author: Gregory Imholte
 */

#ifndef PMA_MULTI_POSVAR_CENSORED_H_
#define PMA_MULTI_POSVAR_CENSORED_H_

struct MH_TUNE{
	int count;
	int total_count;
	double tune;
};

typedef struct MH_TUNE adpt;

double lc_AB(double x, double *argvec, int *arglen);
double lcp_AB(double x, double *argvec, int *arglen);
double lc_alpha(double x, double *argvec, int *arglen);
double lcp_alpha(double x, double *argvec, int *arglen);
double lc_alpha_int(double x, double *argvec, int *arglen);
double lcp_alpha_int(double x, double *argvec, int *arglen);
double logit(double x);

inline double log1m_from_logit(double x);
inline double log_from_logit(double x);
inline double expit(double x);
inline double m1expit(double x);

double truncNorm(double mean, double sigmasqr, RngStream rng);

void update_tuning(adpt *tune_set, int i, int n_burn);

void store_mcmc_output(double *Alpha, double *Mu, double *A, double *B, double *U,
		double *Sig2, double *D, double *Omega_Logit, int* Omega_Ind, double kappa,
		double alpha, double beta, double m, double zeta, double dof,
		int n_peptide, int n_indiv, int n_position, int *cen_num,
		FILE *AFILE, FILE *BFILE, FILE *PFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *OFILE, FILE *ALPHAFILE);

void pepbayes_v1(double *Y, double *hyper_param, int *pstart,
		int *pnum, int *n_position, int *n_peptide, int *n_indiv, int *nP,
		int *cen_ind, int *cen_pep, int *cen_num, int *cen_pos,
		int *n_iter, int *n_sweep, int *n_burn,
		double *OutProbs, int *write);

void update_global_params(double *m, double *zeta, double *alpha, double *beta, double *kappa,
		double *Mu, double *Sig2, double *Alpha_pep,
		double m0, double v0, double alpha0, double beta0,
		int n_indiv, int n_peptide, RngStream rng, adpt *m_tune, adpt *z_tune,
		ARS_workspace *ws, double *xAlpha, int nP);

void update_peptide(double *Exprs, double *Mu, double *Alpha_pep, double *W,
		int *Omega_ind, double *Omega_logit, int *Gamma,
		double *Sig2_pep, double alpha, double beta,
		double u, double a, double b,
		double m, double zeta, double dof, int n_indiv, int n_peptide, int pep,
		RngStream rng, double *RB);

void update_censoring(double *W, double *D, int cen_num, int* cen_ind,
		int* cen_pep, double *Y, double *Exprs,
		int *Gamma, double *Alpha, double *Mu, int n_peptide, double *Sig2,
		RngStream rng);

void initialize_chain(double *ProbSum, double *Exprs, double *Y,
		double *W, int *Omega_Ind, double *Omega_Logit,
		double *Alpha_pep, int *Gamma,
		double *Sig2, double *Mu, double *A, double *B, double *U,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		double **xA, double **xB, double *RB);

void update_dof_integrated(double *dof, double *Exprs, double *W,
		double *Alpha, int *Gamma,
		double *Sig2, double *Mu, double *workspace,
		int n_indiv, int n_peptide,
		RngStream rng);

void update_indiv_mu(double *Exprs, double *W,
		double *Alpha_pep, double *mu_j, int *Gamma,
		double *Sig2, double kappa,
		int n_peptide, int j, RngStream rng);

void update_position_p(int* Omega_Ind, double *Omega_Logit,
		double *A_p, double *B_p, double *U_p,
		double a_0, double b_0, double lambda_a, double lambda_b,
		int p, int n_indiv, int *pnum, int *pstart, RngStream rng,
		double *xA, double *xB, ARS_workspace *workspace);

#endif /* PMA_MULTI_POSVAR_CENSORED_H_ */
