/*
 * pepBayes_v4.h
 *
 *  Created on: Nov 20, 2012
 *      Author: hoblitz
 */

#ifndef PEPBAYES_V4_H_
#define PEPBAYES_V4_H_


struct MH_TUNE{
	int count;
	int total_count;
	double tune;
};

typedef struct MH_TUNE adpt;

double lc_AB(double x, double* restrict argvec, int *arglen);
double lcp_AB(double x, double* restrict argvec, int *arglen);
double lc_alpha(double x, double *argvec, int *arglen);
double lcp_alpha(double x, double *argvec, int *arglen);
double lc_alpha_int(double x, double* restrict argvec, int *arglen);
double lcp_alpha_int(double x, double* restrict argvec, int *arglen);
double logit(double x);

inline double log1m_from_logit(double x);
inline double log_from_logit(double x);
inline double expit(double x);
inline double m1expit(double x);

double truncNorm(double mean, double sigmasqr, RngStream rng);

void update_tuning(adpt *tune_set, int i, int n_burn);

double p1mp(double x);

void store_mcmc_output(double *Alpha, double *Mu, double *Nu,
		double *Sig2, double *D, double *Omega_Logit, double kappa,
		double alpha, double beta, double m, double zeta, double phi,
		double tau, double dof,
		int n_peptide, int n_indiv, int n_position, int *cen_num,
		FILE *NUFILE, FILE *VARFILE, FILE *Sig2FILE, FILE *MUFILE,
		FILE *DFILE, FILE *OFILE, FILE *ALPHAFILE);

void pepbayes_v4(double *Y, double *hyper_param, int *pstart,
		int *pnum, int *n_position, int *n_peptide, int *n_indiv, int *nP,
		int *cen_ind, int *cen_pep, int *cen_num, int *cen_pos,
		int *n_iter, int *n_sweep, int *n_burn,
		double *OutProbs, int *write);

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
		ARS_workspace *ws, double *xAlpha, int nP);

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
		RngStream rng, double* restrict RB, adpt *omega_adpt);

void update_censoring(double *W, double *D, int cen_num, int* cen_ind,
		int* cen_pep, double *Y, double *Exprs,
		int *Gamma, double *Alpha, double *Mu, int n_indiv, double *Sig2,
		RngStream rng);

void initialize_chain(double *ProbSum, double *Exprs, double *Y,
		double *W, double *Omega_Logit,
		double *Alpha_pep, int *Gamma,
		double *Sig2, double *Mu, double *Nu,
		int *n_position, int *pstart, int *pnum, int *n_peptide, int *n_indiv,
		double *RB);

void update_dof_integrated(double* restrict dof,
		double* restrict Exprs,
		double* restrict W,
		double* restrict Alpha,
		int* restrict Gamma,
		double* restrict Sig2,
		double* restrict Mu,
		double* restrict workspace,
		int n_indiv, int n_peptide,
		RngStream rng);

void update_Nu(double *Omega_Logit,
		double *Nu,	double nu0, double tau, double *phi,
		int n_position, double *pnum_double, int *pnum,
		double *Q_work, double *L_cur, double *b_cur, double *mu_cur,
		double *L_prop, double *b_prop, double *mu_prop,
		RngStream rng, adpt *phi_tune);

#endif /* PEPBAYES_V4_H_ */
