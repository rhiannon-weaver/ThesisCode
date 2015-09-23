#ifndef HELPERFUNCTIONS
#define HELPERFUNCTIONS
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include "datparobjects.h"
#include "distributionlist.h"
#define PRINT_PATH_STEP_PROB 0
#define GAMMA_SCALE_MIN 0.00000000001
#define BETA_MAX 0.999999999999999
#define BETA_MIN 0.000000000000001
#define MIN_LOG_P 1e-300
#define MAX_LOG_P 1e300



double trunc_norm_logjump(double val1,
			  double val2,
			  double lim1, 
			  double lim2, 
			  double prop_sd, 
			  DistributionList* distributions);

double trunc_norm_proposal(double current, 
			   double proposalsd, 
			   double lowerlim, 
			   double upperlim, 
			   DistributionList* distributions);

double BetaConstrain(double val);


/* Turning the transition probability parameters into a transition probability */
double sinecurveprob_single(short int pristate, 
			    short int nextstate,
			    short int t,
			    double rho,
			    double gamma,
			    double nu, 
			    curr_par_obj* current_values);


double sum_alpha_k(short int beginval, 
		   short int endval, 
		   curr_par_obj* current_values, //vector<double>& lagstate, 
		   double alpha, 
		   double omega = -1.0);


double log_prod_y_alpha(mach_dat_obj* data, 
			curr_par_obj* current_values, //vector<double>& lagstate, 
			double alpha, 
			double omega,
			short int birthdate,
			short int deathdate);

/* this has been adapted for the model that counts poisson zero states */
double NextStateProb(int pristate, 
		     int which_parameter, 
		     double value, 
		     curr_par_obj* current_values,
		     global_par_obj* global_pars);

void PathToChange(int index, 
		  curr_par_obj* current_values, 
		  mach_dat_obj* data, 
		  int* beginend);


double GetLogPoissonProbabilities(mach_dat_obj* data, 
				  curr_par_obj* current_values, 
				  int startval, 
				  int endval, 
				  vector<double>* indiv_pois);

double CondPoissonLoglike(int x, 
			  double lambda);

double GetLogTransitionProbabilities(mach_dat_obj* data, 
				     curr_par_obj* current_values,
				     global_par_obj* global_pars,
				     int startval, 
				     int endval, 
				     vector<double>* indiv_trans);

double GetPathLogLikelihood(int startval, 
			    int endval, 
			    mach_dat_obj* data,
			    curr_par_obj* current_values,
			    global_par_obj* global_pars,
			    int whichtype, 
			    vector<double>* indiv_trans, 
			    vector<double>* indiv_pois);

double EtaPropH(int index, 
		int startval, 
		int endval, 
		mach_dat_obj* data, 
		curr_par_obj* current_values, 
		global_par_obj* global_pars,
		vector<double>* indiv_trans, 
		vector<double>* indiv_pois);

double LogAplusB(double loga, 
		 double logb);


/*Translate between a "mean, temperature" parameterization of a Beta distribution to the traditional "alpha,beta"  parameterization   */
double BetaDistAlpha(double mean, double temperature);
double BetaDistBeta(double mean, double temperature);

/*Code re-use for OffLambda updater */


/*void ExternalOffStateSummary(int* beginend, 
			     curr_par_obj* current_values, 
			     double* sumprod);*/
double ScaledGammaLogConstant(double lower, double upper, double alpha, double beta);

double LogSpikeProb(double spikemult, short int state1, short int state2);

void StateRunChanges(curr_par_obj* current_values, curr_par_obj* proposed_values, int index, int beginval, int endval, int* add);

double LogGamma(double x);

double LogPriorRatioScaledGamma(double s1, double s2, double s,double minval, double maxval, double pri_k, double pri_theta);
double LogPriorRatioScaledBeta(double s1, double s2, double s, double minval, double maxval, double pri_alpha, double pri_beta);

int FirstLastRun(curr_par_obj* current_values, int which, int use_absolute=0);

int GetBDIIndex(short int birth, short int death, short int immune, int max_death);

void BDIFromIndex(int index, int last_nonoff, int max_death, short int *birth, short int *death, short int *immune);

double LogMoveRatio(vector<double>& moveprobs, int proposed_type, int current_type);

int SimpleSample(vector<double>& prob, DistributionList* dists, double probsum=-1); 

void GetPlusMinusState(int state, curr_par_obj* current_values, double& plusstatecount, double& minusstatecount);

void RhoNuSuccessesAndTotalCounts(short int state, curr_par_obj* current_values, int mod_index, int& successes, int& totals);

void ActiveOverlap(curr_par_obj* par1, curr_par_obj* par2, double& A1, double& A2, double& A1U2, double& A1I2);

double ScaledBetaMode(double l, double u, double alpha, double beta);

void ScaledBetaAlphaBeta(double l, double u, double mode, double samplesize, double& alpha, double& beta);



#endif
