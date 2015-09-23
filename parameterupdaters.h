#ifndef PARAMETERUPDATERS
#define PARAMETERUPDATERS
#define PRINT_PATH_STEP_PROB 0
#define PRINT_NONZERO_COUNT_PATH 0
#include "helperfunctions.h"
#include "distributionlist.h"


class ParameterUpdater{

public:
  string par_name;       /*name of the parameter*/
  double* proposal_value; /*storage for the value that is proposed*/
  double* current_value;  /*current value copied over from a curr_par_obj*/
  short int is_accepted; /*acceptance value (yes/no) */
  double* hyperpar;      /*hyperparameters for the prior */
  int num_hyperpars;     /*number of hyperparameters*/
  double proposal_sd;    /*tuning ratio for proposals */
  int hp_allocated;      /*indicator of whether or not hyperparameters have been allocated*/
  DistributionList* distributions;
  int dim;
  int PRINT_MCMC_STEP;
  double lb, ub;
  int l_is_inf, u_is_inf;
  /* Constructor N = # of hyperparameters */
  ParameterUpdater(int &N, double *&hpar1,DistributionList*& dists,  double &propsd);

  /* for Gibbs steps (no tuning parameter needed */
  ParameterUpdater(int &N, double *&hpar1,DistributionList*& dists);

  /* empty constructor */
  ParameterUpdater(DistributionList*& dists);

  /* for an array with no hyperparameters */
  ParameterUpdater(int &N, DistributionList*& dists);

  void PrintHyperpars(ostream& str);
  /*Virtual functions (different priors and posteriors for each parameter*/
  virtual double Logprior(double value){return(0.0);}
  virtual double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars){return 0.0;}
  virtual void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars){}
  virtual double Logjump(double val1, double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars){return 0.0;}
  virtual void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){ }
  virtual void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){ }
  virtual double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars){return 0.0;}
  void SetValidStartValue(curr_par_obj* current_values, global_par_obj* global_pars);
  /*The MCMC step is (basically) the same for any kind of parameter */
  short int MCMC_Step(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars, int iteration_id, int print_step=0);
  ~ParameterUpdater();
};


class OffLambdaUpdater: public ParameterUpdater
{
 public:
  OffLambdaUpdater(int N, double* hpar1, DistributionList* dists);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){current_values->off_lambda = proposal_value[0];}
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = current_values->off_lambda;}
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  // double Logjump(double val1,double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  //double Logprior(double value);
  //double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
  
};

class QUpdater: public ParameterUpdater
{
public:
  QUpdater(int N, double* hpar1, DistributionList* dists);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){current_values->q = proposal_value[0];}
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = current_values->q;}
   /*Gibbs step; only a proposal value is needed*/
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
 
}; 

class AlphaUpdater: public ParameterUpdater
{
public:
  AlphaUpdater(int N, double* hpar1, DistributionList* dists, double propvar);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){current_values->alpha = proposal_value[0];}
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = current_values->alpha;}
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logprior(double value);
  double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logjump(double val1,double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars); 
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
 
};

class OmegaScaledUpdater: public ParameterUpdater
{
 public:
  OmegaScaledUpdater(int N, double* hpar1, DistributionList* dists, double propvar);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){current_values->omega = proposal_value[0];}
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = current_values->omega;}
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logprior(double value);
  double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logjump(double val1, double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
 
};


class RhoUpdater: public ParameterUpdater
{
public:
  short int state;
  RhoUpdater(short int stateval,int N, double* hpar1, DistributionList* dists, double propvar);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){current_values->rho[state] = proposal_value[0];}
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = current_values->rho[state];} 
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logprior(double value);
  double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logjump(double val1, double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
 
};


class NuUpdater: public ParameterUpdater
{
public:
  short int state;
  NuUpdater(short int stateval, int N, double* hpar1, DistributionList* dists, double propvar);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){current_values->nu[state] = proposal_value[0];}
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = current_values->nu[state];}
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logprior(double value);
  double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logjump(double val1, double val2,  curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
  
};

class GammaUpdater: public ParameterUpdater
{
public:
  short int state;
  GammaUpdater(short int stateval,int N, double* hpar1, DistributionList* dists);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){current_values->gamma[state] = proposal_value[0];}
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = current_values->gamma[state];}
  /*Gibbs step: only a proposal is needed*/
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);

};

class StateUpdater: public ParameterUpdater
{
public:
  curr_par_obj* temp_parameters;
  int index_selected;
  int is_numerator;
  double curr_hfun1, curr_hfun2, prop_hfun1, prop_hfun2;
  //hyperparameters are the bounds for off lambda assuming uniform prior
  StateUpdater(int N, double* hpar1, DistributionList* dists);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars);
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){temp_parameters->SetFromCopy(current_values);}
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars); 
  double SelectionWeight(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
  ~StateUpdater();
};

class BirthDeathImmuneUpdater: public ParameterUpdater
{
 public:
  double likelihoods[6];
  int last_nonoff_transition;
  int first_nonoff_transition;
  double sum_nonoff_transitions;
  curr_par_obj* temp_parameters;
  BirthDeathImmuneUpdater(int N, double* hpar1, DistributionList* dists);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars);
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars);
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
  ~BirthDeathImmuneUpdater();
};


class ImmuneRateUpdater: public ParameterUpdater
{ 
 public:
  ImmuneRateUpdater(int N, double* hpar1, DistributionList* dists);
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars){current_value[0] = global_pars->immune_rate;}
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars){global_pars->immune_rate = proposal_value[0];}
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
  
};

class SurvivalRateUpdater: public ParameterUpdater
{
 public:
  SurvivalRateUpdater(int N, DistributionList* dists);
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars);
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars);
 
};

/*

class SurvivalHyperUpdater: public ParameterUpdater
{
 public:
  int pos;
  int val;
  
 SurvivalHyperUpdater(int N, double* hpar1, DistributionList* dists):ParameterUpdater(N, hpar1, dists){par_name="SurvivalHyperPars"};
  void SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars);
  void SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars);
  void Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logprior(double value);
  double Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);
  double Logjump(double val1, double val2,  curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars);

}
*/


#endif
