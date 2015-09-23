#ifndef MCMCDRIVERS
#define MCMCDRIVERS
#include "parameterupdaters.h"
#include "mergesplitupdaters.h"
#include "networkobjects.h"


class SingleModelUpdater{

 public:
   
  /*Single Host model parameter updaters */
  vector<ParameterUpdater*> parameter_updaters;  /* Estimation of machine-specific parameters */
  vector<ParameterUpdater*>::iterator parup_iter;/* Iterator for above */
  vector<int> repeats;     /*Number of repeats per parameter within the machine-specific MCMC */
  vector<int> estimate_par;/*0-1 indicator of whether or not machine parameter k is estimated */
 
  SingleModelUpdater(const char* parfile, const char* repeatfile, DistributionList* dists, int wp, int ws);
  void UpdateMachine(Machine* machine, global_par_obj* global_parameters, int iterations, int thin, ofstream& par_output, ofstream& state_output, ofstream &acc_output, int export_acc, int iteration_id, int print_steps_to_screen =0);
  ~SingleModelUpdater();
  void PrintRepeatValues(ostream& str=cout);
  void PrintAllHyperpars(ostream& str=cout);
  double GraphLogLikelihood(Machine* machine, global_par_obj* global_pars, int suppress_prior = 0);
  DistributionList* distributions;
 
  int NumPars(){return NUM_TOT;}
 private:
  int Q_REP;  
  int ALPHA_REP;
  int OMEGA_REP;
  int RHO_REP;
  int GAMMA_REP;
  int NU_REP;
  int OFFLAMBDA_REP;
  int BIRTHDEATHIMMUNE_REPEAT;

  int STATE_REPEAT;  
  int OVERALL_PAR_REPEAT;
  int NUM_PAR;
  int NUM_TOT;
  int write_pars;
  int write_state;
};


class CountRearranger
{
 public:
  DistributionList* distributions;

  CountRearranger(double mix,int repeat_iter, int thin_iter, int is_noninfo,int num_hr, int tp, DistributionList* dists);
  int Sample(Network* network, map<machine_id_t, Machine*>& machines, int print_to_screen, int iteration_id, ostream& str, int print_output, int deterministic=0);
  //type = 0: Simple Gibbs Step
  //type = 1: Always do a StateChange step
  //type = 2: Simple Gibbs Step if more than one nonzero rate; StateChange step between existing guy and one other otherwise
  //type = 3: Cointoss?

  void SetMix(double mix);
  int GetRepeat(){return repeat;}
  
  ~CountRearranger();
 private:
  double* pk;
  int* simulated_counts;
  int vectors_allocated;
  vector<Machine*> mapped_machines;
  double mixprob;
  set<machine_id_t>::iterator id_iter;
  map<machine_id_t, Machine*>::iterator machine_iter;
  int min_birth;
  int max_death;
  int repeat;
  int thin;
  int is_noninformative;
  int num_hours;
  int type;
  int SetUpMachines(Network* network, map<machine_id_t, Machine*>& machines);
  int SetUpProbabilities(int t); 
  int GibbsStep(Network* network, int t);
  int StateChangeStep(Network* network, int t);
  void SetNewCount(Network* network, int machine_index, int t);
  
};



#endif
