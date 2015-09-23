#ifndef COUNTUPDATERS
#define COUNTUPDATERS
#include "mergesplitupdaters.h"
#include "mcmcdrivers.h"
#include <gsl/gsl_sf_gamma.h>
#define PRINT_RATIO_STEP 1
class CountMergeSplitUpdater: public MergeSplitUpdater{
 public:
  
  CountMergeSplitUpdater(vector<double>& hpar, DistributionList* dists, CountRearranger* cr);
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
 
 private:
  CountRearranger* count_rearranger;
  map<machine_id_t, Machine*> machine_map;
  double s1_mse;
  double s2_mse;
  double s_mse;
  double s1_tot;
  double s2_tot;
  double s_tot;
  double RateLogLikelihood(unsigned int network_id, int t, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, double mergesplit);
  
};



class MachineMergeSplit{
 public:
  int num_repeats;
  vector<int> estimate_par;
  DistributionList* distributions;
  vector<double> merge_split_prob;
  vector<MergeSplitUpdater*> merge_split_updaters;
  MachineMergeSplit(const char* mergesplitparfile, SingleModelUpdater* single_model, CountRearranger* count_rearranger, DistributionList* dists, int repeats, ostream& info_out);  
  int MCMCStep(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars, int print_to_screen, int iteration_id, ostream& str, int print_output, int record_all_acceptance);
  void PrintAcceptance(ostream& str);
  int NumPars(){return NUM_TOT;}
 
  ~MachineMergeSplit();
 private:
  int num_tries;
  int num_accepted;
  int NUM_TOT;
  map<int, int[2]> Network_Info;
  SingleModelUpdater* machineMCMC;
  double ImmuneOffLambdaRatio(vector<Machine*>& proposed_machines, vector<Machine*>& current_machines, global_par_obj* global_pars, int* immune, double& middle_multiplier, vector<Machine*>& ptr);
  double SinglePathPoissonLoglike(int t, vector<Machine*>& ptr, double middle_multiplier);
  double SinglePathTransLoglike(int t, vector<Machine*>& ptr, global_par_obj* global_pars, int* immunity, double middle_multiplier);
  double BernoulliRatio( int* epsilons, int* mask, double middle_multiplier, double p);
  double LogLikelihoodRatio( vector<Machine*>& proposed_machines, vector<Machine*>& current_machines, global_par_obj* global_pars, int print_to_screen, int iteration_id);
  vector<double> pri_H_par; //p, n for negative binomial prior for H
  
 };



#endif
