#ifndef MERGESPLITUPDATERS
#define MERGESPLITUPDATERS
#include <gsl/gsl_randist.h>
#include <iomanip>
#include "networkobjects.h"
#include "helperfunctions.h"
#include "distributionlist.h"
#include "parameterupdaters.h"
#define LAG_MAX 1000.0

class AsymmetricAverage{
 public:
  AsymmetricAverage(double lb, double ub, double prop_alpha, double prop_beta, int is_u_infinity, int is_l_infinity, double assignment_prob_unif_weight, double assignment_prob_power, DistributionList* distributions);
  void MergePropose(double q1, double q2, double &q);
  void SplitPropose(double q, double &q1, double &q2, double assignment_prob);
  double MergeLogRatio(double q1, double q2, double assignment_prob);
  double SplitLogRatio(double q1, double q2, double assignment_prob);
  
  void gfun(double q, double u1, double& qmin, double& qmax);
  void ginv(double qmin, double qmax, double& q, double& u1);

  double JacobianDqu1Dq1q2(double q1, double q2);
  double JacobianDq1q2Dqu1(double q, double u1);
  double AdjustAssignment(double assignment_prob, double invert);
  ~AsymmetricAverage(){}
private:
  double l, u;
  int l_inf, u_inf;
  //to keep the assignment probability exactly as specified use, assign_alpha = assign_k = 1
  double assign_alpha;
  double assign_k;
  double u1_alpha, u1_beta;
  DistributionList* dists;

};




class BoundedRescaledSampler{

 public:
  BoundedRescaledSampler(double lb, double ub, double alphaplusbeta, DistributionList* dists);
  void MergePropose(double q1, double q2, double &q, vector<double>& proposal_pars, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  void SplitPropose(double q, double &q1, double &q2,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines );
  //default merge: d(1)+d(2) - d(12); use mergesplit=-1 for split, 1, for merge
  double MergeLogRatio(double q1, double q2, double q12,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  double SplitLogRatio(double q1, double q2, double q12,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  void GetAlphaBeta(double q, double & alpha, double & beta);
  virtual void gfun1(double q, double& q1, double& q2,  vector<double>& proposal_pars,   vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  virtual void gfun2(double q1, double q2, double &q,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  virtual double LogRatio(double q1, double q2, double q12,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, int mergesplit);
 protected:
  double l, u, samplesize;
  DistributionList* distributions;

};

class QBoundedRescaledSampler: public BoundedRescaledSampler{
 public:
 QBoundedRescaledSampler(double lb, double ub, double alphaplusbeta, DistributionList* dists):BoundedRescaledSampler(lb, ub, alphaplusbeta, dists){y_stats = new YStateAndCountStats();}
  void gfun1(double q, double& q1, double& q2, vector<double>& proposal_pars,   vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  void gfun2(double q1, double q2, double& q,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  ~QBoundedRescaledSampler(){delete y_stats;}
  YStateAndCountStats* y_stats;
};

class AlphaBoundedRescaledSampler: public BoundedRescaledSampler{
 public:
 AlphaBoundedRescaledSampler(double lb, double ub, double alphaplusbeta, DistributionList* dists):BoundedRescaledSampler(lb, ub, alphaplusbeta, dists){y_stats = new YStateAndCountStats();}
  ~AlphaBoundedRescaledSampler(){delete y_stats;}
  //void gfun1(double q, double& q1, double& q2, vector<double>& pars);
  //void gfun2(double q1, double q2, double& q, vector<double>& pars);
  YStateAndCountStats* y_stats;
};


class OmegaBoundedRescaledSampler: public BoundedRescaledSampler{
 public:
 OmegaBoundedRescaledSampler(double lb, double ub, double alphaplusbeta, DistributionList* dists):BoundedRescaledSampler(lb, ub, alphaplusbeta, dists){y_stats = new YStateAndCountStats();}
  
  void gfun1(double q, double& q1, double& q2, vector<double>& proposal_pars,   vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  void gfun2(double q1, double q2, double& q,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  //double LogRatio(double q1, double q2, double q12,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  ~OmegaBoundedRescaledSampler(){delete y_stats;}
  YStateAndCountStats* y_stats;
};


class MergeSplitUpdater{

 public:
  DistributionList* distributions;
  vector<double> prior_hyperpar;
  vector<double> proposal_par;
  string par_name;
  MergeSplitUpdater(vector<double>& hpar, DistributionList* dists, ParameterUpdater* parameter_updater);
  MergeSplitUpdater(vector<double>& hpar, DistributionList* dists);

  double MergeSplitStep(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars, int print_to_screen); 
  void PrintToScreen(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, double loglike);

  /*These should be designed with the ability only to report -inf or +inf, not nan */
  virtual void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars){}
  virtual double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars){return 0.0;}
 virtual void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars){}
 virtual double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars){return 0.0;}


 virtual double LogPriorRatio(double s, double s1, double s2, int mergesplit){return 0;}
 virtual void PrintValue(vector<Machine*>& machines, ostream& str){}
 virtual void PrintHyperpars(ostream& str=cout);
  ~MergeSplitUpdater(){}
};

class StateMergeSplitUpdater: public MergeSplitUpdater{
 public:
  StateMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par);
   
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
  void PrintHyperpars(ostream& str=cout);
  ~StateMergeSplitUpdater(){delete y_active_stats;}
 private:
  int s_count[3];
  int s1_count[3];
  int s2_count[3];
  YTotalActiveStats* y_active_stats;
  double log_split_prob;
  short int MergeState(short int s1, short int s2);
  double SplitState(short int s, double p_both_on, double pri_lagstate1, double pri_lagstate2, double pri_offrun1, double pri_offrun2, short int& s1, short int& s2);
  double MergeStateProb(short int s1, short int s2, short int s);
  double SplitStateProb(short int s, short int s1, short int s2, double prob_both_on, double pri_lagstate1, double pri_lagstate2, double pri_offrun1, double pri_offrun2);
  void SetLagstateAndOffrun(vector<Machine*>& machines, int i, double& pri_lagstate1, double& pri_lagstate2, double& pri_offrun1, double& pri_offrun2);
  double GetProbabilityBothOn(vector<Machine*>& machines, int i);
  double GetProb1Off(double pri_lagstate1, double pri_lagstate2, double pri_offrun1, double pri_offrun2);
};


class BirthDeathImmuneMergeSplitUpdater: public MergeSplitUpdater{
 public:

  BirthDeathImmuneMergeSplitUpdater(string probfile, vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par);
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
   void PrintHyperpars(ostream& str=cout);
 private:
  //1, 2, 1+2
  double bdi_prob[36][6];
  double row_sum[36];
  double col_sum[6];
  void SetProbsMerge(int bdi1, int bdi2, int sum_only = 0);
  void SetProbsSplit(int bdi, int first_nonoff_1, int first_nonoff_2, int sum_only = 0);	
  void SetFirstLast(vector<Machine*>& two_machines);
  int CheckFirstLast(vector<Machine*>& one_machine);
  vector<int> first_nonoff; //1, 2   (1+2 = min 1,2)
  vector<int> last_nonoff;  //1, 2   (1+2 = max 1,2)
  vector<double> cur_split_prob;
  vector<double> cur_merge_prob;
  double cur_split_sum;
  double cur_merge_sum;
};



class QMergeSplitUpdater:public MergeSplitUpdater{
 public:
  
  QMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, QUpdater* par);
  
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
 void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double LogPriorRatio(double s, double s1, double s2, int mergesplit);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
  ~QMergeSplitUpdater(){delete sampler;}
 private:
  //vector<double> sampler_pars;
  QBoundedRescaledSampler* sampler;
  //void SetGParameters(vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  //samplesize = proposal_par[0]
  //any fixed (not state dependent) parameters for gfun are in the rest of proposal_par.
  //l = prior_hyperpar[0], u = prior_hyperpar[1]
};


class OmegaScaledMergeSplitUpdater: public MergeSplitUpdater{
 public:
  
  OmegaScaledMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par);

  
 void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
 void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
 double LogPriorRatio(double s, double s1, double s2, int mergesplit);
  
   void PrintValue(vector<Machine*>& machines, ostream& str=cout);
   ~OmegaScaledMergeSplitUpdater(); 
 private:
   double l, u, samplesize;
  
  
   //vector<double> sampler_pars;
   OmegaBoundedRescaledSampler* sampler;
   //void SetGParameters(vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
   //samplesize = proposal_par[0]
   //any fixed (not state dependent) parameters for gfun are in the rest of proposal_par.
   //l = prior_hyperpar[0], u = prior_hyperpar[1]
};

class AlphaMergeSplitUpdater: public MergeSplitUpdater{
 public:
  AlphaMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par);
 
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double LogPriorRatio(double s, double s1, double s2, int mergesplit);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
  ~AlphaMergeSplitUpdater(){delete sampler;}
 private:
  //vector<double> sampler_pars;
  AlphaBoundedRescaledSampler* sampler;
  //void SetGParameters(vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  //samplesize = proposal_par[0]
  //any fixed (not state dependent) parameters for gfun are in the rest of proposal_par.
  //l = prior_hyperpar[0], u = prior_hyperpar[1]
  
};

class OffLambdaMergeSplitUpdater: public MergeSplitUpdater{
 public:
 OffLambdaMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par):MergeSplitUpdater(hpar1, dists, par){}

  
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
 void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double LogPriorRatio(double s1, double s2, double s, int mergesplit);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
  
};

class RhoMergeSplitUpdater: public MergeSplitUpdater{
 public:
  short int state;
  RhoMergeSplitUpdater( vector<double>& hpar, DistributionList* dists, RhoUpdater* par);
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double LogPriorRatio(double s1, double s1, double s, int mergesplit);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
 
 private:
  //AsymmetricAverage* avg;
  double LeastSquares(Machine* machine, double& samplesize);
  double rho12_mean, rho1_mean, rho2_mean;
  double rho12_ss, rho1_ss, rho2_ss;
};

class GammaMergeSplitUpdater: public MergeSplitUpdater{
 public:
  short int state;
 GammaMergeSplitUpdater( vector<double>& hpar, DistributionList* dists, GammaUpdater* par):MergeSplitUpdater(hpar, dists, (ParameterUpdater*)(par)){state = par->state;}

  
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double LogPriorRatio(double s1, double s2, double s, int mergesplit);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
 private:
  double plusstatecount_1, plusstatecount_2, plusstatecount_12;
  double minusstatecount_1, minusstatecount_2, minusstatecount_12;
};

class NuMergeSplitUpdater: public MergeSplitUpdater{
 public:
  short int state;
  NuMergeSplitUpdater( vector<double>& hpar, DistributionList* dists, NuUpdater* par);
  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  double LogPriorRatio(double s1, double s2, double s, int mergesplit);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
  ~NuMergeSplitUpdater(){delete sampler;}
 private:
  BoundedRescaledSampler* sampler;
  vector<double> sampler_pars;
  // (1) sample size alpha + beta 
  // (2) exponent for scaling the sample size higher for higher values

};


class InfoMergeSplitUpdater: public MergeSplitUpdater{
 public:
 InfoMergeSplitUpdater(vector<double>& hpar, DistributionList* dists):MergeSplitUpdater(hpar, dists){par_name="MergeSplit_Info";}

  void MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars);
  void PrintValue(vector<Machine*>& machines, ostream& str=cout);
  
};
#endif

