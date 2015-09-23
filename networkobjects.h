#ifndef NETWORKOBJECTS
#define NETWORKOBJECTS
#include "datparobjects.h"
#include "poissonEM.h"
#define BASELINE_LAG 15
extern std::ostream cnull;

int ipdist(unsigned int a, unsigned int b);

class Network{
  /*A Network is a collection of observed data through IP addresses.  It contains the 
    allocated data observed, and points to all the machines that currently are theorized to 
    communicate through it */
public:
  netw_dat_obj* data;
  netw_par_obj* parameters;
  
  
  Network(int N, int rid, netw_par_obj* netw_pars);
  Network(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo);
  Network(int N, int rid, ifstream &spikeinfo , ifstream &flinfo, ifstream &ipinfo, netw_par_obj* netw_pars);
  Network(Network* cpy);
  void Print(ostream& str, int N=-1);
  bool operator<(Network rhs);
  ~Network();
};

class Machine{
  
  /* A Machine is a collection of state parameters and currently allocated data representing a single 
     machine in the network.  It is linked to observed data through network points. */
 public:
  int uid;
  mach_dat_obj* data;
  curr_par_obj* parameters;
  int* par_acc;
  int* par_total;
  int num_par;
  int par_acc_allocated;

  Machine(Network* initial_network, double* psi_params, double q_start, double alpha_start, double omega_start, double* rho_start, double* gamma_start, double * nu_start, double off_min, double off_max, int immune);
  Machine(Network* initial_network, curr_par_obj* initial_pars, hype_par_obj* hyperpars);
  Machine(Machine* cpy, int iteration, int get_new_id);
  void Print(ostream& str, int N=-1);
  int RefineStartValuesUsingEM(ostream& printstream=cnull);
  void ResetIDs();
  double GetRate(int t); //return the poisson count rate at time t
  ~Machine();

 private:
  static machine_id_t id;
  machine_id_t getID();
};

class YTotalActiveStats{

 public:
  YTotalActiveStats(int bl=BASELINE_LAG);
  void Initialize(vector<Machine*>& machines);
  double YCDF(int y);
  double YDensity(int y);
  int y_min;
  int y_max;
  int y_num_active;
  ~YTotalActiveStats(){}
 
 private:
  vector<int> density;
  vector<int> cum_sum;
  int baseline_lag;

};


class YStateAndCountStats{

 public: 
  YStateAndCountStats(int bl = BASELINE_LAG){baseline_lag = bl;}
  void Initialize(vector<Machine*>& current_machines, vector<Machine*>& proposed_machines);
  //order is 1 only, 2 only, both 1 and 2, and 12
  vector<double> spikecountmean;
  vector<double> spikecountssd;
  vector<int> spikecounttotal;
  vector<double> baselineonlymean;
  vector<double> baselineonlyssd;
  vector<int> baselineonlytotal;

  //yes this is redundant.  I don't care at the moment
  vector<double> spikesums;
  vector<double> spikecounts;

  ~YStateAndCountStats(){}
 private: 
  int baseline_lag;
  vector<double> max_lag_s1_s2_s12;
    
};


#endif
