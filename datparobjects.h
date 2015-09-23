#ifndef DATPAROBJECTS
#define DATPAROBJECTS
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include <set> 
#include <gsl/gsl_heapsort.h> 
#include "maintenance_functions.h"

#define KEEP_TEST 1

using namespace std; 
typedef map<short int, int> network_inner_t;
typedef map<unsigned int, network_inner_t> network_outer_t;
typedef long long unsigned int machine_id_t;

const string ippath_name = "/Users/rweaver/Research/data/Conficker/matrix.allipcount.0305-0424.txt";
const string flpath_name = "/Users/rweaver/Research/data/Conficker/matrix.allflcount.0305-0424.txt";
const string inpath_name = "/Users/rweaver/Research/data/Conficker/dat_info.txt";

class hype_par_obj{
 public:
  int off_min;
  int off_max;

  hype_par_obj(ifstream &hyper_par_file);
  void Print(ostream& str){str << "off min, max = " << off_min << "," << off_max << endl;};
  ~hype_par_obj(){};
};

class netw_dat_obj{
 public:
  vector<int> fl;           /* Flow counts associated with data object */
  vector<int> ip;           /* IP address counts associated with object */
  set<machine_id_t> machines; /*Machine IDs associated with this network*/
  int T;              /* Total time of recording */
  short int beginval; /* first nonzero record*/
  short int endval;   /* last nonzero record*/
  short int timezone; /* time zone estimated from prior analyses*/
  short int tstar;    /* "Spike" offset estimated from prior analyses*/
  int row_id;         /* associated with which /24 in the big data set*/ 
  unsigned int block_id; /*Integer representation of the /24 */
  int is_allocated;   /* 1 if fl, ip allocated, 0 otherwise*/
  int total_objects;  /* Total number of objects (IPs) */

  /* constructors */
  /* Copy constructor */
  netw_dat_obj(netw_dat_obj* cpy);
  
  /* File object pointers for processing many lines of input in a row 
     N = size of arrays to make for fl, ip (number of columns in data matrix)
     rid = row id of data to get in the big data matrices (for identification purposes)
     spikeinfo, flinfo, ipinfo = file stream pointers from where the file has already been opened     
  */
  netw_dat_obj(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo);
  
  /*For obtaining data from just one file stream (uses head/tail command line) 
    N = size of arrays to make for fl, ip (number of columns in data matrix)
    rid = row id to get from the main data files 
   */
  netw_dat_obj(int N, int rid);

  /*Printing up to N lines of data (no argument for all)*/
  void Print(ostream &str, int N=-1);

  /*Destructor*/
  ~netw_dat_obj();

 private:
  /*Helper function for reading in data*/
  void GetData(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo);

};


class NetworkBlock{
 public:
 
  NetworkBlock(unsigned int uid, vector<short int>& indexes, vector<int>& totals); //Assume that all observed traffic is coming from the same single network point
  NetworkBlock(unsigned int uid, vector<int>& totals);  //assume that the indexes are equal to 0-->totals.size()
  NetworkBlock(NetworkBlock* cpy);//Copy constructor

  void Add(NetworkBlock* a); 
  void AddNetwork(unsigned int uid);
  int Subtract(NetworkBlock* a);
  int RemoveZeroNetwork(unsigned int uid);
  void Print(ostream &str, int N=-1);
  int isNonZeroNetwork(int netid);
  void ListNetworks(vector<unsigned int>* network_list);
  void SumBlocks(vector<int>& totals);
  int NumNetworks();

  //Both of these check on the existing network list so as to not add a new network without calling "AddNetwork"
  int GetCount(unsigned int uid, short int t);
  void SetCountOnExistingNetwork(unsigned int uid, short int t, int val); 
  ~NetworkBlock(){};  
  
 private:
   network_outer_t block_allocations;

};

class mach_dat_obj{
 public:

  int T;              /* Total time of recording */
  short int beginval; /* first nonzero record*/
  short int endval;   /* last nonzero record*/
  short int timezone; /* time zone estimated from prior analyses*/
  short int tstar;    /* "Spike" offset estimated from prior analyses*/ 
  machine_id_t uid;
  int iteration_born;

  /* constructors */
  /* Copy constructor */
  mach_dat_obj(mach_dat_obj* cpy);
  
  /* File object pointers for processing many lines of input in a row 
     N = size of arrays to make for fl, ip (number of columns in data matrix)
     rid = row id of data to get in the big data matrices (for identification purposes)
     spikeinfo, flinfo, ipinfo = file stream pointers from where the file has already been opened     
  */
  mach_dat_obj(netw_dat_obj* cpy, machine_id_t id, int iteration);
  mach_dat_obj(mach_dat_obj* cpy, machine_id_t id, int iteration);  

  /*Printing up to N lines of data (no argument for all)*/
  void Print(ostream &str, int N=-1);
  
  int GetTotalCount(int t);
  int GetSubCount(unsigned int networkid, int t);
  void SetCountExisting(unsigned int networkid, int hour, int count);
  int SubtractBlocks(mach_dat_obj* machine);
  void AddBlocks(mach_dat_obj* machine);
  void ListNetworks(vector<unsigned int>* network_list);
  void AddNetwork(unsigned int networkid);
  int RemoveZeroNetwork(unsigned int networkid);
  int NumNetworks();

  /*Destructor*/
  ~mach_dat_obj();
 private:
  vector<int> fl;     /* Flow counts associated with entire data object */
  NetworkBlock* networks; /*Sub-allocations when more than one network are involved */  
  int fl_nb_agree;    /* 1 if fl is currently equal to the sum of all network blocks, 0 otherwise */
  /*It will be set up as: the network block structure is what will be changed by the other pieces.  The internal flow vector
    will only get updated to reflect what's in the network block structure, ie the flow vector is the thing that can be 
    out of date */
  void SumBlocks();
  int is_allocated;   /* 1 if fl allocated, 0 otherwise*/
 
};


class curr_par_obj{
public:

  /*local parameters*/
  short int State(int i){return state[i];}  /* [0,1,2] for off, spike, decay */ 
  double Lagstate(int i){return lagstate[i];} /* Representation of state as X from spike (Infinity if following just a baseline) */

  int T;            /* Total length of observations */
  int tstar;        /* "spike offset" of network associated with machine */

  short int Birthdate(){return birthdate;} /* Time "born" (either 0 or 283, default 0) */
  short int Deathdate(){return deathdate;} /* Time "died" (default 1224)*/
  short int Is_immune(){return is_immune;}
  int BDIindex();
  int NumOffOffTrans(){return num_off_off;}
  int NumOffSeq(){return num_off_state_seq;}
  double alpha;       /* decay rate */
  double q;           /* baseline */
  double omega;       /* baseline multiplier for spike */
  double rho[3];      /* amplitude parameter for state->state periodicity */
  double gamma[3];    /* conditional choice for switch */
  double nu[3];       /* Moving baseline up/down */
  double psi_parameterization[3][3][2]; /*parameterization of the rho/gamma/nu products in 3x3 transition matrix*/
  double off_lambda;  /* Poisson parameter for the average number of hours machine stays off given turned off */
  short int is_allocated; /* 0 if state/lagstate aren't allocated, 1 otherwise */

  /* Constructors */
  /* empty constructor */
  curr_par_obj();
  /* data object and starting values for parameters */
  curr_par_obj(mach_dat_obj* data, double* psi_params, double q_start, double alpha_start, double omega_start, double* rho_start, double* gamma_start, double * nu_start, double off_min, double off_max, int immune);
 
  /*Constructor for setting up a par object within a machine that is spawned from a network*/ 
  curr_par_obj(mach_dat_obj* data, curr_par_obj* cpy, hype_par_obj* hyper);
 /* Constructor for reading in inital parameter values only from a file (psi parameterization, q, alpha, omega, rho, gamma, nu */ 
  curr_par_obj(ifstream& init_pars);

  /* Copy constructor */
  curr_par_obj(curr_par_obj* cpy);

  /* Reallocate and set */
  void SetFromCopy(curr_par_obj* cpy);
  
  /* Set the state (and automatically update the lagstate value) */
  void SetState(int index, short int val);
  void SetBirthDeathImmune(short int birth, short int death, short int immune);
  /*3 row-wise 3x2  matrices for off, spike, decay */
  /*rn00, g00, rn01, g01, rn02, g02, rn10, g10,... etc*/
  void SetPsiParameterization(double* pars);	      

 /*Set the starting value of states using a 95% cutoff for spike vs decay, and '3 or more 0s' in a row" for the first off state */
  void SetStateStart(mach_dat_obj* data, double off_min, double off_max);
 

  void Print(ostream& str, int N=-1);
  void PrintState(ostream &fileextension);
  void PrintPars(ostream &fileextension, vector<int>* is_estimated=NULL);
 
  //return the poisson count rate at time t
  double GetRate(int t);


  /* Destructor */
  ~curr_par_obj();

 private:
  void InternalSet(curr_par_obj* cpy);
  vector<short int> state;
  vector<double> lagstate; 
  short int birthdate;
  short int deathdate;
  short int is_immune;
  int num_off_off;
  int num_off_state_seq;
  double SetIndividualLagState(int i);
  void SetLagstate();
  void AdjustSingleOffState(int index, short int oldstate, short int newstate);

};


class netw_par_obj{

 public:
  double churn;
  int type;
  int boundary_id;
  netw_par_obj(ifstream &netw_pars);
  netw_par_obj(netw_par_obj* cpy);
  netw_par_obj();
  void Print(ostream &str);
  ~netw_par_obj(){};
 private:
  static int id;
  int getBoundaryID();
};

class global_par_obj{

 public:
 
  int T;
  double immune_rate;
  double* survival_rate;
  int num_immune;
  int num_not_immune;
  int* num_nonoff_nonimmune_transitions;
  int* num_deaths;
  double immune_hyper[2];
  int survival_hyper_mean_position[5];
  double survival_hyper_mean_height[5];
  double survival_hyper_temperature;
  int num_machines;
  map<int, int*> boundary_info;  //ips and counts for each boundary id
  global_par_obj(ifstream& inputfile);
  void Print(ostream &str, int iteration_id);
  void UpdateGlobalImmunityAndSurvivalCounts(curr_par_obj* parameters, int add_or_delete);
  double SurvivalRatePriorMean(int t);
  
  ~global_par_obj();
};

#endif
