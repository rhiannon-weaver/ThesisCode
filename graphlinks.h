#ifndef GRAPHLINKS
#define GRAPHLINKS
#include <iomanip>
#include "networkobjects.h"
#include "distributionlist.h"
#include "mcmcdrivers.h"
#include "countupdaters.h"
#define NUM_VALID_MOVES 3.0
using namespace std;

const string output_root = "output/";

class ParameterInitializer{

 public:
  hype_par_obj* hyper_par_init;
  netw_par_obj* network_par_init;
  curr_par_obj* machine_par_init;
 
  
  ParameterInitializer(const char* networkfile, const char* machinefile, const char* hyperparfile);
  void Print(ostream& str);
  ~ParameterInitializer();
  /*A container to hold all of the parameter initializations read in from files  */

};

bool NetworkPtrSortPredicate( Network* a, Network* b);

class NetworkGraph{
 public: 
  vector<Network*> network_touchpoints;
  map<machine_id_t, Machine*> machines;
  vector<Network*>::iterator network_iter;
  map<machine_id_t, Machine*>::iterator machine_iter;
  ParameterInitializer* initial_parameters;
  global_par_obj* global_parameters;
  int* global_acc;
  int* global_total;
 
  //Constructor
  NetworkGraph(int N, const char* dataset, const char* netw_init_parfile, const char * mach_init_parfile, const char* hype_init_parfile, const char * globalfile, int SpawnNInitialperNetwork, int print=0);

 //Various print functions
  void PrintNetwork(ostream &str, int id, int N=-1);
  void PrintMachine(ostream &str, machine_id_t id, int N=-1);
  void PrintNetworks(ostream &str, int id_only, int N=-1);
  void PrintMachines(ostream &str, int id_only, int N=-1);
  void PrintNetworkGraph(ostream &str, int id_only, int N=-1);   
  void PrintGraphViz(const char* filename); 
 
  //Move types
  SingleModelUpdater* machineMCMC;
  CountRearranger* countrearrangeMCMC;
  vector<ParameterUpdater*> globalMCMC;
  MachineMergeSplit* machineMergeSplit;

  //Initialize and Run
  void InitializeMCMC(const char* initializationfile);
  void RunMCMC(int iterations);

  void SetAllBirthDeathImmune(short int birth, short int death, short int immune);
  void SetAllDecayRates(double decay_rate=0.0);
  void ReinitializeMCMC(int num_machines_per_network, const char* initializationfile, int print);

  //Insertion/deletion, selection, and helper functions
  void InsertNetwork(int N, int rid,  ifstream& spikeinfo, ifstream &flinfo, ifstream& ipinfo);
  void InsertMachine(Machine* new_machine);
  void DeleteMachine(machine_id_t ex_machine);
  int GetMachineBoundaryID(Machine* machine);
  Machine* SpawnMachine(Network* network);
  void SpawnInitialMachines(Network* network, int num_per_network, int print);  
  Network* RandomlySelectNetwork();
  Machine* RandomlySelectMachine(int pair_provided, Machine* pair = NULL);

  //MCMC output functions (full graph)
  void PrintGraphLogLikelihood(int iteration_id, ostream &str = cnull);
  double GlobalLogLikelihood();
  double GenerationLogLikelihood(Machine* machine, int suppress_prior=0);
  //double AllocationLogLikelihood(ConnectedComponent* component);
  //double ComponentLogLikelihood(ConnectedComponent* component);
  //the above are for the model that we probably won't get to.
  void PrintActivityCounts(int iteration_id, ostream &str = cnull);
  void PrintInfectionCounts(int iteration_id, ostream &str = cnull);
  void SetAllBoundaries(int id);
 
  ~NetworkGraph();
  
 private:

  void Initialize1to1MCMC(const char* directory, const char* repeatfile, const char* hyperparfile, int append_existing, int write_pars, int write_state, int refine_initial_start_values);
  void Write1to1Info();
  void InitializeNoninformativeMCMC( double mix_prob,int count_rearrange_hour_repeat, int count_rearrange_machine_repeat, int count_rearrange_thin, int count_rearrange_switch_type, int merge_rearrange_repeat, int merge_rearrange_thin, const char* merge_parameter_file, const char* boundary_file); 
  void InitializeInformativeMCMC(const char * netw_hyperpar_file);
  void SetOutputDirectory(const char* directory, int append_existing, int refineEM);
  void MCMCCleanUp();  

  //This is the 1-to-1 run option
  
  void RunAllMachinesMCMC(int iterations);
  void RunNonInformativeMCMC(int iterations);
  void RunInformativeMCMC(int iterations);


  //This is the noninformative or informative option
  void RunMachineMCMC(Machine* machine, int iterations,int thin, int export_acc, int iteration_id, int print=1);
  void RandomSweep(int iteration_id);
  void MergeSplitStep(int iteration_id, int print_to_screen, ostream& output, int write_output, int print_all_acceptance);

  
  DistributionList* distributions;
  double* network_weights;
  gsl_ran_discrete_t * network_weight_preproc;
  string hyperparfile;
  int network_weight_preproc_allocated;
  int preproc_up_to_date;  
  void RemoveDuplicateNetworks();
  void PreProcessWeights();			
  void InitializeBoundaryCounts(ostream& str, int print_to_screen);
 
 
  map<unsigned int, int> network_id_to_index;
  ofstream  info_output;
  ofstream  parameter_output;
  ofstream  state_output;
  ofstream  acc_output;
  ofstream  global_output;
  ofstream  parameter_EM_output;
  ofstream  count_rearrange_output;
  ofstream  count_stats_output;
  ofstream  graph_posterior_output;
  ofstream  active_size_output;
  ofstream  infected_size_output;
  ofstream  merge_split_output;
  int print_to_screen;
  int machine_thin_parameter;
  int merge_split_thin_parameter;
  int switch_type;
  int append;
  double logpost_print_prob;
  string indexfilename;
  string estimation_type;
  int split_global;
  int MCMC_1to1_initialized;			
  int MCMC_noninformative_initialized;
  int MCMC_informative_initialized;

};
#endif
