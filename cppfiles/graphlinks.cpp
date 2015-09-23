#include "../graphlinks.h"

/*==========Parameter Initializer and Sorting Functions=============*/

bool NetworkPtrSortPredicate(Network* a, Network* b)
{
  return( a->data->block_id < b->data->block_id);
}

ParameterInitializer::ParameterInitializer(const char* networkfile, const char* machinefile, const char * hyperparfile)
{
  ifstream netw_par;
  ifstream mach_par;
  ifstream hype_par;

  netw_par.open(networkfile, ios::in);
  mach_par.open(machinefile, ios::in);
  hype_par.open(hyperparfile, ios::in);
 
  network_par_init = new netw_par_obj(netw_par);
  hyper_par_init = new hype_par_obj(hype_par);
  machine_par_init = new curr_par_obj(mach_par);

  netw_par.close();
  hype_par.close();
  mach_par.close();

}

void ParameterInitializer::Print(ostream &str)
{
  hyper_par_init->Print(str);
  network_par_init->Print(str);
  machine_par_init->Print(str, 0);
 
}

ParameterInitializer::~ParameterInitializer()
{
  delete network_par_init;
  delete hyper_par_init;
  delete machine_par_init;
  
}

/*==========End Parameter Initializer and Sorting Functions=============*/ 


/*===========Network Graph Constructor and Print functions================*/

NetworkGraph::NetworkGraph(int N, const char* dataset, const char * netw_init_parfile, const char * mach_init_parfile, const char* hype_init_parfile, const char* globalfile, int SpawnNInitialperNetwork, int print)
{
  int i;
  int rid;
  ifstream ind_file;
  ifstream spk_file;
  ifstream flw_file;
  ifstream ipa_file;
  ifstream glb_file;

  string indexfile = string(dataset) + "Index.txt";
  string spikeinfofile = string(dataset) + "Info.txt";
  string flowinfofile = string(dataset) + "FL.txt";
  string ipinfofile = string(dataset) + "IP.txt";

  hyperparfile = hype_init_parfile;
  Machine* current_machine;
  indexfilename = indexfile;
  double hyperpar[2];
  ind_file.open(indexfile.c_str(), ios::in);
  spk_file.open(spikeinfofile.c_str(), ios::in);
  flw_file.open(flowinfofile.c_str(), ios::in);
  ipa_file.open(ipinfofile.c_str(), ios::in);
  glb_file.open(globalfile, ios::in);  
  
  /* Get initial hyperparameters and such */
  initial_parameters = new ParameterInitializer(netw_init_parfile, mach_init_parfile, hype_init_parfile);
  if(print)
    initial_parameters->Print(cout);
  
  distributions = new DistributionList();
  global_parameters = new global_par_obj(glb_file);
  /*Set up all of the networks first. */
  cout <<"Adding Networks" << endl;
  cout.flush();
  while(ind_file >> rid)
    {
      /*cout << rid << endl;
	cout.flush();*/
      /* Read in the index number from the index file */
      InsertNetwork(N, rid, spk_file, flw_file, ipa_file);
    }
  /*now sort the networks and make sure there are no duplicates*/
  sort(network_touchpoints.begin(), network_touchpoints.end(), NetworkPtrSortPredicate);
  RemoveDuplicateNetworks();
  for(i = 0; i < network_touchpoints.size(); i++)
    {
      network_id_to_index[network_touchpoints[i]->data->block_id] = i;
    }
  network_weights = new double[network_touchpoints.size()];
  network_weight_preproc_allocated = 0;
  preproc_up_to_date = 0;
  
 
  /*PrintNetworks(100);*/
  
  /*Now iterate across all the networks and spawn machines, one machine per max non-zero simultaneous IP addresses seen */
  
  cout << "Spawning Initial Machines" << endl;
  cout.flush();
  for(network_iter = network_touchpoints.begin(); network_iter != network_touchpoints.end(); network_iter++)
    {
      SpawnInitialMachines(*network_iter, SpawnNInitialperNetwork, print);
    }
  global_parameters->num_machines=machines.size();

  MCMC_1to1_initialized = 0;
  MCMC_noninformative_initialized = 0;
  MCMC_informative_initialized = 0;
  ind_file.close();
  spk_file.close();
  flw_file.close();
  ipa_file.close();  
}

void NetworkGraph::PrintNetwork(ostream& str, int id,  int N)
{
  map<unsigned int, int>::iterator lookup_iter = network_id_to_index.find(id);
  if(lookup_iter != network_id_to_index.end())
    network_touchpoints[lookup_iter->second]->Print(str, N);
  else
    cout << "Error (NetworkGraph::PrintNetwork): Network " << id << " not found in graph" << endl;
}

void NetworkGraph::PrintMachine(ostream& str, machine_id_t id, int N)
{
  machine_iter = machines.find(id);
  if(machine_iter != machines.end())
    machine_iter->second->Print(str, N);
  else
    cout << "Error (NetworkGraph::PrintNetwork): Machine " << id << " not found in graph" << endl;
}

void NetworkGraph::PrintNetworks(ostream& str, int id_only, int N)
{
  set<machine_id_t>::iterator machine_id_iter;
  int i;
  if(id_only)
    str << "  Network   | Weight | Machine IDs" << endl;

  i = 0;
  for(network_iter = network_touchpoints.begin(); network_iter != network_touchpoints.end(); network_iter++)
    { 
      if(!id_only)
	(*network_iter)->Print(str, N);
      else
	{
	  str << setw(12) << (*network_iter)->data->block_id << "|" << setw(8) << fixed << setprecision(2) << network_weights[i] << "|";
	  for(machine_id_iter = (*network_iter)->data->machines.begin(); machine_id_iter != (*network_iter)->data->machines.end(); machine_id_iter++)
	    str << " " << (*machine_id_iter);
	  str << endl;
	  i++;
	}
      
    }
}

void NetworkGraph::PrintMachines(ostream& str, int id_only, int N)
{
  int i;
  vector<unsigned int> network_ids;
  for(machine_iter = machines.begin(); machine_iter != machines.end(); machine_iter++)
    {
      if(!id_only)
	machine_iter->second->Print(str, N);
      else
	{
	  str << machine_iter->second->uid << ":";
	  machine_iter->second->data->ListNetworks(&network_ids);
	  for(i = 0; i < network_ids.size(); i++)
	    {
	      str << " " << network_ids[i];
	    }
	  str << endl;	    
	}
    }
}

void NetworkGraph::PrintNetworkGraph(ostream& str, int id_only, int N)
{
  str << "%%%%%%%%%%Complete Network Graph%%%%%%%%%" << endl << "Networks" << endl << "--------" << endl;

  PrintNetworks(str, id_only, N);
  str << "Machines" << endl << "--------" << endl;
  PrintMachines(str, id_only, N);
  
  str << "%%%%%%%%End Complete Network Graph%%%%%%%" << endl;
}

void NetworkGraph::PrintGraphViz(const char* filename)
{
  cout << "This will eventually print a graphviz input file" << endl;
}


/*============End Constructor and Print functions===============*/


/*============Public Initialize and Run functions for MCMC======*/

void NetworkGraph::RemoveDuplicateNetworks()
{
  int cur_block_id;
  int num_deleted = 0;
  network_iter = network_touchpoints.begin();
  cur_block_id = (*network_iter)->data->block_id;
  network_iter++;
  while(network_iter != network_touchpoints.end())
    {
      if((*network_iter)->data->block_id == cur_block_id)
	{
	  delete *network_iter;
	  network_iter = network_touchpoints.erase(network_iter);
	  num_deleted++;
	}
      else
	{
	  cur_block_id = (*network_iter)->data->block_id;
	  network_iter++;
	}
    }
  if(num_deleted > 0)
    cout << "Deleted " << num_deleted << " duplicate networks" << endl;
  else
    cout << "All networks unique; no duplicates found " << endl;
}

void NetworkGraph::SpawnInitialMachines(Network* network, int num_per_network, int print)
{
  Network* tmp_network;
  Machine* current_machine;
  vector<unsigned int*> flow_breakdowns;
  double * pk;
  int max_simultaneous = 0;
  int i,j;
  DistributionList distributions(0);
  
  if(num_per_network > 0)
    {
      max_simultaneous = num_per_network;
    }
  else
    {
      for(i = 0; i < network->data->T; i++)
	{
	  if(network->data->ip[i] > max_simultaneous)
	    max_simultaneous = network->data->ip[i];
	}
    }

  if(max_simultaneous == 1)
    {
      current_machine = SpawnMachine(network);
      if(print)
	{
	  current_machine->Print(cout, 100);
	}
      
      if(print)
	{
	  network->Print(cout, 100);
	} 
    }
  else
    {
      /* split up the observed flows amongst the max_simultaneous machines. At this point we don't care if the numbers for
	 simultaneous IP addresses don't 'add up' because there could be NATs (lots of machines going through 1 IP) or there could
	 be some machines switching in that hour (eg, one machine gets all the flows but is actually seen through like 4 IPs).  That
	 second situation at least will be "fixed" in the reallotment phase of the MCMC. 
      */
      if(print)
	{
	  cout <<"max simultaneous = " << max_simultaneous << endl;
	  cout <<"T = " << network->data->T << endl;
	}
      pk = new double[max_simultaneous];
      for(i = 0; i < max_simultaneous; i++)
	pk[i] = 1.0/((double)max_simultaneous);
      for(i = 0; i < network->data->T; i++)
	{
	  flow_breakdowns.push_back(new unsigned int[max_simultaneous]);
	  if(network->data->fl[i] > 0)
	    {
	      distributions.Simulate(8, 1, flow_breakdowns[i], (double)(network->data->fl[i]), (double)(max_simultaneous), pk);
	      /* cout << "Total = " << network->data->fl[i] << ": as ( ";
	      for(j = 0; j < max_simultaneous; j++)
		cout << flow_breakdowns[i][j] << " ";
		cout << ")" << endl;*/
	    }
	  else
	    for(j = 0; j < max_simultaneous; j++)
	      flow_breakdowns[i][j] = 0;
	}
      
      /*  
	  Now copy each of those flow allocations into a new network and use it to spawn a new machine
      */
      tmp_network = new Network(network);
      if(print)
	{
	  tmp_network->Print(cout);
	}
      for(i = 0; i < max_simultaneous; i++)
	{
	  for(j = 0; j < tmp_network->data->T; j++)
	    {
	      network->data->fl[j] = flow_breakdowns[j][i]; 
	    }
	  current_machine = SpawnMachine(network);
	  if(print)
	    {
	      current_machine->Print(cout, 100);
	    }
	}
      if(print)
	{
	  network->Print(cout, 100);
	}

      for(i = 0; i < network->data->T; i++)
	{
	  network->data->fl[i] = tmp_network->data->fl[i];
	}

      for(i = 0; i < network->data->T; i++)
	delete flow_breakdowns[i];
      delete tmp_network; 
      delete pk;
    }
}


void NetworkGraph::RunMCMC(int iterations)
{
  if(!MCMC_1to1_initialized)
    cout << "Call InitializeMCMC before calling RunMCMC" << endl;
  else
    {
      if(estimation_type == "1to1")
	{
	  cout << "Running 1 to 1 model for " << iterations << " iterations" << endl;
	  RunAllMachinesMCMC(iterations);
	}
      if(estimation_type == "Noninformative")
	{
	  cout << "Running Noninformative Network model for " << iterations << " iterations" << endl;
	  RunNonInformativeMCMC(iterations);
	}
      if(estimation_type == "Informative")
	{
	  //RunInformativeMCMC(iterations);
	}
    }
}

void NetworkGraph::RunMachineMCMC(Machine* machine, int iterations,int thin, int export_acc, int iteration_id, int print)
{
  int print_out;
  if(print == 1)
    print_out = print_to_screen;
  else
    print_out = 0;
  machineMCMC->UpdateMachine(machine, global_parameters, iterations, thin, parameter_output, state_output, acc_output, export_acc, iteration_id, print_out);
}

void NetworkGraph::RunAllMachinesMCMC(int iterations)
{ 
  int j;
  int i;
  for(j = 0; j< globalMCMC.size(); j++)
    {
      global_acc[j] = 0;
      global_total[j] = 0;
    }
  info_output << "RunAllMachinesMCMC used, switch_type=" << switch_type << endl;
  for(j = 0; j < iterations; j++)
    {
      if((j % 1000) == 0)
	{
	  cout << j << endl;
	}
      for(machine_iter = machines.begin(); machine_iter != machines.end(); machine_iter++)
	{
	  RunMachineMCMC(machine_iter->second, 1, (j%machine_thin_parameter==0), (j == (iterations -1)), -1);
	}
      for(i = 0; i < globalMCMC.size(); i++)
	{
	  global_total[i]++;
	  global_acc[i]+= (int) globalMCMC[i]->MCMC_Step(NULL, NULL, global_parameters, j, print_to_screen);
	  if( (j % machine_thin_parameter) == 0)
	    {
	      global_parameters->Print(global_output, -1);
	      if(graph_posterior_output.is_open())
		PrintGraphLogLikelihood(j, graph_posterior_output);
	    }
	}
    }
  /*print acceptance values for global totals*/
  for(j = 0; j < globalMCMC.size(); j++)
    {
      acc_output << (double(global_acc[j]))/(double(global_total[j])) << " ";
    }
  acc_output << endl;
}

void NetworkGraph::RunNonInformativeMCMC(int iterations)
{
  int i;
  
  for(i = 0; i < globalMCMC.size(); i++)
    {
      global_acc[i] = 0;
      global_total[i] = 0;
    }
  
  for(i = 0; i < iterations; i++)
    {
      RandomSweep(i);
    }

  //Run each machine one more time and export the acceptance ratios;
  for(machine_iter = machines.begin(); machine_iter != machines.end(); machine_iter++)
    {
      RunMachineMCMC(machine_iter->second, 1, 1, 1, iterations,0);
    }
  //Print the acceptance ratios for the global parameters
  if(acc_output.is_open())
    {
      for(i = 0; i < globalMCMC.size(); i++)
	{
	  acc_output << (double(global_acc[i]))/(double(global_total[i])) << " ";
	}
      acc_output << endl;
      machineMergeSplit->PrintAcceptance(acc_output);
    }
  //Print the log-likelihood values
  if(graph_posterior_output.is_open())
    PrintGraphLogLikelihood(iterations, graph_posterior_output);
}

void NetworkGraph::RunInformativeMCMC(int iterations)
{

}

double NetworkGraph::GenerationLogLikelihood(Machine* machine, int suppress_prior)
{
  double loglike;
  loglike = machineMCMC->GraphLogLikelihood(machine, global_parameters, suppress_prior);
  return loglike;
}

double NetworkGraph::GlobalLogLikelihood()
{
  double loglike = 0.0;
  double val;
  int i;
  for(i = 0; i < globalMCMC.size(); i++)
    {
      val = globalMCMC[i]->GraphLogLikelihood(NULL, global_parameters);
      //cout << globalMCMC[i]->par_name << ":" << val << " ";
      loglike += globalMCMC[i]->GraphLogLikelihood(NULL, global_parameters);
    }
  //cout << endl;
  return loglike;
}

void NetworkGraph::PrintActivityCounts(int iteration_id, ostream &str)
{
  int T;
  int i,j;
  map<int, int*> boundary_counts;
  map<int, int*>::iterator boundary_count_iter;
  set<machine_id_t> ids_visited;
  set<machine_id_t>::iterator id_iterator;
  set<machine_id_t>::iterator ids_visited_iterator;

  T = network_touchpoints[0]->data->T;
  for(i = 0; i < network_touchpoints.size(); i++)
    {
      boundary_count_iter = boundary_counts.find(network_touchpoints[i]->parameters->boundary_id);
      if(boundary_count_iter == boundary_counts.end())
	{
	  boundary_counts[network_touchpoints[i]->parameters->boundary_id] = new int[T];
	  for(j = 0; j < T; j++)
	    boundary_counts[network_touchpoints[i]->parameters->boundary_id][j] = 0;
	}
	  //iterate through machines belonging to this network, skipping the ones you've seen before
      for(id_iterator = network_touchpoints[i]->data->machines.begin();
	  id_iterator != network_touchpoints[i]->data->machines.end();
	  id_iterator++)
	{
	  if(ids_visited.find((*id_iterator)) == ids_visited.end())
	    {
	      ids_visited.insert((*id_iterator));
	      machine_iter = machines.find((*id_iterator)); 
	      for(j = 0; j < T; j++)
		boundary_counts[network_touchpoints[i]->parameters->boundary_id][j] += int(machine_iter->second->parameters->State(j) > 0);
	    }
	}
      
    }
  
  //Print out active counts
  
  for(boundary_count_iter = boundary_counts.begin(); boundary_count_iter != boundary_counts.end(); boundary_count_iter++)
    {
      str << iteration_id << " " << boundary_count_iter->first << " ";
      for(j = 0; j < T; j++)
	str << (boundary_count_iter->second)[j] << " ";
      str << endl;
      delete boundary_count_iter->second;
    }
  return;
}



void NetworkGraph::PrintInfectionCounts(int iteration_id, ostream &str)
{
  int T;
  int i,j;
  map<int, int*> boundary_counts;
  map<int, int*>::iterator boundary_count_iter;
  set<machine_id_t> ids_visited;
  set<machine_id_t>::iterator id_iterator;
  set<machine_id_t>::iterator ids_visited_iterator;

  T = network_touchpoints[0]->data->T;
  for(i = 0; i < network_touchpoints.size(); i++)
    {
      boundary_count_iter = boundary_counts.find(network_touchpoints[i]->parameters->boundary_id);
      if(boundary_count_iter == boundary_counts.end())
	{
	  boundary_counts[network_touchpoints[i]->parameters->boundary_id] = new int[T];
	  for(j = 0; j < T; j++)
	    boundary_counts[network_touchpoints[i]->parameters->boundary_id][j] = 0;
	}
	  //iterate through machines belonging to this network, skipping the ones you've seen before
      for(id_iterator = network_touchpoints[i]->data->machines.begin();
	  id_iterator != network_touchpoints[i]->data->machines.end();
	  id_iterator++)
	{
	  if(ids_visited.find((*id_iterator)) == ids_visited.end())
	    {
	      ids_visited.insert((*id_iterator));
	      machine_iter = machines.find((*id_iterator)); 
	      for(j = machine_iter->second->parameters->Birthdate(); j < machine_iter->second->parameters->Deathdate(); j++)
		boundary_counts[network_touchpoints[i]->parameters->boundary_id][j]++;
	    }
	}
      
    }
  
  //Print out active counts
  
  for(boundary_count_iter = boundary_counts.begin(); boundary_count_iter != boundary_counts.end(); boundary_count_iter++)
    {
      str << iteration_id << " " << boundary_count_iter->first << " ";
      for(j = 0; j < T; j++)
	str << (boundary_count_iter->second)[j] << " ";
      str << endl;
      delete boundary_count_iter->second;
    }
  return;
}


void NetworkGraph::PrintGraphLogLikelihood(int iteration_id, ostream &str)
{
  int i;
  double loglike;
  double val;
  set<machine_id_t> ids_visited;
  set<machine_id_t>::iterator id_iterator;
  set<machine_id_t>::iterator ids_visited_iterator;

  loglike = GlobalLogLikelihood();
  if(iteration_id >= 0)
    str << iteration_id << " " << machines.size() << " "; 
  str << loglike << " ";
  for(i = 0; i < network_touchpoints.size(); i++)
    {
      for(id_iterator = network_touchpoints[i]->data->machines.begin();
	  id_iterator != network_touchpoints[i]->data->machines.end();
	  id_iterator++)
	{
	  if(ids_visited.find((*id_iterator)) == ids_visited.end())
	    {
	      ids_visited.insert(*id_iterator);
	      machine_iter = machines.find((*id_iterator)); 
	      val = GenerationLogLikelihood(machine_iter->second);
	      loglike += val;
	      str << network_touchpoints[i]->parameters->boundary_id << " " << val << " ";
	    }
	}
    }
  str << loglike << endl;
}




void NetworkGraph::InitializeMCMC(const char* initializationfile)
{
  /*
    initialization file should have the following
    Type = "1to1", "Noninformative" or "Informative"
 1to1pars   directory to store the results
            repeatfile for machine parameters
	    split_global
            append_existing (0/1) to append to current directory?
            write_pars (0/1) write machine parameters (y/n)
            write_state (0/1) write state parameters  (y/n)
            refine_initial_start_values (0/1): use Poisson EM to refine starting values
	    machine_thin_parameter (only print every X iterations from the machines)
	    print_to_screen

noninformative pars
            mix_prob: the chance of picking a straight vs. augmented count rearranging step
	    count_rearrange_hour_repeat
	    count_rearrange_machine_repeat
	    count_rearrange_thin
	    count_rearrange_switch_type
	    merge_rearrange_repeat
	    merge_rearrange_thin
	    logpost_print_prob
            merge_parameter_file (filename for the setup parameters for the merge file) (NULL for none)
            boundary file (filename that outlines the boundaries of NATs or DHCP regions for noninformative models (NULL for none, will default to the network names themselves)
	   
	    
informative pars
            network hyperparameter file

   */

  ifstream init_stream;
  string directory;
  string repeatfile;
  int append_existing;
  int write_pars;
  int write_state;
  int refine_initial_start_values;
  double mix_prob;
  int  count_rearrange_hour_repeat;
  int count_rearrange_machine_repeat;
  int  count_rearrange_thin;
  int count_rearrange_switch_type;
  int merge_rearrange_repeat;
  int  merge_rearrange_thin;
  string merge_parameter_file;
  string boundary_file;
  string netw_hyperpar_file;
  
  init_stream.open(initializationfile, ios::in);
  init_stream >> estimation_type;
  //Read in the 1-1 stuff no matter what;
  init_stream >> directory >> repeatfile;
  init_stream >> split_global >> append_existing >> write_pars >> write_state >> refine_initial_start_values >> machine_thin_parameter >> print_to_screen;
  
  //cout << "hyper par file = " << hyperparfile.c_str() << endl;
  Initialize1to1MCMC(directory.c_str(), repeatfile.c_str(), hyperparfile.c_str(), append_existing, write_pars, write_state, refine_initial_start_values);
  if(info_output.is_open())
    info_output << "Print values to screen (0/1): " << print_to_screen << endl;
  if(print_to_screen)
    cout << "Printing values to screen" << endl;
 
  if( (estimation_type == "Noninformative") || (estimation_type == "Informative"))
    {
      /*  mix_prob
            count_rearrange_hour_repeat
	    count_rearrange_machine_repeat
	    count_rearrange_thin
	    count_rearrange_switch_type
	    merge_rearrange_repeat
	    merge_rearrange_thin
            merge_parameter_file (filename for the setup parameters for the merge file) (NULL for none)
            boundary file (filename that outlines the boundaries of NATs or DHCP regions for noninformative models (NULL for none, will default to the network names themselves)
	    
       */
      
      init_stream >> mix_prob
		  >> count_rearrange_hour_repeat
		  >> count_rearrange_machine_repeat
		  >> count_rearrange_thin
		  >> count_rearrange_switch_type
		  >> merge_rearrange_repeat
		  >> merge_rearrange_thin
		  >> logpost_print_prob
		  >> merge_parameter_file
		  >> boundary_file; 
      InitializeNoninformativeMCMC(mix_prob,count_rearrange_hour_repeat, count_rearrange_machine_repeat,count_rearrange_thin,count_rearrange_switch_type, merge_rearrange_repeat,merge_rearrange_thin, merge_parameter_file.c_str(), boundary_file.c_str()); 
     }
  if( (estimation_type == "Informative"))
    {
      init_stream >> netw_hyperpar_file;
      InitializeInformativeMCMC(netw_hyperpar_file.c_str());
    }
  init_stream.close();
  return;
}

void NetworkGraph::Write1to1Info()
{
  int i;
  if(!(info_output.is_open()))
    {
      cout << "Initialize MCMC before calling Write1to1Info()" << endl;
    }
  else
    {
      if(!append)
	{
	  info_output << "Estimation type: " << estimation_type << endl; 
	  info_output << "Data input index file: " << indexfilename << endl;
	  info_output << "Number of network points: " << network_touchpoints.size() << endl;
	  info_output << "Initial number of machines: " << machines.size() << endl;
	  info_output << "Estimation information" << endl;
	  machineMCMC->PrintRepeatValues(info_output);
	  info_output << "Hyperparameters" << endl;
	  machineMCMC->PrintAllHyperpars(info_output);
	  initial_parameters->machine_par_init->Print(info_output,0);
	  for(i=0; i < globalMCMC.size(); i++)
	    globalMCMC[i]->PrintHyperpars(info_output);
	}
    }
 
}

void NetworkGraph::InitializeInformativeMCMC(const char * netw_hyperpar_file)
{

}


void NetworkGraph::RandomSweep(int iteration_id)
{

  /* This is a list of move types, and what to do with each
     When a new move type is added, I just need to update the 
     switch statement */
  machine_iter = machines.begin();
  Machine* selected_machine;
  Network* selected_network;
  int MoveType;
  int i;
  double v;
  if(!print_to_screen)
    {
      if(iteration_id % 1000 == 0)
	cout << endl << iteration_id;
    }

  //
  distributions->Simulate(3,1,&v,0.0,NUM_VALID_MOVES);
  MoveType = int(v);
  distributions->Simulate(3,1,&v,0.0,1.0);
  switch(MoveType)
    {
    case 0: //Update a machine
      if(print_to_screen)
	{
	  cout << endl << "::Machine Parameters, Repeat 1 time::" << endl;
	  cout.flush();
	}
      //selected_machine = machine_iter->second;
      selected_machine = RandomlySelectMachine(0);
      RunMachineMCMC(selected_machine, 1, int(iteration_id%machine_thin_parameter==0), 0, iteration_id);
      break;

    case 1: //Update global machine parameters
      if(print_to_screen)
	{
	  cout << endl << "::Global Parameters, Repeat 1 time::" << endl;
	  cout.flush();
	}
      for(i = 0; i < globalMCMC.size(); i++)
	{
	  global_total[i]++;
	  global_acc[i]+= (int) globalMCMC[i]->MCMC_Step(NULL, NULL, global_parameters, iteration_id, print_to_screen);
	  if(v < logpost_print_prob)
	    {
	      global_output << iteration_id << " ";
	      global_parameters->Print(global_output, -1);
	    }
	}
      break;

    case 2: //Rearrange counts
      if(print_to_screen)
	{
	  cout << endl << "::Count Rearrange, Repeat "<<countrearrangeMCMC->GetRepeat() << " times::" << endl;
	  cout.flush();
	}
      //cout << "skipping for now"<< endl;
           for(i = 0; i < countrearrangeMCMC->GetRepeat(); i++)
      	{
	  selected_network = RandomlySelectNetwork();
	  countrearrangeMCMC->Sample(selected_network, machines, print_to_screen, iteration_id, count_rearrange_output,  int(iteration_id%machine_thin_parameter==0) );
        }
      break;

    case 3: //Merge or Split
      if(print_to_screen)
	{
	  cout << endl << "::Merge/Split, Repeat " << machineMergeSplit->num_repeats << " times::" << endl;
	  cout.flush();
	}
      //cout << "skipping for now"<< endl;

      for(i = 0; i < machineMergeSplit->num_repeats; i++)
      	MergeSplitStep(iteration_id, print_to_screen, merge_split_output, int(iteration_id%merge_split_thin_parameter == 0), 1);
      break;
    }
   if(v < logpost_print_prob) 
     {
       if(graph_posterior_output.is_open())
	 PrintGraphLogLikelihood(iteration_id, graph_posterior_output);
       if(active_size_output.is_open())
	 PrintActivityCounts(iteration_id, active_size_output);
       if(infected_size_output.is_open())
	 PrintInfectionCounts(iteration_id, infected_size_output);
       if(parameter_output.is_open() )
	 {
	   for(machine_iter=machines.begin(); machine_iter != machines.end(); machine_iter++)
	     {
	       if(iteration_id >=0)
		 {
		   parameter_output << iteration_id << " ";
		 }
	       parameter_output << machine_iter->second->uid << " ";
	       machine_iter->second->parameters->PrintPars(parameter_output, &(machineMCMC->estimate_par));
	     }
	 }
     }
   
}


void NetworkGraph::MergeSplitStep(int iteration_id, int print_to_screen, ostream& output, int write_output, int write_all_acceptance)
{
  vector<Network*> networks;
  vector<unsigned int> network_ids;
  set<int> network_id_set;
  set<int>::iterator network_id_iterator;
  set<machine_id_t> deletion_set;
  set<machine_id_t>::iterator deletion_set_iterator;
  vector<Machine*> current_machines;
  vector<Machine*> proposed_machines;
  double val;
  int accept = 1;
  int i,j;
  int NUM_TOT;
  Machine* selected_one;
  Machine* selected_two;
 
  distributions->Simulate(3, 1, &val, 0.0, 1.0);

  //get the first machine randomly from all machines

  
  selected_one = RandomlySelectMachine(0);
  current_machines.push_back(selected_one);
  if(val < machineMergeSplit->merge_split_prob[0])
    {
      //cout << "Merge Step" << endl;
      //cout.flush();
      selected_two = RandomlySelectMachine(1, selected_one);
      if(selected_one == selected_two)
	{
	  accept = 0;
	}
      else
	{
	  current_machines.push_back(selected_two);
	}
    }
  
  if(accept)
    {
      //Now make some proposed machines to hold the data
      for(i = 0; i < (3 - current_machines.size()); i++)
	{
	  //initialize it to the current machine (all of this will change, but might as well use the copy constructor
	  proposed_machines.push_back(new Machine(current_machines[0], iteration_id, 1));
	  //To avoid any difficulties with birth/death/immune vs. min/max 1st/last nonoff state
	  //when we are at a merge step, set the Birth/Death/Immune beforehand to 0T0
	  proposed_machines[i]->parameters->SetBirthDeathImmune(0,current_machines[0]->parameters->T,0);
	}
      
      //Get all associated networks for this set of machines
      for(i = 0; i < current_machines.size(); i++)
	{
	  current_machines[i]->data->ListNetworks(&network_ids);
	  for(j = 0; j < network_ids.size(); j++)
	    network_id_set.insert(network_id_to_index[network_ids[i]]);
	}
      for(network_id_iterator = network_id_set.begin(); network_id_iterator != network_id_set.end(); network_id_iterator++)
	{
	  networks.push_back(network_touchpoints[(*network_id_iterator)]);
	}
      
      
      accept = machineMergeSplit->MCMCStep(networks, current_machines, proposed_machines, global_parameters, print_to_screen, iteration_id, merge_split_output, write_output, write_all_acceptance);
      
      //if an acceptance has been made
      if(accept)
	{
	  if(print_to_screen)
	    {
	      cout << iteration_id <<  " Updating global population count and network boundary count" << endl;
	      cout.flush();
	    }
	  if(proposed_machines.size()==2)
	    {
	      global_parameters->num_machines++;
	      global_parameters->boundary_info[networks[0]->parameters->boundary_id][0] += 1;
	    }
	  else
	    {
	      global_parameters->num_machines--;
	      global_parameters->boundary_info[networks[0]->parameters->boundary_id][0] -= 1;
	    }
	  for(i =0; i< proposed_machines.size(); i++)
	    {
	      //If machine is accepted, first set its parameter acceptance rates to 1 try, 1 acceptance
	      if(print_to_screen)
		{
		  cout << iteration_id << " Adding machine " << proposed_machines[i]->uid << " to the graph" << endl;
		  cout.flush();
		}
	      
	      
	      NUM_TOT = machineMergeSplit->NumPars();
	      proposed_machines[i]->par_acc_allocated = 1;
	      proposed_machines[i]->par_acc = new int[NUM_TOT];
	      proposed_machines[i]->par_total = new int[NUM_TOT];
	      proposed_machines[i]->num_par = NUM_TOT;
	      for(j = 0; j < NUM_TOT; j++)
		{
		  proposed_machines[i]->par_acc[j] = 1;
		  proposed_machines[i]->par_total[j] = 1;
		}
	      if(print_to_screen)
		{
		  cout << iteration_id << " Created acceptance information for machine "<< proposed_machines[i]->uid << endl;
		  cout.flush();
		}
	      //Now add the machine into the graph
	      InsertMachine(proposed_machines[i]);
	      if(print_to_screen)
		{
		  cout << iteration_id << " Inserted machine " << proposed_machines[i]->uid << " into the graph and updated global parameters"<< endl;
		}
	    }
	  proposed_machines.clear();
	  
	  for(i = 0; i < current_machines.size(); i++)
	    {
	      if(print_to_screen)
		{
		  cout << iteration_id << " Removing machine " << current_machines[i]->uid << " from the graph" << endl;
		  cout.flush();
		}
	      //Write the old machine's acceptance information using RunMachineMCMC with write_acc=1 and 0 iterations
	      RunMachineMCMC(current_machines[i], 0, 1, 1, iteration_id);
	      if(print_to_screen)
		{
		  cout << iteration_id << " Wrote acceptance information for machine "<< current_machines[i]->uid << " parameter history"<< endl;
		  cout.flush();
		}
	      //push that machine's uid onto the list of those to delete;
	      deletion_set.insert(current_machines[i]->uid);
	    }
	  current_machines.clear();
	  for(deletion_set_iterator = deletion_set.begin(); deletion_set_iterator != deletion_set.end(); deletion_set_iterator++)
	    {
	      DeleteMachine( (*deletion_set_iterator));
	      if(print_to_screen)
		{
		  cout << iteration_id << " Deleted machine " << (*deletion_set_iterator) << " from the graph and updated global parameters" << endl;
		  cout.flush();
		}
	    }
	}
      else
	{	  
	  for(i = 0; i < proposed_machines.size(); i++)
	    {
	      delete proposed_machines[i];
	      if(print_to_screen)
		{
		  cout << iteration_id << " Deleted proposed machine " << proposed_machines[i]->uid << endl;
		  cout.flush();
		}
	    }
	  proposed_machines.clear();
	  current_machines.clear();
	}
    }
  else
    {
      if(print_to_screen)
	cout << iteration_id << " Merge step with only one valid machine: reject" << endl;
      current_machines.clear();      
    }
} 

void NetworkGraph::InitializeNoninformativeMCMC( double mix_prob,int count_rearrange_hour_repeat, int count_rearrange_machine_repeat, int count_rearrange_thin, int count_rearrange_switch_type, int merge_rearrange_repeat, int merge_rearrange_thin, const char* merge_parameter_file, const char* boundary_file) 
{
  int is_noninformative=1;
  if(estimation_type == "Informative")
    is_noninformative = 0; 
  countrearrangeMCMC = new CountRearranger(mix_prob, count_rearrange_machine_repeat, count_rearrange_thin, is_noninformative, count_rearrange_hour_repeat, count_rearrange_switch_type, distributions);
  info_output << "CR mix_prob = " << mix_prob << endl << "CR hour repeat = " << count_rearrange_hour_repeat << endl << "CR machine repeat = " << count_rearrange_machine_repeat << endl << "CR thin = " << count_rearrange_thin << endl << "CR Switch type = " << count_rearrange_switch_type << endl << "merge Repeat = " << merge_rearrange_repeat << endl << "Merge thin = " << merge_rearrange_thin << endl << "Merge parameter file = " << merge_parameter_file << endl << "boundary file = " << boundary_file << endl;
  machineMergeSplit = new MachineMergeSplit(merge_parameter_file, machineMCMC, countrearrangeMCMC, distributions, merge_rearrange_repeat, info_output);
  merge_split_thin_parameter = merge_rearrange_thin;
  MCMC_noninformative_initialized = 1;		

  InitializeBoundaryCounts(info_output, print_to_screen);

}

void NetworkGraph::SetAllBoundaries(int id)
{
  int i;
  if( MCMC_noninformative_initialized)
    cout << "Set boundaries before calling InitializeMCMC" << endl;
  else
    {
      for(i = 0; i < network_touchpoints.size(); i++)
	network_touchpoints[i]->parameters->boundary_id = id;
    }
}

void NetworkGraph::InitializeBoundaryCounts(ostream& str, int print_to_screen)
{
  set<machine_id_t> machines_visited;
  set<machine_id_t>::iterator mach_iter;
  map<int, int*>::iterator b_iter;
  int i,j;
  int b_id;
  int num_machines;
  int num_ips;
  for(i = 0; i < network_touchpoints.size(); i++)
    {
      num_machines = 0;
      for(mach_iter = network_touchpoints[i]->data->machines.begin(); mach_iter != network_touchpoints[i]->data->machines.end(); mach_iter++)
	if(machines_visited.find((*mach_iter)) == machines_visited.end())
	  {
	    machines_visited.insert( (*mach_iter));
	    num_machines++;
	  }
      b_id = network_touchpoints[i]->parameters->boundary_id;
      num_ips = network_touchpoints[i]->data->total_objects;
      if(global_parameters->boundary_info.find(b_id) == global_parameters->boundary_info.end())
	{
	  global_parameters->boundary_info[b_id] = new int[2];
	  global_parameters->boundary_info[b_id][0] = num_machines;
	  global_parameters->boundary_info[b_id][1] = num_ips;
	}
      else
	{
	  global_parameters->boundary_info[b_id][0] += num_machines;
	  global_parameters->boundary_info[b_id][1] += num_ips;
	}
    }
  str << endl << "boundary information (ID, Machines, IPs)"<< endl;
  if(print_to_screen)
     cout << endl << "boundary information (ID, Machines, IPs)"<< endl;
  for(b_iter = global_parameters->boundary_info.begin(); b_iter != global_parameters->boundary_info.end(); b_iter++)
    {
      str << b_iter->first << " " << b_iter->second[0] << " " << b_iter->second[1] << endl;
      if(print_to_screen)
	 cout << b_iter->first << " " << b_iter->second[0] << " " << b_iter->second[1] << endl;
    }
}


void NetworkGraph::Initialize1to1MCMC(const char* directory, const char* repeatfile, const char* hyperparfile, int append_existing, int write_pars, int write_state, int refine_initial_start_values)
{
  int i;
  SetOutputDirectory(directory, append_existing, refine_initial_start_values);
  machineMCMC = new SingleModelUpdater(hyperparfile, repeatfile, distributions, write_pars, write_state);
  globalMCMC.push_back(new ImmuneRateUpdater(2, global_parameters->immune_hyper, distributions));
  globalMCMC.push_back(new SurvivalRateUpdater(global_parameters->T, distributions));
  global_acc = new int[globalMCMC.size()];
  global_total = new int[globalMCMC.size()];
  if(!append)
    {
      Write1to1Info();
      //write titles for estimated parameters
      acc_output << "Machine ";
      parameter_output << "Machine ";
      for(i = 0; i < machineMCMC->parameter_updaters.size(); i++)
	{
	  if(machineMCMC->estimate_par[i])
	    {
	      acc_output << machineMCMC->parameter_updaters[i]->par_name << " ";
	      if(i < machineMCMC->parameter_updaters.size() - 2)
		parameter_output << machineMCMC->parameter_updaters[i]->par_name << " " ;
	     
	    }
	}
      if(machineMCMC->estimate_par[machineMCMC->parameter_updaters.size()-1] == 1)
	{
	  parameter_output << "birth death immunity";
	}
      acc_output << endl;
      parameter_output << endl;
      for(i = 0; i < globalMCMC.size(); i++)
	global_output << globalMCMC[i]->par_name << " ";
      global_output << endl;
    }
  if(refine_initial_start_values)
    {
      //Run through all the machines and do a Data Driven EM starting update, then close the file
      parameter_EM_output << "Machine Nonzero InitialSpike InitialDecay Par1Lambda Par1BIC Par2Lambda1 Par2Lambda2 Par2Pi1 Par2Pi2 Par2BIC Converged" << endl;
      for(machine_iter = machines.begin(); machine_iter != machines.end(); machine_iter++)
	{
	  parameter_EM_output << machine_iter->second->uid << " ";
	  machine_iter->second->RefineStartValuesUsingEM(parameter_EM_output);
	}
      parameter_EM_output.close();
    }
  MCMC_1to1_initialized = 1;
}

void NetworkGraph::SetOutputDirectory(const char* directory, int append_existing, int refineEM)
{
  int ind;
  stringstream str;
  string str_ind;
  int found;

  int cmdtest;
  int cmdtestall;
  string command;
  string test_dir = "test -d ";
  string test_file = "test -f ";
  string make_dir = "mkdir ";
  string directory_par(directory);
  string directory_loc(directory);
  directory_par = output_root + directory_par;

  string acc_file = "acceptance.txt";
  string par_file = "parameters.txt";
  string state_file = "states.txt";
  string info_file = "info.txt";
  string global_file = "global.txt";
  string em_file = "EM_starting_values.txt";
  string cr_file = "countsrearranged.txt";
  string ct_file = "countstats.txt";
  string lp_file = "graph_log_posterior.txt";
  string as_file = "active_size_counts.txt";
  string in_file = "infection_size_counts.txt";
  string mg_file = "merge_split.txt";
  string pad0to9 = "00";
  string pad10to99 = "0";
  string sep = "-";
  string slash = "/";

  time_t T;
  struct tm* timeinfo;
  char fmt[11];
  char fmt2[21];

  time(&T);
  timeinfo = localtime(&T);
  strftime(fmt, 11, "%Y-%m-%d", timeinfo);
  string dateval(fmt, 10);
  strftime(fmt2, 20, "%Y-%m-%d:%h:%m:%s", timeinfo);
  string datetimeval(fmt2,19);
  directory_loc = directory_loc + dateval;

  append = append_existing;
  if(append == 1)
    {
      cmdtestall = 1;
      //find the exact directory needed, test to be sure all files exist
      //then open output stream iterators with ios::append
      command = test_dir + directory_par;
      cout << command.c_str() << endl;
      cmdtest = system(command.c_str());
      cmdtestall = (cmdtestall && (cmdtest == 0 ));
      if(!cmdtestall)
	cout << "Directory " << directory_par << " not found for appending output" << endl;

      command = test_file + directory_par + slash + par_file;
      cout << command.c_str() << endl;
      cmdtest = system(command.c_str());
      cmdtestall = (cmdtestall && (cmdtest == 0));
      if(!cmdtestall)
	cout << "File " << directory_par << "/" << par_file << "not found for appending output " <<endl;
      command = test_file + directory_par + slash + acc_file;
      cout << command.c_str() << endl;
      cmdtest = system(command.c_str());
      cmdtestall = (cmdtestall && (cmdtest == 0));
      if(!cmdtestall)
	cout << "File " << directory_par << "/" << acc_file << "not found for appending output " <<endl;
      
      command = test_file + directory_par + slash + state_file;
      cout << command.c_str() << endl;
      cmdtest = system(command.c_str());
      cmdtestall = (cmdtestall && (cmdtest == 0));
      if(!cmdtestall)
	cout << "File " << directory_par << "/" << state_file << "not found for appending output " <<endl;
   
      command = test_file + directory_par + slash + info_file;
      cout << command.c_str() << endl;
      cmdtest = system(command.c_str());
      cmdtestall = (cmdtestall && (cmdtest == 0));
      if(!cmdtestall)
	cout << "File " << directory_par << "/" << info_file << "not found for appending output " <<endl;

      command = test_file + directory_par + slash + global_file;
      cout << command.c_str() << endl;
      cmdtest = system(command.c_str());
      cmdtestall = (cmdtestall && (cmdtest == 0));
      if(!cmdtestall)
	cout << "File " << directory_par << "/" << global_file << "not found for appending output " <<endl;

      command = test_file + directory_par + slash + lp_file;
      cout << command.c_str() << endl;
      cmdtest = system(command.c_str());
      cmdtestall = (cmdtestall && (cmdtest == 0));
      if(!cmdtestall)
	cout << "File " << directory_par << "/" << lp_file << "not found for appending output " <<endl;
      
      if(estimation_type != "1to1")
	{
	  command = test_file + directory_par + slash + cr_file;
	  cout << command.c_str() << endl;
	  cmdtest = system(command.c_str());
	  cmdtestall = (cmdtestall && (cmdtest == 0));
	  if(!cmdtestall)
	    cout << "File " << directory_par << "/" << cr_file << "not found for appending output " <<endl;

	  command = test_file + directory_par + slash + ct_file;
	  cout << command.c_str() << endl;
	  cmdtest = system(command.c_str());
	  cmdtestall = (cmdtestall && (cmdtest == 0));
	  if(!cmdtestall)
	    cout << "File " << directory_par << "/" << ct_file << "not found for appending output " <<endl;
	  	  
	  command = test_file + directory_par + slash + as_file;
	  cout << command.c_str() << endl;
	  cmdtest = system(command.c_str());
	  cmdtestall = (cmdtestall && (cmdtest == 0));
	  if(!cmdtestall)
	    cout << "File " << directory_par << "/" << as_file << "not found for appending output " <<endl;
	  
	  command = test_file + directory_par + slash + in_file;
	  cout << command.c_str() << endl;
	  cmdtest = system(command.c_str());
	  cmdtestall = (cmdtestall && (cmdtest == 0));
	  if(!cmdtestall)
	    cout << "File " << directory_par << "/" << in_file << "not found for appending output " <<endl;
	  
	  command = test_file + directory_par + slash + mg_file;
	  cout << command.c_str() << endl;
	  cmdtest = system(command.c_str());
	  cmdtestall = (cmdtestall && (cmdtest == 0));
	  if(!cmdtestall)
	    cout << "File " << directory_par << "/" << mg_file << "not found for appending output " <<endl;
	  
	  
	}
      if(refineEM)
	{
	  command = test_file + directory_par + slash + em_file;
	  cout << command.c_str() << endl;
	  cmdtest = system(command.c_str());
	  cmdtestall = (cmdtestall && (cmdtest == 0));
	  if(!cmdtestall)
	    cout << "File " << directory_par << "/" << em_file << "not found for appending output " <<endl;
	}

      if(cmdtestall)
	{
	  cout << "Found files:" << endl;
	  cout << (directory_par +slash + par_file).c_str() << endl 
	       << (directory_par +slash + state_file).c_str() << endl
	       << (directory_par +slash + acc_file).c_str() << endl
	       << (directory_par +slash + info_file).c_str() << endl
	       << (directory_par +slash + global_file).c_str() << endl
	       << (directory_par +slash + lp_file).c_str() << endl;
	  if(estimation_type != "1to1")
	    {
	      cout << (directory_par +slash + cr_file).c_str() << endl; 
	      cout << (directory_par +slash + ct_file).c_str() << endl;
	      cout << (directory_par +slash + as_file).c_str() << endl;
	      cout << (directory_par +slash + in_file).c_str() << endl;
	      cout << (directory_par +slash + mg_file).c_str() << endl;
	    }
	  parameter_output.open( (directory_par +slash + par_file).c_str(), ios::app);
	  state_output.open( (directory_par + slash + state_file).c_str(), ios::app);
	  acc_output.open( (directory_par + slash + acc_file).c_str(), ios::app);
	  info_output.open( (directory_par + slash + info_file).c_str(), ios::app);
	  global_output.open( (directory_par + slash + global_file).c_str(), ios::app);
	  graph_posterior_output.open((directory_par + slash + lp_file).c_str(), ios::app);
	  if(refineEM)
	    parameter_EM_output.open( (directory_par + slash + acc_file).c_str(), ios::app);
	  if(estimation_type != "1to1")
	    {
	      count_rearrange_output.open((directory_par + slash + cr_file).c_str(), ios::app);
	      count_stats_output.open((directory_par + slash + ct_file).c_str(), ios::app);

	      active_size_output.open((directory_par + slash + as_file).c_str(), ios::app);
	      infected_size_output.open((directory_par + slash + in_file).c_str(), ios::app);
	      merge_split_output.open((directory_par + slash + mg_file).c_str(), ios::app);
	    }
	  info_output << "Appended on " << datetimeval << endl;
	  Write1to1Info();
	}
    }
  else
    {
      //create a new directory with the given name, also the date and a unique ID,
      //Open output stream iterators with ios::out 
      ind = 0;
      found = 1;
      while(found)
	{
	  str << ind;
	  str_ind = str.str();
	  if(ind < 10)
	    str_ind = pad0to9 + str_ind;
	  if( (ind >= 10) && (ind < 100))
	    str_ind = pad10to99 + str_ind;
	  command = test_dir + directory_par + dateval + sep + str_ind;
	  cout << command << endl;
	  found = !(system(command.c_str()));
	  if(found)
	    {
	      ind++;
	      str.str("");
	      str_ind.clear();
	    }
	}
      command = make_dir + directory_par + dateval + sep + str_ind;
      cout << command << endl;
      cout << command.c_str() << endl;
      system(command.c_str());
      parameter_output.open( (directory_par + dateval + sep + str_ind + slash + par_file).c_str(), ios::out);
      state_output.open( (directory_par + dateval + sep + str_ind + slash + state_file).c_str(), ios::out);
      acc_output.open( (directory_par + dateval + sep + str_ind + slash + acc_file).c_str(), ios::out);
      info_output.open( (directory_par  + dateval + sep + str_ind + slash + info_file).c_str(), ios::out);
      global_output.open( (directory_par + dateval + sep + str_ind + slash + global_file).c_str(), ios::out);
      graph_posterior_output.open( (directory_par + dateval + sep + str_ind + slash + lp_file).c_str(), ios::out);
      if(refineEM)
	parameter_EM_output.open( (directory_par + dateval + sep + str_ind + slash + em_file).c_str(), ios::out);
      if(estimation_type != "1to1")
	{
	  count_rearrange_output.open( (directory_par + dateval + sep + str_ind + slash + cr_file).c_str(), ios::out);
	  count_rearrange_output << "Iteration Network Time Total #Machines [ID Mean NewCount]" << endl;
	  count_stats_output.open( (directory_par + dateval + sep + str_ind + slash + ct_file).c_str(), ios::out);	  
	  active_size_output.open( (directory_par + dateval + sep + str_ind + slash + as_file).c_str(), ios::out);
	  infected_size_output.open( (directory_par + dateval + sep + str_ind + slash + in_file).c_str(), ios::out);
	  merge_split_output.open( (directory_par + dateval + sep + str_ind + slash + mg_file).c_str(), ios::out);
	}
      info_output << "Created on " << datetimeval << endl;
    
    }
  
}

void NetworkGraph::InsertNetwork(int N, int rid,  ifstream& spikeinfo, ifstream &flinfo, ifstream& ipinfo)
{
  Network* tmp_network = new Network(N, rid, spikeinfo, flinfo, ipinfo, initial_parameters->network_par_init);
  //tmp_network->Print(100);
  /*if(network_touchpoints.find(tmp_network->data->block_id) == network_touchpoints.end())*/
  network_touchpoints.push_back(tmp_network);
  /*else
    cout << "Error: Attempting to insert Network " << tmp_network->data->block_id << " (row " << rid << ") when it already exists in the graph" << endl;*/
  return;// tmp_network;
}

Machine* NetworkGraph::SpawnMachine(Network* network)
{
  /*first create a new machine*/
  machine_id_t uid;
  double immune_value;
  distributions->Simulate(3,1,&(immune_value),0.0,1.0);
  initial_parameters->machine_par_init->SetBirthDeathImmune(0, (short int)network->data->T, (short int)(immune_value < global_parameters->immune_rate));
  Machine* tmp_machine = new Machine(network, initial_parameters->machine_par_init, initial_parameters->hyper_par_init);
  uid = tmp_machine->uid;
  /*Now insert that machine into the machine map*/
  if(machines.find(uid) == machines.end())  /*if uid does not already exist in the map (which it shouldn't)*/
    {
      global_parameters->UpdateGlobalImmunityAndSurvivalCounts(tmp_machine->parameters, 1);
      machines[uid] = tmp_machine;
      /*now update the network's machine list to include the machine ID*/
      network->data->machines.insert(network->data->machines.end(), uid);
      /* the machine adds weight 1 to the network */
      network_weights[network_id_to_index.find(network->data->block_id)->second] += 1.0;
      preproc_up_to_date=0;
    }
  else
    cout << "Error: Machine with id " << uid << " created but already exists in the graph" << endl;
  return tmp_machine;
  
}

void NetworkGraph::InsertMachine(Machine* new_machine)
{
  map<unsigned int,int>::iterator index_iter;
  int i;
  vector<unsigned int> network_ids;
  double num_networks;
  /*Add the new machine to the set of machines*/
  if(machines.find(new_machine->uid) == machines.end())  /*if uid does not already exist in the map */
    {
      global_parameters->UpdateGlobalImmunityAndSurvivalCounts(new_machine->parameters, 1);
      new_machine->data->ListNetworks(&network_ids);
      num_networks = (double) network_ids.size();
      machines[new_machine->uid] = new_machine;
      
      /*Add the machine's uid to all its corresponding networks */
      /*Adjust the network sampling weights to account for the new machine */
      for(i =0; i < network_ids.size(); i++)
	{
	  index_iter = network_id_to_index.find(network_ids[i]);
	  if(index_iter != network_id_to_index.end())
	    {
	      network_touchpoints[index_iter->second]->data->machines.insert(network_touchpoints[index_iter->second]->data->machines.end(), new_machine->uid);
	      network_weights[index_iter->second] += 1.0/num_networks;
	    }
	  else
	    cout << "Network " << network_ids[i] << " accessed by machine " << new_machine->uid << " but does not exist in network graph"<< endl;
	}
      
      /*Tell the rng that the optimized discrete weights per machine are no longer up to date */
      preproc_up_to_date = 0;
    }
  else
    {
      cout << "Warning in InsertMachine: Machine " << new_machine->uid << " requesting insertion but already exists in the Network Graph, it will not be inserted again" << endl;
    }
}

void NetworkGraph::DeleteMachine(machine_id_t ex_machine)
{
  vector<unsigned int> network_ids;
  int i;
  map<unsigned int,int>::iterator index_iter;
  int cur_netw_index;
  set<machine_id_t>::iterator uid_iter;
  double num_networks;
  /*Find the machine in question*/
  machine_iter = machines.find(ex_machine);
  if(machine_iter != machines.end())
    {									
      machine_iter->second->data->ListNetworks(&(network_ids));
      global_parameters->UpdateGlobalImmunityAndSurvivalCounts(machine_iter->second->parameters, -1);
      num_networks = (double)network_ids.size();
      for(i =0; i < network_ids.size(); i++)
	{	  
	  index_iter = network_id_to_index.find(network_ids[i]);
	  cur_netw_index = index_iter->second;
	  uid_iter = network_touchpoints[cur_netw_index]->data->machines.find(ex_machine);
	  if(uid_iter != network_touchpoints[cur_netw_index]->data->machines.end())
	    {
	      /*delete the machine's uid from each corresponding network sets*/
	      network_touchpoints[cur_netw_index]->data->machines.erase(uid_iter);
	      /*Adjust the network sampling weights to account for the deletion of the machine*/
	      network_weights[cur_netw_index] -= 1.0/num_networks;
	      if(network_weights[cur_netw_index] < 0)
		{
		  cout << "Warning in DeleteMachine: network weight = " << network_weights[cur_netw_index] << " : setting to 0.0" << endl;
		  network_weights[cur_netw_index] = 0;
		}
	    }
	}
      /*Tell the rng that the discrete weights per machine are no longer up to date */
      preproc_up_to_date = 0;
      /*Remove the machine from the set of machines*/
      machines.erase(machine_iter);
    }
  else
    {
      cout << "Error in DeleteMachine: Machine " << ex_machine << " requested for deletion but does not exist " << endl;
    }
}

void NetworkGraph::PreProcessWeights()
{
  if(!preproc_up_to_date)
    {
      if(network_weight_preproc_allocated)
	{
	  gsl_ran_discrete_free(network_weight_preproc);
	}
      network_weight_preproc = gsl_ran_discrete_preproc(network_touchpoints.size(), network_weights);
      network_weight_preproc_allocated = 1;
      preproc_up_to_date = 1;
    }  
}


Network* NetworkGraph::RandomlySelectNetwork()
{
  /*This one is easy; since we have a vector of networks, just select a uniform value from 1 to network size and return the pointer*/
  double val;
  int ind;
  distributions->Simulate(3, 1, &val, 0.0, (double)network_touchpoints.size());
  ind = (int)val;
  return network_touchpoints[ind];
}

int NetworkGraph::GetMachineBoundaryID(Machine* machine)
{
  vector<unsigned int> network_ids;
  machine->data->ListNetworks(&network_ids);
  return(network_touchpoints[network_id_to_index[network_ids[0]]]->parameters->boundary_id);
}

Machine* NetworkGraph::RandomlySelectMachine(int pair_provided, Machine* pair)
{
  set<machine_id_t>::iterator machine_id_iter;
  unsigned int ind;
  vector<unsigned int> all_inds;
  int ind2;
  machine_id_t uid;
  Network* selected_network;
  Machine* selected_machine;
  double sum_total = 0.0;
  double sum_aggregate;
  machine_id_t pair_uid = 0;

  int pair_boundary = -1;
  vector<Machine*> candidates;
  vector<double> candidate_prob;
  double selected;
  int i;

  /*first select a network according to its weight*/
  if(pair_provided == 1)
    {
      pair_uid = pair->uid;
      pair_boundary = GetMachineBoundaryID(pair);
      for(ind = 0; ind < network_touchpoints.size(); ind++)
	{
	  if(network_touchpoints[ind]->parameters->boundary_id == pair_boundary)
	    all_inds.push_back(ind);
	}
    }
  else
    {
      PreProcessWeights();
      distributions->Simulate(9, 1, &ind, network_weight_preproc);
      all_inds.push_back(ind);
    }

  /*Now get the network associated with that guy and simulate a machine according to 1/networks pointed to for each machine*/

  for(ind = 0; ind < all_inds.size(); ind++)
    {
      selected_network = network_touchpoints[all_inds[ind]];
      
      if((selected_network->data->machines.size() == 1) && (all_inds.size() == 1))
	{
	  uid = *(selected_network->data->machines.begin());
	  machine_iter = machines.find(uid);
	  if(machine_iter != machines.end())
	    {
	      selected_machine = machine_iter->second;
	    }
	  else
	    {
	      cout << "Error in Select Machine: Network " << selected_network->data->block_id << " (index " << ind << ") points to machine id " << uid << " but that machine ID is not found in the Network Graph" << endl;
	      selected_machine = NULL;
	    }
	}
      else
	{
	  //either the selected network has more than one machine, or we are looking at a pair selection already with more than one network candidate
	  sum_total = 0.0;
	  for(machine_id_iter = selected_network->data->machines.begin(); machine_id_iter != selected_network->data->machines.end();  machine_id_iter++)
	    {
	      uid = *machine_id_iter;
	      machine_iter = machines.find(uid);
	      if(machine_iter != machines.end())
		{
		  if(  (!pair_provided) || ((pair_provided) &&(machine_iter->second->uid != pair_uid)))
		    {
		      candidates.push_back(machine_iter->second);
		      if(pair_boundary == -1)
			candidate_prob.push_back(1.0/(double)(machine_iter->second->data->NumNetworks()));
		      else
			candidate_prob.push_back(1.0);
		      sum_total += candidate_prob[candidate_prob.size()-1];
		    }
		}
	      else
		{
		  cout << "Error: Network " << selected_network->data->block_id << " (index " << ind <<") references machine " << uid << " which does not exist in the graph" << endl;
		}	  
	    }
	  
	  
	  if(sum_total > 0)
	    {
	      distributions->Simulate(3, 1, &selected, 0.0, sum_total);
	      //cout << "total probability weight is " << sum_total << " : selection is " << selected << endl;
	      
	      ind2 = -1;
	      sum_aggregate = 0.0;
	      while((selected > sum_aggregate) && (ind2 < (int)(candidates.size())))
		{
		  ind2++;
		  sum_aggregate += candidate_prob[ind2];
		  //cout << "index= " << ind2 << "; sum_weight = " << sum_aggregate << endl;
		  //cout.flush();
		}
	      //cout << "Index chosen is " << ind2 << endl;
	      selected_machine = candidates[ind2];
	    }
	  else
	    {
	      selected_machine=pair;
	    }
	  
	}
    }
  candidates.clear();
  candidate_prob.clear();
  all_inds.clear();
  return selected_machine;
}
 
void NetworkGraph::SetAllBirthDeathImmune(short int birth, short int death, short int immune)
{
  for(machine_iter = machines.begin(); machine_iter != machines.end(); machine_iter++)
    {
      global_parameters->UpdateGlobalImmunityAndSurvivalCounts(machine_iter->second->parameters, -1);
      machine_iter->second->parameters->SetBirthDeathImmune(birth, death, immune);
      global_parameters->UpdateGlobalImmunityAndSurvivalCounts(machine_iter->second->parameters, 1);
    }
  if(info_output.is_open())
    {
      info_output << endl << "Setting all initial birth/death/immune to birth="<< birth <<", death="<<death<<", immune="<< immune<< endl;
    }
}

void NetworkGraph::SetAllDecayRates(double decay_rate)
{
  for(machine_iter = machines.begin(); machine_iter != machines.end(); machine_iter++)
    {
      machine_iter->second->parameters->alpha = decay_rate;
    }
}

void NetworkGraph::MCMCCleanUp()
{
  int i;
  vector<machine_id_t> ex_machines;
  //Undo the Global Parameter Updaters
  if(globalMCMC.size() > 0)
    {
      delete global_acc;
      delete global_total;
    }
  for(i = 0; i < globalMCMC.size(); i++)
    delete globalMCMC[i];
  while(globalMCMC.size() > 0)
    globalMCMC.pop_back();

  for(machine_iter = machines.begin();
      machine_iter != machines.end();
      machine_iter++)
    {
      ex_machines.push_back(machine_iter->first);
    }

  //Undo all of the machines 
  for(i = 0; i < ex_machines.size(); i++)
    {
      if(i == ex_machines.size() -1)
	{
	  machines[i]->ResetIDs();
	}
      DeleteMachine(ex_machines[i]);
    } 
  ex_machines.clear();
  
  //Undo the single model updater
  if(MCMC_1to1_initialized)
    delete machineMCMC;
  if(MCMC_noninformative_initialized)
    {
      delete countrearrangeMCMC;
      delete machineMergeSplit;
    }
  if(parameter_output.is_open())
      parameter_output.close();
  if(state_output.is_open())
    state_output.close();
  if(acc_output.is_open())
    acc_output.close();
  if(info_output.is_open())
    info_output.close();
  if(global_output.is_open())
    global_output.close();
  if(parameter_EM_output.is_open())
    parameter_EM_output.close();
  if(count_rearrange_output.is_open())
    count_rearrange_output.close();
  if(count_stats_output.is_open())
    count_stats_output.close();
  if(graph_posterior_output.is_open())
    graph_posterior_output.close();
  if(active_size_output.is_open())
    active_size_output.close();			
  if(infected_size_output.is_open())
    infected_size_output.close();
  if(merge_split_output.is_open())
    merge_split_output.close();
  MCMC_1to1_initialized = 0;
  MCMC_noninformative_initialized = 0;
  MCMC_informative_initialized = 0;
}

void NetworkGraph::ReinitializeMCMC(int num_machines_per_network, const char* initializationfile, int print)
{
  MCMCCleanUp();
  for(network_iter = network_touchpoints.begin(); network_iter != network_touchpoints.end(); network_iter++)
    {
      SpawnInitialMachines(*network_iter, num_machines_per_network, print);
    }
  InitializeMCMC(initializationfile);
}

NetworkGraph::~NetworkGraph()
{
  int i;
  MCMCCleanUp();
  for(network_iter = network_touchpoints.begin(); network_iter != network_touchpoints.end(); network_iter++)
    delete *network_iter;
  network_touchpoints.clear();
  delete network_weights;
  delete distributions;
  delete initial_parameters;
  delete global_parameters;
  if(network_weight_preproc_allocated)
    {
      gsl_ran_discrete_free(network_weight_preproc);
    }

}
