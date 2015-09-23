#include "../mcmcdrivers.h"

SingleModelUpdater::SingleModelUpdater(const char* parfile, const char* repeatfile, DistributionList* dists, int wp, int ws)
{
  double hpar1[4];
  double tunevar;
  int i,m;
  ifstream hyperpars;
  ifstream repeatstream;
  double lambdapars[2];
  write_pars = wp;
  write_state = ws;
  distributions = dists;
  hyperpars.open(parfile, ios::in);
  repeatstream.open(repeatfile, ios::in);
  repeatstream >> OVERALL_PAR_REPEAT >> STATE_REPEAT >> ALPHA_REP >> Q_REP >> OMEGA_REP >> RHO_REP >> GAMMA_REP >> NU_REP >> OFFLAMBDA_REP >> BIRTHDEATHIMMUNE_REPEAT >> NUM_PAR >> NUM_TOT;
  cout << "Setting up machine specific parameter updaters" << endl;
  cout.flush();
  hyperpars >> lambdapars[0] >> lambdapars[1];
  parameter_updaters.push_back(new OffLambdaUpdater(2, lambdapars, distributions));
  repeats.push_back(OFFLAMBDA_REP);
  
  hyperpars >> hpar1[0] >> hpar1[1] >> hpar1[2] >> hpar1[3];
  parameter_updaters.push_back(new QUpdater(4, hpar1, distributions));
  repeats.push_back(Q_REP);
  
  hyperpars >> hpar1[0] >> hpar1[1] >> tunevar;
  parameter_updaters.push_back(new AlphaUpdater(2, hpar1, distributions, tunevar));
  repeats.push_back(ALPHA_REP);
  
  hyperpars >> hpar1[0] >> hpar1[1] >> hpar1[2] >> hpar1[3] >> tunevar;
  parameter_updaters.push_back(new OmegaScaledUpdater(4, hpar1, distributions, tunevar));
  repeats.push_back(OMEGA_REP);
  
  for( i = 0; i < 3; i++)
    {
      hyperpars >> hpar1[0] >> hpar1[1] >> tunevar;
      parameter_updaters.push_back(new RhoUpdater(i, 2, hpar1, distributions, tunevar));
      repeats.push_back(RHO_REP);
      
      hyperpars >> hpar1[0] >> hpar1[1];
      parameter_updaters.push_back(new GammaUpdater(i, 2, hpar1, distributions));
      repeats.push_back(GAMMA_REP);
      
      hyperpars >> hpar1[0] >> hpar1[1] >> hpar1[2] >> hpar1[3] >> tunevar;
      parameter_updaters.push_back(new NuUpdater(i, 4, hpar1, distributions, tunevar));	
      repeats.push_back(NU_REP);
    }
  
  parameter_updaters.push_back(new StateUpdater(2, lambdapars, distributions));
  repeats.push_back(STATE_REPEAT);
  
  parameter_updaters.push_back(new BirthDeathImmuneUpdater(2, lambdapars, distributions));
  repeats.push_back(BIRTHDEATHIMMUNE_REPEAT);

  cout << "Getting Estimation information" << endl;
  for(i = 0; i < NUM_TOT; i++)
    {
      hyperpars >> m;
      estimate_par.push_back(m);
    }
  hyperpars.close();
  repeatstream.close();
}

void SingleModelUpdater::PrintRepeatValues(ostream& str)
{
  str << "Name: Repeat IsEstimated" << endl;
  str << "OffLambda : " << OFFLAMBDA_REP << " " << estimate_par[0] << endl
      << "Q : " << Q_REP << " " << estimate_par[1] << endl 
      << "Alpha : " << ALPHA_REP << " " << estimate_par[2] << endl 
      << "Omega : " << OMEGA_REP << " " << estimate_par[3] << endl 
      << "Rho : " << RHO_REP << " " << estimate_par[4] << estimate_par[7] << estimate_par[10] << endl 
      << "Gamma : " << GAMMA_REP << " " << estimate_par[5] << estimate_par[8] << estimate_par[11] << endl
      << "Nu : " << NU_REP << " " << estimate_par[6] << estimate_par[9] << estimate_par[12] << endl
      << "States : " << STATE_REPEAT << " " << estimate_par[13] << endl
      << "Birth/Death/Immune: " << BIRTHDEATHIMMUNE_REPEAT << " " << estimate_par[14] << endl;

  str <<"Overall Parameter repeats : " << OVERALL_PAR_REPEAT << endl
      << "Number of parameters : " << NUM_PAR << endl 
      << "Number of updaters (including states) : " << NUM_TOT << endl;
}

void SingleModelUpdater::UpdateMachine(Machine* machine, global_par_obj* global_parameters, int iterations, int thin, ofstream &par_output, ofstream &state_output, ofstream& acc_output, int export_acc, int iteration_id, int print_steps_to_screen)
{
  if(iterations > 100)
    {
      cout << "Updating Machine "<< machine->uid << " (" << iterations << " iterations, thinned by " << thin << ")" << endl;
    }
  int j,k,l,m;
 
  if(!machine->par_acc_allocated)
    {
      machine->par_acc = new int[NUM_TOT];
      machine->par_total = new int[NUM_TOT];
      for(j = 0; j < NUM_TOT; j++)
	{
	  machine->par_acc[j] = 0;
	  machine->par_total[j] = 0;
	  parameter_updaters[j]->SetValidStartValue(machine->parameters, global_parameters);
	}      
      machine->num_par = NUM_TOT;
      machine->par_acc_allocated = 1;
    }
  
  curr_par_obj* current_parameter_values = machine->parameters;
  mach_dat_obj* data_values = machine->data;
  
  //PrintRepeatValues();
  //for(parup_iter = parameter_updaters.begin(); parup_iter != parameter_updaters.end(); parup_iter++)
  // (*parup_iter)->PrintHyperpars();

  //cout << "Summing machine blocks...";
  //cout.flush();
 
  for(j = 0; j < iterations; j++)
    {
      /* parameter runs*/
      for(m = 0; m < OVERALL_PAR_REPEAT; m++)
	{
	  for(k = 0; k < NUM_PAR; k++)
	    {
	      for(l = 0; l < repeats[k]; l++)
		{
		  if(estimate_par[k])
		    {
		      //cout << "updating: ";
		      //parameter_updaters[k]->PrintHyperpars();
		      machine->par_total[k]++;
		      machine->par_acc[k] += (int) parameter_updaters[k]->MCMC_Step(current_parameter_values, data_values, global_parameters, iteration_id, print_steps_to_screen);
		    }
		}
	    }
	  if(write_pars)
	    {
	      if(thin != 0)
		{
		  if( (j%thin) == 0)
		    {
		      if(iteration_id >=0)
			{
			  par_output << iteration_id << " ";
			}
		      par_output << machine->uid << " ";
		      current_parameter_values->PrintPars(par_output, &(estimate_par));
		    }
		}
	    }
	}
      /*state runs and birth/death/immune updates*/
      if(estimate_par[NUM_TOT-2])
	{
	  for(k = 0; k < repeats[NUM_TOT-2]; k++)
	    {
	      machine->par_total[NUM_TOT-2]++;
	      machine->par_acc[NUM_TOT-2] += (int) parameter_updaters[NUM_TOT-2]->MCMC_Step(current_parameter_values, data_values, global_parameters, iteration_id, print_steps_to_screen);
	    }
	  if(estimate_par[NUM_TOT-1])
	    {
	      /*update birth/death/immunity*/
	      for(k = 0; k < repeats[NUM_TOT-1]; k++)
		{
		  machine->par_total[NUM_TOT-1]++;
		  machine->par_acc[NUM_TOT-1] += (int) parameter_updaters[NUM_TOT-1]->MCMC_Step(current_parameter_values, data_values, global_parameters, iteration_id, print_steps_to_screen);
		}
	    }
	  //current_parameter_values->PrintPars(outfile_par);
	  if(write_state)
	    {
	      if(thin != 0)
		{
		  if( (j%thin) == 0)
		    {
		      //cout << "writing states"<< endl;
		      if(iteration_id >= 0)
			{
			  state_output << iteration_id << " ";
			}
		      state_output << machine->uid << " " << machine->parameters->Birthdate() << " " << machine->parameters->Deathdate() << " " << machine->parameters->Is_immune() << " ";
		      current_parameter_values->PrintState(state_output);
		    }
		}
	    }
	}
    }
  if(export_acc)
    {
      acc_output << machine->uid << " " ;
      if(machine->num_par > 0)
	{
	  for(j=0; j< machine->num_par; j++)
	    {
	      if(estimate_par[j])
		{
		  if(machine->par_total[j] > 0)
		    acc_output << (double)(machine->par_acc[j])/(double)(machine->par_total[j]) << " ";
		  else
		    acc_output << 0 << " ";
		}
	    }
	  acc_output << endl << machine->uid << " ";
	  for(j = 0; j < machine->num_par; j++)
	    {
	      if(estimate_par[j])
		{
		  acc_output << machine->par_total[j] << " ";	     	    
		}
	    }
	  acc_output << endl;
	}
      else
	{
	  acc_output << machine->uid << " ";
	  for(j=0; j< NUM_TOT; j++)
	    {
	      if(estimate_par[j])
		acc_output << 0 << " ";
	    }
	  acc_output << endl << machine->uid << " " ;
	  for(j = 0; j < NUM_TOT; j++)
	    {
	      if(estimate_par[j])
		acc_output << 0 << " ";
	    }
	  acc_output << endl;
	}
    }
}

void SingleModelUpdater::PrintAllHyperpars(ostream& str)
{
  for(parup_iter = parameter_updaters.begin(); parup_iter != parameter_updaters.end(); parup_iter++)
    {
      (*parup_iter)->PrintHyperpars(str);
    }
}

double SingleModelUpdater::GraphLogLikelihood(Machine* machine, global_par_obj* global_pars, int suppress_prior)
{
  double loglike;
  loglike = 0.0;
  vector<double> indiv_pois;
  vector<double> indiv_trans;
  int i;
  if(!suppress_prior)
    {
      for(i=0; i < parameter_updaters.size(); i++)
	loglike += parameter_updaters[i]->GraphLogLikelihood(machine->parameters, global_pars);
    }
  loglike += GetPathLogLikelihood(0, machine->parameters->T, machine->data, machine->parameters, global_pars, 0, &indiv_pois, &indiv_trans);
  indiv_pois.clear();
  indiv_trans.clear();
  
  return(loglike);
}

SingleModelUpdater::~SingleModelUpdater()
{
 for(parup_iter = parameter_updaters.begin(); parup_iter != parameter_updaters.end(); parup_iter++)
   {
     delete *parup_iter;
   }
 parameter_updaters.clear();
 //delete distributions;
}

CountRearranger::CountRearranger(double mix, int repeat_iter, int thin_iter, int is_noninfo, int num_hr, int tp, DistributionList* dists)
{
  distributions = dists;
  vectors_allocated = 0;
  mixprob=mix;
  repeat = repeat_iter;
  is_noninformative = is_noninfo;
  thin = thin_iter;
  num_hours = num_hr;
  type = tp;
}

/*void CountRearranger::Print(ostream& str)
{
  
}*/

int CountRearranger::SetUpMachines(Network* network, map<machine_id_t, Machine*>& machines)
{
  int i;
  int multiple_machines = 1;
  machine_id_t cur_machine;
  min_birth = network->data->T;
  max_death = -1;
  int num_machines = network->data->machines.size();
  if(num_machines > 1)
    {
      mapped_machines.clear();
      //First get the pointers to the machines associated with this network
      for(id_iter = network->data->machines.begin(); id_iter != network->data->machines.end(); id_iter++)
	{
	  //pop the machine pointers onto the vector of mapped machines
	  cur_machine = *id_iter;
	  //we could push all back in one step
	  //but it is a good idea to check to make sure that the machine exists 
	  machine_iter = machines.find(cur_machine);
	  if(machine_iter != machines.end())
	    {
	      mapped_machines.push_back(machine_iter->second);
	      if(machine_iter->second->parameters->Birthdate() < min_birth)
		min_birth = machine_iter->second->parameters->Birthdate();
	      if(machine_iter->second->parameters->Deathdate() > max_death)
		max_death = machine_iter->second->parameters->Deathdate();
	      if(max_death >= machine_iter->second->data->T)
		max_death = machine_iter->second->data->T-1;
	    }
	  else
	    {
	      cout << "Error: Network " << network->data->block_id << " (row " << network->data->row_id << ") is linked to machine "<< cur_machine << " which does not exist in the graph" << endl;
	      num_machines--;
	    }
	}
      if(vectors_allocated)
	{
	  delete pk;
	  delete simulated_counts;
	}
      pk = new double[num_machines];
      simulated_counts = new int[num_machines]; 
      vectors_allocated = 1;
      for(i = 0 ; i < num_machines; i++)
	{
	  pk[i] = 0.0;
	  simulated_counts[i] = 0;
	}
    }
  else
    {
      multiple_machines = 0;
    }
  return(multiple_machines);
}

int CountRearranger::SetUpProbabilities(int t)
{
  int i;
  int num_nonzero_rate = 0;
  if(!vectors_allocated)
    {
      cout << "Error in SetUpProbabilities: No machines set up for setting probabilities" << endl;
    }
  else
    {
      for(i = 0; i < mapped_machines.size(); i++)
	{
	  if(is_noninformative)
	    pk[i] = mapped_machines[i]->GetRate(t);
	  else
	    cout << "Haven't done informative yet" << endl;
	  if(pk[i] > 0.0)
	    num_nonzero_rate++;
	}
    }
  
  return(num_nonzero_rate);
}

int CountRearranger::Sample(Network* network, map<machine_id_t, Machine*>& machines, int print_to_screen, int iteration_id, ostream& str, int print_output, int deterministic)
{
  //type = 0: Simple Gibbs Step
  //type = 1: Always do a StateChange step
  //type = 2: Simple Gibbs Step if more than one nonzero rate; StateChange step between existing guy and one other otherwise
  //type = 3: Cointoss?
  int i;
  int j;
  double t;
  int t_int;
  int non_trivial;
  int num_nonzero_rate_t;
  int num_hours_loop;
  int num_acc = 0;
  non_trivial = SetUpMachines(network, machines);
  if(non_trivial)
    {
      if(deterministic)
	{
	  num_hours_loop = network->data->T;
	}
      else
	{
	  num_hours_loop = num_hours;
	}
      for(i = 0; i < num_hours_loop; i++)
	{  
	  if(!deterministic)
	    {
	      distributions->Simulate(3, 1, &t, (double)(min_birth), (double)(max_death)+1.0);
	      t_int = (int)floor(t);
	    }
	  else
	    {
	      t_int = i;
	    }
	  //Check the count at this time, if zero, skip to the next hour for obvious reasons
	  if( network->data->fl[t_int] > 0)
	    {
	      num_nonzero_rate_t = SetUpProbabilities(t_int);
	      if(!deterministic)
		{
		  switch(type){
		  case 0:
		    if(num_nonzero_rate_t > 1)
		      num_acc += GibbsStep(network,t_int); 		  
		    break;
		  case 1:
		    num_acc += StateChangeStep(network, t_int);
		    break;
		  case 2:
		    if(num_nonzero_rate_t > 1)
		      {
			num_acc += GibbsStep(network, t_int);
		      }
		    else
		      {
			num_acc += StateChangeStep(network, t_int);
		      }		    
		    break;
		  case 3:
		    distributions->Simulate(3, 1, &t, 0.0, 1.0);
		    if(t < mixprob)
		      {  
			if(num_nonzero_rate_t > 1)
			  num_acc += GibbsStep(network, t_int);
			else
			  num_acc +=  StateChangeStep(network, t_int);
		      }
		    else
		      {
			num_acc += StateChangeStep(network, t_int);
		      }
		    break;
		  }
		}
	      else
		{
		  //do an assignment even if one probability is equal to 0.  
		    num_acc += GibbsStep(network,t_int);
		}
	    }
	  if(!deterministic)
	    {
	      if( (i % thin) == 0)
		{
		  if(print_output)
		    {
		      if(iteration_id >= 0)
			{
			  str << iteration_id << " ";
			  
			}
		      str << network->data->block_id << " " << t_int << " " << network->data->fl[t_int] << " " << mapped_machines.size() << " ";		   
		      for(j = 0; j < mapped_machines.size(); j++)
			{
			  str << mapped_machines[j]->uid << " " << pk[j] << " ";
			  if(network->data->fl[t_int] > 0)
			    str << simulated_counts[j] << " ";
			  else
			    str << 0 << " ";
			}
		      str << endl;
		    }
		}
	      if(print_to_screen)
		{
		  cout << iteration_id << " ";
		  cout <<  network->data->block_id << " " << t_int << " " << network->data->fl[t_int] << " " << mapped_machines.size() << " ";
		  for(j = 0; j < mapped_machines.size(); j++)
		    {
		      cout << mapped_machines[j]->uid << " " << pk[j] << " ";
		      if(network->data->fl[t_int] > 0)
			cout << simulated_counts[j] << " ";
		      else
			cout << 0 << " ";
		      
		    }
		  cout << endl;
		}
	    }
	}
    }
  else
    {
      /*
      if(!deterministic)
	{
	  if(print_output)
	    {
	      if(iteration_id >= 0)
		{
		  str << iteration_id << " ";	
		}
	      str << network->data->block_id << " -1 -1 1 ";	   
	      for(id_iter = network->data->machines.begin(); id_iter != network->data->machines.end(); id_iter++)
		{
		  str << *(id_iter);
		  str << " 1.0 -1" << endl;
		  
		}
	      str << " 1.0 -1" << endl;
	      
	      
	    }
	  if(print_to_screen)
	    {
	      for(id_iter = network->data->machines.begin(); id_iter != network->data->machines.end(); id_iter++)
		{
		  cout << iteration_id << " ";
		  cout << network->data->block_id << " -1 -1 1 ";
		  cout << *(id_iter);
		  cout << " 1.0 -1" << endl;
		}
	      
	    }
	}
      */
    }
  
  return num_acc;
}

int CountRearranger::GibbsStep(Network* network, int t)
{
  int i;
  double flowcount;
  double K;
  flowcount = (double)network->data->fl[t];
  K = (double)mapped_machines.size();
  
  distributions->Simulate(8, 1, simulated_counts, flowcount, K, pk);
  for(i = 0; i < mapped_machines.size(); i++)
    {
      //Reset the counts on the machines to the ones in the simulated counts
      SetNewCount(network, i, t);
    }
  return 1;
}

void CountRearranger::SetNewCount(Network* network, int machine_index,  int t)
{
  //Set the counts for the machine in mapped_networks[machine_index] at hour t
  mapped_machines[machine_index]->data->SetCountExisting(network->data->block_id, t, simulated_counts[machine_index]);
}


void CountRearranger::SetMix(double mix)
{
  mixprob = mix;
}

int CountRearranger::StateChangeStep(Network* network, int  t)
{
  return(0);

}

CountRearranger::~CountRearranger()
{
  if(vectors_allocated)
    {
      delete pk;
      delete simulated_counts;
    }
  mapped_machines.clear();
}

/* ===================Graph Object Updaters===================== 

int GraphObjectUpdater::MCMC_Step(vector<Network*>& Networks, vector<Machine*>& Machines, global_par_obj* global_pars, int iterations, int print_step=0)
{
  int i;
  double alpha;
  double logratio;
  int is_accepted = 0;
  for(i = 0; i < iterations; i++)
    {
      SetCurrentValue(Networks, Machines, global_pars);
      Propose(Networks, Machines);
      if(!ProposeEqCurrent())
	{
	  logratio = LogRatio(Networks, Machines, global_pars);
	  distributions->Simulate(3, 1, &alpha, 0.0, 1.0);
	  if(log(alpha) < logratio)
	    {
	      is_accepted++;
	      SetCurrentProposed(Networks, Machines, global_pars);
	    }
	}
      if(print_step == 1)
	PrintToScreen();
    }
  return is_accepted;
}
*/
