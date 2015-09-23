#include "../countupdaters.h"


/*======================================Counts========================================*/

CountMergeSplitUpdater::CountMergeSplitUpdater(vector<double>& hpar, DistributionList* dists, CountRearranger* cr):MergeSplitUpdater(hpar, dists)
{
  par_name="MergeSplit_Counts";
  count_rearranger = cr;
  machine_map.clear();
}

void CountMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  int i;

  //the proposal is currently set to Machine 0; all we need to do is add the blocks for Machine 1 to it
  proposed_machines[0]->data->AddBlocks(current_machines[1]->data);
  // cout << "Counts:";
  //for(i = 0; i < proposed_machines[0]->parameters->T; i++)
  //  {
  //    cout << endl << i << " " << current_machines[0]->data->GetTotalCount(i) << " " << current_machines[1]->data->GetTotalCount(i) << " -> " << proposed_machines[0]->data->GetTotalCount(i);
  //  }
  //cout << endl;
}

double CountMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{

  //log(1->2) - log(2->1)
  double loglike =0;
  vector<unsigned int> network_list;
  int i, j;
  s1_mse=0;
  s2_mse=0;
  s_mse =0;
  s1_tot = 0;
  s2_tot = 0;
  s_tot = 0;
  current_machines[0]->data->ListNetworks(&network_list);
  for(i = 0; i < network_list.size(); i++)
    {
      for(j = 0; j < current_machines[0]->parameters->T; j++)
	{
	  loglike += RateLogLikelihood(network_list[i], j, current_machines, proposed_machines, 1);
	  if(PRINT_RATIO_STEP == 1)
	    {
	      if(j < current_machines[0]->parameters->T -1 )
		cout << " partial sum  = " << loglike;
	    }
	}
    }

  return loglike;
}

void CountMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  int i;
  int t;
  Network* cur_network;
  machine_map.clear();
  machine_map[proposed_machines[0]->uid] = proposed_machines[0];
  machine_map[proposed_machines[1]->uid] = proposed_machines[1];
  //this may be a bit overkill for big DHCP networks but that shall be reserved for a later date
  //when I know that there is a chance that networks.size() > 1 (eg, not quite yet)
  for(i = 0; i < networks.size(); i++)
    {
      //Set up a new "fake" network with all of the properties of the selected network 
      cur_network = new Network(networks[i]);
      //Except the only machines it has are the two proposed
      cur_network->data->machines.clear();
      cur_network->data->machines.insert(proposed_machines[0]->uid);
      cur_network->data->machines.insert(proposed_machines[1]->uid);
      //And its counts are equal to the sum of the two proposed networks (eg, the current machine counts thru that network)
      for(t = 0; t < cur_network->data->T; t++)
	{
	  cur_network->data->fl[t] = current_machines[0]->data->GetSubCount(networks[i]->data->block_id, t);
	}
      //Now rearrange these counts according to the multinomial (binomial) probability for each of the time periods, using the count rearranger
      count_rearranger->Sample(cur_network, machine_map, 0,0, cnull, 0, 1);
      delete cur_network;
    }
  //cout << "Counts:";
  //for(i = 0; i < proposed_machines[0]->parameters->T; i++)
  //  {
  //    cout << endl << i << " " << current_machines[0]->data->GetTotalCount(i) << " -> " << proposed_machines[1]->data->GetTotalCount(i) << " " << proposed_machines[0]->data->GetTotalCount(i);
  //  }
  //cout << endl;
}

double CountMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  //log(2->1) - log(1->2)
  double loglike =0;
  vector<unsigned int> network_list;
  int i, j;
  s1_mse=0;
  s2_mse=0;
  s_mse = 0;
  s1_tot = 0;
  s2_tot = 0;
  s_tot = 0;
  current_machines[0]->data->ListNetworks(&network_list);
  for(i = 0; i < network_list.size(); i++)
    {
      for(j = 0; j < current_machines[0]->parameters->T; j++)
	{
	   loglike += RateLogLikelihood(network_list[i], j, current_machines, proposed_machines, -1);
	   if(PRINT_RATIO_STEP == 1)
	     {
	       if(j < current_machines[0]->parameters->T -1 )
		 cout << " partial sum  = " << loglike;
	     }
	}
    }

  return loglike;
}

void CountMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  if(s_tot > 0)
    {
      s1_mse = sqrt(s1_mse/s1_tot);
      s2_mse = sqrt(s2_mse/s2_tot);
      s_mse = sqrt(s_mse/s_tot);
    }
   if(machines.size() == 1)
    {
      str << machines[0]->uid << ":(" << s_mse << " Chi-sq "<< s_tot << " df) " ;
    }
  else
    {
      str << machines[0]->uid << ":(" << s1_mse << " Chi-sq, "<< s1_tot << " df) " ;
      str << machines[1]->uid << ":(" << s2_mse << " Chi-sq, "<< s2_tot << " df) " ;
    }
}


double CountMergeSplitUpdater::RateLogLikelihood(unsigned int network_id, int t, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, double mergesplit)
{
  double rate1, rate2, rate_12;
  double count1, count2, count_12;
  short int s1, s2, s12;
  double p;
  double loglike = 0.0;
  //1 = merge, -1 = split (multiply by logratio)
  if(mergesplit == 1)
    {
      //merge assignments
      rate1 = current_machines[0]->parameters->GetRate(t);
      count1 = (double)(current_machines[0]->data->GetSubCount(network_id, t));
      s1 = current_machines[0]->parameters->State(t);

      rate2 = current_machines[1]->parameters->GetRate(t);
      count2 = (double)(current_machines[1]->data->GetSubCount(network_id, t));
      s2 = current_machines[1]->parameters->State(t);

      rate_12 = proposed_machines[0]->parameters->GetRate(t);
      count_12 = (double)(proposed_machines[0]->data->GetSubCount(network_id, t));
      s12 = proposed_machines[0]->parameters->State(t);
    }
  else
    {
      //split assignments
      rate1 = proposed_machines[0]->parameters->GetRate(t);
      count1 = (double)(proposed_machines[0]->data->GetSubCount(network_id, t));
      s1 = proposed_machines[0]->parameters->State(t);

      rate2 = proposed_machines[1]->parameters->GetRate(t);
      count2 = (double)(proposed_machines[1]->data->GetSubCount(network_id, t));
      s2 = proposed_machines[1]->parameters->State(t);

      rate_12 = current_machines[0]->parameters->GetRate(t);
      count_12 =(double)(current_machines[0]->data->GetSubCount(network_id, t));
      s12 = current_machines[0]->parameters->State(t);
    }
  
  if(rate_12 > 0)
    {
      s_mse += (rate_12-count_12)*(rate_12-count_12)/(rate_12);
      s_tot++;
      if(rate1 > 0)
	{
	  s1_mse += ((rate1-count1)*(rate1-count1))/(rate1);
	  s1_tot++;
	}
      if(rate2 > 0)
	{
	  s2_mse += (rate2-count2)*(rate2-count2)/(rate2);
	  s2_tot++;
	}
      //Factoring out the Poisson likelihood ratio at each step (and not calculating it in the likelihood) 
      //gives a simple form
      loglike = rate1 + rate2 - rate_12 + count_12*( log(rate_12) - log(rate1 + rate2));
      
      

      if(PRINT_RATIO_STEP ==1)
	{
      
	  cout << endl << t <<  ": Counts " << count_12;
	  if(mergesplit == -1)
	    cout << " -> ";
	  else
	    cout << " <- ";
	  cout << count1 << " " << count2;
	  
	  cout <<  " : States " << s12;
	  if(mergesplit == -1)
	    cout << " -> ";
	  else
	    cout << " <- ";
	  cout << s1 << " " << s2;
	  
	  cout << " : Rates " << rate_12 ;
	  if(mergesplit == -1)
	    cout << " -> ";
	  else
	    cout << " <- ";
	  cout << rate1 << " " << rate2;
	  cout << " : loglikelihood = " << (double)(mergesplit)*loglike;
	  if(t == current_machines[0]->parameters->T - 1)
	    cout << endl;
	}
	/*
	p = rate1/(rate1 + rate2);
	loglike = log(gsl_ran_binomial_pdf((unsigned int)(count1), p, (unsigned int)(count_12)));
      */     
    }
  return(((double)(mergesplit))*loglike);
}

/*====================================End Counts======================================*/

/*==============================Overall Machine Merge Split==========================*/

MachineMergeSplit::MachineMergeSplit(const char* mergesplitparfile, SingleModelUpdater* single_model, CountRearranger* count_rearranger, DistributionList* dists, int repeats, ostream& info_out)
{
  ifstream infile;
  int i,j;
  string bdi_pars;
  num_accepted = 0;
  num_tries = 0;
  vector<double> hyperpars;
  double val;
  const int update_order[] = {13,14,1,2,3,0,5,6,4,8,9,7,11,12,10};
  info_out << endl << "Merge-Split parameters:" << endl;

  distributions = dists;
  machineMCMC = single_model;
  infile.open(mergesplitparfile, ios::in);
  num_repeats = repeats;
  info_out << "Number of merge-split repeats per iteration when chosen:" << num_repeats<<endl;
  infile >> val;
  merge_split_prob.push_back(val);
  info_out << "Prob(Merge chosen | Merge-Split chosen):" << val <<endl; 
  merge_split_prob.push_back(1.0 - val);
  
  infile >> val;
  pri_H_par.push_back(val);
  infile >> val;
  pri_H_par.push_back(val);
  //prior parameters for the number of machines (n, p)
  info_out << "Prior(H) parameters (n, p) for negative binomial " << pri_H_par[0] << " " << pri_H_par[1] << endl;


  infile >> bdi_pars;
   //get estimation info
  
  for(i = 0; i < single_model->estimate_par.size(); i++)
    {
      j = update_order[i];
      estimate_par.push_back(single_model->estimate_par[j]);
    }
  
  NUM_TOT = single_model->NumPars();
  
  //State needs to be done first, so that the things that depend on it will make sense
  merge_split_updaters.push_back(new StateMergeSplitUpdater(hyperpars, distributions, single_model->parameter_updaters[13]));
  merge_split_updaters.push_back(new BirthDeathImmuneMergeSplitUpdater(bdi_pars, hyperpars, distributions, single_model->parameter_updaters[14]));

  /*===========Q============*/
  infile >> val;
  hyperpars.push_back(val); // proposal_sample_size for bounded average
  merge_split_updaters.push_back(new QMergeSplitUpdater(hyperpars, distributions, (QUpdater*)single_model->parameter_updaters[1]));
  hyperpars.clear();

  /*==========Alpha==========*/
  infile >> val;
  hyperpars.push_back(val); // proposal_sample_size for bounded average
  merge_split_updaters.push_back(new AlphaMergeSplitUpdater(hyperpars, distributions, single_model->parameter_updaters[2]));
  hyperpars.clear();

  /*========Omega==========*/
  infile >> val;
  hyperpars.push_back(val); //proposal_sample size for bounded average
  merge_split_updaters.push_back(new OmegaScaledMergeSplitUpdater(hyperpars, distributions, single_model->parameter_updaters[3]));
  hyperpars.clear();
 
  merge_split_updaters.push_back(new OffLambdaMergeSplitUpdater(hyperpars, distributions, single_model->parameter_updaters[0]));
  merge_split_updaters.push_back(new GammaMergeSplitUpdater(hyperpars, distributions, (GammaUpdater*)(single_model->parameter_updaters[5])));

  infile >> val;
  hyperpars.push_back(val); //sample size (a + b)
  merge_split_updaters.push_back(new NuMergeSplitUpdater(hyperpars, distributions, (NuUpdater*)(single_model->parameter_updaters[6])));
  hyperpars.clear();


 /*========Rho-Off==========*/
  merge_split_updaters.push_back(new RhoMergeSplitUpdater(hyperpars, distributions, (RhoUpdater*)(single_model->parameter_updaters[4])));


  merge_split_updaters.push_back(new GammaMergeSplitUpdater(hyperpars, distributions,(GammaUpdater*)(single_model->parameter_updaters[8])));

  infile >> val;
  hyperpars.push_back(val); //sample size (a + b)
  merge_split_updaters.push_back(new NuMergeSplitUpdater(hyperpars, distributions, (NuUpdater*)(single_model->parameter_updaters[9])));
  hyperpars.clear();

  /*========Rho-Spike==========*/
   merge_split_updaters.push_back(new RhoMergeSplitUpdater(hyperpars, distributions, (RhoUpdater*)(single_model->parameter_updaters[7])));
   merge_split_updaters.push_back(new GammaMergeSplitUpdater(hyperpars, distributions, (GammaUpdater*)(single_model->parameter_updaters[11])));
  infile >> val;
  hyperpars.push_back(val); //sample size (a + b)
  merge_split_updaters.push_back(new NuMergeSplitUpdater(hyperpars, distributions, (NuUpdater*)(single_model->parameter_updaters[12])));
  hyperpars.clear();
  
 /*========Rho-Decay==========*/
  merge_split_updaters.push_back(new RhoMergeSplitUpdater(hyperpars, distributions, (RhoUpdater*)(single_model->parameter_updaters[10])));
  
 
  merge_split_updaters.push_back(new InfoMergeSplitUpdater(hyperpars, distributions));
  estimate_par.push_back(1);
  info_out << "Info:1" << endl;
 
  //Counts need to be done last for the split, because they depend on the weighted means (they are just a sum for the merge)
  merge_split_updaters.push_back(new CountMergeSplitUpdater(hyperpars, distributions, count_rearranger));
  info_out << "Counts:1"<< endl;
  estimate_par.push_back(1);

  info_out << "estimation information [IsEstimated: Name (Prior Hyperparameters) (MergeSplit parameters)]"<<endl;
  for(i = 0; i < merge_split_updaters.size(); i++)
    {
      info_out << estimate_par[i] << ": ";
      merge_split_updaters[i]->PrintHyperpars(info_out);
    }

  infile.close();
}

MachineMergeSplit::~MachineMergeSplit()
{
  vector<MergeSplitUpdater*>::iterator mergesplit_iterator;
  for(mergesplit_iterator = merge_split_updaters.begin(); mergesplit_iterator != merge_split_updaters.end(); mergesplit_iterator++)
    {
      delete (*mergesplit_iterator);
    }
  merge_split_updaters.clear();
}

int MachineMergeSplit::MCMCStep(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars, int print_to_screen, int iteration_id, ostream& str, int print_output, int record_all_acceptance)
{
  double loglike;
  double temp_inf_check;
  int i = 0;
  int accepted = 1;
  double logalpha;
  double curr_H;
  double num_machines; 
  loglike = 0.0;
  double machine_likelihood;
  double machine_middle_mult;
  double prior_weights;
  if(current_machines.size() == 2)  
    {
      if(current_machines[0]->uid == current_machines[1]->uid)
	{
	  loglike = -1.0*INFINITY;
	  accepted = 0;
	}
    }
  if(accepted)
    {    
      num_machines = (double)(global_pars->num_machines);
      curr_H = (double)(global_pars->boundary_info[networks[0]->parameters->boundary_id][0]);
      
      //taking care of disparate probabilities of choosing merge v. split
      prior_weights = LogMoveRatio(merge_split_prob, 0, 1);
      if(print_to_screen)
	{
	  cout << iteration_id << " "<<  "MergeSplit Prior: ";
	}
      if(current_machines.size() == 1)
	{	
	  //split step
	  prior_weights = -1.0*prior_weights;
	  //taking care of the prior on H which is negative binomial (n, p)   
	  prior_weights += log(curr_H+pri_H_par[0]) + log(1.0 - pri_H_par[1]) -  log(curr_H+1.0);
	  //multiplicity of the allocations in the split step: symmetry in assigning state configurations to machine1, machine2 vs. machine2, machine1
	  prior_weights -= log(2.0);

	  //taking care of sampling two vs. one machine
	  prior_weights +=  log(2.0) + log(num_machines) - log(num_machines + 1) - log(curr_H);
	  
	  if(print_to_screen)
	    cout << "Split Step: ";
	}
      else
	{
	  //taking care of the prior on H which is negative binomial (n, p)   
	  prior_weights += log(curr_H) - log(1.0 - pri_H_par[1]) -  log(curr_H-1.0+pri_H_par[0]);
	  //multiplicity of the allocations in the split step: symmetry in assigning state configurations to machine1, machine2 vs. machine2, machine1
	  prior_weights += log(2.0);

	  //taking care of sampling one vs. two machine
	  prior_weights += log(num_machines) + log(curr_H - 1.0) - log(2.0) - log(num_machines - 1.0);
	  if(print_to_screen)
	    cout << "Merge Step: ";
	}
      if(print_to_screen)
	cout << "current H = " << num_machines << " Pair boundary = " << curr_H << " prior and selection log-likelihood = " << prior_weights << endl;

      loglike = 0.0;
      //Proposing a new set of states and counts and getting back the jumping ratios
      while( (i < merge_split_updaters.size()) && accepted)
	{
	  if(estimate_par[i])
	    {
	      if(print_to_screen)
		{
		  cout << iteration_id << " ";
		}
	      temp_inf_check = merge_split_updaters[i]->MergeSplitStep(networks, current_machines, proposed_machines, global_pars, print_to_screen);
	      if(temp_inf_check == -1.0*INFINITY)
		{
		  loglike = -1.0*INFINITY;
		  //accepted = 0;
		}
	      else
		{
		  loglike += temp_inf_check;	      
		}
	    }
	  i++;
	}
      if(print_to_screen)
	cout <<iteration_id << " Total log-jump ratio = " << loglike << ": log-jump + prior = " << loglike + prior_weights << endl;
  
      loglike = loglike + prior_weights;
      
      machine_likelihood = LogLikelihoodRatio(proposed_machines, current_machines, global_pars, print_to_screen, iteration_id);
      loglike += machine_likelihood;

      if(print_to_screen)
	{
	  cout << iteration_id << " Total log-likelihood ratio for the proposal = " << loglike ; 
	}

      //Now do the likelihoods of the machines.  Note that this method should protect against infinities or ratios thereof
      //but we will do a check anyway
      if(loglike == -1.0*INFINITY)
	{
	  accepted = 0;
	} 
      

      if(accepted)
	{
	 
	  distributions->Simulate(3, 1, &logalpha, 0.0, 1.0);
	   if(print_to_screen)
	    {
	      cout << " : log(alpha) = "<< log(logalpha); 
	    }
 
	  if(log(logalpha) < loglike)
	    {
	      if(print_to_screen)
		cout << " : accepted" << endl;
	      accepted = 1;
	    }
	  else
	    {
	      if(print_to_screen)
		cout << " : rejected" << endl;
	      accepted = 0;
	    }
	}
    }
  if(print_output)
    {
      str << iteration_id << " Sampled ";
      if(current_machines.size() == 2)
	{
	  str << "Merge ";
	}
      else
	{
	  str << "Split ";
	}
      for(i = 0; i < current_machines.size(); i++)
	{
	  str << current_machines[i]->uid << " " << current_machines[i]->parameters->q << " "<< current_machines[i]->parameters->alpha << " " << current_machines[i]->parameters->omega << " ";
	}
      for(i = 0; i < proposed_machines.size(); i++)
	{
	  str << proposed_machines[i]->uid << " " << proposed_machines[i]->parameters->q << " "<< proposed_machines[i]->parameters->alpha << " " << proposed_machines[i]->parameters->omega << " ";
	}
      if(proposed_machines.size() == 0)
	{
	  str << "0 0 0 0 ";
	}
      str << loglike << " " << accepted << endl;
      
    }

  if( (accepted) && (record_all_acceptance))
    {
      str << iteration_id << " Trace ";
      if(current_machines.size() == 2)
	{
	  str << "Merge ";
	}
      else
	{
	  str << "Split ";
	}
      for(i = 0; i < current_machines.size(); i++)
	{
	  str << current_machines[i]->uid << " " << current_machines[i]->parameters->q << " "<< current_machines[i]->parameters->alpha << " " << current_machines[i]->parameters->omega << " ";
	}
      for(i = 0; i < proposed_machines.size(); i++)
	{
	  str << proposed_machines[i]->uid << " " << proposed_machines[i]->parameters->q << " "<< proposed_machines[i]->parameters->alpha << " " << proposed_machines[i]->parameters->omega << " ";
	}
      str << loglike << " " << accepted << endl;
      
    }
  num_tries++;
  num_accepted += accepted;
  return accepted;
}

double MachineMergeSplit::ImmuneOffLambdaRatio(vector<Machine*>& proposed_machines, vector<Machine*>& current_machines, global_par_obj* global_pars, int* immune, double& middle_multiplier, vector<Machine*>& ptr)
{
  double loglike = 0.0;
  double offlambda[3];
  double offofftrans[3];
  double offseq[3];
  double logprodgamma[3];
  int mask[3] = {1,1,1};
  int i, j;
  int end_trunc;
  int T;
  double offstatecount;
  ptr.push_back(proposed_machines[0]);
  if(current_machines.size() == 2)
    {
      ptr.push_back(current_machines[1]);
      middle_multiplier = -1.0;
    }
  else
    {
      ptr.push_back(proposed_machines[1]);
      middle_multiplier = 1.0;
    }
  ptr.push_back(current_machines[0]);

  for(j = 0; j < 3; j++)
    {
      offlambda[j] = ptr[j]->parameters->off_lambda;

      if(offlambda[j] < MIN_LOG_P)
	offlambda[j] = MIN_LOG_P;

      offofftrans[j] = ptr[j]->parameters->NumOffOffTrans();
      offseq[j] = ptr[j]->parameters->NumOffSeq();
      immune[j] = ptr[j]->parameters->Is_immune();     
    }

  //log(exp(-N_i*lambda_i))
  loglike = (offseq[2]*offlambda[2]) - (middle_multiplier*offseq[1]*offlambda[1]) - (offseq[0]*offlambda[0]);
  // Sum_i log(lambda_i)
  loglike += offofftrans[0]*log(offlambda[0]) + (middle_multiplier)*offofftrans[1]*log(offlambda[1]) - offofftrans[2]*offlambda[2];

  //now get the runs of consecutive states
  for(j = 0; j < 3; j++)
    { 
      offstatecount = 0;
      logprodgamma[j] = 0.0;
      T=ptr[j]->parameters->T;
      end_trunc = ptr[j]->parameters->Deathdate();
      if(end_trunc >= T)
	end_trunc = T-1;
      for(i = ptr[j]->parameters->Birthdate(); i < end_trunc; i++)
	{
	  if(ptr[j]->parameters->State(i) != 0)
	    {
	      logprodgamma[j] += LogGamma(offstatecount + 1.0);
	      offstatecount = 0;
	    }
	  else
	    {
	      if(ptr[j]->parameters->State(i+1) == 0)
		offstatecount += 1.0;
	    }
	}
      logprodgamma[j] += LogGamma(offstatecount + 1.0);
    }
  loglike += logprodgamma[2] - (middle_multiplier)*logprodgamma[1] - logprodgamma[0];

  //adding in the immunity stuff;
  loglike += BernoulliRatio(immune, mask, middle_multiplier, global_pars->immune_rate);

  return loglike;
}

double MachineMergeSplit::SinglePathPoissonLoglike(int t,  vector<Machine*>& ptr, double middle_multiplier )
{
  //gsl_sf_lnchoose (n,k)
  double rate[3];
  double y[3];
  double loglike;
  int i;
 
  for(i = 0; i < 3; i++)
    {
      //Poisson Counts
      y[i] = ptr[i]->data->GetTotalCount(t);
      //rates
      rate[i] = ptr[i]->parameters->GetRate(t);
    }

  loglike = rate[2] - (middle_multiplier)*(rate[1]) - rate[0];
  loglike += LogGamma(y[2]+1.0) - (middle_multiplier)*LogGamma(y[1]+1.0) - LogGamma(y[0]+1.0);
  if(rate[0] > 0)
    {
      loglike += y[0]*log(rate[0]);
    }
  else
    {
      if(y[0] != 0)
	loglike = -1.0*INFINITY;
    }
  if(rate[1] > 0)
    {
      loglike += (middle_multiplier)*y[1]*log(rate[1]);
    }
  else
    {
      if(y[1] != 0)
	loglike = -1.0*middle_multiplier*INFINITY;
    }
  if(rate[2] > 0)
    {
      loglike -= y[2]*log(rate[2]);
    }
  else
    {
      if(y[2] != 0)
	loglike = -1.0*INFINITY;
    }

  return loglike;

}

double MachineMergeSplit::SinglePathTransLoglike(int t, vector<Machine*>& ptr, global_par_obj* global_pars, int* immunity, double middle_multiplier)
{

  double loglike;
  int i;
  int survival_mask[3];
  int is_deathdate[3];
  double prob[3];
  double transition_ratio = 1.0;
  int state[2];
  double rho, gamma, nu;
  int T = ptr[0]->parameters->T;

  for(i = 0; i < 3; i++)
    {
      state[0] = ptr[i]->parameters->State(t);     
      is_deathdate[i] = (ptr[i]->parameters->Deathdate() != t);
      //any states after the deathdate will be masked out because they are equal to off
      survival_mask[i] = (state[0] != 0) * (immunity[i] != 0);
      //if we are not at either the deathdate or the end of the array
      //get the appropriate transition probability
      if( !( is_deathdate[i] || (t == T) ))
	{
	  state[1] = ptr[i]->parameters->State(t+1);
	  rho = ptr[i]->parameters->rho[state[0]];
	  gamma= ptr[i]->parameters->gamma[state[0]];
	  nu = ptr[i]->parameters->nu[state[0]];
	  prob[i] = sinecurveprob_single(state[0], state[1], t, rho, gamma, nu, ptr[i]->parameters); 
	  if(prob[i] < MIN_LOG_P)
	    prob[i] = MIN_LOG_P;
	}
      else
	{
	  prob[i] = 1.0;
	}
    }
  
  loglike = BernoulliRatio(is_deathdate, survival_mask, middle_multiplier, global_pars->survival_rate[t]);
  loglike += log(prob[0]) + middle_multiplier*log(prob[1]) - log(prob[2]);

  return loglike;

}

double MachineMergeSplit::BernoulliRatio(int* epsilons, int* mask, double middle_multiplier, double p)
{
  double loglike = 0.0;
  if(p < MIN_LOG_P)
    p = MIN_LOG_P;
 
  loglike += (double)(((epsilons[0])*(mask[0]) +(middle_multiplier)*(epsilons[1])*(mask[1]) - (epsilons[2])*(mask[2])))*log(p);
  loglike += (double)(((1-epsilons[0])*(mask[0]) + (middle_multiplier)*(1-epsilons[1])*(mask[1]) - (1-epsilons[2])*(mask[2])))*log(1.0-p);

  return loglike;
}

double MachineMergeSplit::LogLikelihoodRatio( vector<Machine*>& proposed_machines, vector<Machine*>& current_machines, global_par_obj* global_pars, int print_to_screen, int iteration_id)
{
  vector<Machine*> ptr;
  int immunity[3];
  double middle_multiplier;
  double loglike;
  double offstate_ratio;
  double transition_ratio;
  double poisson_ratio;
  int i;
  ptr.clear();

  offstate_ratio = ImmuneOffLambdaRatio( proposed_machines, current_machines, global_pars, immunity, middle_multiplier, ptr);
  if(print_to_screen)
    {
      cout << iteration_id<< " Off state and immunity log-ratio = " << offstate_ratio << endl;
    }
  poisson_ratio = 0.0;
  transition_ratio = 0.0;

  for(i = 0; i < ptr[0]->parameters->T; i++)
    {
      //we will take care of this one in the count updater.
      //poisson_ratio += SinglePathPoissonLoglike(i, ptr, middle_multiplier);
      transition_ratio += SinglePathTransLoglike(i, ptr, global_pars, immunity, middle_multiplier);
    }
  if(print_to_screen)
    {
      cout << iteration_id << " Transition probability log-ratio = " << transition_ratio << endl;
      //cout << iteration_id << " Poisson probability log-ratio = " << poisson_ratio << endl;
    }
  loglike = offstate_ratio + poisson_ratio + transition_ratio;
  if(print_to_screen)
    {
      cout << iteration_id << " Full Log-likeihood ratio for proposal vs. current = " << loglike << endl;
    }

  return loglike;

}


void MachineMergeSplit::PrintAcceptance(ostream& str)
{
  str << num_tries << " ";
  if(num_tries > 0)
    str << (double)(num_accepted)/(double)(num_tries);
  else
    str << 0;
  str << endl;
}
