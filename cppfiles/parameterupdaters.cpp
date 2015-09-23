#include "../parameterupdaters.h"

/*=============Main ParameterUpdater Class=====================*/
ParameterUpdater::ParameterUpdater(int &N,
				   double *&hpar1, 
				   DistributionList*& dists, 
				   double &propsd)
{
  int i;
  hyperpar = new double[N];
  current_value = new double[1];
  proposal_value = new double[1];
  dim = 1;
  hp_allocated = 1;
  num_hyperpars = N;
  for(i = 0; i < N; i++)
    hyperpar[i] = hpar1[i];
  distributions = dists;  
  proposal_sd = propsd;
}



ParameterUpdater::ParameterUpdater(int &N,
				   double *&hpar1,  
				   DistributionList*& dists)
{
  int i;
  current_value = new double[1];
  proposal_value = new double[1];
  dim = 1;
  hyperpar = new double[N];
  num_hyperpars = N;
  hp_allocated = 1;
  for(i = 0; i < N; i++)
    hyperpar[i] = hpar1[i];
  distributions = dists; 
  proposal_sd = 0;
}

ParameterUpdater::ParameterUpdater(DistributionList*& dists)
{
  current_value = new double[1];
  proposal_value = new double[1];
  dim = 1;
  hp_allocated = 0;
  num_hyperpars = 0;
  proposal_sd = 0;
  distributions = dists;
}

ParameterUpdater::ParameterUpdater(int &N, DistributionList*& dists)
{
  current_value = new double[N];
  proposal_value = new double[N];
  dim=N;
  hp_allocated= 0;
  num_hyperpars = 0;
  proposal_sd = 0;
  distributions = dists;
}

void ParameterUpdater::PrintHyperpars(ostream& str)
{
  int i;
  str << par_name << ":";
  if(num_hyperpars > 0)
  {
  for(i = 0 ; i < num_hyperpars ; i++)
	str << hyperpar[i] << " ";
  }
  else
  {
   str << "no hyperparameters";
  }
  str << " Proposal StdDev: " << proposal_sd; 
  if(proposal_sd == 0)
    {
      str << "(Gibbs Step)";
    }
  str << endl;
}

void ParameterUpdater::SetValidStartValue(curr_par_obj* current_values, global_par_obj* global_pars)
{
  int i;
  if(!(l_is_inf && u_is_inf))
    {
      SetCurrentValue(current_values, global_pars);
      for(i = 0; i < dim; i++)
	{
	  proposal_value[i] = current_value[i];
	  if(!(l_is_inf))
	    {
	      if(current_value[i] < lb)
		proposal_value[i] = lb + BETA_MIN;		
	    }
	  if(!(u_is_inf))
	    {
	      if(current_value[i] > ub)
		proposal_value[i] = ub - BETA_MIN;
	    }	  
	}
      SetCurrentProposed(current_values, global_pars);
    }
}

short int ParameterUpdater::MCMC_Step(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars, int iteration_id, int print_step )
{
  double mcmc_ratio;
  double logalpha;
  is_accepted = 0;
  int i;
  PRINT_MCMC_STEP=print_step;
  if(PRINT_MCMC_STEP)
    {
      if(iteration_id >= 0)
	cout << iteration_id << " ";
      cout << par_name;
    }
  SetCurrentValue(current_values, global_pars);
  Propose(current_values, data, global_pars);
  if(PRINT_MCMC_STEP)
    {
      cout << ": current = " << current_value[0] << ", proposed = " << proposal_value[0] << " ";
    }
  if(proposal_value[0] != current_value[0])
    {
      mcmc_ratio = Logposterior(proposal_value[0], current_values, data, global_pars) 
	- Logposterior(current_value[0], current_values, data, global_pars)
	+ Logjump(current_value[0], proposal_value[0], current_values, data, global_pars) 
	- Logjump(proposal_value[0], current_value[0], current_values, data, global_pars);
      
      distributions->Simulate(3, 1, &logalpha, 0.0, 1.0);
      if(PRINT_MCMC_STEP)
	{
	  cout << " ratio = " << mcmc_ratio << ", logalpha = " << log(logalpha);
	}
      if(log(logalpha) < mcmc_ratio)
	{
	  is_accepted += 1;
	  SetCurrentProposed(current_values, global_pars);
	  if(PRINT_MCMC_STEP)
	    cout << " accepted " << endl;
	}
      else
	{
	  if(PRINT_MCMC_STEP)
	    cout << " rejected " << endl;
	}      
    }
  else
    {
      if(PRINT_MCMC_STEP)
	{
	  cout << " no move " << endl;
	}
    }
  
  return(is_accepted);
}

ParameterUpdater::~ParameterUpdater()
 {
   if(hp_allocated)
     delete hyperpar;
   delete proposal_value;
   delete current_value;
 }


/*======================Individual Parameters=======================*/

/*-----------*/
OffLambdaUpdater::OffLambdaUpdater(int N, double* hpar1, DistributionList* dists):ParameterUpdater(N, hpar1, dists)
{
  par_name="OffLambda";
  lb = 0;
  ub =0;
  l_is_inf = 0;
  u_is_inf = 1;
}

void OffLambdaUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{   
  double kappa, theta;
  kappa = hyperpar[0];
  theta = hyperpar[1];
  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid;
    }
  SetCurrentValue(current_values, global_pars);

  distributions->Simulate(4,1,proposal_value,(double)(current_values->NumOffOffTrans())+kappa,(double)1.0/(((double)(current_values->NumOffSeq()))+ (1.0/theta) ));
      // if( (proposal_value[0] < hyperpar[0]) || (proposal_value[0] > hyperpar[1]))
      //proposal_value[0] = current_value[0];
 
}


/*double OffLambdaUpdater::Logjump(double val1,double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  return( trunc_norm_logjump(val1, val2, hyperpar[0], hyperpar[1],proposal_sd, distributions) );
}
double OffLambdaUpdater::Logprior(double value)
{
  return( distributions->LogDensity(3,value, hyperpar[0], hyperpar[1]));
}

double OffLambdaUpdater::Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double logpost;
  int num_off_states;
  int num_off_sequences;
  int in_off_sequence;
  int i, end_trunc;
  end_trunc = current_values->Deathdate();
  if(end_trunc >= current_values->T)
    {
      end_trunc = current_values->T -1;
    }
  logpost = Logprior(value);
  in_off_sequence = 0;
  num_off_sequences = 0;
  num_off_states = 0;
  for(i = current_values->Birthdate(); i <= end_trunc; i++)
    {
      if(current_values->State(i) > 0)
	{
	  if(in_off_sequence)
	    {
	      in_off_sequence = 0;
	      num_off_sequences++;
	    }
	}
      else
	{
	  num_off_states++;
	  if(!in_off_sequence)
	    in_off_sequence = 1;
	}
    }
  if(in_off_sequence)
    num_off_sequences++;
  
  logpost += (float)(num_off_states)*log(value) - value*(float)(num_off_sequences) - (float)(num_off_sequences)*log(1.0-exp(-1.0*value));
  
  return(logpost);
}*/

double OffLambdaUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  SetCurrentValue(current_values, global_pars);
  return(Logprior(current_value[0]));
}

/*-----------*/
QUpdater::QUpdater(int N, double* hpar1, DistributionList* dists):ParameterUpdater(N, hpar1, dists)
{
 par_name = "q";
 lb = hyperpar[0];
 ub = hyperpar[1];
 l_is_inf = 0;
 u_is_inf = 0;
}

void QUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
 {
   int i, end_trunc;
   int num_obs;
   int num_on;
   int sum_y;
   double sumalpha;
   double newpsi;
   double newphi;
   int iter_try=0;
   double l, u, kappa, theta;
   l = hyperpar[0];
   u = hyperpar[1];
   kappa = hyperpar[2];
   theta = hyperpar[3];

   if(PRINT_MCMC_STEP)
     {
       cout << "_" << data->uid;
     }
   

   if(current_values->Deathdate() >= current_values->T)
     end_trunc = current_values->T -1;
   else
     end_trunc = current_values->Deathdate();
   sum_y = 0;
   for(i = current_values->Birthdate(); i <= end_trunc; i++)
     sum_y += data->GetTotalCount(i);
   sumalpha = sum_alpha_k(current_values->Birthdate(), 
			  current_values->Deathdate(), 
			  current_values, 
			  current_values->alpha, 
			  current_values->omega);
   num_obs = end_trunc - current_values->Birthdate() + 1;
   num_on = 0;
   for(i = current_values->Birthdate(); i <= end_trunc; i++)
     {
       if(current_values->State(i) != 0)
	 num_on++;
     }
   // cout << "q_pripsi = " << hyperpar1 << " : q_priphi = " << hyperpar2;
   newpsi = kappa + (double)sum_y;
   newphi = 1.0/((1.0/theta) +(double) num_on + sumalpha);
   //cout << " : q_newpsi = " << newpsi << " : q_newphi = " <<newphi << endl;   
   distributions->Simulate(4,1,proposal_value, newpsi, newphi);    
   if((proposal_value[0] < l) || (proposal_value[0] > u))
     {
       proposal_value[0] = current_value[0];
     }
 }

double QUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  double loglike; 
  double l, u, kappa, theta;
  lb = hyperpar[0];
  ub = hyperpar[1];
  kappa = hyperpar[2];
  theta = hyperpar[3];

  loglike = distributions->LogDensity(4,current_values->q, kappa, theta) - log(gsl_cdf_gamma_Q(l, kappa, theta));
  if( (current_values->q < l) ||  (current_values->q > u) )
    loglike = -1.0*INFINITY;
  return loglike;
}

/*-----------*/



AlphaUpdater::AlphaUpdater(int N, double* hpar1, DistributionList* dists, double propvar):ParameterUpdater(N, hpar1, dists, propvar)
{
  par_name = "alpha";
  lb = 0.0;
  ub = 1.0;
  l_is_inf= 0;
  u_is_inf= 0;
}

double AlphaUpdater::Logprior(double value)
{
  double alpha, beta;
  alpha= hyperpar[0];
  beta = hyperpar[1];
  return(distributions->LogDensity(2, value, alpha, beta));
}

double AlphaUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  SetCurrentValue(current_values, global_pars);
  return(Logprior(current_value[0]));
}

double AlphaUpdater::Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double sumalpha = sum_alpha_k(current_values->Birthdate(), 
				current_values->Deathdate(), 
				current_values, 
				value);
  double logpost;
  logpost =  Logprior(value);
  if(logpost != (-1,0*INFINITY))
    {
      logpost += log_prod_y_alpha(data, current_values, value, current_values->omega, current_values->Birthdate(), current_values->Deathdate());
      logpost += -1.0*current_values->q*(current_values->omega -1.0)*sumalpha;
    }
  return(logpost);
}

void AlphaUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid;
    }
  

  SetCurrentValue(current_values, global_pars);
  proposal_value[0] = trunc_norm_proposal(current_value[0], 
				       proposal_sd, 
				       BETA_MIN, 
				       BETA_MAX,
				       distributions);
}

double AlphaUpdater::Logjump(double val1,double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  return( trunc_norm_logjump(val1, val2, 0.0, 1.0,proposal_sd, distributions) );
}


/*-------------*/
OmegaScaledUpdater::OmegaScaledUpdater(int N, double* hpar1, DistributionList* dists, double propvar):ParameterUpdater( N, hpar1, dists, propvar)
{
 par_name = "omega";
 lb = hyperpar[0];
 ub = hyperpar[1];
 l_is_inf = 0;
 u_is_inf = 0;
}

void OmegaScaledUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double l, u;
  l = hyperpar[0];
  u = hyperpar[1];
  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid;
    }
  
  SetCurrentValue(current_values, global_pars);
  proposal_value[0] = trunc_norm_proposal(current_value[0], 
				       proposal_sd, 
				       l, 
				       u,
				       distributions);
}

double OmegaScaledUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  SetCurrentValue(current_values, global_pars);
  return(Logprior(current_value[0]));
}

double OmegaScaledUpdater::Logprior(double value)
{
  double l, u, alpha, beta;
  l = hyperpar[0];
  u = hyperpar[1];
  alpha = hyperpar[2];
  beta = hyperpar[3];
  return(distributions->LogDensity(7,value,l, u, alpha, beta)); 
}

double OmegaScaledUpdater::Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double logpost;
  if(value < 1.0)
    logpost = -1*INFINITY;
  else
    {
      logpost = Logprior(value);
      if(logpost != (-1.0*INFINITY))
	{
	  logpost = -1*current_values->q * value;
	  logpost *= sum_alpha_k(current_values->Birthdate(), current_values->Deathdate(), current_values, current_values->alpha);
	  logpost += log_prod_y_alpha(data, current_values, current_values->alpha,value, current_values->Birthdate(), current_values->Deathdate());
	}
    }
  return(logpost);
}

double OmegaScaledUpdater::Logjump(double val1, double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double l, u;
  l = hyperpar[0];
  u = hyperpar[1];
  return( trunc_norm_logjump(val1, val2, l, u, proposal_sd, distributions) );
}

/*-------------*/

RhoUpdater::RhoUpdater(short int stateval,int N, double* hpar1, DistributionList* dists, double propvar):ParameterUpdater(N, hpar1, dists, propvar)
{
  state = stateval;
  if(state == 0)
    par_name = "rho.off";
  if(state == 1)
    par_name = "rho.spike";
  if(state == 2)
    par_name = "rho.decay";

  lb = 0.0;
  ub = 1.0;
  l_is_inf = 0;
  u_is_inf = 0;
}

double RhoUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  SetCurrentValue(current_values, global_pars);
  return(Logprior(current_value[0]));
}

void RhoUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid;
    }

  SetCurrentValue(current_values, global_pars);
  proposal_value[0] = trunc_norm_proposal(current_value[0], 
				       proposal_sd, 
				       BETA_MIN, 
				       BETA_MAX,
				       distributions);
}

double RhoUpdater::Logprior(double value)
{
  double alpha, beta;
  alpha = hyperpar[0];
  beta = hyperpar[1];
  return(distributions->LogDensity(2, value, alpha, beta));
}

double RhoUpdater::Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double logpost;
  logpost = Logprior(value);
  if(logpost != (-1.0*INFINITY))
    logpost += NextStateProb(state, 0, value, current_values, global_pars);
  return(logpost);
}

double RhoUpdater::Logjump(double val1, double val2, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  return( trunc_norm_logjump(val1, val2, 0.0, 1.0, proposal_sd, distributions));
}

/*-------------*/
NuUpdater::NuUpdater(short int stateval, int N, double* hpar1, DistributionList* dists, double propvar):ParameterUpdater(N, hpar1, dists, propvar)
  {
    state = stateval;
    if(state == 0)
      par_name = "nu.off";
    if(state == 1)
      par_name = "nu.spike";
    if(state == 2)
      par_name = "nu.decay";
    lb = hyperpar[0];
    ub = hyperpar[1];
    l_is_inf =0;
    u_is_inf =0;

  }

double NuUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  SetCurrentValue(current_values, global_pars);
  return(Logprior(current_value[0]));
}

void NuUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double l, u;
  l = hyperpar[0];
  u = hyperpar[1];
  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid;
    }

  SetCurrentValue(current_values, global_pars);
  proposal_value[0] = trunc_norm_proposal(current_value[0], 
				       proposal_sd, 
				       l, 
				       u,
				       distributions);
}

double NuUpdater::Logprior(double value)
{
  double l, u, alpha, beta;
  l = hyperpar[0];
  u = hyperpar[1];
  alpha = hyperpar[2];
  beta = hyperpar[3];
  return(distributions->LogDensity(7, value,l, u, alpha, beta));
}

double NuUpdater::Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double logpost;
  logpost = Logprior(value);
  if(logpost != (-1.0*INFINITY))
    logpost += NextStateProb(state, 2, value, current_values, global_pars);
  return(logpost);
}

double NuUpdater::Logjump(double val1, double val2,  curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double l, u;
  l = hyperpar[0];
  u = hyperpar[1];
  return( trunc_norm_logjump(val1, val2, l, u,proposal_sd, distributions));
}
 
/*-------------*/

GammaUpdater::GammaUpdater(short int stateval,int N, double* hpar1, DistributionList* dists):ParameterUpdater(N, hpar1, dists)
{
  state = stateval;
  if(state == 0)
    par_name = "gamma.off";
  if(state == 1)
    par_name = "gamma.spike";
  if(state == 2)
    par_name = "gamma.decay";
  lb = 0.0;
  ub = 1.0;
  l_is_inf = 0;
  u_is_inf = 0;

}

double GammaUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  double loglike;
  double alpha, beta;
  alpha = hyperpar[0];
  beta = hyperpar[1];
  loglike = distributions->LogDensity(2, current_values->gamma[state], alpha, beta);
  return loglike;
}

void GammaUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  
  double plusstatecount = 0.0;
  double minusstatecount = 0.0;
  double alpha, beta;
  alpha = hyperpar[0];
  beta = hyperpar[1];
  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid;
    }

  GetPlusMinusState(state, current_values, plusstatecount, minusstatecount);
  
  distributions->Simulate(2,1,proposal_value, 
			  alpha + plusstatecount , 
			  beta + minusstatecount);
  proposal_value[0] = BetaConstrain(proposal_value[0]);
  
}

/*-------------------*/

StateUpdater::StateUpdater(int N, double* hpar1, DistributionList* dists):ParameterUpdater(N, hpar1, dists)
{
  temp_parameters = new curr_par_obj();
  is_numerator = 1;
  par_name = "state";
  l_is_inf = 1;
  u_is_inf = 1;
}

void StateUpdater::SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars)
{
  if(!(current_values->Is_immune()))
    {
      if( (current_values->State(index_selected) != 0) && ((short int)(proposal_value[0]) == 0)) 
	{
	  global_pars->num_nonoff_nonimmune_transitions[index_selected]--;
	}
      if( (((short int)(proposal_value[0])) != 0)  && ( current_values->State(index_selected) == 0  ))
	{
	  global_pars->num_nonoff_nonimmune_transitions[index_selected]++;
	}
    } 
  current_values->SetState(index_selected, (short int)(proposal_value[0]));
  if(temp_parameters->off_lambda != current_values->off_lambda)
    current_values->off_lambda = temp_parameters->off_lambda;
}

double StateUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  return(0.0);
}

double StateUpdater::Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double logpost;
  if(is_numerator == 1)
    {
      logpost = LogAplusB(prop_hfun1, prop_hfun2);
    }
  else
    {
      logpost = LogAplusB(curr_hfun1, curr_hfun2);
    }
  is_numerator = 1 - is_numerator;
  return(logpost);
}

void StateUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double val;
  int range;
  int end_range;
  short int original_state;
  short int states_to_choose[2];
  int pre_index_nondecay;
  int post_index_nondecay;
  int path_off_off_transitions;
  int path_num_off_runs;
  int j;
  double alphaval;
  double lambdaval;
  double spikemult;
  int beginend[4];
  int addsub[2];
  vector<double> prop1_trans;
  vector<double> prop1_pois;
  vector<double> prop2_trans;
  vector<double> prop2_pois;
  vector<double> curr_trans;
  vector<double> curr_pois;
  double s1, s2, s3;
  double kappa, theta;
  kappa = hyperpar[0];
  theta = hyperpar[1];
  curr_hfun1 = curr_hfun2 = prop_hfun1 = prop_hfun2 = 0.0;
  

  /*uniform assumption of weights */
  if(current_values->Deathdate() >= current_values->T)
    end_range = current_values->T-1;
  else
    end_range = current_values->Deathdate();
  //range = end_range - current_values->birthdate + 1;
  distributions->Simulate(3, 1, &val, (double)current_values->Birthdate(), (double)(end_range)+1.0);
  index_selected = (int)floor(val);
  if( (index_selected > end_range) || (index_selected < current_values->Birthdate()))
    {
      cout << "In StateUpdater::Logposterior : index selected is out of range" << endl;
    }
  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid << "_" << index_selected << "_" << "ct_" << data->GetTotalCount(index_selected);
    }
  //  if((index_selected != current_values->beginval) && (index_selected != current_values->endval))
  // {
  //	if((current_values->state[index_selected - 1] == 1) || (current_values->state[index_selected + 1] == 1))
  //cout << endl << " spike vicinty : "<< current_values->state[index_selected - 1] << " " << current_values->state[index_selected] << " " << current_values->state[index_selected + 1] << endl;
  //      }
  
  temp_parameters->SetFromCopy(current_values);
  //current_values->Print(0);
  //temp_parameters->Print(0);
  
  original_state = current_values->State(index_selected);
  current_value[0] = (double) original_state;
  
  PathToChange(index_selected, current_values, data, beginend);
  pre_index_nondecay = beginend[0];
  post_index_nondecay = beginend[1];
  path_off_off_transitions = beginend[2];
  path_num_off_runs = beginend[3];

  curr_hfun1 = EtaPropH(index_selected,pre_index_nondecay, post_index_nondecay, data, temp_parameters, global_pars, &(curr_trans), &(curr_pois));
  
  if(data->GetTotalCount(index_selected) > 0)
    {
      /* This may need to be adapted for transition states with prob 0 (eg spike->spike)  */
      if(temp_parameters->State(index_selected) == 1)
	proposal_value[0] = 2;
      else
	proposal_value[0] = 1;
      states_to_choose[0] = (short int)(proposal_value[0]);
      states_to_choose[1] = -1;
      temp_parameters->SetState(index_selected, (short int)(proposal_value[0]  ));
      prop_hfun1 = EtaPropH(index_selected, pre_index_nondecay, post_index_nondecay, data, temp_parameters, global_pars, &(prop1_trans), &(prop1_pois));
      if(prop_hfun1 == (-1.0*INFINITY))
	{
	  //cout << endl << prop_hfun1 << endl << (short int)prop_hfun1 ;
	  //cout << endl << "Cannot change state: " << endl << "st:"<< current_values->state[max(0,index_selected - 1)] << " " << original_state << " " << current_values->state[min(current_values->T, index_selected + 1)] << endl << "fl:"<< data->fl[max(0,index_selected - 1)] << " " << data->fl[index_selected] << " " << data->fl[min(current_values->T, index_selected + 1)] << endl;
	  proposal_value[0] = original_state;
	}
      prop_hfun2 = 0.0;
      curr_hfun2 = 0.0;
      temp_parameters->SetState(index_selected, original_state);
      
    }
  else
    {
      /* y = 0; choose between two other states */
      /*And do a lambda update as well*/
      prop_hfun2 = 0.0;
      curr_hfun2 = 0.0;
      j = 0;
      spikemult = 0.5;
      if(original_state != 0)
	{
	  states_to_choose[j] = 0;
	  j++;
	}
      if(original_state != 2)
	{
	  states_to_choose[j] = 2;
	  j++;
	}
      if(original_state != 1)
	{
	  states_to_choose[j] = 1;
	  j++;
          spikemult = 0.8;
	}
      distributions->Simulate(3, 1, &alphaval, 0.0, 1.0);
      if(alphaval < (spikemult))
	proposal_value[0] = states_to_choose[0];
      else
	proposal_value[0] = states_to_choose[1];
      temp_parameters->SetState(index_selected, (short int)(proposal_value[0]));
      
      //cout << "current state = " << (int)(original_state) << " num-off-off, num-off-seq = " << current_values->NumOffOffTrans() << ", " << current_values->NumOffSeq() << endl;
      //cout << "proposed state= " << (int)(proposal_value[0]) << " num-off-off, num-off-seq = " << temp_parameters->NumOffOffTrans() << ", " << temp_parameters->NumOffSeq() << endl;

      //if(temp_parameters->NumOffSeq() > 0)
      distributions->Simulate(4,1, &lambdaval,(double)(temp_parameters->NumOffOffTrans())+kappa,(double)1.0/(((double)(temp_parameters->NumOffSeq()))+ (1.0/theta) ));


	//else
	//lambdaval = current_values->off_lambda;
      //cout << "Lambda simulated = " << lambdaval << " current = " << current_values->off_lambda << endl;
      /*  if(  (lambdaval < hyperpar[0]) || (lambdaval > hyperpar[1]))
	{
	  lambdaval = current_values->off_lambda;
	  //proposal_value[0] = (double)(current_values->State(index_selected));
	}
      else
	{
	  temp_parameters->off_lambda = lambdaval;	  
	  }*/

      prop_hfun1 = EtaPropH(index_selected, pre_index_nondecay, post_index_nondecay, data, temp_parameters, global_pars, &(prop1_trans), &(prop1_pois));

      //add in log jump probabilities
      prop_hfun1 += LogSpikeProb(0.8, (short int)original_state, (short int)(proposal_value[0]));
      curr_hfun1 += LogSpikeProb(0.8, (short int)(proposal_value[0]), (short int)original_state);
      //add in the extra bits from the gamma proposal constants
      prop_hfun1 += LogGamma((double)((double)(temp_parameters->NumOffOffTrans()) + hyperpar[0]));
      curr_hfun1 += LogGamma((double)((double)(current_values->NumOffOffTrans()) + hyperpar[0])); 
      prop_hfun1 += (double)((double)(current_values->NumOffOffTrans())+kappa)*log((double)(current_values->NumOffSeq()) + (1.0/theta)   );
	  //curr_hfun1 += ScaledGammaLogConstant(hyperpar[0], hyperpar[1], (double)(current_values->NumOffOffTrans())+1.0 , 1.0/double(current_values->NumOffSeq()));
      
      curr_hfun1 += (double)((double)(temp_parameters->NumOffOffTrans())+kappa)*log((double)(temp_parameters->NumOffSeq()) + (1.0/theta) );
      //prop_hfun1 += ScaledGammaLogConstant(hyperpar[0], hyperpar[1], (double)(temp_parameters->NumOffOffTrans())+1.0 , 1.0/double(temp_parameters->NumOffSeq()));

      //last bit: add in the exp for the additional path variables      
      //get addition/subtraction values for the proposed path as opposed to original
      StateRunChanges(current_values, temp_parameters, index_selected, pre_index_nondecay, post_index_nondecay, addsub);
      prop_hfun1 += temp_parameters->off_lambda*( (double)path_num_off_runs + (double)(addsub[1]) );
      curr_hfun1 += current_values->off_lambda*((double)path_num_off_runs);
      
      prop_hfun1 += -1.0*((double)path_off_off_transitions + (double)(addsub[0]))* log(temp_parameters->off_lambda);
      curr_hfun1 += -1.0*((double)path_off_off_transitions)*log(current_values->off_lambda);	  
    }      
      
      
      /*
      for(j = 0; j < 2; j++)
	{
	  temp_parameters->SetState(index_selected, states_to_choose[j]);
	  if(j == 0)
	    prop_hfun1 = EtaPropH(index_selected, pre_index_nondecay, post_index_nondecay, data, temp_parameters, global_pars, &(prop1_trans), &(prop1_pois));
	  if(j == 1)
	    prop_hfun2 = EtaPropH(index_selected, pre_index_nondecay, post_index_nondecay, data, temp_parameters, global_pars, &(prop2_trans), &(prop2_pois));
	}
      if( ( prop_hfun1 == (-1.0*INFINITY)) || (prop_hfun2 == (-1.0*INFINITY)))
	{
	  //cout << endl<<"_ct_0s impossible jump " << endl;
	  //cout << endl << "_ct_0s Cannot change state at index: "<< index_selected << endl << "_ct_0s st:"<< current_values->state[max(0,index_selected - 1)] << " " << original_state << " " << current_values->state[min(current_values->T, index_selected + 1)] << endl << "_ct_0s fl:"<< data->fl[max(0,index_selected - 1)] << " " << data->fl[index_selected] << " " << data->fl[min(current_values->T, index_selected + 1)] << endl;
	  if(prop_hfun1 != (-1.0*INFINITY))
	    proposal_value[0] = states_to_choose[0];
	  else
	    if(prop_hfun2 != (-1.0*INFINITY))
	      proposal_value[0] = states_to_choose[1];
	    else
	      proposal_value[0] = original_state;
	}
      else
	{
	  distributions->Simulate(3, 1, &alphaval, 0.0, 1.0 + exp(prop_hfun2 - prop_hfun1));
	  if(alphaval < 1.0)
	    {
	      proposal_value[0] = states_to_choose[0];
	      curr_hfun2 = prop_hfun1;
	    }
	  else
	    {
	      proposal_value[0] = states_to_choose[1];
	      curr_hfun2 = prop_hfun2;
	    }
	}
      */
    
  return;
}

/*------------Update Birth, Death and Immunity Parameters for a single Machine ----------------*/

BirthDeathImmuneUpdater::BirthDeathImmuneUpdater(int N, double* hpar1, DistributionList* dists):ParameterUpdater(N, hpar1, dists)
{
  temp_parameters = new curr_par_obj();
  par_name = "Birth/Death/Immune";
  l_is_inf = 1;
  u_is_inf = 1;
}


void BirthDeathImmuneUpdater::SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars)
{
  int i;
  int found;
  int subind;
  temp_parameters->SetFromCopy(current_values);
  first_nonoff_transition = FirstLastRun(current_values, 1, 0);
  last_nonoff_transition = FirstLastRun(current_values, -1, 0);  
  sum_nonoff_transitions = 0.0;
  for(i = current_values->Birthdate(); i < current_values->Deathdate(); i++)
    if(current_values->State(i) !=0)
      {
	sum_nonoff_transitions += log(global_pars->survival_rate[i]);
      }
  current_value[0] = (double)(GetBDIIndex(current_values->Birthdate(), current_values->Deathdate(), current_values->Is_immune(), current_values->T));
}

double BirthDeathImmuneUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  return(0.0);
}

void BirthDeathImmuneUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{  /* proposal value is one of 0 through 5
        BDI 
    0 = 000
    1 = 010
    2 = 011
    3 = 100
    4 = 110
    5 = 111
  where B=0 means birth at the start of the observation window, B=1 means birth at 284
  D = 0 means death at the last non-zero transition, D = 1 means no death, eg death after time T
  I = 0 means not immune, I = 1 means immune
*/
  int i,x;
  vector<int> proposal_states;
  int state_selected;
  short int prop_birth, prop_death, prop_immune;
  double lambdaval;
  double kappa, theta;
  kappa = hyperpar[0];
  theta = hyperpar[1];

  if(PRINT_MCMC_STEP)
    {
      cout << "_" << data->uid;
    }



  for(i = 0; i < 3; i++)
    {
      if(((int)(current_value[0])) != i)
	proposal_states.push_back(i);
    }
  if(first_nonoff_transition > 284)
    {
      for(i = 3; i < 6; i++)
	{
	  if(current_value[0] != i)
	    proposal_states.push_back(i);
	}
    }
  distributions->Simulate(3,1,proposal_value,0.0, (double)(proposal_states.size()));
  state_selected = proposal_states[(int)(floor(proposal_value[0]))];
  proposal_value[0] = (double)(state_selected);
  BDIFromIndex(state_selected, last_nonoff_transition, current_values->T, &prop_birth, &prop_death, &prop_immune);
  temp_parameters->SetBirthDeathImmune(prop_birth, prop_death, prop_immune);

  
  /*Simulate a new lambda value*/
 
  lambdaval = current_values->off_lambda;

  distributions->Simulate( 4, 1, &lambdaval, (double)( (double)(temp_parameters->NumOffOffTrans())+kappa), (double)(1.0/((double)(temp_parameters->NumOffSeq()) + (1.0/theta )  )));
   

  temp_parameters->off_lambda = lambdaval;

  /*
    for(i = 0; i < 6; i++)
    likelihoods[i] = 0.0;
  
   //add in the immunity states
  likelihoods[0] += sum_nonoff_transitions + log(1-global_pars->immune_rate);
  likelihoods[1] += sum_nonoff_transitions + log(1-global_pars->immune_rate); 
  likelihoods[3] += sum_nonoff_transitions + log(1-global_pars->immune_rate);
  likelihoods[4] += sum_nonoff_transitions + log(1-global_pars->immune_rate);

  likelihoods[2] += log(global_pars->immune_rate);
  likelihoods[5] += log(global_pars->immune_rate);

  //add in the death if it occurs
  likelihoods[0] += log(1.0 - global_pars->survival_rate[last_nonoff_transition]);
  likelihoods[3] += log(1.0 - global_pars->survival_rate[last_nonoff_transition]);

  if(first_nonoff_transition < 284)
    {
      //birthday cannot be equal to 284
      prob_284_multiplier = 0.0;
    }
  else
    {
      //For birth=284 we do not include the prior off states
      
      likelihoods[0] += CondPoissonLoglike(first_nonoff_transition, current_values->off_lambda);
      likelihoods[1] += CondPoissonLoglike(first_nonoff_transition, current_values->off_lambda);
      likelihoods[2] += CondPoissonLoglike(first_nonoff_transition, current_values->off_lambda);
      if(first_nonoff_transition > 284)
	{
	  likelihoods[3] += CondPoissonLoglike(first_nonoff_transition-284, current_values->off_lambda);
	  likelihoods[4] += CondPoissonLoglike(first_nonoff_transition-284, current_values->off_lambda);
	  likelihoods[5] += CondPoissonLoglike(first_nonoff_transition-284, current_values->off_lambda);
	}
    }
  
  likelihood_sum = 0;
  for(i = 0; i < 6; i++)
    {
      likelihoods[i] = exp(likelihoods[i]);
      if(i > 2)
	likelihoods[i] *= prob_284_multiplier;
      likelihood_sum += likelihoods[i];
    }

  distributions->Simulate(3, 1, proposal_value, 0.0, likelihood_sum);
  i = 0;
  likelihood_sum = likelihoods[0];
  while(proposal_value[0] > likelihood_sum)
    {
      i++;
      likelihood_sum += likelihoods[i];
    }
  proposal_value[0] = (double)i;
*/
}

double BirthDeathImmuneUpdater::Logposterior(double value, curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  curr_par_obj* ptr;
  double deathdatemultiplier=0.0;
  double runcount;
  double loglike;
  double kappa = hyperpar[0];
  double theta = hyperpar[1];
  if(value == current_value[0])
    {
      ptr = current_values;
    }
  else
    {
      if(value == proposal_value[0])
	{
	  ptr = temp_parameters;
	}
      else
	{
	  cout << "BirthDeathImmuneUpdater Error: neither proposal nor current value entered into log posterior" << endl;
	}
    }
  if(ptr->Deathdate() < ptr->T)
    {
      deathdatemultiplier = global_pars->survival_rate[ptr->Deathdate()];
    }
  loglike = (ptr->Is_immune())*log(global_pars->immune_rate) + (1-ptr->Is_immune())*( log(1.0-global_pars->immune_rate) + log(1.0-deathdatemultiplier) + sum_nonoff_transitions );

  loglike += LogGamma((double)(ptr->NumOffOffTrans())+kappa)  - ((double)(ptr->NumOffOffTrans())+kappa)*log((double)(ptr->NumOffSeq()) + (1.0/theta) );
  
  runcount = (double)first_nonoff_transition - (double)(ptr->Birthdate());
  if( runcount > 0)
    loglike -= LogGamma(runcount + 1.0);
  runcount = (double)(ptr->Deathdate()) - (double)last_nonoff_transition;
  if( runcount > 0)
    loglike -= LogGamma(runcount + 1.0);

  return(loglike);
}

void BirthDeathImmuneUpdater::SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars)
{
  short int newbirth;
  short int newdeath;
  short int newimmune;
  
  /* proposal value is one of 0 through 5
        BDI 
    0 = 000
    1 = 010
    2 = 011
    3 = 100
    4 = 110
    5 = 111
  where B=0 means birth at the start of the observation window, B=1 means birth at 284
  D = 0 means death at the last non-zero transition, D = 1 means no death, eg death after time T
  I = 0 means not immune, I = 1 means immune
*/
  /*First undo the existing counts for this machine*/
  global_pars->UpdateGlobalImmunityAndSurvivalCounts(current_values, -1);
  BDIFromIndex((int)(proposal_value[0]), last_nonoff_transition, temp_parameters->T, &newbirth, &newdeath, &newimmune);
  current_values->SetBirthDeathImmune(newbirth, newdeath, newimmune);
  /*Now add in the new values for the current machine */
  global_pars->UpdateGlobalImmunityAndSurvivalCounts(current_values, 1);
  if(temp_parameters->off_lambda != current_values->off_lambda)
    current_values->off_lambda = temp_parameters->off_lambda;
}

BirthDeathImmuneUpdater::~BirthDeathImmuneUpdater()
{
  delete temp_parameters;
}
/*Gibbs steps for global parameter values: survival rate and immunity rate*/

ImmuneRateUpdater::ImmuneRateUpdater(int N, double* hpar1, DistributionList* dists):ParameterUpdater(N, hpar1, dists)
{
  par_name = "ImmuneRate";
  lb = 0.0;
  ub = 1.0;
  l_is_inf = 0;
  u_is_inf = 0;
}

void ImmuneRateUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  double alpha=hyperpar[0];
  double beta =hyperpar[1];
  distributions->Simulate(2,1,proposal_value,alpha + global_pars->num_immune, beta + global_pars->num_not_immune); 
  proposal_value[0] = BetaConstrain(proposal_value[0]);
}

double ImmuneRateUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  double loglike;
  double alpha=hyperpar[0];
 double beta =hyperpar[1];
  loglike = distributions->LogDensity(2, global_pars->immune_rate, alpha, beta);
  return loglike;
}

SurvivalRateUpdater::SurvivalRateUpdater(int N, DistributionList* dists):ParameterUpdater(N, dists)
{
  par_name="SurvivalRate";
  lb = 0.0;
  ub = 1.0;
  l_is_inf = 0;
  u_is_inf = 0;

}

void SurvivalRateUpdater::SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars)
{
  int i;
  for(i = 0; i < dim; i++)
    current_value[i] = global_pars->survival_rate[i];
}

void SurvivalRateUpdater::SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars)
{
  int i;
  for(i = 0; i < dim; i++)
    global_pars->survival_rate[i] = proposal_value[i];
}

double SurvivalRateUpdater::GraphLogLikelihood(curr_par_obj* current_values, global_par_obj* global_pars)
{
  int i;
  double loglike = 0.0;
  double val;
  for(i = 0; i < dim; i++)
    {
      val = distributions->LogDensity(2, global_pars->survival_rate[i], BetaDistAlpha(global_pars->SurvivalRatePriorMean(i), global_pars->survival_hyper_temperature), BetaDistBeta(global_pars->SurvivalRatePriorMean(i), global_pars->survival_hyper_temperature));
      //cout << i << " " <<  global_pars->survival_rate[i] << " " << global_pars->SurvivalRatePriorMean(i) << " " << global_pars->survival_hyper_temperature << " " << BetaDistAlpha(global_pars->SurvivalRatePriorMean(i), global_pars->survival_hyper_temperature) << " " << BetaDistBeta(global_pars->SurvivalRatePriorMean(i), global_pars->survival_hyper_temperature) << " " << val << endl;
      loglike += val;
    }
  return(loglike);
}

void SurvivalRateUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  //Will update all of them at once in a line.
  int i;
  for(i = 0; i < global_pars->T; i++)
    {
      distributions->Simulate(2,1,&(proposal_value[i]), BetaDistAlpha(global_pars->SurvivalRatePriorMean(i),global_pars->survival_hyper_temperature) + global_pars->num_nonoff_nonimmune_transitions[i], BetaDistBeta(global_pars->SurvivalRatePriorMean(i),global_pars->survival_hyper_temperature) + global_pars->num_deaths[i]);
      proposal_value[i] = BetaConstrain(proposal_value[i]);
    } 
}


/*
void SurvivalHyperUpdater::SetCurrentValue(curr_par_obj* current_values, global_par_obj* global_pars)
{
  double v;
  distributions->Simulate(3, 1, &v, 0.0, 5.0);
  pos = (int)(floor(v));
  current_value[0] = global_pars->survival_hyper_mean_height[pos];
  val = global_pars->survival_hyper_mean_position[pos];
}

void SurvivalHyperUpdater::SetCurrentProposed(curr_par_obj* current_values, global_par_obj* global_pars)
{
  global_pars->survival_hyper_mean_height[pos] = current_value[0];
  global_pars->survival_hyper_mean_position[pos] = val;
}

void SurvivalHyperUpdater::Propose(curr_par_obj* current_values, mach_dat_obj* data, global_par_obj* global_pars)
{
  
}

*/
